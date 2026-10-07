# MCCCS Towhee (High-Performance Optimization Fork)

An optimized, high-performance fork of **[MCCCS Towhee 8.2.3](https://towhee.sourceforge.net)** (Monte Carlo for Complex Chemical Systems), a multi-ensemble Monte Carlo molecular simulation engine written in Fortran 77 and C.

This fork focuses on single-node acceleration, memory accessor elimination, vectorization, and multi-node MPI scaling while preserving numerical correctness within standard simulation tolerances (`1e-6` relative).

---

## Table of Contents

- [Overview](#overview)
- [Key Optimizations](#key-optimizations)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Build Instructions](#build-instructions)
  - [1. HPI Serial (Recommended Production Build)](#1-hpi-serial-recommended-production-build)
  - [2. High-Performance Vectorized + Ewald](#2-high-performance-vectorized--ewald)
  - [3. MPI Parallel Build](#3-mpi-parallel-build)
  - [4. Plain Serial Baseline](#4-plain-serial-baseline)
  - [5. Safe Compare (Strict Bitwise Mode)](#5-safe-compare-strict-bitwise-mode)
- [Running Simulations](#running-simulations)
- [Testing & Benchmarking](#testing--benchmarking)
  - [Regression Testing](#regression-testing)
  - [A/B Comparison Benchmarks](#ab-comparison-benchmarks)
- [Performance Summary](#performance-summary)
- [Roadmap & Future Improvements](#roadmap--future-improvements)
- [Documentation & References](#documentation--references)
- [License](#license)

---

## Overview

[MCCCS Towhee](https://towhee.sourceforge.net) simulates classical molecular systems across various ensembles (NVT, NPT, $\mu$VT, Gibbs, CBMC, etc.) using detailed force fields.

In original Towhee 8.2.3, simulation bottlenecks are dominated by:
1. Millions of accessor function calls (`twh_*`) across simulation steps.
2. Huge lookup tables (e.g., the 266 MB CBMC bend info table `twh_cbmcbendinfo_`).
3. Overhead in differential energy evaluations and minimum image calculations.
4. Redundant MPI barriers in replica exchange and broadcasting.

This optimization fork introduces **Hot Path Inlining (HPI)**, **zero-copy COMMON block bindings**, **vectorized SIMD inner loops**, and **streamlined MPI communication** to significantly accelerate simulation throughput.

---

## Key Optimizations

### 1. Hot Path Inlining (HPI) & Zero-Copy Accessors (Phases 1 & 2)
- Replaced hundreds of millions of scalar and array accessor dispatches (`twh_hmatrix_`, `twh_halfbox_`, `twh_mimage_`, `twh_nbcoeff_`, `twh_rcut_`, `twh_qqscale_`, etc.) with zero-overhead `COMMON /hot_*/` shared memory mappings between `globaldata.F` and compute kernels (`engtotal_globals.F`, `vtwobody.F`, `vcoulomb.F`, `mimage.F`).
- Eliminated 5-index arithmetic overhead in the 266 MB CBMC table (`twh_cbmcbendinfo_`) via direct memory access in `getcbangle.F` and `angle_b_distribution.F`, reducing lookup function calls from ~540k to only 5 during initialization.

### 2. MPI Modernization (Phase 3)
- Removed redundant `MPI_Barrier` calls before `MPI_Bcast` and `MPI_Finalize`.
- Replaced sequential `MPI_Send` fan-outs in replica exchange (`rex.c`) with non-blocking `MPI_Isend` + `MPI_Waitall` loops.
- Integrated thread-safe MPI initialization (`MPI_Init_thread(MPI_THREAD_FUNNELED)`).

### 3. Vectorization & SIMD Acceleration (Phase 4)
- Optimized compiler configurations (`-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize`).
- Added OpenMP SIMD directives (`!$OMP DECLARE SIMD`, `!$OMP SIMD REDUCTION`) to pair energy routines, Lennard-Jones evaluation (`vtwobody.F`), and Ewald reciprocal summations.

### 4. Ewald Summation Enhancements (Phase 5)
- Direct zero-copy access for Ewald parameters (`calp`).
- Vectorized inner atom loops with SIMD reduction across reciprocal space wavevectors in `recipsum.F`.

### 5. OpenMP Root Cause Fix
- Isolated and resolved a longstanding F77 `-fopenmp` $\rightarrow$ `-frecursive` issue in `globaldata.F` where un-`SAVE`d Ewald partial-sum arrays (`ssumi`, `ssuminew`, `ssumr`, `ssumrnew`) became stack-allocated and clobbered state across simulation steps.

---

## Repository Structure

```plaintext
MCCSTowhee_opt/
├── AGENTS.md                  # Detailed developer agent instructions and project overview
├── SUPERVISOR_REPORT.md       # Executive summary of milestones, performance, and roadmap
├── README.md                  # Repository documentation (this file)
├── openspec/
│   └── live_specs/
│       ├── optimization_plan.md           # Detailed Phases 1–5 spec and verification
│       └── phase6_additional_improvements.md # Phase 6+ roadmap (differential zero-copy, RNG, cell list)
├── tests/
│   ├── comparison_benchmark.py # Automated A/B performance and tolerance benchmark harness
│   └── SPC_E_WATER/           # Standard 500-step and 50,000-step validation system
└── src/
    └── towhee-8.2.3/          # Main build worktree
        ├── configure.ac       # Autotools build configuration
        ├── Makefile.am        # Top-level automake template
        ├── Examples/          # ~50 self-contained simulation test cases & validation suites
        ├── ForceFields/       # 87+ classical force field parameter definitions
        ├── Manual/            # Official HTML documentation and input reference
        ├── Utils/             # Conversion tools (e.g. pdb2towhee)
        └── Source/            # Fortran 77 / C engine source code
            ├── towhee.c       # Main entrypoint
            ├── mainloop.F     # Core Monte Carlo driver
            ├── energy.F       # Differential energy evaluations
            ├── engtotal.F     # Full O(N²) box energy calculations
            ├── engtotal_globals.F # HPI zero-copy COMMON mappings
            ├── hot_common.h   # COMMON block memory definitions
            ├── globaldata.F   # Engine data storage and accessors
            ├── vtwobody.F     # Non-bonded pairwise potentials
            ├── vcoulomb.F     # Coulombic and real-space electrostatics
            └── recipsum.F     # Reciprocal-space Ewald summation
```

---

## Prerequisites

- **Fortran Compiler**: GNU Fortran (`gfortran` 11+ recommended, GCC 11–13 tested).
- **C Compiler**: `gcc` (matches `gfortran`).
- **Build Utilities**: `GNU make` (3.82+), `autotools` (`autoconf`, `automake` only if modifying `configure.ac` / `Makefile.am`).
- **MPI** (Optional for parallel runs): OpenMPI (`mpicc`, `mpifort`) or MPICH.
- **Python**: Python 3.8+ (for benchmark and tolerance testing scripts).

---

## Build Instructions

All builds take place inside `src/towhee-8.2.3/`.

### 1. HPI Serial (Recommended Production Build)

The Hot Path Inlining build delivers the best serial performance with strict regression parity (`0 diff` vs reference outputs).

```bash
cd src/towhee-8.2.3
./configure
make -j$(nproc) -C Source towhee HOT_PATH_CPPFLAGS=-DHOT_PATH_INLINE
```

### 2. High-Performance Vectorized + Ewald

Enables AVX/SIMD vectorization and `-ffast-math` along with HPI for maximum single-node throughput:

```bash
cd src/towhee-8.2.3
./configure
make -j$(nproc) -C Source towhee HOT_PATH_CPPFLAGS=-DHOT_PATH_INLINE
```
*(Vector flags `-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize` are automatically set by `./configure` unless `--enable-safe-compare` is selected).*

### 3. MPI Parallel Build

Enables multi-node MPI replica exchange and parallel execution:

```bash
cd src/towhee-8.2.3
./configure --enable-mpi
make -j$(nproc) -C Source towhee HOT_PATH_CPPFLAGS=-DHOT_PATH_INLINE
```

### 4. Plain Serial Baseline

Standard baseline build (without HPI optimizations):

```bash
cd src/towhee-8.2.3
./configure FFLAGS="-g -O2 -ffixed-line-length-none"
make -j$(nproc) -C Source towhee
```

### 5. Safe Compare (Strict Bitwise Mode)

Disables vectorization and fast-math optimizations for strict bit-for-bit reproducible cross-platform comparisons and bounds checking:

```bash
cd src/towhee-8.2.3
./configure --enable-safe-compare
make -j$(nproc) -C Source towhee
```

---

## Running Simulations

Towhee requires a configuration file named `towhee_input` in the working directory (and optional restart files `towhee_initial`).

### Serial Execution
```bash
cd /path/to/simulation/directory
/path/to/src/towhee-8.2.3/Source/towhee
```

### MPI Execution
```bash
cd /path/to/simulation/directory
mpirun -np 4 /path/to/src/towhee-8.2.3/Source/towhee
```

### Primary Output Files
- Standard Output (`towhee.out` or terminal): Energetics, step acceptance rates, block averages.
- `towhee_backup` / `towhee_final`: Checkpoint and restart files.
- `box_*_step_*.pdb` / `towhee_movie`: Coordinate trajectory files (if enabled).

---

## Testing & Benchmarking

### Regression Testing

Towhee includes ~50 regression tests in `src/towhee-8.2.3/Examples/`. Each example contains a `towhee_input`, reference `answer_current`, and a `do_test` test runner:

```bash
cd src/towhee-8.2.3/Examples/SPC_E_Water
../../do_test
```

To run all serial tests in batch:
```bash
cd src/towhee-8.2.3
./run_serial_test
```

### A/B Comparison Benchmarks

Use `tests/comparison_benchmark.py` to compare execution speed and verify numerical tolerance (`1e-6` relative error) between an original binary and the optimized binary:

```bash
python3 tests/comparison_benchmark.py \
  --original /path/to/original/towhee \
  --optimized src/towhee-8.2.3/Source/towhee \
  --run-dir tests/SPC_E_WATER \
  --threads 1 \
  --repeats 5 \
  --warmup 2 \
  --check-tolerance 1e-6
```

---

## Performance Summary

Measured performance improvements on native Linux (`gfortran` 11, `-O3 -march=native -ffast-math`):

| Benchmark System | Baseline Runtime | Optimized (HPI + SIMD) | Speedup | Numerical Parity |
| :--- | :---: | :---: | :---: | :---: |
| `tests/SPC_E_WATER` | 6.648 s | 5.799 s | **1.15×** | Strict `0 diff` |
| `src/Gibbs/EPM_VLCC` | 0.894 s | 0.809 s | **1.11×** | `1e-6` Pass |
| `tests/TraPPE_Pentane` | 0.554 s | 0.529 s | **1.05×** | `1e-6` Pass |
| `src/TMMC/SPC_E_WATER` | 1.069 s | 1.048 s | **1.02×** | `1e-6` Pass |

*Note: On small test systems (<1,000 atoms), runtimes are heavily dominated by CBMC angle tables and random number generation. Larger systems with heavy Ewald reciprocal sums or long trajectories exhibit higher speedups.*

---

## Roadmap & Future Improvements

Active plans are detailed in [`openspec/live_specs/phase6_additional_improvements.md`](file:///home/sharir/MCCSTowhee_opt/openspec/live_specs/phase6_additional_improvements.md):

- **P6 Differential Zero-Copy**: Direct flat coordinate indexing (`c_coord_real`) and precomputed atom arrays (`atomnum_arr`) in `energy.F`, bypassing `twh_coordinates_` (est. 5–12% gain).
- **P7 High-Speed RNG**: Upgrading default RNG from `RANLUX` to `DX-1597-2-7` in `initialize.F` (est. 5–8% gain on CBMC moves).
- **P8 Coefficient Caching**: Caching bend and torsion parameters in `vangle.F` (est. 1–2% gain).
- **P9 Ewald Hoisting**: Hoisting `twh_hinverse` out of reciprocal space loops (est. 10–18% on large Ewald systems).
- **P10 Branchless Minimum Image**: Vectorized `mimage.F` using `anint`.
- **P11 LTO / PGO**: Link-Time Optimization and Profile-Guided Optimization pipelines.
- **P13 Cell Lists**: Linked-cell list for large systems (>50,000 atoms) to achieve $O(N)$ pair scaling (projected **5–10×** speedup on large systems).

---

## Documentation & References

- **Project Specs**:
  - [`openspec/live_specs/optimization_plan.md`](file:///home/sharir/MCCSTowhee_opt/openspec/live_specs/optimization_plan.md) — Optimization phases 1–5 report and profile data.
  - [`openspec/live_specs/phase6_additional_improvements.md`](file:///home/sharir/MCCSTowhee_opt/openspec/live_specs/phase6_additional_improvements.md) — Phase 6 technical roadmap.
  - [`SUPERVISOR_REPORT.md`](file:///home/sharir/MCCSTowhee_opt/SUPERVISOR_REPORT.md) — Summary milestone report.
- **Upstream Documentation**:
  - Official Web Documentation: [towhee.sourceforge.net](https://towhee.sourceforge.net)
  - Input Manual: [`src/towhee-8.2.3/Manual/projects/towhee/input/towhee_input_v8_2_x.html`](file:///home/sharir/MCCSTowhee_opt/src/towhee-8.2.3/Manual/projects/towhee/input/towhee_input_v8_2_x.html)

---

## License

MCCCS Towhee is distributed under the GNU General Public License (GPL). See individual source files and `src/towhee-8.2.3/` for full licensing notices.
