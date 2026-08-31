# AGENTS.md — Towhee 8.2.3 Optimization

## Overview

[Towhee](https://towhee.sourceforge.net) is a Monte Carlo molecular simulation engine (Fortran 77 + C). This repo is an **optimization fork** of Towhee 8.2.3.

Active work is documented in `openspec/live_specs/`: `optimization_plan.md` (Phases 1–5 complete, HPI serial deliverable) and `phase6_additional_improvements.md` (Phases 6+ proposal). The only build worktree is `src/towhee-8.2.3/`; everything else (`.git`, `openspec/`, `tests/`, root reports) is meta.

**Current status (2026-09-01):** Phase 5 complete — Ewald `calp` COMMON zero-copy + `recipsum` SIMD, Phase 4 vector flags + SIMD + `1e-6` tolerance, Phase 3 MPI sync + `MPI_THREAD_FUNNELED`, Phases 1–2 accessor/CBMC zero-copy. **HPI serial (`HOT_PATH_INLINE`) is the production build** — `0 diff` vs `answer_current` on `tests/SPC_E_WATER`. OpenMP (`--enable-openmp`, `engtotal.F:335` `twh_engtotal_twobody`) is **out-of-scope**: pre-existing deterministic corruption from move 1 (see `optimization_plan.md:25`).

---

## Prerequisites

- `gfortran` (GCC 11+, tested 11–13)
- `GNU make`
- `autotools` (`autoconf`, `automake`) only for full rebuild from `Source/Makefile.am` / `Source/configure.ac`
- `OpenMPI` (`mpicc`/`mpifort`) for `--enable-mpi` builds
- `python3` for `tests/comparison_benchmark.py` + `src/towhee-8.2.3/Examples/tolerance_diff.py`

---

## Build

### Quick Start (Recommended — HPI Serial, Strictly Correct)

```bash
cd src/towhee-8.2.3
./configure
make -j16 -C Source towhee HOT_PATH_CPPFLAGS=-DHOT_PATH_INLINE
# or: make -j16 -C Source towhee   (plain serial, also correct)
```

Produces `Source/towhee` + `Source/forcefield`. Plain serial and HPI serial both pass `0 diff`; HPI adds `COMMON /hot_*/` zero-copy (`Source/engtotal_globals.F`, `Source/hot_common.h`, `Source/globaldata.F`).

### Reproducible Variants

| Variant | Configure | Flags | When to use |
|---------|-----------|-------|-------------|
| **Plain serial baseline** | `./configure FFLAGS="-g -O2 -ffixed-line-length-none"` | `make -C Source towhee` | Regression gold, `diff -w` strict |
| **HPI serial (production)** | `./configure` | `make -C Source towhee HOT_PATH_CPPFLAGS=-DHOT_PATH_INLINE` | `1.02–1.15×` on small systems (`6.648s→5.799s` `tests/SPC_E_WATER`, `optimization_plan.md:369`) |
| **HPI + Vector + Ewald** | `./configure` (auto `VEC_FFLAGS=-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize` `Source/Makefile.am:38` unless `--enable-safe-compare`) | `make -C Source towhee HOT_PATH_CPPFLAGS=-DHOT_PATH_INLINE` | Best single-node; needs `1e-6` tolerance for `Wall_Water`/`EPM_VLCC` fast-math diverge |
| **MPI** | `./configure --enable-mpi` | `make -C Source towhee` (add `HOT_PATH_CPPFLAGS` for HPI+MPI) | `mpirun -np N ./towhee`; `control.c`/`rex.c` `Isend` + `MPI_Init_thread(FUNNELED)` `Source/towhee.c:395` |
| **Safe compare** | `./configure --enable-safe-compare` | clears `VEC` to `-fbounds-check -DSAFE_COMPARE` `Source/Makefile.am:19` | Cross-platform bitwise, slower |

Notes:
- `engtotal_globals.F` + `hot_common.h` are already in `Source/Makefile.am`/`Makefile.in`; no manual `OBJEXT` edits needed (historical `BUILD_RECIPE.md:34` fix is now upstream).
- `VEC_FFLAGS` is appended *after* user `FFLAGS` in generated `Source/Makefile` so `-O3` wins over autoconf default `-g -O2` (`Source/configure.ac:95` clears defaults when unset).
- `--enable-openmp` is supported but **not recommended** — see OpenMP scope decision below.

### Source Layout (`src/towhee-8.2.3/`)

| Path | Purpose |
|------|---------|
| `Source/towhee.c` | C main entrypoint; argc/argv dispatch |
| `Source/mainloop.F` | Fortran main MC loop (`mainloop.F:384` move driver) |
| `Source/energy.F` | Differential two-body `twh_eng_two_body` (`energy.F:164`, dominant hot path) |
| `Source/energy_change.F` | Per-trial differential driver (`energy_change.F:268`) |
| `Source/engtotal.F` | Full `O(N²)` `twh_engtotal_twobody` (`engtotal.F:335`, `0.25%`, OMP target) |
| `Source/engtotal_globals.F` + `Source/hot_common.h` | HPI zero-copy `COMMON /hot_*/` (`hmatrix`, `nbcoeff`, `cbmcbendinfo` `266 MB`, `calp`) |
| `Source/vtwobody.F` | Pairwise LJ `twh_vtwobody` (`vtwobody.F:149`, `inline_lj_12_6`) |
| `Source/vcoulomb.F` / `Source/recipsum.F` | Coulomb + Ewald `twh_vcoulomb` / `twh_recipsum` (`recipsum.F:116` `hinverse` hoist target) |
| `Source/getcbangle.F` / `Source/angle_b_distribution.F` | CBMC angle generation (`getcbangle.F:138` `/hot_cbmcbendinfo/`, `13.5%`) |
| `Source/control.c` / `Source/rex.c` | MPI sync (`control.c:49` `Bcast` barrier removal) / replica exchange (`rex.c:115` `Isend`) |
| `Source/globaldata.F` | `506` `twh_*` accessors (`globaldata.F:16771`, `GLB_GET 2002` `preproc.h:360`) |
| `Source/globalinfo.c` | `52` C wrappers (`globalinfo.c:1034` `twh_coordstorage_` `PNT_COORDSTORAGE`) |
| `Source/towhee` | Compiled executable |
| `Source/*.F` / `Source/*.c` | `~150` Fortran / `7` C sources |
| `Source/forcefield` | Standalone force field writer |
| `ForceFields/` | `87+` parameter files |
| `Examples/` | `~50` self-contained simulations by ensemble (`towhee_input`, `header`, `answer_current`) |
| `Utils/` / `Manual/` | Utilities (`pdb2towhee`) / HTML docs (`projects/towhee/input/towhee_input_v8_2_x.html`) |

---

## Tests

Each `Examples/*/` is a self-contained regression test:

```bash
cd Examples/Some_Example     # contains towhee_input, header, answer_current
../../do_test                # runs ../../Source/towhee, diff -w strict + tolerance_diff.py 1e-6
cat results
```

Other helpers:

```bash
./reset_files                # rebuild gold files after verifying output
./run_serial_test            # all serial tests from Examples/, aggregates test_report
python3 ../../tests/comparison_benchmark.py --original /path/plain/towhee --optimized ./Source/towhee --run-dir ../../tests/SPC_E_WATER --threads 1 --repeats 10 --check-tolerance 1e-6
python3 Examples/tolerance_diff.py answer_current new_output 1e-6  # standalone 1e-6 relative
rm -f gmon.out && ./Source/towhee && gprof ./Source/towhee gmon.out | head -n 30  # profile
```

Regression matrix (`optimization_plan.md:358`): HPI serial `SPC_E_WATER` strict `0 diff`; `TraPPE_Pentane` `1e-6` pass; `Wall_Water` HPI strict `4-line -0.314` CBMC diverge (baseline+Vec `1.03×` passes); `EPM_VLCC`/`TraPPE_Molecules` fast-math `1.2e-4` needs `1e-6` relaxed.

---

## Optimization Project

### Completed (Phases 1–5) — See `openspec/live_specs/optimization_plan.md`

* **Phase 1** Accessor zero-copy (`hmatrix`/`halfbox` `150M` calls `8.1%/6.4%`, `mimage` `6.1%`, `rcut`/`qqscale`/`calp` `COMMON /hot_*/`, `twh_nbcoeff_` `440K→352K`) — `~9%` CPU accounting `0.87s→0.79s`.
* **Phase 2** CBMC `twh_cbmcbendinfo_ 13.5% 540K` → `COMMON /hot_cbmcbendinfo/` `getcbangle.F:138` (`2986→5` calls, `45%` flat is one-time `266 MB` page-fault).
* **Phase 3** MPI — remove redundant `MPI_Barrier` (`control.c:49` `towhee.c:136`), `rex.c:115` `MPI_Send → Isend+Waitall` (distinct `&(J[k])` buffers), `MPI_Init_thread(FUNNELED)` (`towhee.c:395`); domain decomposition deferred (`<10K` atoms replicate, spatial only `>50K`).
* **Phase 4** Vector — `VEC_FFLAGS` `Source/Makefile.am:38` + `configure.ac:95` fix, `!$OMP DECLARE SIMD`/`SIMD` (`engtotal.F:272`, `vtwobody.F:1046`, `recipsum.F:144`), tolerance harness.
* **Phase 5** Ewald — `twh_ewald_calp` `COMMON /hot_ewald_calp/` (`25M` calls bypassed), `recipsum.F:144` `SIMD REDUCTION`.

### Proposed (Phase 6+) — See `openspec/live_specs/phase6_additional_improvements.md`

P6 differential zero-copy (`energy.F:164` `twh_coordinates_` `65M/2.5%`, `twh_globalpos` `engtotal.F:256` → `atomnum_arr`), P7 RNG `initialize.F:42` `RANLUX 0 → DX_1597_2_7 5` (`11% → ~2%`), P10 branchless `mimage.F:199`, P8 coeff caches (`vangle.F:126`), P9 Ewald hoist (`recipsum.F:116` `15k→9`), P11 LTO/PGO, P12 I/O + `cbmcbendinfo` layout, P13 cell list (`>50K` only). Projected `1.35–1.55×` small / `1.80–2.50×` large `/ 5–10× (>50K with cell list)` vs plain.

### OpenMP Scope Decision

**Dropped.** OMP on `engtotal.F:335` ceiling `1.03×` (`0.25%`). Build corrupts from move 1 even with `PARALLEL DO` disabled / `-fno-recursive`; plain `--enable-openmp` without HPI fails `121 diff` at `1T`. Suspected `gfortran -fopenmp → -frecursive` vs F77 `SAVE`/`COMMON`. HPI serial is the deliverable.

### Historical Fixes Consolidated (from `BUILD_RECIPE.md`)

* `engtotal.F` OMP region (when enabled): `lfailure` `PRIVATE → FIRSTPRIVATE` (dummy arg uninitialized), `lovrlap` explicit `.false.` before `twh_vtwobody`/`twh_vcoulomb` (beyond-cutoff leaves unset), `lonefour_global` `PRIVATE → FIRSTPRIVATE`, debug `lomp_failure` write removed.
* `engtotal_globals.F` `engtotal_globals_init()` sequencing — now correctly after `setmixterms`/`setclassical` in `readtowhee.F`; `nbcoeff` `COMMON` copies verified stable.
* Benchmark note: earlier `BUILD_RECIPE.md:54` `8.5s` baseline / `~17%` slower OMP was WSL `2×` inflated (`bench_results.json:13`); native Linux is `6.648s → 5.799s` `1.15×`.

---

## Input/Output Convention

- Input: `towhee_input` (free-format keyword file)
- Restart in: `towhee_initial` (`linit .false.`)
- Restart out: `towhee_backup`, `towhee_final`
- Trajectory: `box_*_step_*.pdb` (`writepdb.F:22`, `pdboutfreq`), `towhee_movie` (`writemovie.F:22`, `moviefreq`)
- `towhee_parallel` (`mainloop.F:346`) for MPI
- Standard output contains averages, block averages, move statistics

---

## Key Gotchas

* **Fortran/C mixed code**: `.F` = Fortran 77 + CPP. Headers `preproc.h`, `functions.h`, `globalc.h`; cross-language `wrapperF.h`/`wrapperC.h`. Generated `Makefile.in` tracks `Makefile.am`; `hot_common.h` is now part of build.
* **RNG**: `preproc.h:669` `RNG_RANLUX 0`, `RNG_DX_1597_2_7 5`, `RNG_KISS99 7`, `RNG_MRG32K3A 8` — runtime `initialize.F:42` (`RNG_DX_1597` proposed default for speed `phase6: P7`).
* **Regression tests** are `diff -w` + `1e-6` relative `tolerance_diff.py` — `Wall_Water`/`EPM` fast-math diverges are known and tolerated (`optimization_plan.md:369`).
* **Force field names** must match `ForceFields/towhee_ff_*` exactly.
* **No build artifacts** tracked in git (only `Makefile.in`); `Source/*.o`, `Source/towhee`, `autom4te.cache/` are ignored.
* **Array sizes**: `preproc.h:27` `MAXBOX 3`, `NUMAX 2202`, `NTMAX 7`, `MAXKMAX 12` → `VECTORMAX 5184`, `CROSSTYPEMAX 70500`; `cbmcbendinfo` `266 MB` — do not bump without dynamic allocation.
* **Box-global scalars**: `rcelect`/`rcut`/`box_volume` change on volume moves (`mainloop.F:713` `updatevolmaxdisp`) — cache invalidation required.

---

## Commands Reference

```bash
# Benchmark (uses bench_config.json schema, results → bench_results.json)
python3 tests/comparison_benchmark.py --original /home/sharir/towhee-8.2.3/Source/towhee --optimized src/towhee-8.2.3/Source/towhee --run-dir tests/SPC_E_WATER --threads 1 --repeats 10 --warmup 2

# MPI
mpirun -np 4 --bind-to core --map-by socket ./Source/towhee
export OMPI_MCA_btl=^openib  # if no IB

# Supervisor report
cat SUPERVISOR_REPORT.md
```

## References

* `openspec/live_specs/optimization_plan.md` — full Phase 1–5 plan, profiling, benchmarks, file map
* `openspec/live_specs/phase6_additional_improvements.md` — detailed Phase 6+ spec with `file:line` pins and speedup model
* `SUPERVISOR_REPORT.md` — 2-page summary for review
* Towhee manual `src/towhee-8.2.3/Manual/projects/towhee/input/towhee_input_v8_2_x.html`
