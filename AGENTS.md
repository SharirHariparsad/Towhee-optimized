# AGENTS.md — Towhee 8.2.3 Optimization

## Overview

[Towhee](https://towhee.sourceforge.net) is a Monte Carlo molecular simulation engine (Fortran 77 + C). This repo is an **optimization fork** of Towhee 8.2.3 with active MPI sync reduction and OpenMP parallelization work under `openspec/live_specs/`.

## Build

Only the `src/towhee-8.2.3/` subtree is the build worktree. Everything else (`.git`, `openspec/`, `tests/`, root guides) is meta.

```bash
cd src/towhee-8.2.3
./configure                # or ./configure --enable-mpi
make                       # produces Source/towhee + Source/forcefield
```

Use `./configure --enable-mpi` for MPI support, `--enable-safe-compare` for cross-platform consistency.

## Tests

Each example in `Examples/*/` is a self-contained regression test:

```bash
cd Examples/Some_Example     # contains towhee_input, header, answer_current
../../do_test                # runs ../../Source/towhee, diffs output vs answer_current
```

To rebuild gold files: `./reset_files` after verifying new output is correct.

To run all serial tests: `./run_serial_test` from `Examples/`. Results aggregated into `test_report`.

## Source Layout (`src/towhee-8.2.3/`)

| Path | Purpose |
|------|---------|
| `Source/towhee.c` | C main entrypoint; argc/argv dispatch |
| `Source/mainloop.F` | Fortran main MC loop |
| `Source/control.c` | MPI communication (focus of MPI sync work) |
| `Source/rex.c` | Replica exchange (MPI async target) |
| `Source/energy.F` | Energy calculation (OpenMP target) |
| `Source/vtwobody.F` | Pairwise nonbonded (OpenMP target) |
| `Source/towhee` | Compiled executable |
| `Source/*.F` | Fortran source (~150 files) |
| `Source/*.c` | C source (7 files: control, RNG, jobfarm, rex, etc.) |
| `Source/forcefield` | Standalone force field writer binary |
| `ForceFields/` | 87+ force field parameter files |
| `Utils/` | Utility programs (pdb2towhee, analyse_movie, fitcoex, etc.) |
| `Examples/` | ~50 self-contained simulation examples by ensemble |
| `Manual/` | HTML docs to build; key file: `projects/towhee/input/towhee_input_v8_2_x.html` |

## Optimization Project

Active specs in `openspec/live_specs/`:

- **spec1_MPI_sync.md** — Remove redundant MPI_Barrier calls in `control.c`/`towhee.c`; convert serial MPI_Send to non-blocking MPI_Isend in `rex.c`
- **spec3_openmp_parallelization_v2.md** — Parallelize the `ichain`/`jchain` pair loop in `energy.F`/`vtwobody.F` with carefully scoped PRIVATE/REDUCTION/SHARED clauses and SCHEDULE(dynamic) for load balance

## Input/Output Convention

- Input: `towhee_input` (free-format keyword file)
- Restart in: `towhee_initial` (set `linit .false.`)
- Restart out: `towhee_backup`, `towhee_final`
- Trajectory: `box_*_step_*.pdb`
- Standard output contains all averages, block averages, move statistics

## Key Gotchas

- **Fortran/C mixed code**: `.F` files are Fortran 77 with C preprocessor. Headers: `preproc.h`, `functions.h`, `globalc.h`. `wrapperF.h`/`wrapperC.h` handle cross-language linkage.
- **RNG**: Multiple generators (dx-1597-2-7, RANLUX, mrg32k3a, kiss) — selected at compile time via `preproc.h` defines.
- **Regression tests** are simple output diffs — small numerics differences (platform/compiler) can cause failures.
- **Force field atom type names** must match exactly what appears in `ForceFields/towhee_ff_*` files.
- **No build artifacts** are tracked in git (only generated `Makefile.in` from autotools).
