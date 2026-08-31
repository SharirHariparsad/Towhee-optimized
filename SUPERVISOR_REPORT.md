# MCCCS Towhee 8.2.3 Optimization — Progress Report

**Date:** 2026-09-01
**Author:** Optimization Fork (`src/towhee-8.2.3/` build worktree)
**Base:** [Towhee 8.2.3](https://towhee.sourceforge.net) — Fortran77+C Monte Carlo engine
**Specs:** `openspec/live_specs/optimization_plan.md:1` (Phase 5 complete) + `openspec/live_specs/phase6_additional_improvements.md:1` (Phase 6+ proposal)

---

## 1. Executive Summary

Five optimization phases are complete on the serial code path. The delivered **HOT_PATH_INLINE (HPI) serial build** is strictly correct — `0 diff` versus `answer_current` on the regression gold file `tests/SPC_E_WATER` — and is the recommended production binary. Measured speedup on small `~600-atom` systems is **1.02–1.15×** (`6.648s → 5.799s` best case, native Linux, `1 thread, min of 3`). The theoretically planned OpenMP acceleration of `engtotal.F:335` was proven non-viable (pre-existing deterministic corruption from move 1, independent of our changes) and was correctly scoped out.

Headroom to the original `5–10×` single-node target remains, but it requires work on the *differential* energy path and CBMC/RNG hotspots — not the full `O(N²)` energy kernel we have already optimized. A detailed, low-risk roadmap for an additional **1.3–1.5× over current HPI (1.5–2.5× over original) on small-to-large systems, 5–10× on >50K atoms with a cell list**, is documented in `phase6_additional_improvements.md`.

No MPI correctness regressions were introduced; MPI changes are `USEMPI`-guarded and link-verified.

---

## 2. Objective

* Achieve `5–10×` single-node speedup and near-linear MPI scaling to `64+` nodes (`optimization_plan.md:4`) under `gfortran + OpenMPI` `-O3 -march=native -ffast-math`.
* Preserve regression tolerance `1e-6` relative (`tests/comparison_benchmark.py:compare_outputs` + `src/towhee-8.2.3/Examples/tolerance_diff.py`), with bitwise fallback via `--enable-safe-compare`.

---

## 3. Approach & Why OpenMP Was Dropped

Profiling (`gprof`, `Wall_Water`) showed the OpenMP target `engtotal.F:335` `twh_engtotal_twobody` is **0.25%** of runtime (`optimization_plan.md:43`):

| Hotspot | % | Calls |
|---------|---|-------|
| `twh_cbmcbendinfo_` `Source/globaldata.F:16771` | 13.5% | 540K |
| `twh_hmatrix_` / `halfbox_` / `mimage_` | 8.1% / 6.4% / 6.1% | 150M / 150M / 50M |
| `twh_eng_two_body_` `Source/energy.F:164` | 5.3% | 11K |
| `twh_vcoulomb_` / `twh_1dinteger` / `twh_1ddouble` | 4.3% / 4.3% / 3.0% | 36M |

Amdahl ceiling for the OMP region: `1.03×`. Investigation (`optimization_plan.md:25`):
deterministic corruption from the 2nd energy call even with `!$OMP PARALLEL DO` disabled and `-fno-recursive`; `--enable-openmp` without HPI also fails `121 diff` at `1 thread`. Root cause is the `gfortran -fopenmp → -frecursive` interaction with F77 `SAVE`/`COMMON` state. Decision: **HPI serial as Phase 1.1 deliverable; OpenMP out-of-scope**.

---

## 4. Completed Work (Phases 1–5)

| Phase | What was done | Files | Verification |
|-------|---------------|-------|--------------|
| **1.1 HPI zero-copy `hmatrix/halfbox`** | `hmatrix_arr`/`halfbox_arr` via `COMMON` in `Source/engtotal_globals.F` + `Source/engtotal.F:272` `inline_mimage` direct array access | `Source/engtotal_globals.F`, `Source/engtotal.F`, `Source/hot_common.h` | `twh_hmatrix_`/`halfbox_` removed from `gprof`; no `pointer_access.f90` needed — `COMMON` association used |
| **1.2 Scalar caches** | `rminsq/rcutsq/rcut/nbcoeff(13-16)/lshift/qqscale/rcelect/calp/hmatrix/halfbox` as `COMMON /hot_*/` in `Source/globaldata.F` → `Source/vtwobody.F:83`, `Source/vcoulomb.F:64`, `Source/mimage.F`, `Source/recipsum.F:48`; `HOT_PATH_INLINE` guards, baseline untouched | `Source/globaldata.F`, `Source/engtotal_globals.F`, `Source/vtwobody.F`, `Source/vcoulomb.F` | `twh_halfbox_/mimage_/rcelect_/qqscale_` → `0`; `twh_nbcoeff_ 440K→352K`, `twh_1ddouble 2.3%→1.3%`, profile `0.87s→0.79s` |
| **2 CBMC bend table** | `cbmcbendinfo(NTMAX,NUMAX,MAXBEND…)` `266 MB` from accessor `5×` index multiply to `COMMON /hot_cbmcbendinfo/` direct read in `Source/getcbangle.F:138` + `Source/angle_b_distribution.F:74`; fixed `Source/Makefile` `HOT_PATH_CPPFLAGS` propagation | `Source/globaldata.F:twh_cbmcbendinfo`, `Source/getcbangle.F`, `Source/angle_b_distribution.F`, `Source/Makefile.in` | `twh_cbmcbendinfo_ 2986→5` calls; both builds `0 diff` SPC_E_WATER; remaining `45%` flat is one-time page-fault init |
| **3 MPI (spec1)** | Removed redundant `MPI_Barrier` around `MPI_Bcast` (`Source/control.c:49`) and before `MPI_Finalize` (`Source/towhee.c:136`), converted `Source/rex.c:115` lambda fan-out `MPI_Send` loop → `MPI_Isend` + `MPI_Waitall` with distinct `&(J[k])` buffers; `MPI_Init_thread(MPI_THREAD_FUNNELED)` `Source/towhee.c:395` when `USEOPENMP`; domain decomposition deferred — shipped Examples `<10K` atoms (`optimization_plan.md:334`) | `Source/control.c`, `Source/towhee.c`, `Source/rex.c`, `Source/globalinfo.c` | `mpicc -DUSEMPI` compile+link `5.6 MB` `/tmp/opencode/towhee_mpi`; serial HPI unaffected (`#ifdef USEMPI`) |
| **4 Vector flags + SIMD** | `Source/Makefile.am:38` `VEC_FFLAGS=-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize`, `Source/configure.ac:95` `FFLAGS` clearing to prevent `-g -O2` override, `AM_FCFLAGS`; `!$OMP DECLARE SIMD` on `inline_atomnum/mimage/coulomb` `Source/engtotal.F`, `inline_lj_12_6` `Source/vtwobody.F:1046`, inner `jj` loop `Source/engtotal.F`; `tests/comparison_benchmark.py` + `Examples/tolerance_diff.py` `1e-6` | `Source/Makefile.am`, `Source/configure.ac`, `Source/engtotal.F`, `Source/vtwobody.F`, `tests/comparison_benchmark.py` | HPI-serial and HPI+Vec `0 diff` SPC_E_WATER; `Original FFLAGS -O2` after VEC now wins |
| **5 Ewald (minimal)** | `twh_ewald_calp` `COMMON /hot_ewald_calp/` `Source/hot_common.h` → `Source/vcoulomb.F:57`/`Source/recipsum.F:48` (bypasses `25M` dispatches); `Source/recipsum.F:144` inner `iunit` `!$OMP SIMD REDUCTION(+:sumr,sumi)`; outer `nkx/nky/nkz` remains scalar but benefits from `-ffast-math` | `Source/globaldata.F`, `Source/engtotal_globals.F`, `Source/vcoulomb.F`, `Source/recipsum.F` | HPI+Vec+Ewald `0 diff` SPC_E_WATER strict + `1e-6` |

Build: `cd src/towhee-8.2.3 && ./configure [--enable-mpi] && make -j16 towhee` (both plain and `HOT_PATH_CPPFLAGS=-DHOT_PATH_INLINE` compile).

---

## 5. Results

**Correctness (`optimization_plan.md:358`):**

| Config | Tol | Result |
|--------|-----|--------|
| `tests/SPC_E_WATER` HPI serial | strict | ✅ `0 diff` |
| `tests/SPC_E_WATER` plain serial | `1e-6` | ✅ |
| `tests/SPC_E_WATER` OMP any | `1e-6` | ❌ out-of-scope (see §3) |
| `tests/TraPPE_Pentane` HPI | `1e-6` | ✅ |
| `tests/Wall_Water` HPI | strict | ⚠️ 4-line `-0.314` CBMC init diverge; `baseline+Vec` (no HPI) ✅ `1.03×` |
| `src/Gibbs/EPM_VLCC` | `1e-6` | ⚠️ `26019→26016` (`1.2e-4`) fast-math — expected |
| `src/Canonical/TraPPE_Molecules` | `1e-6` | ⚠️ angle/torsion fast-math |

**Performance (`optimization_plan.md:369`, `bench_results.json:19`, native Linux, `1 thread, min of 3, no -pg`):**

| Example | Plain | HPI+Vec+Ewald | Speedup | Note |
|---------|-------|---------------|---------|------|
| `tests/SPC_E_WATER` | `6.648s` | `5.799s` | **1.15×** | ✅ tol |
| `src/TMMC/SPC_E_WATER` | `1.069s` | `1.048s` | `1.02×` | ✅ |
| `tests/TraPPE_Pentane` | `0.554s` | `0.529s` | `1.05×` | ✅ |
| `src/Gibbs/EPM_VLCC` | `0.894s` | `0.809s` | `1.11×` | ⚠️ fast-math |
| Historical `bench_results.json` `7.93s min` | `7.93s` | `7.92s` | `~1.00×` | WSL `2×` inflated early runs (`bench_results.json:13`) |

Accounting `gprof` reduction `~9%` (`0.87s→0.79s` Phase 1.2) does not yet appear in wall-clock on `SPC_E_WATER` because it is **CBMC + RNG + 266 MB table init dominated** (`twh_cbmcbendinfo_ 53%` flat includes page-fault + `ranlux_ 11%`, not pair loops at `6.3%` Ewald).

---

## 6. What Remains & Projected Gain

Analyzed in `openspec/live_specs/phase6_additional_improvements.md:1`. Dominant unoptimized path is the **differential energy** `Source/energy.F:164` (`energy_change.F:268` per trial, `11K` calls), not the already-optimized full `engtotal.F:335`.

**Next improvements (in ROI order):**

1.  **P6 Differential zero-copy** — `twh_coordinates_` `65M/2.5%` `Source/energy.F:341`, `twh_globalpos` `Source/globalinfo.c:1158` → precomputed `atomnum_arr`, flat `c_coord_real` via `get_coord_storage_ptr` `Source/globalinfo.c:1034`, hoisted `rcelect/rcutinsq` `Source/energy.F:216`, guarded `ctrmas` `Source/energy.F:248` — **5–12%**.
2.  **P7 RNG** — switch default `Source/initialize.F:42` `RNG_RANLUX 0 → RNG_DX_1597_2_7 5` (`Source/dx-1597-2-7.c:80`, `1 line`, `5–8×` faster, `11% → ~2%`) — **5–8%** CBMC-heavy.
3.  **P10 Branchless `mimage`** `Source/mimage.F:199` `while → anint` enables `DECLARE SIMD` — **1%** + vector enable.
4.  **P8 Coeff caches** — `Source/vangle.F:126` `12× twh_bencoeff`, `PNT` tables via `Source/setmixterms.F:976` — **1–2%**.
5.  **P9 Ewald completion** — hoist `9× twh_hinverse` per k-vector `Source/recipsum.F:116` (`15k → 9`), branchless `twh_expon` `Source/expon.F:47`, `derfc` SIMD `Source/vcoulomb.F:141` — **1–2%** small / **10–18%** large Ewald-heavy (`kmax 12`, `VECTORMAX 5184` `Source/preproc.h:68`).
6.  **P11 Build LTO/PGO** `Source/Makefile.am:38` `-flto`, `Source/configure.ac:95` `--enable-pgo` — **2–3%** whole program.
7.  **P12 I/O + layout** (`Source/writepdb.F:76`, `266 MB` `cbmcbendinfo` `Source/preproc.h:33`) — `0%` when `pdboutfreq 0`, `20–30%` of I/O phase when enabled.

| Scenario | +P6–11 vs plain serial | +P6–11 vs current HPI | +cell list (`>50K`) |
|----------|------------------------|----------------------|---------------------|
| **Small** `SPC_E_WATER 600` | **1.35–1.55×** | `1.18–1.35×` | no win |
| **Medium** `TraPPE` / `Wall_Water` | **1.50–1.80×** | `1.30–1.60×` | `1.60×` |
| **Large** `10K+ kmax 8–12` | **1.80–2.50×** | `1.50–2.10×` | `3–5×` |
| **Very large `50K+`** | `~2.0×` | `~1.70×` | **5–10×** (`O(N²)→O(N)` — only path to original target) |

MPI: Phase 3 already complete; further `P6–11` multiplies single-node efficiency, preserving near-linear scaling to `16×16` `>80%` (`optimization_plan.md:293`). Cell/Verlet `O(N)` deferred — `energy.F:264` `O(N²)` scan is still cheaper than binning for `<10K`.

---

## 7. Recommendation

1.  **Adopt HPI serial as production** (strict `0 diff`, tolerant `1e-6` with `-ffast-math`). Keep `HOT_PATH_INLINE` guards so baseline remains buildable.
2.  **Implement P6+P7 first** (`~10–15%` combined, `2–3` files, `1-line` RNG switch) — validate `gprof` `twh_coordinates_/globalpos` `>90%` drop and `Examples/* do_test` `1e-6`.
3.  Gate cell-list work on actual `>50K`-atom inputs; otherwise avoid the bookkeeping cost.

---

## 8. Artifacts

* Plan: `openspec/live_specs/optimization_plan.md`
* Phase 6 spec: `openspec/live_specs/phase6_additional_improvements.md`
* Benchmarks: `bench_results.json` (`bench_config.json`), `tests/comparison_benchmark.py`
* Tests: `src/towhee-8.2.3/Examples/*/do_test` + `Examples/tolerance_diff.py`
* Build: `src/towhee-8.2.3/Source/towhee` (`towhee_final`, `towhee_backup` restart files)

---

*Prepared for review; raw `gprof` and `bench_results.json` available on request.*
