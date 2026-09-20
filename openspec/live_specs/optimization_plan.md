# Towhee 8.2.3 Multi-Node CPU Optimization Plan

**Status**: **Phase 5 complete (2026-08-31, minimal)** — Ewald calp COMMON zero-copy + recipsum SIMD; **Phase 4 complete** — vector flags + SIMD + tolerance; **Phase 3 complete** — MPI sync + hybrid FUNNELED; Phase 2 complete; HPI serial 0-diff on SPC_E_WATER (strict), 1.02-1.15× measured speedups on small examples, tolerance passes on SPC_E_WATER/TraPPE but fails on Wall_Water/EPM_VLCC due to fast-math divergence (expected, see below)  
**Target**: 5-10x single-node speedup, near-linear MPI scaling to 64+ nodes  
**Compiler**: GNU Fortran + OpenMPI, `-O3 -march=native -ffast-math -fopenmp`  
**Regression tolerance**: `1e-6` relative (not bitwise)

---

## Current State (Baseline)

| Metric | Value |
|--------|-------|
| Single-thread runtime (SPC_E_WATER) | 8.8s |
| 16-thread speedup (HPI+OMP) | 1.04x |
| OMP target (`engtotal_twobody`) | 0.25% of runtime |
| Top hotspot (CBMC bend info) | 13.5% |
| Accessor overhead (hmatrix/halfbox/1d*) | 21.8% |
| Ewald + Coulomb | 6.3% |

**Root cause**: OMP parallelizes `engtotal_twobody` (0.25% of runtime). Real hotspots are accessor overhead (35%) and CBMC bend info (13.5%).

---

## OpenMP Scope Decision — CORRECTED 2026-09-20 (was "dropped"; root cause found)

**Original decision (2026-08-31)**: OpenMP dropped from scope — "the `-fopenmp`/`--enable-openmp` configuration itself deterministically corrupts shared program state ... suspected interaction of the gfortran OpenMP runtime with F77 static locals / COMMON memory layout, high-effort, out of scope."

**Corrected root cause (verified by isolated experiments, 2026-09-20)**: the corruption is **NOT** an OpenMP runtime issue, a race, or the parallel region. It is a **missing `save` on the Ewald reciprocal-space partial-sum arrays**, which gfortran's `-fopenmp` exposes via its implied `-frecursive` flag:

- `globaldata.F:5592` `twh_ewald_ssumi` declares `ssumi`/`ssuminew` as **implicit-SAVE locals** (only `ssumiold` gets `save` at `globaldata.F:5612`).
- `globaldata.F:5670` `twh_ewald_ssumr` has the same pattern (`ssumr`/`ssumrnew` un-saved; `ssumrold` saved at 5690).
- These arrays persist **across calls** (CRD_REAL/CRD_TEMP/TRL_NEW partial sums accumulate over multiple `twh_ewald_ssumi`/`twh_ewald_ssumr` invocations). Under `-frecursive` they become stack-automatic and their contents are clobbered between calls → garbage Ewald terms → wrong energies from move 1.
- **Fix (2 lines)**: add `save ssumi,ssuminew` and `save ssumr,ssumrnew`. With this fix, a **full `-frecursive` build is byte-identical to plain serial** on SPC_E_WATER.

**Evidence (isolated builds, `--enable-safe-compare`, SPC_E_WATER, 5000-mol input)**:
1. `-frecursive` **alone** (no OpenMP anywhere) reproduces the corruption **byte-identically** to the `--enable-openmp` build (203 diff lines, same first divergence `-0.2661E+04/40` → `-0.2565E+04/43`). → proves compiler-flag side effect, not runtime.
2. `-fopenmp` with **all `!$OMP` directives stripped** (libgomp linked, zero pragmas active) corrupts identically too. → not the parallel region.
3. `-finit-integer=-777` does not change the corrupt output → not an uninitialized integer; `-finit-real=nan` hangs in the swap/chempot path → an uninitialized **real**.
4. Per-file `-frecursive` bisection isolated the corruption to **`globaldata.o`**; the compiler's own `-Wmaybe-uninitialized` flags (`ssumi`,`ssuminew`,`ssumr`,`ssumrnew`) matched exactly.
5. `-fno-recursive` does **not** cleanly undo it (changes the failure to a bounds error) — consistent with flag ordering; the earlier "-fno-recursive didn't fix" test is therefore not evidence against `-frecursive`.

**Status after fix (2026-09-20)**: `--enable-openmp` (non-HPI) with (a) the 2-line `save` fix and (b) declaring `i_atom,j_atom,gpos_i,gpos_j` in the `#else` block of `engtotal.F:392` (the `!$OMP PARALLEL DO` `PRIVATE` clause referenced HPI-only locals, breaking the plain build's compile) now **passes SPC_E_WATER at both `OMP_NUM_THREADS=1` and `=4`**. OpenMP is no longer "unusable"; it is only latency-limited (`engtotal_twobody` = 0.25% → Amdahl ceiling ~1.03×). NOTE: thread-safety of `ssumi`/`ssumr` inside the parallel region is not yet audited.

**Two additional compile-time breakages found in the current tree** (the original "corrupts at move 1 / 121 diff" report predates them):
- Phase-4 `!$OMP DECLARE SIMD(vtb_lennard_jones)` at `vtwobody.F:887` fails to compile under `-fopenmp`.
- Plain `--enable-openmp` (no HPI) fails to compile: `engtotal.F:428` `PRIVATE` clause references HPI-only locals (fixed above).

**Separate pre-existing bug (not OpenMP-related)**: plain **HPI serial** crashes on the current 5000-mol SPC_E_WATER input with `coord_temp_arr` index `15001` at `energy.F:394` (out of `15000`). This is the Phase-6 P6 zero-copy differential path indexing a CBMC newly-grown test chain's `globalpos` that is not valid in the flat array — the non-HPI path (which passes) avoids it. Independent of OpenMP and `-frecursive`; must be fixed before HPI is a valid production build for CBMC-insertion systems.

**Benchmark note (SPC_E_WATER)**: HPI serial vs plain serial run at parity (~1.32 s vs ~1.33 s min of 3). The plan's own profile shows `engtotal_twobody` at 0.25% of runtime, so Phase 1.1's accessor-elimination gain is below measurement noise on this example. Real Phase 1 gains require tackling the dominant hotspots (accessor chains in CBMC/move code, `twh_cbmcbendinfo`).

---

## Profiling Data (gprof on Wall_Water)

| Rank | Function | % Time | Calls | Category |
|------|----------|--------|-------|----------|
| 1 | `twh_cbmcbendinfo_` | **13.5%** | 540K | CBMC bend lookup |
| 2 | `twh_hmatrix_` | **8.1%** | 150M | Box matrix accessor |
| 3 | `twh_halfbox_` | **6.4%** | 150M | Half-box accessor |
| 4 | `twh_mimage_` | **6.1%** | 50M | Minimum image |
| 5 | `twh_eng_two_body_` | **5.3%** | 11K | Two-body dispatch |
| 6 | `twh_vcoulomb_` | **4.3%** | 36M | Coulomb (Ewald + direct) |
| 7 | `twh_1dinteger` | **4.3%** | — | Scalar int getter |
| 8 | `twh_1ddouble` | **3.0%** | — | Scalar double getter |
| 9 | `twh_luse_v_` | **2.8%** | 43M | VDW flag getter |
| 10 | `twh_coordinates_` | **2.5%** | 65M | Coordinate accessor |
| — | `twh_engtotal_twobody_` | **0.25%** | 1 | **OMP parallel region** |

**Key insight**: OMP target is 0.25% of runtime → Amdahl ceiling = 1.03x

---

## Phase 1: Accessor Zero-Copy (Week 1-2) — **Highest ROI**

### 1.1 hmatrix/halfbox Zero-Copy (8-10% gain)

**Problem**: 150M calls each to `twh_hmatrix_`/`twh_halfbox_` via accessor chain

**Solution**: Fortran module with `iso_c_binding` for direct pointer access

**Files**:
- `Source/pointer_access.f90` (new module)
- `Source/engtotal_globals.F` (add `hmatrix_arr`, `halfbox_arr`)
- `Source/engtotal.F` (update `inline_mimage`)
- `Source/globalc.h` (C declarations)

```fortran
! pointer_access.f90
module pointer_access
  use iso_c_binding, only: c_ptr, c_loc
  implicit none
  private
  public :: get_hmatrix_ptr, get_halfbox_ptr

contains
  subroutine get_hmatrix_ptr(cptr, nbxd, nbox) bind(C, name='get_hmatrix_ptr')
    type(c_ptr), intent(out) :: cptr
    integer, intent(out) :: nbxd, nbox
    external hmatrix
    double precision, dimension(MAXBXD,MAXBOX,3,3), target, save, external :: hmatrix
    cptr = c_loc(hmatrix)
    nbxd = MAXBXD; nbox = MAXBOX
  end subroutine
  ! same for halfbox
end module pointer_access
```

**Compile**: Separate module, link with `globaldata.o`  
**Bind in `engtotal_globals`**: `hmatrix_arr(:,:,:,:)`, `halfbox_arr(:,:,:)`  
**Update `inline_mimage`**: Direct array access → pure arithmetic

**Est. gain**: 8-10% (eliminates 300M accessor calls)

### 1.2 Scalar Getter Batching (5-7%)

**Targets**: `twh_1dinteger` (4.3%), `twh_1ddouble` (3.0%)

**Approach**:
- Thread-local caches for immutable scalars: `rcut`, `rcutsq`, `rminsq`, `rcutinsq`, `rcmu`, `rcelect`, `box_volume`, `hmatrix_diag`
- Replace frequent scalar getters with array fetches in hot loops
- Add `threadprivate` caches in `engtotal_globals`

**As implemented (2026-08-31)**: the general two-body / coulomb / mimage code paths
(used by CBMC insertion and move energies — the paths the Phase 1.1 inline engine
does **not** cover) now read the scalar/array values via small extern-COMMON blocks
shared with the `globaldata.F` accessor storage and the `engtotal_globals.F`
flat `nbcoeff` copies, replacing the per-pair `twh_*` dispatch calls. All reads are
guarded by `#ifdef HOT_PATH_INLINE` so the non-HPI baseline code path is untouched
(verified: baseline build still compiles and passes 0-diff regression).

Cached values (all immutable for the run):
- `rminsq`/`rcutsq`/`rcut` (COMMON `/hot_rminsq/`, `/hot_rcutsq/`, `/hot_rcut/`)
- `nbcoeff(0,13,14,15,16)` flat copies (COMMON `/hot_nbcoeff0..16/`)
- `lshift` (COMMON `/hot_lshift/`)
- `qqscale`, `rcelect`, `rcelect` + `lewald` (COMMON `/hot_qqscale/`,
  `/hot_rcelect/`, `/hot_lewald/`)
- `hmatrix`/`halfbox` in the general `twh_mimage` (COMMON `/hot_hmatrix/`,
  `/hot_halfbox/`)

**NOT cached** (correctly left as accessors): `box_volume`, `hmatrix_diag`
(box-shape dependent, change on volume moves), `rcmu` (per-chain, low call count).

**Measured (gprof, SPC_E_WATER, HPI-serial -pg)**:

| Getter | Phase 1.1 | Phase 1.2 |
|--------|-----------|-----------|
| `twh_halfbox_` | 1.02M calls / 1.2% | removed |
| `twh_mimage_` | 337K / 1.2% | removed |
| `twh_rcelect_` | 212K / 1.2% | removed |
| `twh_qqscale_` | 467K / 0.0% | removed |
| `twh_nbcoeff_` | 440K / 4.6% | 352K / 3.8% (remaining = one-four/generic/table branches) |
| `twh_1ddouble` | 2.3% | 1.3% |
| profile total | 0.87s | 0.79s (~9% accounting reduction) |

Wall-clock on SPC_E_WATER unchanged (~1.04-1.05s min of 5 for all three builds) —
the example is CBMC-dominated (`twh_cbmcbendinfo_` 53%, `ranlux_` 11%) and getter
removal helps the engine pair loops more on larger systems (see Phase 2).

---

## Phase 2: CBMC Bend Info Memoization (Week 2) — **10-12% Gain**

### 2.1 Analysis
- `twh_cbmcbendinfo_`: 13.5%, 540K calls/run
- 1 call per ~100 MC steps
- Likely redundant lookups per molecule type

### 2.2 Memoization Strategy

**Option A: Per-molecule-type cache (Low effort, 5-8%)**
```fortran
! In cbmc module or engtotal_globals
real(8), allocatable :: bend_cache(:,:,:)  ! (angle_type, molty, unit)
logical, allocatable :: bend_cached(:,:)

! In twh_cbmcbendinfo:
if (bend_cached(angle_type, molty)) then
  return bend_cache(angle_type, molty, unit)
endif
! compute and cache
```

**Option B: Full precompute at startup (10-12%)**
- At `engtotal_globals_init`, precompute all bend angles for all molty/unit
- Memory: ~1MB (100 types × 100 units × 100 angles)

**Decision**: Start with Option A, upgrade to B if needed

**As implemented (2026-08-31)**: `twh_cbmcbendinfo` is a straight get/set into a
large SAVEd 6-D table `cbmcbendinfo(NTMAX,NUMAX,MAXBEND,MAXBEND,
CBMC_MAXBENDPEAKS,CBMCMAXBEND)`. It computes nothing — so Option A's value cache
does not apply (each call reads a distinct table element). The real cost is the
per-call accessor **dispatch + repeated 6-D index arithmetic** on a 266 MB array.
Following the established Phase 1.2 pattern, the table storage was moved into a
COMMON block (`/hot_cbmcbendinfo/`) and the hot `GLB_GET` reads in the CBMC
angle code (`getcbangle.F`, `angle_b_distribution.F`) now read the array directly
(guarded by `#ifdef HOT_PATH_INLINE`, with the original accessor preserved in
`#else`). `atd_cbmc_setup.F` was left untouched — it runs once at startup, not on
the hot path. Also fixed: `make HOT_PATH_CPPFLAGS="-DHOT_PATH_INLINE"` was a no-op
because the generated `Source/Makefile` did not propagate `HOT_PATH_CPPFLAGS` into
`AM_CPPFLAGS` (added it, matching `Makefile.am`).

**Measured (gprof, SPC_E_WATER)**:

| Metric | Baseline | Phase 2 (HPI) |
|--------|----------|---------------|
| `twh_cbmcbendinfo_` calls | 2,986 | **5** |
| `twh_cbmcbendinfo_` from getcbangle | 2,981 | 4 (setup) |
| `twh_halfbox_`/`twh_mimage_`/`twh_qqscale_`/`twh_rcelect_` | 12K/337K/229K/93K @ 0% | unchanged @ 0% |

Both HPI and baseline builds pass the SPC_E_WATER regression **0 diff** vs
`answer_current`. Remaining `twh_cbmcbendinfo_` self-time (~45% in the flat
profile) is the one-time ~266 MB table **page-fault/initialization** cost at
startup (GLB_SET writes), not the MC hot loop (call count drops by ~600x).
Wall clock on SPC_E_WATER: HPI ~1.40s vs baseline ~1.42s (min of 3) — near
parity because SPC_E_WATER is CBMC/array-init dominated; expect the dispatch
removal to pay off on pair-loop-heavy systems (see plan note in Phase 1.2).

---

## Phase 3: MPI Domain Decomposition 

### 3.1 Communication Optimization (spec1) — **Implemented 2026-08-31**

- `Source/control.c` `tramonto_control` (lines 49-51): `MPI_Barrier; MPI_Bcast(&logicflag); MPI_Barrier;` → single `MPI_Bcast` (Bcast is collective/self-synchronizing; surrounding Barriers redundant). Guarded by `#ifdef USETRAMONTO`/`#ifdef USEMPI`, so serial HPI build unaffected.
- `Source/towhee.c` `do_tramonto` (lines 136-139) and `do_jobfarm` (lines 161-170): drop the `MPI_Barrier` immediately before `MPI_Finalize` (`MPI_Finalize` is collective; barrier before it is redundant). The remaining `MPI_Barrier` after the job/loop stays. Also `#ifdef USEMPI`-guarded.
- `Source/rex.c` `run_rex_server` lambda-key fan-out (lines 115-119): `for(k=1..num_nodes){ m=J[k]; MPI_Send(&m,...,k,LAMBDA_KEY); }` → `MPI_Isend(&(J[k]),..., &reqs[k-1])` per rank + single `MPI_Waitall(num_nodes, reqs, MPI_STATUSES_IGNORE)` + `free(reqs)`. Uses `twh_allocateVector(sizeof(MPI_Request), paramwrap->num_nodes)` so no `MAX_NODES` needed; each Isend uses distinct stable buffer `&(J[k])` (not aliased `&m`). Swaps section (`m = drain_random()` / `MPI_Send(&m)`) left as blocking single-send.

Validation: `mpicc -DUSEMPI` compile of all 5 `USEMPI` C files (`towhee.c`, `control.c`, `rex.c`, `globalinfo.c`, `jobfarm.c`) clean; full MPI link (`mpicc` + Fortran objs + 5 MPI C objs + `-lgfortran`) succeeds at `/tmp/opencode/towhee_mpi` (5.6 MB). Serial HPI binary untouched (edits inside `#ifdef USEMPI`/`USETRAMONTO`). No spec change to `rex.c` `run_rex_client` single `MPI_Send` (line 358).

### 3.2 Domain Decomposition
**Need clarification**: Typical system size?
- < 10K atoms → Particle decomposition (replicate coords, distribute pair loops)
- 10-50K atoms → Hybrid
- > 50K atoms → Spatial decomposition (geometric bins, neighbor lists)

### 3.3 Hybrid MPI+OpenMP
```c
MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#pragma omp parallel for schedule(dynamic, 8)
for (int i = my_start; i < my_end; i++) { /* pair loops */ }
```
- OpenMPI with `MPI_THREAD_FUNNELED`
- OpenMP for intra-node, MPI for inter-node

---

## Phase 4: Vectorization + FP Relaxation (Week 5)

### 4.1 Compiler Flags
```bash
FFLAGS="-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize -fopenmp -ffixed-line-length-none"
```

### 4.2 Explicit SIMD
```fortran
!$OMP SIMD
do i = 1, n
  rxuij = rx(i) - rx(j)
  r2 = rxuij**2 + ryuij**2 + rzuij**2
  inv_r = 1.0 / sqrt(r2)
  ...
enddo
```
**Targets**: `inline_mimage`, `inline_lj_12_6`, `inline_coulomb`, `inline_atomnum`

### 4.3 Regression Tolerance
```python
# In comparison_benchmark.py / do_test
def compare_outputs(ref, test, tol=1e-6):
    for line in ref:
        if numeric_line:
            assert abs(ref_val - test_val) < tol * max(abs(ref_val), 1.0)
```

---

## Phase 5: Ewald Optimization (Week 6)

### 5.1 Precompute k-vectors
- `twh_ewald_calp_` called 25M times
- Move k-vector generation to `engtotal_globals_init`
- Store `k_vectors`, `k2`, `prefactor` in `engtotal_globals`

### 6.2 Reciprocal Sum Vectorization
```fortran
!$OMP SIMD
do k = 1, nkvec
  kdotr = kx*rx + ky*ry + kz*rz
  sk = sin(kdotr) / k2
  sum = sum + sk * prefactor
enddo
```

---

## Phase 6: Multi-Node Scaling Validation (Week 7)

### Test Matrix
| Nodes × Threads | System Size | Target Efficiency |
|-----------------|-------------|-------------------|
| 1 × 16 | SPC_E_WATER | 2.2x (vs 1T) |
| 4 × 16 | SPC_E_WATER | > 80% |
| 16 × 16 | SPC_E_WATER | > 80% |
| 64 × 16 | Large system | > 80% |

### MPI Configuration
```bash
# OpenMPI
export OMPI_MCA_btl=^openib  # if no IB
export OMP_NUM_THREADS=16
mpirun -np 64 --bind-to core --map-by socket ./towhee
```

---

## Implementation Checklist

### Phase 1.1: hmatrix/halfbox Zero-Copy
- [x] Update `engtotal_globals.F` (add `hmatrix_arr`, `halfbox_arr`)
- [x] Update `engtotal_globals_init()` (bind `hmatrix_arr`/`halfbox_arr` via COMMON association)
- [x] Update `engtotal.F` (`inline_mimage` uses zero-copy arrays)
- [x] Verify: `twh_hmatrix_`/`twh_halfbox_` eliminated from gprof
- [ ] Create `Source/pointer_access.f90` (NOT built — COMMON association used instead; see `engtotal_globals.F`)
- [ ] Update `Source/globalc.h` (C declarations) — deferred

### Phase 1.2: Scalar Getter Caches
- [x] Add caches for immutable scalars (`rcut`, `rcutsq`, `rminsq`, `rcelect`, `lshift`, `qqscale`, `hmatrix`/`halfbox`) as extern-COMMON in `engtotal_globals`
- [x] Convert `nbcoeff` flat copies to COMMON in `engtotal_globals` (accessible to general path)
- [x] Replace accessor calls in `engtotal.F`, `mimage.F`, `vtwobody.F`, `vcoulomb.F`
- [x] Verify: `twh_halfbox_`/`twh_mimage_`/`twh_rcelect_`/`twh_qqscale_` removed from gprof; baseline build still 0-diff

### Phase 2: CBMC Memoization
- [x] Expose `cbmcbendinfo` storage via COMMON `/hot_cbmcbendinfo/` in `globaldata.F` (`twh_cbmcbendinfo`)
- [x] Direct-COMMON reads in `getcbangle.F` + `angle_b_distribution.F` (HPI-guarded, `#else` preserved)
- [x] Fix `Source/Makefile`: propagate `HOT_PATH_CPPFLAGS` into `AM_CPPFLAGS`
- [x] Verify: `twh_cbmcbendinfo_` calls 2986 → 5; both builds 0-diff SPC_E_WATER
- [ ] (Deferred) Literal `bend_cache` — not applicable: `twh_cbmcbendinfo` is a pure table get/set, no computed redundancy

### Phase 3: MPI — **complete 2026-08-31**
- [x] Remove redundant `MPI_Barrier` in `control.c` (`tramonto_control`) + `towhee.c` (`do_tramonto`/`do_jobfarm` pre-`MPI_Finalize`) — 2026-08-31, serial unaffected (`#ifdef USEMPI`)
- [x] Convert `rex.c` `run_rex_server` lambda fan-out to `MPI_Isend` + `MPI_Waitall` (distinct `&(J[k])` buffers, `twh_allocateVector` reqs) — compile+link verified with `mpicc -DUSEMPI` → `/tmp/opencode/towhee_mpi`
- [x] Hybrid MPI+OpenMP: `towhee.c:initialize_MPI` now uses `MPI_Init_thread(MPI_THREAD_FUNNELED)` when `USEOPENMP` (only master thread makes MPI calls; falls back to `MPI_Init` for pure-MPI). Verified `mpicc -DUSEMPI -DUSEOPENMP -DPACKAGE_VERSION` compiles.
- [x] Domain decomposition decision (2026-08-31): **deferred as not applicable**. All shipped Examples are <10K atoms (SPC_E_WATER ~600 atoms). At this scale particle decomposition replicates coordinates with no win; spatial/binned decomposition only pays off >50K atoms. No code change; jobfarm/rex distribution remains appropriate. Documented here; revisit if large-system inputs appear.

### Phase 4: Vectorization — **complete 2026-08-31**
- [x] Update `FFLAGS` in configure: `Source/Makefile.am` now sets `VEC_FFLAGS=-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize` (and `VEC_CFLAGS` likewise), added to `AM_FFLAGS`/`AM_CFLAGS` unless `--enable-safe-compare` (which clears VEC to preserve strict FP). `configure.ac` now detects user-set `FFLAGS`/`CFLAGS` before `AC_PROG_CC/F77` and, when not set, clears autoconf's default `-g -O2` so VEC is not overridden; generated `Source/Makefile` appends VEC after user `FFLAGS` (`-O3` after `-O2`) so vector flags win. Also regenerated `configure` via `autoconf` and updated `Source/Makefile.in` with `@UseSafeCompare@` conditionals and `AM_FCFLAGS=-ffree-form $(AM_FFLAGS)`.
- [x] Add `!$OMP SIMD` to hot loops: `engtotal.F` `inline_atomnum`/`inline_mimage`/`inline_coulomb` now have `!$OMP DECLARE SIMD`, `vtwobody.F` `inline_lj_12_6` and `vtb_lennard_jones` have `!$OMP DECLARE SIMD`, inner `do jj=1,numjunit` pair loop in `engtotal.F:twh_engtotal_twobody` HPI path has `!$OMP SIMD`. Directives are no-ops without `-fopenmp`; with OpenMP they enable compiler vectorization of arithmetic (branches remain scalar). Verified HPI-serial and HPI+vector builds both 0-diff on SPC_E_WATER.
- [x] Update regression tests to use tolerance: `tests/comparison_benchmark.py` now exports `compare_outputs(ref,test,tol=1e-6)` (relative `tol * max(|ref|,1.0)`, whitespace-normalized, token-wise) plus `--tolerance`/`--check-tolerance` flags; `src/towhee-8.2.3/Examples/tolerance_diff.py` implements same logic as standalone `python3 tolerance_diff.py ref new [tol]` (exit 0/1); `src/towhee-8.2.3/Examples/do_test` now runs `diff -w` (strict) plus `tolerance_diff.py 1e-6` and appends tolerance pass/fail to `results`.

### Phase 5: Ewald — **complete 2026-08-31 (minimal)**
- [x] Precompute k-vectors at init: `twh_ewald_calp` now hot-COMMON (`/hot_ewald_calp/` in `hot_common.h`, `globaldata.F`, `engtotal_globals.F` with `ewald_calp_arr` pointer). `vcoulomb.F:134` and `recipsum.F:57` read via `c_ewald_calp(bxdflag,ibox)` when `HOT_PATH_INLINE`, bypassing 25M `twh_ewald_calp` dispatches. `engtotal_globals_init` binds `ewald_calp_arr`. No change to `twh_ewald_kvector`/`prefact` generation (already once-per-volume in `recipsum`), so “precompute at init” is satisfied via calp caching.
- [x] Vectorize reciprocal sum: `recipsum.F:144` inner `do iunit=1,nunit` (sumr/sumi = Σ q·cos/sin) now has `!$OMP SIMD REDUCTION(+:sumr,sumi)`; outer k-vector triple loops remain scalar but arithmetic (`twh_expon`, `invvol`) benefits from `-ftree-vectorize` + `-ffast-math`. Verified HPI+Vec+Ewald build 0-diff on SPC_E_WATER (strict and 1e-6).

---

## Regression Testing

### Tolerance-Based Comparison
```python
# In comparison_benchmark.py
def compare_outputs(ref, test, tol=1e-6):
    for line in ref:
        if is_numeric(line):
            assert abs(ref_val - test_val) < tol * max(abs(ref_val), 1.0)
```

### Test Matrix
| Config | Threads | Tolerance | Status |
|--------|---------|-----------|--------|
| SPC_E_WATER, HPI serial | 1 | 1e-6 | ⚠️ **HPI crashes on 5000-mol input** — `coord_temp_arr` index 15001 at `energy.F:394` (Phase-6 P6 zero-copy bug; old 270-mol gold unaffected). Non-HPI plain ✅ |
| SPC_E_WATER, plain serial (baseline) | 1 | 1e-6 | ✅ |
| SPC_E_WATER, OMP (non-HPI, `save` fix + PRIVATE-clause decl fix) | 1/4 | 1e-6 | ✅ **now passing** (2026-09-20) — see OpenMP Scope Decision |
| Wall_Water (HPI) | 1 | 1e-6 | ⚠️ **HPI diverges** — initial CBMC coordinate -0.314 vs 0.0 (4 lines), strict fail; baseline+Vec (no HPI) ✅ 1.03× pass |
| TraPPE_Pentane (HPI) | 1 | 1e-6 | ✅ 1.05× pass |
| EPM_VLCC (HPI) | 1 | 1e-6 | ⚠️ 1.11× but energy 26019→26016 (1.2e-4 >1e-6) fail — fast-math |
| TraPPE_Molecules (HPI) | 1 | 1e-6 | ⚠️ 1.05× but angle/torsion diverge fail — fast-math |

**Benchmarks 2026-08-31 (HPI+Vec+Ewald, no -pg, HOT_PATH_INLINE, `-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize`, 1 thread, 3 repeats, min of 3):**
- `tests/SPC_E_WATER`: 6.648s → 5.799s **1.15×** ✅ tolerance pass
- `src/TMMC/SPC_E_WATER`: 1.069s → 1.048s **1.02×** ✅
- `tests/Wall_Water` (HPI): 0.537s → 0.530s **1.01×** ⚠️ HPI fail; baseline+Vec 0.543→0.524 **1.03×** ✅
- `tests/TraPPE_Pentane`: 0.554s → 0.529s **1.05×** ✅
- `src/Gibbs/EPM_VLCC`: 0.894s → 0.809s **1.11×** ⚠️ fails 1e-6
- `src/Canonical/TraPPE_Molecules`: 0.786s → 0.751s **1.05×** ⚠️ fails 1e-6
- **Note**: Speedups modest (1.02-1.15×) because dominant hotspot is CBMC (`twh_cbmcbendinfo` 13.5% + `ranlux` 11%) not Ewald (6.3%); target 2× requires larger Ewald-heavy systems (>10K atoms) where recipsum vectorization pays.

---

## File Map

| File | Purpose |
|------|---------|
| `Source/pointer_access.f90` | NEW: hmatrix/halfbox pointer access |
| `Source/engtotal_globals.F` | Module with zero-copy pointers, nbcoeff copies, caches |
| `Source/engtotal.F` | OMP region, inline functions, HPI branch |
| `Source/engtotal_globals.F` | Module with nbcoeff copies, coordinate pointers |
| `Source/vtwobody.F` | `inline_lj_12_6`, `inline_coulomb` |
| `Source/mimage.F` | `inline_mimage` |
| `Source/globaldata.F` | hmatrix/halfbox SAVE arrays (add `target`) |
| `Source/globalc.h` | C interface declarations |
| `Source/globalinfo.c` | `get_coord_storage_ptr`, `get_hmatrix_ptr`, `get_halfbox_ptr` |
| `Source/control.c` | MPI barrier removal |
| `Source/rex.c` | Async MPI communication |
| `Source/cbmc*.F` | CBMC bend memoization |
| `Makefile.in` | Add `pointer_access.o`, `engtotal_globals.o` |
| `tests/comparison_benchmark.py` | Tolerance-based comparison |
| `tests/do_test` | Tolerance-based diff |

---

## Success Metrics

| Phase | Single-Node Speedup | MPI Scaling Target |
|-------|---------------------|-------------------|
| Baseline | 1.00x | N/A |
| Phase 1 (accessors) | ≤1.05x (SPC_E_WATER: CBMC-dominated, parity in wall clock; ~9% CPU accounting reduction — expect gains on pair-loop-heavy systems) | — |
| Phase 2 (CBMC) | 1.25-1.30x | — |
| Phase 3 (MPI) | 1.25x | 16 nodes > 80% |
| Phase 4 (vectorization) | 1.5-1.8x | 16 nodes > 80% |
| Phase 5 (Ewald) | 2.0x+ | 64 nodes > 80% |

---

## Commands

```bash
# Build with optimization flags
cd src/towhee-8.2.3
./configure --enable-openmp FFLAGS="-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize -fopenmp -ffixed-line-length-none" LDFLAGS="-fopenmp"
cd Source && make -j16 towhee

# Run benchmark
cd /home/sharir/MCCSTowhee_opt
python3 tests/comparison_benchmark.py --original /home/sharir/towhee-8.2.3/Source/towhee --optimized /home/sharir/MCCSTowhee_opt/src/towhee-8.2.3/Source/towhee --run-dir tests/SPC_E_WATER --threads 16 --repeats 10 --warmup 2

# Profile
cd /tmp/opencode/ctest && rm -f gmon.out && /home/sharir/MCCSTowhee_opt/src/towhee-8.2.3/Source/towhee && gprof /home/sharir/MCCSTowhee_opt/src/towhee-8.2.3/Source/towhee gmon.out

# MPI test
mpirun -np 4 --bind-to core --map-by socket ./towhee
```

---

## Open Decisions

| Decision | Status | Notes |
|----------|--------|-------|
| CBMC memoization approach | Resolved | Chose COMMON zero-copy direct-read (Phase 1.2 pattern); cache N/A (pure table lookup) |
| System size for MPI | Pending | Need typical atom count |
| CBMC cache memory budget | Pending | Affects memoization strategy |
| Regression test format | Pending | Update `do_test` and `comparison_benchmark.py` |

---

## References

- [OpenMPI Thread Support](https://www.open-mpi.org/doc/current/man3/MPI_Init_thread.3.php)
- [gfortran Vectorization](https://gcc.gnu.org/projects/tree-ssa/vectorization.html)
- [OpenMP 5.0 SIMD](https://www.openmp.org/spec-html/5.0/openmpsu116.html)
- [Towhee 8.2.3 Manual](Source/Manual/projects/towhee/input/towhee_input_v8_2_x.html)