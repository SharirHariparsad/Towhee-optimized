# Towhee 8.2.3 Multi-Node CPU Optimization Plan

**Status**: Phase 1.2 complete (Phase 1 done — HPI serial deliverable verified, getter caches in place); OpenMP dropped from scope (see below)  
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

## OpenMP Scope Decision (2026-08-31)

**Decision**: OpenMP is dropped from scope for Phase 1.1. The deliverable is the **HOT_PATH_INLINE (HPI) serial** build, which passes the SPC_E_WATER regression exactly (0 diff lines vs `answer_current`). The OpenMP build is **not usable** — see below.

**Findings (isolated by experiment, all deterministic)**:

1. **The OpenMP build corrupts results from the very first MC move.** The 2nd energy call reads the same CRD_REAL coordinate buffer but at atom 241 instead of 85 (wrong `inline_atomnum`/molecule-to-atom mapping), producing an impossible negative VDW energy. The trajectory then diverges from step 0.
2. **Not the parallel loop.** With the `!$OMP PARALLEL DO` directive disabled (plain serial loop), the `-fopenmp` build still reproduces the identical wrong energy on call 2.
3. **Not `-frecursive`.** Adding `-fno-recursive` (gfortran OpenMP implies `-frecursive`, breaking F77 implicit-SAVE statics) does not change the corruption.
4. **Not an OpenMP-gated code path elsewhere.** `engtotal.F` is the only file with `#ifdef USEOPENMP`; the RNG has no OpenMP branches.
5. **Pre-existing, independent of HPI.** `--enable-openmp` WITHOUT `HOT_PATH_INLINE` also fails: 121 diff lines vs gold at `OMP_NUM_THREADS=1`.

**Conclusion**: the `-fopenmp`/`--enable-openmp` configuration itself deterministically corrupts shared program state (atom indexing) from the first move, in code untouched by the HPI path. Root-causing it (suspected interaction of the gfortran OpenMP runtime with F77 static locals / COMMON memory layout) is a separate, high-effort effort and is explicitly **out of scope** for Phase 1.1.

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

---

## Phase 3: MPI Domain Decomposition (Week 3-4)

### 3.1 Communication Optimization (spec1)
- `Source/control.c`: Remove redundant `MPI_Barrier`
- `Source/rex.c`: `MPI_Send` → `MPI_Isend`/`MPI_Irecv` + `MPI_Waitall`

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
- [ ] Add `bend_cache` / `bend_cached` to cbmc module
- [ ] Modify `twh_cbmcbendinfo` to check cache first
- [ ] Profile: verify `twh_cbmcbendinfo_` calls reduced

### Phase 3: MPI
- [ ] Remove redundant `MPI_Barrier` in `control.c`
- [ ] Convert `rex.c` to `MPI_Isend`/`MPI_Irecv`
- [ ] Implement domain decomposition (await system size clarification)

### Phase 4: Vectorization
- [ ] Update `FFLAGS` in configure
- [ ] Add `!$OMP SIMD` to hot loops
- [ ] Update regression tests to use tolerance

### Phase 5: Ewald
- [ ] Precompute k-vectors at init
- [ ] Vectorize reciprocal sum

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
| SPC_E_WATER, HPI serial | 1 | 1e-6 | ✅ (0 diff, byte-identical) |
| SPC_E_WATER, plain serial (baseline) | 1 | 1e-6 | ✅ |
| SPC_E_WATER, OMP (any) | 1/16 | 1e-6 | ❌ **out of scope** — pre-existing corruption from move 1 (see above) |
| Wall_Water | 1 | 1e-6 | 🕓 pending |

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
| CBMC memoization approach | Pending | Option A (cache) vs B (precompute) |
| System size for MPI | Pending | Need typical atom count |
| CBMC cache memory budget | Pending | Affects memoization strategy |
| Regression test format | Pending | Update `do_test` and `comparison_benchmark.py` |

---

## References

- [OpenMPI Thread Support](https://www.open-mpi.org/doc/current/man3/MPI_Init_thread.3.php)
- [gfortran Vectorization](https://gcc.gnu.org/projects/tree-ssa/vectorization.html)
- [OpenMP 5.0 SIMD](https://www.openmp.org/spec-html/5.0/openmpsu116.html)
- [Towhee 8.2.3 Manual](Source/Manual/projects/towhee/input/towhee_input_v8_2_x.html)