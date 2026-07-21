## OBJECTIVE
Eliminate the Fortran-C dispatch-call-chain overhead in the O(N²) nonbonded energy evaluation loop by bypassing the `twh_*(GLB_GET, ...)` flag-dispatch accessors in favor of direct module-POINTER array reads. Target: reduce dispatch-layer operations from ~8-15 per atom pair to 0, while retaining the physics computation calls (erfc, energy functions) as module-internal pure functions. The optimization adds a bypass path — the original `twh_*` accessors remain intact for non-hot-path callers.

## CONTEXT
The hot path in `engtotal.F:334-458` evaluates all intermolecular atom pairs. Each pair currently traverses:
```
twh_coordinates(GLB_GET, ...) → twh_atomnum → twh_globalpos → C twh_coordstorage_ → C twh_1ddouble (switch dispatch)
twh_vtwobody(...)              → twh_nbcoeff(GLB_GET, ...) called 3-15× (each: globaldata.F switch on flag)
twh_vcoulomb(...)              → twh_qqscale(GLB_GET, ...)   (same dispatch pattern)
twh_ntype, twh_luse_c/v        → globaldata.F switch on flag
```

For a simulation with 1000 molecules × 100 atoms each, this is ~5×10⁹ atom pairs × ~8-15 call overheads = **40-75 billion subroutine call + switch-dispatch operations** spent on data access, not on physics.

The fix: create a Fortran module with `POINTER` aliases to the specific `SAVE`d arrays used in the hot path, and access them directly by index.

## TARGET FILES
- **NEW:** `Source/engtotal_globals.F` — Fortran module with `POINTER` aliases (data only)
- **NEW:** `Source/engtotal_inlines.F` — Module with inline pure functions for LJ and Coulomb (computation only)
- `Source/engtotal.F:334-458` — Hot loop body (modify ~50 lines)
- `Source/globaldata.F` — Add `TARGET` attribute to ~8 key SAVEd arrays
- `Source/vtwobody.F:849-882` — Replace `twh_nbcoeff(GLB_GET, ...)` with module array access
- `Source/vcoulomb.F:100-140` — Replace `twh_qqscale(GLB_GET, ...)` with module array access
- `Source/allocate.F` — Add module pointer initialization after allocations
- `Source/Makefile.am:44` — Add `engtotal_globals.F` and `engtotal_inlines.F` to `towhee_SOURCES`
- `Source/energy_change.F` — Apply same inline pattern when this function dispatches to engtotal sub-loops

## ARCHITECTURE

### New Module: `engtotal_globals.F`
```
module engtotal_globals
  implicit none
  save

  ! Coordinate arrays (SoA layout for cache-friendly access)
  real(8), pointer, contiguous :: crd_x(:), crd_y(:), crd_z(:)
  real(8), pointer, contiguous :: crd_temp_x(:), crd_temp_y(:), crd_temp_z(:)

  ! Nonbonded coefficients (flattened 1D for stride-1 access)
  real(8), pointer, contiguous :: nbcoeff_0(:)    ! shift coeffs
  real(8), pointer, contiguous :: nbcoeff_13(:)   ! sigma^(1/6)
  real(8), pointer, contiguous :: nbcoeff_14(:)   ! epsilon
  real(8), pointer, contiguous :: nbcoeff_15(:)   ! sigma^(1/6) for 1-4
  real(8), pointer, contiguous :: nbcoeff_16(:)   ! epsilon for 1-4

  ! Atom type data
  integer, pointer, contiguous :: ntype_arr(:,:)
  integer, pointer, contiguous :: nboxi_arr(:)
  integer, pointer, contiguous :: moltyp_arr(:)
  integer, pointer, contiguous :: nunit_arr(:)
  integer, pointer, contiguous :: globalpos_arr(:)

  ! Interaction site flags
  logical, pointer, contiguous :: luse_c_arr(:,:)
  logical, pointer, contiguous :: luse_v_arr(:,:)

  ! Charge and scaling data
  real(8), pointer, contiguous :: qqatom_arr(:,:)
  real(8), pointer, contiguous :: qqscale_arr(:,:)

  ! Cutoff data
  real(8), pointer :: rcut_val, rcelect_val
  logical, pointer :: lewald_val, lelect_val, lshift_val

  ! Box matrix (for minimum image)
  real(8), pointer, contiguous :: hmatrix_arr(:,:,:)
  real(8), pointer, contiguous :: hinverse_arr(:,:,:)

  ! Box pair content
  integer, pointer, contiguous :: chainlist_arr(:,:,:)
  integer, pointer, contiguous :: nboxpair_arr(:)
  integer, pointer, contiguous :: pairbox_arr(:,:)

  ! Lambda scaling
  real(8), pointer :: lambda_lj_val, lambda_c_val
  integer, pointer :: scp_num_foreign_lambda_val
  real(8), pointer, contiguous :: foreign_lambda_lj_arr(:)
  real(8), pointer, contiguous :: foreign_lambda_c_arr(:)
  logical, pointer :: scp_scale_twobody_val

end module engtotal_globals
```

### New Module: `engtotal_inlines.F`
Separate computation functions from the data module. This module contains pure functions for LJ and Coulomb energy evaluation that read from `engtotal_globals` pointers:

```fortran
module engtotal_inlines
  use engtotal_globals
  implicit none
contains

  real(8) function inline_lj_12_6(rijsq, ntij, lonefour) result(energy)
    real(8), intent(in) :: rijsq
    integer, intent(in) :: ntij
    logical, intent(in) :: lonefour
    real(8) :: sixterm, prefactor, srsix

    if (lonefour) then
      sixterm  = nbcoeff_15(ntij)
      prefactor = nbcoeff_16(ntij)
    else
      sixterm  = nbcoeff_13(ntij)
      prefactor = nbcoeff_14(ntij)
    endif
    srsix = sixterm / (rijsq * rijsq * rijsq)
    energy = srsix * (srsix - 1.0d0) * prefactor
    if (lshift_val) energy = energy - nbcoeff_0(ntij)
  end function inline_lj_12_6

  real(8) function inline_coulomb(rij, qqij, calp, lewald) result(energy)
    real(8), intent(in) :: rij, qqij, calp
    logical, intent(in) :: lewald
    real(8) :: erfc_val

    if (lewald) then
      erfc_val = erfc(calp * rij)
      energy = qqij * erfc_val / rij
    else
      energy = qqij / rij
    endif
  end function inline_coulomb

end module engtotal_inlines
```

Note: The functions are regular (not `PURE`) because they access module-level pointers. Fortran `PURE` requires all data access through `INTENT(IN)` arguments, which would defeat the purpose of eliminating argument passing. The module pointers are read-only during simulation, so the functions are side-effect-free in practice.

### Pointer Wiring in `allocate.F`
After all `twh_nbcoeff`, `twh_coordinates`, etc. arrays are allocated, add a call to `engtotal_globals_init()` which associates module pointers to the TARGET arrays in `globaldata.F`:

```fortran
subroutine engtotal_globals_init()
  use engtotal_globals
  use iso_c_binding, only: c_double, c_f_pointer, c_ptr
  use globaldata, only: natoms
  type(c_ptr) :: cptr
  real(c_double), pointer :: coord_flat(:)

  ! Coordinate arrays: retrieve C pointer and bind via c_f_pointer
  call get_coord_storage_ptr(CRD_REAL, cptr)
  call c_f_pointer(cptr, coord_flat, [3 * natoms])
  ! SoA decomposition: flat array stores x0,y0,z0, x1,y1,z1, ...
  ! Access pattern: crd_x(i) = coord_flat(3*(i-1)+1), etc.
  crd_x_ptr => coord_flat(1::3)   ! Fortran pointer alias with stride
  crd_y_ptr => coord_flat(2::3)
  crd_z_ptr => coord_flat(3::3)

  ! ... same for crd_temp_x/y/z using CRD_TEMP ...
  ! ... same for nbcoeff, ntype, luse_c, etc. using TARGET arrays ...
end subroutine
```

NOTE: Fortran `POINTER` cannot directly stride a C pointer via `=>` with stride syntax on a section. Implementation approach: use a shim array `coord_flat(0:3*natoms-1)` with `TARGET`, set via a single `memcpy` from the C heap on initialization and on each accept/reject that modifies coordinates. This adds one memcpy per accept, but eliminates per-call overhead.

### Expected Memory Layout (Hot Loop Access Pattern)

Per atom pair (i, j), the ideal access pattern after optimization:
```
L1 cache line 1: crd_x(i), crd_x(j)     — 16 bytes (2 doubles)
L1 cache line 2: crd_y(i), crd_y(j)     — 16 bytes
L1 cache line 3: crd_z(i), crd_z(j)     — 16 bytes
L1 cache line 4: ntype_arr(molty_i, ii), ntype_arr(molty_j, jj) — 8 bytes
L1 cache line 5: nbcoeff_13(ntij)        — 8 bytes (single read, triangular index)
```
Total: ~64 bytes per pair ≈ 1 cache line. Original path: ~64 bytes + 8-15× subroutine stack frames + switch branches.

## STEP-BY-STEP IMPLEMENTATION PLAN

### Phase 1: Create Module and Wire Arrays (Estimated: 3 hours)

#### Step 1.1 — Create `Source/engtotal_globals.F` and `Source/engtotal_inlines.F`
Write the Fortran modules as specified above in ARCHITECTURE. `engtotal_globals.F` is the data-only module with `POINTER` aliases (every `POINTER` uses `CONTIGUOUS` attribute). `engtotal_inlines.F` contains the inline pure functions that read from the data module.

Constraint: Do NOT introduce any USE statements in `globaldata.F`. The modules are consumed only by the hot-path functions (`engtotal.F`, `vtwobody.F`, `vcoulomb.F`), not by the global data layer.

#### Step 1.2 — Add `TARGET` to Key Arrays in `globaldata.F`
Locate each SAVEd array accessed from the hot path and add `TARGET`:

| Subroutine | Array | Line (approx.) | Change |
|-----------|-------|---------------|--------|
| `twh_nbcoeff` | `nbcoeff(MINNBCOEFF:MAXNBCOEFF, CROSSTYPEMAX)` | `globaldata.F:10685` | `dimension, target :: nbcoeff(...)` |
| `twh_ntype` | `ntype(NTMAX, NUMAX)` | `globaldata.F:5162` | `integer, save, target :: ntype(...)` |
| `twh_nboxi` | `nboxi(NTMAX, MAXCHAIN)` | `globaldata.F:6604` | Add TARGET |
| `twh_moltyp` | `moltyp(MAXCHAIN)` | `globaldata.F:6990` | Add TARGET |
| `twh_nunit` | `nunit(NTMAX)` | `globaldata.F:7641` | Add TARGET |
| `twh_luse_c` | `luse_c(NTMAX, NUMAX)` | `globaldata.F:6461` | Add TARGET |
| `twh_luse_v` | `luse_v(NTMAX, NUMAX)` | `globaldata.F:6491` | Add TARGET |
| `twh_rcut` | `rcut` | `globaldata.F:9177` | Add TARGET |
| `twh_lshift` | `lshift` | `globaldata.F:6417` | Add TARGET |
| `twh_hmatrix` | `hmatrix(3,3,MAXBOX)` | `globaldata.F:5342` | Add TARGET |

Constraint: Only add TARGET to the variable declaration line. Do NOT change any logic in the flag-dispatch body. The existing `twh_*` subroutines remain unchanged for use by other parts of the code.

#### Step 1.3 — Initialize Module Pointers in `allocate.F`
Create subroutine `engtotal_globals_init()` called once after all allocations complete. This associates every module POINTER with the TARGET array it shadows:

```fortran
subroutine engtotal_globals_init()
  use engtotal_globals
  ! nbcoeff: extract slices for hot-path coefficients
  ! Method: use TARGET arrays and pointer assignment to sections
  ! nbcoeff is declared as: real(8), save, target :: nbcoeff(MINNBCOEFF:MAXNBCOEFF, CROSSTYPEMAX)
  ! Pointers alias specific coefficient rows for stride-1 access
  nbcoeff_0  => nbcoeff(0,  1:CROSSTYPEMAX)
  nbcoeff_13 => nbcoeff(13, 1:CROSSTYPEMAX)
  nbcoeff_14 => nbcoeff(14, 1:CROSSTYPEMAX)
  nbcoeff_15 => nbcoeff(15, 1:CROSSTYPEMAX)
  nbcoeff_16 => nbcoeff(16, 1:CROSSTYPEMAX)
  ...
end subroutine
```

For the coordinate arrays (which live as C heap pointers in `double_pointer[]`), the implementation depends on C-Fortran interop. Create a C function `get_coord_storage_ptr(int crdflag, double **ptr)` returning the base address, then use `c_f_pointer` in Fortran to bind:

```fortran
use iso_c_binding, only: c_double, c_f_pointer, c_ptr
type(c_ptr) :: cptr
real(c_double), pointer :: coord_flat(:)
call get_coord_storage_ptr(CRD_REAL, cptr)
call c_f_pointer(cptr, coord_flat, [3 * natoms])
! coord_flat stores interleaved: x0,y0,z0, x1,y1,z1, ...
! Crd arrays use strided pointer aliasing (see MEMORY LAYOUT section)
```

Constraint: The coordinate pointer must be re-queried after any `CRD_T_TO_R` swap (accept step), because the underlying C pointer changes. Add a hook in the accept path.

### Phase 2: Rewrite the Hot Loop Body (Estimated: 4 hours)

#### Step 2.1 — Rewrite `engtotal.F:385-452` (inner atom-pair loop)

**Before:**
```fortran
do ii = 1, numiunit
  call twh_ntype(GLB_GET, imolty, ii, ntii)
  call twh_coordinates(GLB_GET, crdflag, ichain, ii, rxui, ryui, rzui)
  call twh_luse_c(GLB_GET, imolty, ii, luse_ci)
  call twh_luse_v(GLB_GET, imolty, ii, luse_vi)
  do jj = 1, numjunit
    call twh_ntype(GLB_GET, jmolty, jj, ntjj)
    call twh_coordinates(GLB_GET, crdflag, jchain, jj, rxuj, ryuj, rzuj)
    call twh_luse_c(GLB_GET, jmolty, jj, luse_cj)
    call twh_luse_v(GLB_GET, jmolty, jj, luse_vj)
    ...
    rijsq = ...
    vtwointer = vtwointer + twh_vtwobody(lfailure, rijsq, ntii, ntjj, ...)
    velectinter = velectinter + vcoulomb
```

**After:**
```fortran
use engtotal_globals, only: crd_x, crd_y, crd_z, ntype_arr, luse_c_arr, luse_v_arr

! Select coordinate set based on crdflag
if (crdflag == CRD_REAL) then
  x_ptr => crd_x; y_ptr => crd_y; z_ptr => crd_z
else
  x_ptr => crd_temp_x; y_ptr => crd_temp_y; z_ptr => crd_temp_z
endif

do ii = 1, numiunit
  ntii = ntype_arr(imolty, ii)
  idx_i = globalpos_arr(ichain) + (ii - 1)
  rxui = x_ptr(idx_i); ryui = y_ptr(idx_i); rzui = z_ptr(idx_i)
  luse_ci = luse_c_arr(imolty, ii)
  luse_vi = luse_v_arr(imolty, ii)
  do jj = 1, numjunit
    ntjj = ntype_arr(jmolty, jj)
    idx_j = globalpos_arr(jchain) + (jj - 1)
    rxuj = x_ptr(idx_j); ryuj = y_ptr(idx_j); rzuj = z_ptr(idx_j)
    luse_cj = luse_c_arr(jmolty, jj)
    luse_vj = luse_v_arr(jmolty, jj)
    ...
```

#### Step 2.2 — Inline `twh_vtwobody` in `vtwobody.F` for LJ 12-6 (the dominant case)

Create a new function `inline_lj_12_6(rijsq, ntij, lonefour)` in `engtotal_inlines.F`:

```fortran
real(8) function inline_lj_12_6(rijsq, ntij, lonefour) result(energy)
  use engtotal_globals, only: nbcoeff_13, nbcoeff_14, nbcoeff_15, nbcoeff_16, &
                               nbcoeff_0, lshift_val
  real(8), intent(in) :: rijsq
  integer, intent(in) :: ntij
  logical, intent(in) :: lonefour
  real(8) :: sixterm, prefactor, srsix

  if (lonefour) then
    sixterm  = nbcoeff_15(ntij)
    prefactor = nbcoeff_16(ntij)
  else
    sixterm  = nbcoeff_13(ntij)
    prefactor = nbcoeff_14(ntij)
  endif
  srsix = sixterm / (rijsq * rijsq * rijsq)
  energy = srsix * (srsix - 1.0d0) * prefactor
  if (lshift_val) energy = energy - nbcoeff_0(ntij)
end function
```

Constraint: This function must have no side effects — no I/O, no writes to module pointers (which are read-only during simulation). The function is not declared `PURE` because it reads module-level pointers, but it is side-effect-free in practice. Verify `lshift_val` is read-only during simulation.

#### Step 2.3 — Inline `twh_vcoulomb` in `vcoulomb.F`

Create `inline_coulomb(rij, imolty, ii, jmolty, jj)` that reads `qqatom_arr` and `qqscale_arr` directly. For the Ewald short-range case, use the Fortran intrinsic `erfc()` (available in F2008) instead of the custom `twh_derfc()`.

Constraint: The Coulomb inline function must check `lelect_val` and `lewald_val` (module pointers to read-only SAVEd logicals). If `lewald_val` is true, apply `erfc(calp*rij)` scaling.

#### Step 2.4 — Replace pure function calls in `engtotal.F` hot loop

```fortran
if (leval_v) then
  ntij = getnbtype(ntii, ntjj)
  energy_vdw = inline_lj_12_6(rijsq, ntij, lonefour)
  vtwointer = vtwointer + energy_vdw
endif
if (lelect_val .and. leval_c) then
  rij = dsqrt(rijsq)
  vcoulomb = inline_coulomb(rij, imolty, ii, jmolty, jj, ntii, ntjj)
  velectinter = velectinter + vcoulomb
endif
```

Constraint: `twh_mimage` (minimum image) remains a regular subroutine call — it is a pure geometric computation with no data access overhead. Inlining it is not necessary.

Constraint: `twh_getnbtype` (ntii,ntjj → ntij triangular index map) is pure integer arithmetic and is efficient as-is. Keep the existing `getnbtype.F` interface but call it directly without the `twh_` dispatch layer if possible.

### Phase 3: Coordinate Refresh on Accept (Estimated: 1 hour)

#### Step 3.1 — Hook in Accept Path
After every MC accept, the C coordinate pointer swaps (`CRD_REAL` and `CRD_TEMP` exchange via `twh_coordstorage_` pointer swap). After this swap, the module pointers must be updated:

```fortran
! In engtotal_inlines.F or a new subroutine engtotal_refresh_coords()
subroutine engtotal_refresh_coords()
  use engtotal_globals
  if (associated(crd_x)) call refresh_coord_ptr(CRD_REAL, crd_x, crd_y, crd_z)
  if (associated(crd_temp_x)) call refresh_coord_ptr(CRD_TEMP, crd_temp_x, crd_temp_y, crd_temp_z)
end subroutine
```

Call this in every accept path immediately after `CRD_T_TO_R`.

Constraint: Verify that the pointer refresh does not introduce race conditions in the OpenMP parallel case (Optimization 1). The refresh should happen outside the parallel region, between MC moves.

#### Coordinate Refresh Audit Checklist
Every file that calls `CRD_T_TO_R` must have a subsequent `engtotal_refresh_coords()` call:

| File | Function | Accept path location | Refresh call added? |
|------|----------|---------------------|-------------------|
| `trancom.F` | translation move | ~line 250 (after CRD_T_TO_R) | ☐ |
| `rotate.F` | rotation move | ~line 218-263 | ☐ |
| `tranatom.F` | atom translation | (find accept block) | ☐ |
| `swapmoves.F` | identity swap | (find accept block) | ☐ |
| `volumemoves.F` | volume change | (find accept block) | ☐ |
| `insert.F` | insertion | (find accept block) | ☐ |
| `delete.F` | deletion | (find accept block) | ☐ |

Add a grep-based CI check:
```bash
# Verify every CRD_T_TO_R has a subsequent refresh call
grep -rn "CRD_T_TO_R" Source/ | while read line; do
  file=$(echo "$line" | cut -d: -f1)
  linenum=$(echo "$line" | cut -d: -f2)
  if ! grep -q "engtotal_refresh_coords" "$file"; then
    echo "ERROR: $file calls CRD_T_TO_R but not engtotal_refresh_coords"
    exit 1
  fi
done
```

### Phase 4: Update energy_change.F (Estimated: 1.5 hours)

#### Step 4.1 — Find equivalent hot loops in `energy_change.F`
The file `energy_change.F` contains `ec_two_body_dv` which evaluates pair energy changes for moved molecules vs all others. It has the same nested structure as `engtotal.F` but over (moved molecule, all others) instead of all pairs. Apply identical inlining.

Search for the pattern:
```
do ichain = 1, nchainmove
  ...
  call twh_coordinates(GLB_GET, ...)  ! target
```

Replace with direct module array access.

Constraint: `energy_change.F` uses `chaininfo` (array of flags per moved molecule) to determine whether to use OLD (CRD_REAL) or NEW (CRD_TEMP) coordinates. The coordinate pointer selection must be done per-moved-molecule, not once for the entire loop. Implement this with a pointer array:
```fortran
x_ptr_chain(1:nchainmove)  ! array of pointers to x coords for each chain
y_ptr_chain(1:nchainmove)
z_ptr_chain(1:nchainmove)
```
Set per-chain before the outer loop based on `chaininfo(:,2:3)` flags.

### Phase 5: Build System (Estimated: 0.5 hours)

#### Step 5.1 — Update `Makefile.am`
Add `engtotal_globals.F` and `engtotal_inlines.F` to `towhee_SOURCES` in alphabetical position (after `energy_change.F`, before `engatom.F`):
```
engtotal_globals.F \
engtotal_inlines.F \
```

Constraint: Ensure the modules are compiled before `engtotal.F`, `vtwobody.F`, and `vcoulomb.F` which USE them. GNU `make` with `.mod` dependency tracking handles this automatically if the module files are listed first, but add explicit dependencies:
```
engtotal.o vtwobody.o vcoulomb.o: engtotal_globals.o engtotal_inlines.o
```

#### Step 5.2 — Add Preprocessor Guard

In `engtotal_globals.F`, `engtotal_inlines.F`, and all files that USE them, wrap the inline path with a preprocessor guard:

```fortran
#ifndef HOT_PATH_INLINE
#define HOT_PATH_INLINE 0
#endif

#if HOT_PATH_INLINE
  use engtotal_globals, only: crd_x, crd_y, crd_z, ...
  ! inline path — direct array access
#else
  ! original dispatch path (unchanged, uses twh_* subroutines)
#endif
```

In `Makefile.am`, add a configure-time option:
```makefile
AC_ARG_ENABLE([hot-path-inline],
  [AS_HELP_STRING([--enable-hot-path-inline],
    [Enable inline hot-path optimization @<:@default=no@:>@])],
  [hot_path_inline=$enableval],
  [hot_path_inline=no])

if test "x$hot_path_inline" = "xyes"; then
  AM_CPPFLAGS += -DHOT_PATH_INLINE=1
else
  AM_CPPFLAGS += -DHOT_PATH_INLINE=0
endif
```

This allows A/B comparison without code changes:
```bash
./configure --enable-hot-path-inline   # optimized build
./configure                            # baseline build
```

## RISK & ROLLBACK

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| C coordinate pointer stride layout incompatible with Fortran POINTER | Medium | High — blocks Phase 1 | Fallback: copy coordinates to a Fortran-allocated flat array (`coord_flat(1:3*natoms)`) with `TARGET`. Accept one memcpy per step; still eliminates dispatch overhead. |
| Pointer aliasing violates Fortran aliasing rules under aggressive optimization (ifort -O3) | Low | Medium — silent wrong results | Build both gfortran and ifort; compare energies. Add `-fcheck=all` in debug builds. |
| `TARGET` attribute on SAVEd arrays causes unexpected codegen in `globaldata.F` | Low | Low — compile error only | The compiler will reject invalid TARGET placement; fix is mechanical. |
| Coordinate refresh missed in an accept path → stale coordinates | Medium | High — wrong energies | Add a debug-mode assertion: after each energy eval, verify `crd_x(i) == twh_coordinates(GLB_GET, ...)` for a random sample of atoms. |
| Performance gain < 1.5× | Medium | Medium — effort wasted | Profile before starting. If dispatch overhead < 30% of hot-loop time, pivot to loop tiling or vectorization instead. |

## TESTING & VALIDATION

### Regression Test
1. Build without the optimization (`HOT_PATH_INLINE=0`), run a reference simulation (e.g. NPT N₂, 500 molecules, 2000 MC steps). Record all per-step energies.
2. Build with the optimization (`HOT_PATH_INLINE=1`), run the same simulation with same random seed.
3. Compare per-step energies. Must match within 1×10⁻¹² relative error (machine epsilon for double precision).

### Performance Test
1. Benchmark a large system (e.g., 10,000-atom SPC/E water, 10,000 MC steps) with and without the optimization.
2. Profile with `perf stat` or `gprof` — measure reduction in function call counts in the hot path.
3. Expected result: 1.5-3× speedup on the nonbonded energy evaluation, 20-50% reduction in total simulation time.

### Memory Test
1. Run with `valgrind --tool=memcheck` to detect any dangling pointer references introduced by the POINTER aliasing.
2. Run with `-fbounds-check` (gfortran) to detect any array index out-of-bounds from incorrect index math.

### Index Map Validation Test
Write a unit test `test_getnbtype.F` that:
1. Iterates all `(ntii, ntjj)` pairs in `[1, NNTYPE]`
2. Computes `ntij` via the inline formula from Constraint #5
3. Asserts `ntij == getnbtype(ntii, ntjj)` (existing function)
4. Asserts `ntij` is in `[1, CROSSTYPEMAX]`
5. Runs as part of `make check` before and after optimization

## CONSTRAINTS & GUARDRAILS

1. **Do NOT remove or alter the original `twh_*` accessor subroutines** in `globaldata.F`. Other parts of the code (not in the hot path) still use them. The optimization adds a bypass, not a replacement.

2. **Do NOT add `USE engtotal_globals` to `globaldata.F`** or to any file outside the direct hot path. The module is a private optimization detail.

3. **Maintain bitwise reproducibility.** The inline functions must compute identical results to the original dispatch-based calls. Verify with `SAFE_COMPARE` mode.

4. **The coordinate pointer refresh on accept** must be called before every energy evaluation that follows a coordinate change. Audit all accept paths: `trancom.F`, `rotate.F`, `tranatom.F`, `swapmoves.F`, `volumemoves.F`.

5. **The `nbcoeff` extraction must use the correct `ntij` index** — the triangular map from `getnbtype.F:22-69` computes `ntij` from `(ntii, ntjj)` as:
   ```
   low = min(ntii, ntjj)
   high = max(ntii, ntjj)
   istep = (low-1) * NNTYPE - ((low-1)*(low-2))/2
   ntij = istep + high - low + 1
   ```
   This is the same index used by the original `twh_nbcoeff` array. Verify with a precomputed lookup table test.

6. **The `lonefour` flag** (1-4 interactions) is computed elsewhere (from bonding topology). The inline function must accept it as a parameter — same as the original `vtb_lennard_jones`.

7. **The `lovrlap` check** (rijsq < rminsq → overlap) and `twh_engtotal_overlap` call must remain. Move the overlap detection inline:
   ```fortran
   if (rijsq < rminsq_val) then
     lovrlap = .true.
     ! overlap recovery logic (preserve existing behavior)
   endif
   ```

8. **Do NOT use Fortran 2008+ features** beyond `CONTIGUOUS` and `c_f_pointer`. The codebase uses F77 + extensions. The module approach is compatible with gfortran -std=gnu and ifort -standard-semantics.

9. **Thread safety (pending Optimization 1).** All module pointers in `engtotal_globals` are read-only during the parallel region. The coordinate refresh (`engtotal_refresh_coords`) must be called outside any `!$OMP PARALLEL` region. Verify with ThreadSanitizer or `valgrind --tool=helgrind` on a small test case.

## SUCCESS CRITERIA

| Metric | Target | Measurement |
|--------|--------|------------|
| Function calls per atom pair in hot loop | ≤2 (was 5-20) | `perf stat` on engtotal function |
| Wall-clock time for nonbonded evaluation | 33-67% reduction | Runtime comparison |
| Per-step energy drift vs baseline | ≤1×10⁻¹² (relative) | Compare final energies |
| Memory increase | ≤1% | `valgrind` heap usage |
| Build time increase | ≤10 seconds | `make` timing |

## POST-DEPLOYMENT MONITORING

1. **Energy conservation drift**: In NVE runs, monitor total energy drift over 10⁶ steps. The optimization should not change drift characteristics.
2. **Regression trigger**: If any CI test shows >1×10⁻¹⁰ relative energy deviation, bisect to identify the offending change.
3. **Performance regression**: Add a `perf stat` step to CI that records branch-miss rate and L1-dcache-load-misses for the hot loop. Flag regressions >10%.
4. **Pointer validity**: In debug builds (`-fbounds-check -fcheck=all`), run periodically to detect dangling pointer references or out-of-bounds access from the aliased arrays.
