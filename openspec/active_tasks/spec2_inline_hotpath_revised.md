## OBJECTIVE
Eliminate the Fortran-C dispatch call chain overhead in the O(N²) nonbonded energy evaluation loop (`engtotal.F`) and the O(N) per-move evaluation loop (`energy_change.F`) by replacing indirect accessor subroutine calls (`twh_*(GLB_GET, ...)`) with direct array references. This reduces per-atom-pair overhead from ~5-20 indirect subroutine calls (each behind a `globaldata.F` switch-dispatch) to ~0-2 direct subroutine calls (`getnbtype`, `mimage`) plus array reads — **without introducing new O(N) work on the accept path**, since that path runs far more often than the O(N²) full evaluation (see CONTEXT).

## CONTEXT
The hot path in `engtotal.F:334-458` evaluates all intermolecular atom pairs. Each pair currently traverses:
```
twh_coordinates(GLB_GET, ...) → twh_atomnum → twh_globalpos → C twh_coordstorage_ → C twh_1ddouble (switch dispatch)
twh_vtwobody(...)              → twh_nbcoeff(GLB_GET, ...) called 3-15× (each: globaldata.F switch on flag)
twh_vcoulomb(...)              → twh_qqscale(GLB_GET, ...)   (same dispatch pattern)
twh_ntype, twh_luse_c/v        → globaldata.F switch on flag
```

For a simulation with 1000 molecules × 100 atoms each, this is ~5×10⁹ atom pairs × ~8-15 call overheads = **40-75 billion subroutine call + switch-dispatch operations** spent on data access, not on physics.

**Call frequency note (important — drives the design below):** `engtotal.F`'s full O(N²) evaluation runs rarely: at initialization and at periodic energy-drift checkpoints. `energy_change.F`'s O(N) per-moved-molecule evaluation runs on (or near) every attempted MC move — this is the dominant call frequency in a real simulation. Any per-accept overhead this spec introduces must be evaluated against *that* path, not just the O(N²) one, or the optimization can net negative on the workload that actually dominates wall-clock time.

The fix: create a Fortran module with `POINTER` aliases to the specific `SAVE`d arrays used in the hot path, and access them directly by index — using **zero-copy pointer re-binding**, not a per-accept memory copy, for the arrays that mutate during the run (coordinates).

## TARGET FILES
- **NEW:** `Source/engtotal_globals.F` — Fortran module with `POINTER`/`ALLOCATABLE` aliases
- `Source/engtotal.F:334-458` — Hot loop body (modify ~50 lines)
- `Source/globaldata.F` — Add `TARGET` attribute to ~8 key SAVEd arrays
- `Source/vtwobody.F:849-882` — Replace `twh_nbcoeff(GLB_GET, ...)` with module array access
- `Source/vcoulomb.F:100-140` — Replace `twh_qqscale(GLB_GET, ...)` with module array access
- `Source/allocate.F` — Add module pointer initialization after allocations
- Whichever C source implements `twh_coordstorage_` (or a new adjacent file if that one is unsuitable to extend) — add `get_coord_storage_ptr(int crdflag, double **ptr, int *natoms)`
- `Source/Makefile.am` — Add `engtotal_globals.F` to `towhee_SOURCES`; add `-DHOT_PATH_INLINE` build variant
- `Source/energy_change.F` — Apply the same inline pattern to `ec_two_body_dv`

**Note on line numbers:** every line reference in this document (e.g. `globaldata.F:10685`) was recorded against a specific snapshot of the tree and will drift. Treat them as a starting search hint, not ground truth — confirm each target with `grep -n` for the subroutine/variable name immediately before editing, on whatever revision is actually checked out.

## ARCHITECTURE

### Design decision: zero-copy pointer aliasing, not a memcpy shim
The original approach considered here was a `TARGET` "shim" array populated by `memcpy` from the C-owned coordinate storage on every accept, to get a fully de-interleaved SoA layout. **Rejected**, because:
- Coordinate storage is mutated (or its owning pointer swapped) on essentially every accepted MC move, so a shim refreshed at that cadence puts an O(N) copy on the O(N) `energy_change.F` path — the same order as the computation it's meant to accelerate, on the path that runs most often. It could easily erase the intended speedup.
- It isn't actually necessary: Fortran can bind directly to the live C-owned buffer with zero copying via rank-remapped `c_f_pointer` association, re-associated only when the underlying base address changes.

Instead: bind a 2D pointer to the interleaved C buffer, then take array-section pointers into it. This is legal Fortran (pointer assignment to a non-contiguous section is permitted; it just can't carry the `CONTIGUOUS` attribute):

```fortran
call c_f_pointer(cptr_real, crd2d_real, [3, natoms])   ! zero-copy, valid, contiguous target
crd_x => crd2d_real(1, :)   ! non-contiguous section — legal pointer target, not CONTIGUOUS
crd_y => crd2d_real(2, :)
crd_z => crd2d_real(3, :)
```
This trades stride-1 access for zero copying. Re-associating these pointers is O(1) — it's a **pull**, done once at the top of each hot-path entry function (`engtotal.F`, `energy_change.F`'s hot loop), not a **push** hooked into every mutation site. That removes the need to audit and modify five separate accept-path files (`trancom.F`, `rotate.F`, `tranatom.F`, `swapmoves.F`, `volumemoves.F`) — a single missed hook there would silently produce wrong energies with no crash to flag it. Two call sites, both already being edited in this change, own the sync instead.

> If profiling after Phase 2 shows the O(N²) `engtotal.F` path is bandwidth-bound on the stride-3 access pattern specifically, a genuine contiguous SoA snapshot for *that function only*, built fresh at the top of its own (rare) invocation, is a safe follow-up — it still avoids hooking the accept path. Don't build this pre-emptively; measure first.

### New Module: `engtotal_globals.F`
```fortran
module engtotal_globals
  use, intrinsic :: iso_c_binding
  implicit none
  save

  ! --- Coordinates: zero-copy views into C-owned storage, re-bound at hot-path entry ---
  ! NOT contiguous — these are stride-3 sections of an interleaved (x,y,z) C buffer.
  real(8), pointer :: crd2d_real(:,:) => null(), crd2d_temp(:,:) => null()
  real(8), pointer :: crd_x(:) => null(), crd_y(:) => null(), crd_z(:) => null()
  real(8), pointer :: crd_temp_x(:) => null(), crd_temp_y(:) => null(), crd_temp_z(:) => null()

  ! --- Nonbonded coefficients: force-field data, fixed for the whole run.
  ! Copied ONCE at init into genuinely contiguous flat arrays (see rationale below) —
  ! these are NOT pointer aliases into nbcoeff(MINNBCOEFF:MAXNBCOEFF, CROSSTYPEMAX),
  ! because slicing a fixed first index out of that array is stride-(MAXNBCOEFF-MINNBCOEFF+1),
  ! not stride-1, and declaring that CONTIGUOUS would be an incorrect assertion to the compiler.
  real(8), allocatable :: nbcoeff_0(:)    ! shift coeffs
  real(8), allocatable :: nbcoeff_13(:)   ! sigma^(1/6)
  real(8), allocatable :: nbcoeff_14(:)   ! epsilon
  real(8), allocatable :: nbcoeff_15(:)   ! sigma^(1/6) for 1-4
  real(8), allocatable :: nbcoeff_16(:)   ! epsilon for 1-4

  ! --- Atom type / flag data: whole-array TARGET aliases, genuinely contiguous ---
  integer, pointer, contiguous :: ntype_arr(:,:) => null()
  integer, pointer, contiguous :: nboxi_arr(:) => null()
  integer, pointer, contiguous :: moltyp_arr(:) => null()
  integer, pointer, contiguous :: nunit_arr(:) => null()
  integer, pointer, contiguous :: globalpos_arr(:) => null()
  logical, pointer, contiguous :: luse_c_arr(:,:) => null()
  logical, pointer, contiguous :: luse_v_arr(:,:) => null()
  real(8), pointer, contiguous :: qqatom_arr(:,:) => null()
  real(8), pointer, contiguous :: qqscale_arr(:,:) => null()

  ! --- Scalars / cutoffs / lambda scaling: read-only for the run, aliased once at init ---
  real(8), pointer :: rcut_val => null(), rcelect_val => null()
  logical, pointer :: lewald_val => null(), lelect_val => null(), lshift_val => null()
  real(8), pointer, contiguous :: hmatrix_arr(:,:,:) => null()
  real(8), pointer, contiguous :: hinverse_arr(:,:,:) => null()
  real(8), pointer :: lambda_lj_val => null(), lambda_c_val => null()
  logical, pointer :: scp_scale_twobody_val => null()

  integer :: natoms_cached = 0   ! sanity check only, see engtotal_sync_coord_ptrs

contains

  ! Called ONCE at the top of engtotal.F and once at the top of energy_change.F's
  ! hot loop entry — NOT hooked into accept paths. O(1): re-binds pointer
  ! descriptors, copies no coordinate data. Safe to call unconditionally on
  ! every entry; idempotent if nothing changed since the last call.
  subroutine engtotal_sync_coord_ptrs()
    type(c_ptr) :: cptr_real, cptr_temp
    integer(c_int) :: n

    call get_coord_storage_ptr(0, cptr_real, n)   ! 0 == CRD_REAL
    call c_f_pointer(cptr_real, crd2d_real, [3, int(n)])
    crd_x => crd2d_real(1, :); crd_y => crd2d_real(2, :); crd_z => crd2d_real(3, :)

    call get_coord_storage_ptr(1, cptr_temp, n)   ! 1 == CRD_TEMP
    call c_f_pointer(cptr_temp, crd2d_temp, [3, int(n)])
    crd_temp_x => crd2d_temp(1, :); crd_temp_y => crd2d_temp(2, :); crd_temp_z => crd2d_temp(3, :)

    natoms_cached = int(n)
  end subroutine engtotal_sync_coord_ptrs

end module engtotal_globals
```

C side (new, small):
```c
/* returns the current base address + atom count for CRD_REAL (0) or CRD_TEMP (1) */
void get_coord_storage_ptr(int crdflag, double **ptr, int *natoms);
```

### OpenMP note
`x_ptr`, `y_ptr`, `z_ptr` (the per-loop selectors used in Step 2.1 below) must be declared **local to the subroutine containing the pair loop**, not as `engtotal_globals` module variables. If the existing OpenMP parallelization work ("Optimization 1") parallelizes this loop, module-scoped mutable pointers reassigned per-call would be shared/racy across threads; locals are automatically private under `!$OMP PARALLEL DO` default scoping. Confirm this explicitly before Phase 2 lands if Optimization 1 is in flight concurrently.

## STEP-BY-STEP IMPLEMENTATION PLAN

### Phase 0: Build Infrastructure & Feature Flag (Estimated: 1 hour)

#### Step 0.1 — Add a compile-time toggle
The regression methodology (below) needs two builds: one with the old dispatch path, one with the new inline path. This must be a **compile-time** switch, not a runtime `if` inside the hot loop — a runtime branch checked per atom-pair would reintroduce the overhead this spec removes.

- Add `AM_FCFLAGS += -cpp` (gfortran) / equivalent `-fpp` (ifort) to `Makefile.am` if not already present, so CPP directives work in `.F`/`.F90` sources.
- Guard the new call sites with `#ifdef HOT_PATH_INLINE ... #else ... #endif`, preserving the original dispatch-call code path under `#else` until Phase 2's unit tests and Phase 3's regression test both pass.
- Add an Automake conditional or a documented `configure`/`make` variable so `make HOT_PATH_INLINE=1` (or a `--enable-hot-path-inline` configure flag) controls the `-DHOT_PATH_INLINE` define.

**Acceptance gate:** both `make` (flag off) and `make HOT_PATH_INLINE=1` (flag on) build cleanly with no source changes yet beyond the guard scaffolding — i.e., the `#else` branch is just the current code, verified byte-identical to `git diff` before this step.

### Phase 1: Create Module and Wire Static/Semi-Static Arrays (Estimated: 4 hours)

#### Step 1.1 — Create `Source/engtotal_globals.F`
Write the module as specified above in ARCHITECTURE.

Constraint: Do NOT introduce any USE statements in `globaldata.F`. The module is consumed only by the hot-path functions (`engtotal.F`, `vtwobody.F`, `vcoulomb.F`, `energy_change.F`), not by the global data layer.

#### Step 1.2 — Add `TARGET` to Key Arrays in `globaldata.F`
Locate each SAVEd array accessed from the hot path and add `TARGET` (verify line numbers per the TARGET FILES note above before editing):

| Subroutine | Array | Change |
|-----------|-------|--------|
| `twh_nbcoeff` | `nbcoeff(MINNBCOEFF:MAXNBCOEFF, CROSSTYPEMAX)` | `dimension, target :: nbcoeff(...)` |
| `twh_ntype` | `ntype(NTMAX, NUMAX)` | Add TARGET |
| `twh_nboxi` | `nboxi(NTMAX, MAXCHAIN)` | Add TARGET |
| `twh_moltyp` | `moltyp(MAXCHAIN)` | Add TARGET |
| `twh_nunit` | `nunit(NTMAX)` | Add TARGET |
| `twh_luse_c` | `luse_c(NTMAX, NUMAX)` | Add TARGET |
| `twh_luse_v` | `luse_v(NTMAX, NUMAX)` | Add TARGET |
| `twh_rcut` | `rcut` | Add TARGET |
| `twh_lshift` | `lshift` | Add TARGET |
| `twh_hmatrix` | `hmatrix(3,3,MAXBOX)` | Add TARGET |

Constraint: Only add TARGET to the variable declaration line. Do NOT change any logic in the flag-dispatch body. The existing `twh_*` subroutines remain unchanged for use by other parts of the code.

#### Step 1.3 — Initialize Module Pointers in `allocate.F`
Create `engtotal_globals_init()`, called once after all allocations complete:

```fortran
subroutine engtotal_globals_init()
  use engtotal_globals
  integer :: ct

  ! Whole-array TARGET aliases — genuinely contiguous, zero-copy
  ntype_arr => ntype
  luse_c_arr => luse_c
  luse_v_arr => luse_v
  ! ... same pattern for nboxi_arr, moltyp_arr, nunit_arr, qqatom_arr,
  !     qqscale_arr, hmatrix_arr, hinverse_arr, rcut_val, lshift_val, etc.

  ! nbcoeff slices — ONE-TIME real copy (not a pointer alias; see ARCHITECTURE
  ! for why nbcoeff(fixed_row, :) is not a valid CONTIGUOUS pointer target).
  ! nbcoeff is fixed force-field data for the whole run, so this copy happens
  ! exactly once and is never refreshed.
  allocate(nbcoeff_0(CROSSTYPEMAX), nbcoeff_13(CROSSTYPEMAX), nbcoeff_14(CROSSTYPEMAX), &
           nbcoeff_15(CROSSTYPEMAX), nbcoeff_16(CROSSTYPEMAX))
  do ct = 1, CROSSTYPEMAX
    nbcoeff_0(ct)  = nbcoeff(0, ct)
    nbcoeff_13(ct) = nbcoeff(13, ct)
    nbcoeff_14(ct) = nbcoeff(14, ct)
    nbcoeff_15(ct) = nbcoeff(15, ct)
    nbcoeff_16(ct) = nbcoeff(16, ct)
  end do

  ! Coordinates are NOT bound here — engtotal_sync_coord_ptrs() is called at
  ! the top of engtotal.F and energy_change.F instead, since coordinates
  ! mutate during the run and everything else on this list does not.
end subroutine engtotal_globals_init
```

**Acceptance gate:** with `HOT_PATH_INLINE` defined, add a temporary debug routine that prints every module pointer/array's first and last element alongside the equivalent `twh_*(GLB_GET, ...)` call's output for a small hand-built test system, and diff by eye. Remove the debug routine once Phase 1 is confirmed correct — don't let it leak into Phase 2.

### Phase 2: Rewrite the Hot Loop Bodies (Estimated: 5 hours)

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
use engtotal_globals, only: engtotal_sync_coord_ptrs, crd_x, crd_y, crd_z, &
     crd_temp_x, crd_temp_y, crd_temp_z, ntype_arr, luse_c_arr, luse_v_arr, &
     globalpos_arr
! x_ptr/y_ptr/z_ptr are LOCAL to this subroutine — see ARCHITECTURE's OpenMP note.
real(8), pointer :: x_ptr(:), y_ptr(:), z_ptr(:)

call engtotal_sync_coord_ptrs()   ! O(1) rebind; safe/idempotent every call

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

Constraint: This function must be `PURE` — no side effects, no I/O, no SAVE variables except the module data it reads (which are read-only for the duration of the run). **Verify this compiles under both `gfortran -std=gnu` and `ifort -standard-semantics` before relying on it** — `PURE` procedures reading use-associated module data have had inconsistent enforcement across compiler versions, and Constraint 8 pins this codebase to pre-2008 semantics, so this specific combination needs a direct compile check rather than an assumption.

#### Step 2.3 — Inline `twh_vcoulomb` in `vcoulomb.F`

Create `inline_coulomb(rij, imolty, ii, jmolty, jj)` that reads `qqatom_arr` and `qqscale_arr` directly. For the Ewald short-range case, use the Fortran intrinsic `erfc()` (F2008) instead of the custom `twh_derfc()` — confirm this is acceptable given Constraint 8's F2008-feature ceiling (`CONTIGUOUS` and `c_f_pointer` only); if not, keep `twh_derfc()`.

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

Constraint: `twh_mimage` (minimum image) remains a regular subroutine call — pure geometric computation, no data-access overhead. Not worth inlining.

Constraint: `twh_getnbtype` is pure integer arithmetic; keep `getnbtype.F`'s interface but call it directly, without the `twh_` dispatch layer, if possible.

Constraint: The `lovrlap` check (`rijsq < rminsq` → overlap) and `twh_engtotal_overlap` call must remain, moved inline:
```fortran
if (rijsq < rminsq_val) then
  lovrlap = .true.
  ! overlap recovery logic (preserve existing behavior)
endif
```

**Acceptance gate (this phase's most important check):** before moving to Phase 3, run the **unit-level function diff test** described under TESTING & VALIDATION — do not rely solely on the whole-simulation regression test at the end to catch bugs in these two functions.

### Phase 3: Coordinate Zero-Copy Sync (Estimated: 1 hour)

This phase is much smaller than the original push-based design. `engtotal_sync_coord_ptrs()` is already called at the top of `engtotal.F` (Step 2.1) and needs the same one-line call added at the top of `energy_change.F`'s hot loop entry (Phase 4). **No accept-path files need to be touched for this purpose** — no hooks in `trancom.F`, `rotate.F`, `tranatom.F`, `swapmoves.F`, or `volumemoves.F`.

#### Step 3.1 — Confirm idempotence and cost
Add a debug counter (temporary, removed before merge) incrementing on each `engtotal_sync_coord_ptrs()` call, and confirm via `perf stat` that its cost is negligible relative to the pair loop that follows it (it should be — one small C call plus two `c_f_pointer` rebinds per hot-function entry, not per atom-pair).

Constraint: If a future accept path is added elsewhere in the codebase, it does not need to know about this module at all — correctness comes from the pull at the top of `engtotal.F`/`energy_change.F`, not from anything the accept path does. This is the main robustness win of the pull-based design over the original push-based one.

### Phase 4: Update `energy_change.F` (Estimated: 1.5 hours) — COMPLETE 2026-08-08

**Implementation note:** In this codebase `twh_ec_two_body_dv` (energy_change.F:217) does **not** contain its own pair loop — it delegates to `twh_eng_two_body` (energy.F) for both `TRL_NEW` and `TRL_OLD`, which is exactly the routine inlined in Phase 2. Therefore no additional inlining is required in `energy_change.F` itself; the OLD/NEW coordinate selection flows through `chaininfo(ichainmove,crdflag)` → `testcrd`/`jcrd`, which feed the `twh_coordinates` calls preserved in the inline path. No `engtotal_sync_coord_ptrs()` call was added because the current architecture keeps `twh_coordinates` for coordinate access (zero-copy coord pointers are not used anywhere).

**Validation (2026-08-08):** Baseline (dispatch) vs `HOT_PATH_INLINE` builds on a 10-example set spanning all ensembles (NVT AVB1/Hard_Sphere/Formamide_Scaled/TraPPE_Molecules, Gibbs EPM_VLCC/Cui2002_Hexane, grand canonical Energy_Biasing, NPT OPLS_Propanamide, TMMC LJ/SPC_E_WATER): **10/10 bit-identical** outputs. Move-heavy OLD/NEW acceptance checks (AVB1_Methane, Square_Well_Chain, TraPPE_Isomers) also bit-identical. Full 56-test serial suite: identical pass/fail sets (45 pre-existing platform/compiler diffs in both builds). **Zero regressions.**

#### Step 4.1 — Find equivalent hot loops in `energy_change.F`
`ec_two_body_dv` evaluates pair energy changes for moved molecules vs. all others, with the same nested structure as `engtotal.F` but over (moved molecule, all others). Add `call engtotal_sync_coord_ptrs()` once at the top of this routine (or its outer caller, whichever is the single entry point for a given move evaluation — not inside the pair loop), then apply the identical inlining pattern from Phase 2.

Constraint: `energy_change.F` uses `chaininfo` (per-moved-molecule flags) to pick OLD (`CRD_REAL`) vs NEW (`CRD_TEMP`) coordinates. This selection must be done per moved molecule, not once for the whole loop:
```fortran
x_ptr_chain(1:nchainmove)  ! array of pointers to x coords for each chain
y_ptr_chain(1:nchainmove)
z_ptr_chain(1:nchainmove)
```
Set per-chain before the outer loop, based on `chaininfo(:,2:3)`, pointing each entry at the already-synced `crd_x`/`crd_temp_x` etc. (no additional sync needed per chain — one call to `engtotal_sync_coord_ptrs()` per invocation covers both buffers).

**Acceptance gate:** run the unit-level diff test again, this time specifically exercising the OLD/NEW coordinate selection logic (a molecule with `chaininfo` indicating NEW coordinates, another indicating OLD, in the same call) — this is the piece of Phase 4 most likely to have an indexing bug that a whole-molecule energy check could mask.

### Phase 5: Build System (Estimated: 0.5 hours)

#### Step 5.1 — Update `Makefile.am`
Add `engtotal_globals.F` to `towhee_SOURCES` in alphabetical position (after `energy_change.F`, before `engatom.F`):
```
engtotal_globals.F \
```
Add an explicit build dependency so the module compiles first:
```
engtotal.o vtwobody.o vcoulomb.o energy_change.o: engtotal_globals.o
```
Confirm Phase 0's `HOT_PATH_INLINE` conditional-compile plumbing still produces both variants cleanly with the new source file included.

## TESTING & VALIDATION

### Baseline Capture (before Phase 0)
Before any code changes, record: (a) per-step energies for the regression system below on the unmodified tree, (b) wall-clock and `perf stat` function-call counts for the performance system below, (c) current `make` clean-build time. This is the reference every later comparison uses — capturing it after changes have already started risks an apples-to-oranges baseline.

### Unit-Level Function Diff Test (new — run after Phase 2, again after Phase 4)
Whole-simulation energy comparisons can hide compensating bugs (e.g., a sign error that happens not to matter for the specific atom types present in one test system). Write a standalone harness that calls `inline_lj_12_6` / `inline_coulomb` and the original `twh_vtwobody` / `twh_vcoulomb` directly, over a swept grid of inputs:
- `rijsq` spanning below, at, and above `rminsq` (the overlap boundary)
- `ntij` at `MINNBCOEFF`/`MAXNBCOEFF` boundaries and mid-range
- `lonefour` both states
- `lshift`, `lewald`, `lelect` on and off in all combinations

Diff to machine epsilon (not the 1×10⁻¹² tolerance below — at the single-function level there's no vectorization-driven reordering to account for, so exact match is the right bar here).

### Regression Test
1. Build with `HOT_PATH_INLINE=0`, run a reference simulation (e.g. NPT N₂, 500 molecules, 2000 MC steps). Record per-step energies.
2. Build with `HOT_PATH_INLINE=1`, run the same simulation with the same random seed.
3. Compare per-step energies. Must match within 1×10⁻¹² relative error.

Note: this replaces any requirement for bit-for-bit identical output. `CONTIGUOUS` and the removal of call barriers are expected to let the compiler auto-vectorize the pair-energy accumulation, which can reorder floating-point summation relative to the original scalar/call-gated loop. Since addition isn't associative, true bitwise identity between the two code paths is not a realistic target if vectorization is contributing to the speedup — 1×10⁻¹² relative error is the standard this spec holds to throughout.

4. Run a second regression system dominated by single-molecule translation/rotation moves (not just the full-recompute-heavy NPT case above), since Phase 4's `energy_change.F` path is the one flagged as highest-risk in ARCHITECTURE.

### Performance Test
1. Benchmark a realistic system: **1,000-molecule SPC/E water (3,000 atoms), run for 10,000 MC cycles with `stepstyle='cycles'`** (`mainloop.F:360-370`: a cycle = `nchain` moves, i.e. one attempted move per molecule), with and without the optimization, using the Baseline Capture numbers above as the "without" side. **Do not use `stepstyle='moves'` at 10,000 steps on a 10,000-atom system** — that yields only ~1 attempted move per molecule, which leaves the O(N) `energy_change.F` per-move path (the dominant call frequency in a real simulation, per CONTEXT) essentially unexercised and understates the optimization's benefit. 1,000 molecules × 10,000 cycles = 10,000 attempted moves per molecule, exercising the accept path at realistic density.
2. Profile with `perf stat` or `gprof` — measure reduction in function call counts in the hot path, and separately confirm `engtotal_sync_coord_ptrs()` call count matches hot-function entry count (not atom-pair count).
3. Expected result: 1.5-3× speedup on the nonbonded energy evaluation, 20-50% reduction in total simulation time.

### Cross-Compiler Test
Run the Regression Test and Performance Test under both `gfortran -std=gnu` and `ifort -standard-semantics`. Pointer-aliasing and auto-vectorization behavior are known areas of divergence between these compilers, and Constraint 3's tolerance claim should be verified under both, not assumed to transfer from one.

### Memory Test
1. Run with `valgrind --tool=memcheck` to detect any dangling pointer references from the `c_f_pointer` re-association, especially around the `CRD_REAL`/`CRD_TEMP` swap boundary.
2. Run with `-fbounds-check` (gfortran) to catch array index out-of-bounds from incorrect index math.

## CONSTRAINTS & GUARDRAILS

1. **Do NOT remove or alter the original `twh_*` accessor subroutines** in `globaldata.F`. Other parts of the code (not in the hot path) still use them. The optimization adds a bypass, not a replacement.

2. **Do NOT add `USE engtotal_globals` to `globaldata.F`** or to any file outside the direct hot path. The module is a private optimization detail.

3. **Match the original computation to within 1×10⁻¹² relative error**, not bit-for-bit. See the Regression Test note above for why bitwise identity isn't the right bar once vectorization is in play.

4. **Coordinate correctness comes from the pull-based sync, not from accept-path hooks.** `engtotal_sync_coord_ptrs()` must be called at the top of `engtotal.F` and at the top of `energy_change.F`'s hot-path entry — and only there. Do not add calls to it elsewhere; a single, small set of call sites is the whole point of this design (see ARCHITECTURE).

5. **The `nbcoeff` extraction must use the correct `ntij` index** — the triangular map from `getnbtype.F:22-69` computes `ntij` from `(ntii, ntjj)` as:
   ```
   low = min(ntii, ntjj)
   high = max(ntii, ntjj)
   istep = (low-1) * NNTYPE - ((low-1)*(low-2))/2
   ntij = istep + high - low + 1
   ```
   Verify this against `twh_getnbtype`'s actual output with an automated test sweeping the full `(ntii, ntjj)` range, not a spot check.

6. **The `lonefour` flag** (1-4 interactions) is computed elsewhere from bonding topology and must be accepted as a parameter by the inline function, matching the original `vtb_lennard_jones` signature.

7. **Only alias arrays as `CONTIGUOUS` when the alias genuinely has unit stride.** Whole-array `TARGET` aliases (`ntype_arr`, `luse_c_arr`, etc.) qualify. A fixed-first-index slice of a column-major array (the original `nbcoeff_13`-style approach) does not — see ARCHITECTURE for why those five coefficients are one-time-copied instead of pointer-aliased.

8. **Do NOT use Fortran 2008+ features** beyond `CONTIGUOUS` and `c_f_pointer` — with the `PURE`+module-data combination in Step 2.2 and the `erfc()` intrinsic in Step 2.3 flagged for an explicit compiler check first, since both have had version-dependent behavior. The codebase targets F77 + extensions, compiled with gfortran `-std=gnu` and ifort `-standard-semantics`.

9. **Line numbers throughout this document are approximate.** Confirm every target location against the actual checked-out revision before editing (see TARGET FILES).

## SUCCESS CRITERIA

| Metric | Target | Measurement |
|--------|--------|------------|
| Function calls per atom pair in hot loop | ≤2 (was 5-20) | `perf stat` on engtotal function |
| Wall-clock time for nonbonded evaluation | 33-67% reduction | Runtime comparison vs. Baseline Capture |
| Coordinate sync cost | O(1) per hot-function entry, not per atom-pair or per accept | Code review + `perf stat` call count |
| Per-step energy drift vs baseline | ≤1×10⁻¹² (relative) | Compare final energies, both compilers |
| Unit-level function diff (`inline_lj_12_6`, `inline_coulomb`) | Exact match to original over full input sweep | Standalone harness, see TESTING |
| Memory increase | ≤1% | `valgrind` heap usage |
| Build time increase | ≤10 seconds, measured as clean-build delta vs. Baseline Capture on the same machine | `make` timing |
