## OBJECTIVE
Implement OpenMP multi-threading over the primary molecule-pair loop (`do ichain = 1, nchain - 1`) in `twh_engtotal_twobody` (`engtotal.F:334-458`) to parallelize the O(N²) intermolecular nonbonded energy evaluation across CPU cores. Every variable updated inside the parallel region must be either included in a `REDUCTION` clause or declared `PRIVATE`. No `GLB_SET` operations are allowed within the parallel region.

## TARGET FILES
- `Source/engtotal.F:334-458` — The `ichain` outer loop body (`twh_engtotal_twobody`)
- `Source/globaldata.F` — Audit all `twh_*` subroutines accessed from the loop to confirm safe read-only access
- `Source/Makefile.am` — Add OpenMP compiler flag (`-fopenmp` / `-qopenmp`)

## HOT LOOP CONTEXT (engtotal.F:334-458)

The target loop is in subroutine `twh_engtotal_twobody`:

```fortran
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1, nchain - 1
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,numiunit)
            if ( lmolecmove ) then
               call twh_com(GLB_GET,crdflag,ichain,xcmi,ycmi,zcmi)
               call twh_rcmu(GLB_GET,ichain,rcmi)
            endif
            do jchain = ichain + 1, nchain
               call twh_nboxi(GLB_GET,jchain,testboxj)
               if ( testboxj .eq. ibox ) then
                  call twh_moltyp(GLB_GET,jchain,jmolty)
                  call twh_nunit(GLB_GET,jmolty,numjunit)
                  lcompute = .true.
                  if ( lmovelist ) then
                     call twh_lmoveme(GLB_GET,ichain,lmovemei)
                     call twh_lmoveme(GLB_GET,jchain,lmovemej)
                     if ( .not. (lmovemei .or. lmovemej) ) then
                        lcompute = .false.
                     endif
                  endif
                  if ( lmolecmove .and. lusecom ) then
                     call twh_com(GLB_GET,crdflag,jchain,xcmj,ycmj,zcmj)
                     call twh_mimage(bxdflag,ibox,xcmi,ycmi,zcmi,xcmj,ycmj,zcmj,rxuij,ryuij,rzuij)
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     call twh_rcmu(GLB_GET,jchain,rcmj)
                     call twh_rcut(GLB_GET,rcm)
                     if ( lelect ) then
                        if ( rcelect .gt. rcm ) rcm = rcelect
                     endif
                     rcm = rcm + rcmi + rcmj
                     rcmsq = rcm*rcm
                     if ( rijsq .gt. rcmsq ) lcompute = .false.
                  endif
                  if ( lcompute ) then
                     do ii = 1,numiunit
                        call twh_ntype(GLB_GET,imolty,ii,ntii)
                        call twh_coordinates(GLB_GET,crdflag,ichain,ii,rxui,ryui,rzui)
                        call twh_luse_c(GLB_GET,imolty,ii,luse_ci)
                        call twh_luse_v(GLB_GET,imolty,ii,luse_vi)
                        do jj = 1, numjunit
                           call twh_ntype(GLB_GET,jmolty,jj,ntjj)
                           call twh_coordinates(GLB_GET,crdflag,jchain,jj,rxuj,ryuj,rzuj)
                           call twh_luse_c(GLB_GET,jmolty,jj,luse_cj)
                           call twh_luse_v(GLB_GET,jmolty,jj,luse_vj)
                           leval_v = luse_vi .and. luse_vj
                           leval_c = luse_ci .and. luse_cj
                           if ( leval_v .or. leval_c ) then
                              call twh_mimage(bxdflag,ibox,rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
                              rijsq = rxuj*rxuj+ryuj*ryuj+rzuj*rzuj
                              if (leval_v) then
                                 vtwointer = vtwointer + twh_vtwobody(...)
                                 dl_twointer = dl_twointer + dvdl
                                 do l = 1, scp_num_foreign_lambda
                                    vftwo(l)=vftwo(l) + vf(l)
                                 enddo
                                 if ( lovrlap ) call twh_engtotal_overlap(...)
                              endif
                              if ( lelect .and. leval_c) then
                                 rij = dsqrt(rijsq)
                                 call twh_vcoulomb(...)
                                 velectinter = velectinter + vcoulomb
                                 if ( lovrlap ) call twh_engtotal_overlap(...)
                                 dl_elect = dl_elect + dvdl
                                 do l = 1, scp_num_foreign_lambda
                                    vfelec(l) = vfelec(l) + vf(l)
                                 enddo
                                 if ( lmadelung ) then
                                    if ( rij .lt. maddist ) maddist = rij
                                 endif
                              endif
                           endif
                        enddo
                     enddo
                  endif
               endif
            enddo
         endif
      enddo
```

## VARIABLE AUDIT & THREAD-SAFETY CLASSIFICATION

### Scalars updated inside the loop body (must be REDUCTION or PRIVATE)

| Variable | Nature | Classification | Rationale |
|----------|--------|---------------|-----------|
| `vtwointer` | Energy accumulator | `REDUCTION(+:vtwointer)` | Summed with `+` across pairs |
| `velectinter` | Energy accumulator | `REDUCTION(+:velectinter)` | Summed with `+` across pairs |
| `dl_twointer` | Derivative accumulator | `REDUCTION(+:dl_twointer)` | Summed with `+` across pairs |
| `dl_elect` | Derivative accumulator | `REDUCTION(+:dl_elect)` | Summed with `+` across pairs |
| `maddist` | `min` reduction (Madelung) | `REDUCTION(min:maddist)` | Updated via `min()`; requires `.false.` initial |
| `vftwo(1:scp_num_foreign_lambda)` | Array accumulator | `REDUCTION(+:vftwo)` | Summed element-wise |
| `vfelec(1:scp_num_foreign_lambda)` | Array accumulator | `REDUCTION(+:vfelec)` | Summed element-wise |
| `ichain` | Outer loop variable | Automatically PRIVATE | OMP loop index |
| `jchain` | Inner loop variable | PRIVATE | Varies per thread |
| `testboxi` | Per-ichain cache | PRIVATE | Local to iteration |
| `testboxj` | Per-jchain cache | PRIVATE | Local to iteration |
| `imolty` | Per-ichain cache | PRIVATE | Varies per iteration |
| `jmolty` | Per-jchain cache | PRIVATE | Varies per iteration |
| `numiunit` | Per-ichain cache | PRIVATE | Varies per iteration |
| `numjunit` | Per-jchain cache | PRIVATE | Varies per iteration |
| `xcmi, ycmi, zcmi` | COM i | PRIVATE | Varies per iteration |
| `rcmi` | COM radius i | PRIVATE | Varies per iteration |
| `xcmj, ycmj, zcmj` | COM j | PRIVATE | Varies per iteration |
| `rcmj` | COM radius j | PRIVATE | Varies per iteration |
| `rcm` | Cutoff variable | PRIVATE | Local to iteration |
| `rcmsq` | Cutoff² variable | PRIVATE | Local to iteration |
| `rijsq` | Pair distance² | PRIVATE | Local to pair |
| `rij` | Pair distance | PRIVATE | Local to pair |
| `rxuij, ryuij, rzuij` | Min-image vector | PRIVATE | Local to pair |
| `rxui, ryui, rzui` | Coords i | PRIVATE | Local to (i,ii) |
| `rxuj, ryuj, rzuj` | Coords j | PRIVATE | Local to (j,jj) |
| `ii` | Atom loop i | PRIVATE | Nested loop |
| `jj` | Atom loop j | PRIVATE | Nested loop |
| `ntii` | Atom type i | PRIVATE | Local to (i,ii) |
| `ntjj` | Atom type j | PRIVATE | Local to (j,jj) |
| `luse_ci, luse_vi` | Flag i | PRIVATE | Local to (i,ii) |
| `luse_cj, luse_vj` | Flag j | PRIVATE | Local to (j,jj) |
| `leval_v` | Boolean | PRIVATE | Local to pair |
| `leval_c` | Boolean | PRIVATE | Local to pair |
| `vcoulomb` | Coulomb energy | PRIVATE | Local to pair |
| `dvdl` | Derivative | PRIVATE | Returned by vtwobody/vcoulomb |
| `vf(1:MAX_FOREIGN_LAMBDA)` | Temp array | PRIVATE | Passed to vtwobody/vcoulomb |
| `lcompute` | Boolean | PRIVATE | Local to (i,j) |
| `lmovemei, lmovemej` | Boolean | PRIVATE | Local to (i,j) |
| `l` | Loop index | PRIVATE | Implicit in DO, but declare for safety |

### Scalars and arrays READ-ONLY (shared, no conflict)

| Variable | Rationale |
|----------|-----------|
| `nchain` | Set once before the parallel region, read-only during loop |
| `ibox` | Function argument, constant |
| `crdflag` | Function argument, constant |
| `bxdflag` | Function argument, constant |
| `lmolecmove` | Function argument, constant |
| `lmovelist` | Function argument, constant |
| `cutstyle` | Local scalar set before loop |
| `lelect` | Set before loop via `twh_lelect(GLB_GET,lelect)` |
| `lusecom` | Set before loop via `twh_lusecom(GLB_GET,lusecom)` |
| `lovrlap` | Logical read, only set to `.true.` on overlap within iteration |
| `lfailure` | Logical read, only set to `.true.` on error within iteration |
| `rcelect` | Read via `twh_rcelect(GLB_GET,bxdflag,ibox,rcelect)` |
| `scp_num_foreign_lambda` | Read via `twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)` |
| `lmadelung` | Set before loop, read-only during |
| `rcmi` (zero-init) | Set before loop to 0.0d0 |

### GLB_GET calls (all read-only, safe for shared access)

Each `twh_*(GLB_GET, ...)` call in the loop body reads from a SAVEd array or scalar in `globaldata.F` without writing. Since OpenMP with `PARALLEL DO` and `DEFAULT(NONE)` with explicit PRIVATE/REDUCTION/SHARED ensures no unintended writes, all GLB_GET access is thread-safe.

Specifically verified safe:
- `twh_nboxi(GLB_GET, ...)` — reads `nboxi()` array, read-only
- `twh_moltyp(GLB_GET, ...)` — reads `moltyp()` array, read-only
- `twh_nunit(GLB_GET, ...)` — reads `nunit()` array, read-only
- `twh_ntype(GLB_GET, ...)` — reads `ntype()` array, read-only
- `twh_coordinates(GLB_GET, ...)` — reads coordinate storage, read-only
- `twh_luse_c/luse_v(GLB_GET, ...)` — reads `luse_c/luse_v` arrays, read-only
- `twh_com(GLB_GET, ...)` — reads COM storage, read-only
- `twh_rcmu(GLB_GET, ...)` — reads per-molecule COM radius, read-only
- `twh_lmoveme(GLB_GET, ...)` — reads movelist flags, read-only
- `twh_rcut(GLB_GET, ...)` — reads scalar cutoff, read-only
- `twh_mimage(...)` — pure geometric subroutine, no global state

### GLB_SET audit — NONE in the loop body

Confirmed: There are **zero** `GLB_SET` operations inside `do ichain = 1, nchain - 1` in `twh_engtotal_twobody`. The subroutines called during the loop (`twh_vtwobody`, `twh_vcoulomb`, `twh_mimage`, `twh_engtotal_overlap`) also perform only GLB_GET or local side-effects (`lfailure`, `lovrlap`).

**Important exception:** `twh_engtotal_overlap` inside the loop may set `lfailure = .true.` — this is thread-safe as a write to PRIVATE or per-iteration state, but `lfailure` is passed into the subroutine as an argument and is thread-private by virtue of being local. No global data is modified.

## STEP-BY-STEP IMPLEMENTATION PLAN

### Phase 1: Add OpenMP Compiler Flag to Build System (Estimated: 0.5 hours)

#### Step 1.1 — Modify `configure.ac`
Add a new Automake conditional `UseOpenMP` that accepts `--enable-openmp`:

```
AC_ARG_ENABLE([openmp],
  [AS_HELP_STRING([--enable-openmp], [enable OpenMP parallelization])],
  [],
  [enable_openmp=no])
AM_CONDITIONAL([UseOpenMP], [test "x$enable_openmp" = xyes])
AS_IF([test "x$enable_openmp" = xyes], [
  AC_DEFINE([USEOPENMP], [1], [Define if OpenMP is enabled])
])
```

#### Step 1.2 — Modify `Makefile.am`
Add OpenMP compiler flags and add `engtotal.o` dependency:

```
if UseOpenMP
OMP_FFLAGS = -fopenmp
OMP_CPPFLAGS = -DUSEOPENMP
endif

AM_FFLAGS += $(OMP_FFLAGS)
AM_CPPFLAGS += $(OMP_CPPFLAGS)
```

**Note:** `-fopenmp` is for gfortran. For ifort, use `-qopenmp`. The `Makefile.am` can detect the compiler and set appropriately, or this can be a configure-time decision. Recommend a configure-time check adding `-fopenmp` to `FCFLAGS` when `--enable-openmp` is used.

### Phase 2: Insert OpenMP Directives in engtotal.F (Estimated: 2 hours)

#### Step 2.1 — Add the PARALLEL DO directive

Insert before the `do ichain` loop at line 335 in `twh_engtotal_twobody`:

```fortran
#ifdef USEOPENMP
!$OMP PARALLEL DO DEFAULT(NONE)                                         &
!$OMP& PRIVATE(ichain, jchain, testboxi, testboxj, imolty, jmolty,     &
!$OMP&         numiunit, numjunit, lcompute, lmovemei, lmovemej,        &
!$OMP&         xcmi, ycmi, zcmi, rcmi, xcmj, ycmj, zcmj, rcmj,         &
!$OMP&         rcm, rcmsq, rijsq, rij,                                  &
!$OMP&         rxuij, ryuij, rzuij,                                     &
!$OMP&         rxui, ryui, rzui, rxuj, ryuj, rzuj,                     &
!$OMP&         ii, jj, ntii, ntjj, luse_ci, luse_vi, luse_cj, luse_vj, &
!$OMP&         leval_v, leval_c, vcoulomb, dvdl, l,                     &
!$OMP&         vf,                                                       &
!$OMP&         lovrlap, lfailure)                                        &
!$OMP& SHARED(nchain, ibox, crdflag, bxdflag, lmolecmove, lmovelist,    &
!$OMP&        cutstyle, lelect, lusecom, lmadelung,                      &
!$OMP&        scp_num_foreign_lambda, rcelect,                           &
!$OMP&        vftwo, vfelec)                                             &
!$OMP& REDUCTION(+:vtwointer, velectinter, dl_twointer, dl_elect)        &
!$OMP& REDUCTION(+:vftwo, vfelec)                                       &
!$OMP& REDUCTION(min:maddist)
#endif
      do ichain = 1, nchain - 1
```

#### Step 2.2 — Add OMP END PARALLEL DO

Insert after the `enddo` at line 458 (before the `do l=1,MAX_FOREIGN_LAMBDA` summation):

```fortran
#ifdef USEOPENMP
!$OMP END PARALLEL DO
#endif
```

#### Step 2.3 — Handle early exit (lfailure / lovrlap)

**Critical constraint:** `RETURN` and `STOP` cannot be executed from inside a parallel region in standard OpenMP. The existing early-exit logic (`if ( lfailure ) return` and `if ( lovrlap ) ... return`) must be replaced with a flag-based approach:

1. Declare a shared flag variable `lomp_failure` (add to `twh_engtotal_twobody` locals):
   ```fortran
   logical lomp_failure
   ```

2. Initialize to `.false.` before the parallel region:
   ```fortran
   lomp_failure = .false.
   ```

3. Inside the parallel region, replace all `if ( lfailure ) return` with:
   ```fortran
   if ( lfailure ) then
      lomp_failure = .true.
   endif
   ```

4. Replace all `if ( lovrlap ) ... return` with:
   ```fortran
   if ( lovrlap ) then
      ! Overlap recovery already done by twh_engtotal_overlap
   endif
   ```

5. After the parallel region, check and propagate:
   ```fortran
   if ( lomp_failure ) then
      lfailure = .true.
      return
   endif
   ```

**Design note on `lovrlap`:** The `lovrlap` variable is logically OR'd across threads — if *any* thread detects an overlap, the overall result should be `.true.`. Since `lovrlap` is currently per-iteration, an OpenMP reduction `REDUCTION(.OR.:lovrlap)` is the cleanest approach:
```fortran
!$OMP& REDUCTION(.OR.:lovrlap)
```
However, Fortran OpenMP requires the `.TRUE.` initial value rule for .OR. reduction to work correctly. Since `lovrlap` is initialized to `.false.` before the parallel region, the reduction is safe.

### Phase 3: Variable Audit and globaldata.F Confirmation (Estimated: 1 hour)

#### Step 3.1 — Confirm all GLB_GET-only access

For each `twh_*` subroutine called from the loop body, verify the code path guarded by `flag .eq. GLB_GET` only reads from SAVEd arrays:

| Subroutine | globaldata.F line | Safe? | 
|-----------|-------------------|-------|
| `twh_nboxi` | ~8388 | Yes — GLB_GET reads `nboxi(array)` |
| `twh_moltyp` | ~11464 | Yes — GLB_GET reads `moltyp(array)` |
| `twh_nunit` | ~11629 | Yes — GLB_GET reads `nunit(array)` |
| `twh_ntype` | ~11568 | Yes — GLB_GET reads `ntype(array)` |
| `twh_luse_c` | ~10053 | Yes — GLB_GET reads `luse_c(array)` |
| `twh_luse_v` | ~10090 | Yes — GLB_GET reads `luse_v(array)` |
| `twh_coordinates` | ~4153 | Yes — GLB_GET reads coordinate arrays, read-only |
| `twh_com` | ~3756 | Yes — GLB_GET reads COM arrays, read-only |
| `twh_rcmu` | ~14847 | Yes — GLB_GET reads rcmu array, read-only |
| `twh_lmoveme` | ~9603 | Yes — GLB_GET reads moveme array, read-only |
| `twh_rcut` | ~14447 | Yes — GLB_GET reads scalar rcut |
| `twh_rcelect` | ~14407 | Yes — GLB_GET reads per-box scalar |
| `twh_scp_num_foreign_lambda` | ~14975 | Yes — GLB_GET reads scalar |

#### Step 3.2 — Verify called subroutines (vtwobody, vcoulomb, mimage)

- **`twh_vtwobody`** (`vtwobody.F`): Reads nonbonded coefficients via `twh_nbcoeff(GLB_GET, ...)` — read-only. No global state mutation. Returns `vf` array and `dvdl` scalar as output arguments (local to iteration).
- **`twh_vcoulomb`** (`vcoulomb.F`): Reads charge, qqscale via `twh_qqatom(GLB_GET, ...)`, `twh_qqscale(GLB_GET, ...)`. No global state mutation.
- **`twh_mimage`** (`mimage.F`): Pure geometric minimum-image computation — no global state access.
- **`twh_engtotal_overlap`** (`engtotal.F:1573`): Reads `twh_rminsq(GLB_GET)` (read-only), writes `lfailure` (intent out — thread-private). No GLB_SET.

**All clear — the entire call tree is read-only on global data.**

### Phase 4: Testing & Validation (Estimated: 2 hours)

#### Step 4.1 — Build test
```bash
./configure --enable-openmp
make -j
```
Confirm the binary links and `engtotal.o` was compiled with `-fopenmp`.

#### Step 4.2 — Correctness test (bitwise comparison with single thread)
```bash
export OMP_NUM_THREADS=1
./towhee < test_input > output_1t.txt
export OMP_NUM_THREADS=4
./towhee < test_input > output_4t.txt
```
Compare per-step energies — must match to within 1x10⁻¹² relative error (floating-point reassociation across threads can cause non-bitwise-identical results even for the same arithmetic).

#### Step 4.3 — Performance test
Run a large system (10,000+ atoms, e.g., SPC/E water) and measure wall time at OMP_NUM_THREADS=1,2,4,8. Expected scaling: near-linear up to available cores, with diminishing returns when memory-bandwidth-bound.

#### Step 4.4 — Race condition / thread-safety stress test
Run the correctness test 5 times with 4 threads — verify bitwise-identical output across all 5 runs (within the 1e-12 tolerance). Non-deterministic results indicate a missing PRIVATE clause or a data race.

#### Step 4.5 — valgrind DRD/Helgrind test
```bash
export OMP_NUM_THREADS=2
valgrind --tool=drd --check-stack-var=yes ./towhee < test_input
```
Check for any data race reports. A clean run confirms thread safety of the parallel region.

## CONSTRAINTS & GUARDRAILS

1. **No `GLB_SET` allowed inside the parallel region.** If future code changes introduce a `GLB_SET` call in the loop body, the parallelization must be re-audited.

2. **The `PRIVATE` list must be exhaustive.** Omitting a variable from `PRIVATE` that is written-to per-iteration causes undefined behavior (data race). Use `DEFAULT(NONE)` to force the compiler to error on any unlisted variable.

3. **`lfailure` and `lovrlap` carry special meaning across threads.** `lovrlap` should use `REDUCTION(.OR.:lovrlap)`. `lfailure` should be handled via a shared flag that threads set but do not read for control flow within the parallel region (reading it for early exit is illegal — OpenMP forbids `RETURN`/early exit from parallel regions).

4. **`STOP` and `RETURN` are forbidden in OpenMP parallel regions.** The existing early-return pattern must be replaced with flag-setting. Specifically:
   - `if ( lfailure ) return` → `if ( lfailure ) lomp_failure = .true.`
   - `if ( lovrlap ) ... return` → set flag, skip inner work, continue loop

5. **`maddist` min-reduction initialization:** `maddist` is initialized to `1d10` before the parallel region. The `REDUCTION(min:maddist)` clause handles the per-thread initialization automatically (each thread's private copy starts at the initial value). This is correct.

6. **Array reductions (`vftwo`, `vfelec`):** OpenMP 3.1+ supports array reductions. Fortran compilers (gfortran 5+, ifort 16+) support `REDUCTION(+:vftwo)` for explicit-shape arrays. Confirm with `gfortran --version` before relying on this. For older compilers, manually sum into thread-local arrays and combine after the parallel region.

7. **`vforeign` recombination** happens AFTER the parallel region (line 460-462): `vforeign(l) = vfelec(l) + vftwo(l)`. Since `vftwo` and `vfelec` use REDUCTION, their values are correctly combined across threads before this summation. **No change needed to this code.**

8. **Coordinate data consistency:** The coordinates read via `twh_coordinates(GLB_GET, crdflag, ...)` are set before `twh_engtotal_twobody` is called and are not mutated during the loop. The caller guarantees that the coordinate arrays are stable for the duration of this call.

9. **OpenMP version:** Target OpenMP 3.1 minimum (for array reductions). Fortran 2008 is not required.

## SUCCESS CRITERIA

| Metric | Target | Measurement |
|--------|--------|-------------|
| Wall-clock speedup at 4 threads | 2.5-3.5x vs 1 thread | `time ./towhee` with `OMP_NUM_THREADS=1` vs `4` |
| Per-step energy deviation | ≤1×10⁻¹² relative | Diff energies across OMP runs |
| Reproducibility (same thread count) | Bitwise identical | 5 runs at 4 threads, diff output |
| Data races | 0 | `valgrind --tool=drd` clean |
| Compile with `--enable-openmp` | Clean build | `make` with flag enabled |
| Compile without `--enable-openmp` | No regression | Default build unchanged |
