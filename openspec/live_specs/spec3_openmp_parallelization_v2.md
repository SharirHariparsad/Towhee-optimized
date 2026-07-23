## REVISION NOTES (v2)

This revises the original spec after an implementation-readiness review. Nothing about
the overall approach changed — parallelize the `ichain` loop, reduce the accumulators —
but the following were fixed or added because they would have caused a compile failure,
a silent correctness bug, or an untested behavioral change:

1. **Compile-blocking bug:** `vftwo`/`vfelec` were listed in both `SHARED` and
   `REDUCTION` in the same directive. Fixed — removed from `SHARED`.
2. **Contradiction:** `lovrlap` was `PRIVATE` in the directive but `REDUCTION(.OR.:...)`
   in the guardrails section. Resolved in favor of `PRIVATE` (see Phase 2, Step 2.1)
   because the visible code only ever reads `lovrlap` within the same iteration it's
   set — it doesn't need to escape the region itself; only its downstream effect on
   `lfailure` does.
3. **Contradiction:** `lfailure` was "shared/read-only" in one table and `PRIVATE` in
   the directive. Resolved as `PRIVATE`, with a new `lomp_failure` shared flag set via
   `ATOMIC WRITE` for clean cross-thread propagation (see Phase 2, Step 2.3).
4. **Untested behavior change:** the flag-based early-exit replaces an immediate abort
   with "keep going until the region ends." Added named-loop `EXIT` logic to stop the
   *current thread's* remaining pair evaluations as soon as overlap is confirmed,
   instead of grinding through the rest of the O(N²) sweep on bad geometry. Added a
   dedicated overlap-path test (Phase 4) — the original test plan never exercised this
   code path at all.
5. **Audit gap:** the original GLB_GET safety table used approximate (`~`) line numbers
   and only checked for `GLB_SET`. Added an explicit requirement to re-verify against
   the current source and to check for `SAVE`/`COMMON` state, which a GLB_SET-only
   search would miss entirely (see Phase 3, Step 3.0).
6. **Missing load-balancing directive:** the `ichain`/`jchain` workload is triangular
   *and* filtered by box membership (multi-box GEMC runs), so default static scheduling
   will imbalance badly. Added an explicit `SCHEDULE` clause requirement.
7. **Build gap:** OpenMP flags were added to `FCFLAGS`/`CPPFLAGS` but not to link flags —
   added `AM_LDFLAGS` so linking against `libgomp` doesn't silently fail.
8. **Testing gap:** added a pre-change baseline capture so correctness is validated
   against the *original* serial behavior, not just OMP-1-thread vs OMP-4-thread
   self-consistency (which would hide a bug the refactor itself introduces).
9. Minor: corrected the mechanism description of the `maddist` min-reduction, and moved
   the race-detector run to a small/fast system instead of the 10,000+ atom performance
   system.

An "OPEN QUESTIONS" checklist is at the end — these need to be answered by reading the
actual current source before writing any `!$OMP` directive, not assumed from this spec.

---

## OBJECTIVE
Implement OpenMP multi-threading over the primary molecule-pair loop (`do ichain = 1, nchain - 1`) in `twh_engtotal_twobody` (`engtotal.F:334-458`) to parallelize the O(N²) intermolecular nonbonded energy evaluation across CPU cores. Every variable updated inside the parallel region must be either included in a `REDUCTION` clause or declared `PRIVATE`. No `GLB_SET` operations are allowed within the parallel region.

## TARGET FILES
- `Source/engtotal.F:334-458` — The `ichain` outer loop body (`twh_engtotal_twobody`)
- `Source/globaldata.F` — Audit all `twh_*` subroutines accessed from the loop to confirm safe read-only access
- `Source/vtwobody.F`, `Source/vcoulomb.F`, `Source/mimage.F` — Audit for hidden shared state (see Phase 3, Step 3.0)
- `Source/configure.ac`, `Source/Makefile.am` — Add OpenMP compiler and linker flags

## HOT LOOP CONTEXT (engtotal.F:334-458)

**Caveat before relying on this:** the block below is a representative excerpt, not a
confirmed complete transcript of lines 334-458. The 125-line target range is longer than
what's reproduced here, and the guardrails section references a `vforeign` summation at
lines 460-462 that falls outside even the stated range — meaning there is code adjacent
to (and possibly inside) the target loop that hasn't been walked line-by-line. **Before
writing the `!$OMP` directive, view the actual current file and confirm this excerpt
matches it**, including anything between the last visible `enddo` below and line 458.

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
            douter: do jchain = ichain + 1, nchain
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
                                 if ( lovrlap ) then
                                    call twh_engtotal_overlap(...)
                                    if ( lfailure ) then
                                       !$OMP ATOMIC WRITE
                                       lomp_failure = .true.
                                       exit douter
                                    endif
                                 endif
                              endif
                              if ( lelect .and. leval_c) then
                                 rij = dsqrt(rijsq)
                                 call twh_vcoulomb(...)
                                 velectinter = velectinter + vcoulomb
                                 if ( lovrlap ) then
                                    call twh_engtotal_overlap(...)
                                    if ( lfailure ) then
                                       !$OMP ATOMIC WRITE
                                       lomp_failure = .true.
                                       exit douter
                                    endif
                                 endif
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
            enddo douter
         endif
      enddo
```

**What changed from the original excerpt:** the `jchain` loop is now a named construct
(`douter`) and each `if (lovrlap) call twh_engtotal_overlap(...)` site now checks the
resulting `lfailure` and, if set, atomically raises the shared `lomp_failure` flag and
exits the named loop — abandoning only the *remaining pairs for the current `ichain`*
rather than either (a) illegally returning from inside a parallel region or (b) silently
computing the rest of the entire O(N²) sweep on top of a known-bad configuration. See
Phase 2, Step 2.3 for the reasoning.

## VARIABLE AUDIT & THREAD-SAFETY CLASSIFICATION

### Scalars updated inside the loop body (must be REDUCTION or PRIVATE)

| Variable | Nature | Classification | Rationale |
|----------|--------|---------------|-----------|
| `vtwointer` | Energy accumulator | `REDUCTION(+:vtwointer)` | Summed with `+` across pairs |
| `velectinter` | Energy accumulator | `REDUCTION(+:velectinter)` | Summed with `+` across pairs |
| `dl_twointer` | Derivative accumulator | `REDUCTION(+:dl_twointer)` | Summed with `+` across pairs |
| `dl_elect` | Derivative accumulator | `REDUCTION(+:dl_elect)` | Summed with `+` across pairs |
| `maddist` | `min` reduction (Madelung) | `REDUCTION(min:maddist)` | See note below on how this actually combines |
| `vftwo(1:scp_num_foreign_lambda)` | Array accumulator | `REDUCTION(+:vftwo)` | Summed element-wise. **Verify `vftwo` is a fixed-size (explicit-shape) array, not allocatable** — see Open Questions |
| `vfelec(1:scp_num_foreign_lambda)` | Array accumulator | `REDUCTION(+:vfelec)` | Summed element-wise, same caveat |
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
| `rijsq` | Pair distance² (reused for COM check and atom-atom check) | PRIVATE | Local to pair |
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
| `lovrlap` | Overlap flag | PRIVATE | **Resolved (was contradictory in v1).** Set as an output of `twh_vtwobody`/`twh_vcoulomb` and consumed within the same iteration — does not need to escape the region itself. Only `lomp_failure` needs to escape. *Confirm this against the real call signature — see Open Questions.* |
| `lfailure` | Failure flag | PRIVATE | **Resolved (was contradictory in v1).** Set by `twh_engtotal_overlap` within the current iteration; immediately checked and, if true, propagated to the shared `lomp_failure` via `ATOMIC WRITE`. *Confirm `lfailure` isn't actually a dummy argument of `twh_engtotal_twobody` itself — see Open Questions.* |

### Scalars and arrays READ-ONLY (shared, no conflict)

| Variable | Rationale |
|----------|-----------|
| `nchain` | Set once before the parallel region, read-only during loop |
| `ibox` | Function argument, constant |
| `crdflag` | Function argument, constant |
| `bxdflag` | Function argument, constant |
| `lmolecmove` | Function argument, constant |
| `lmovelist` | Function argument, constant |
| `cutstyle` | Local scalar set before loop. *Not visible in the excerpt above — confirm it's actually referenced in the real 334-458 range before including it* |
| `lelect` | Set before loop via `twh_lelect(GLB_GET,lelect)` |
| `lusecom` | Set before loop via `twh_lusecom(GLB_GET,lusecom)` |
| `rcelect` | Read via `twh_rcelect(GLB_GET,bxdflag,ibox,rcelect)` |
| `scp_num_foreign_lambda` | Read via `twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)` |
| `lmadelung` | Set before loop, read-only during |
| `lomp_failure` | **New.** Shared flag, initialized `.false.` before the region, written only via `ATOMIC WRITE`, read only after the region ends |

### GLB_GET calls — claimed read-only, needs re-verification

Each `twh_*(GLB_GET, ...)` call in the loop body is claimed to read from a SAVEd array or
scalar in `globaldata.F` without writing. **This has not been independently confirmed
line-by-line** (see Phase 3, Step 3.0/3.1) — treat the line numbers below as a starting
point for verification, not a confirmed audit:

- `twh_nboxi`, `twh_moltyp`, `twh_nunit`, `twh_ntype` — read arrays by index
- `twh_coordinates`, `twh_com`, `twh_rcmu`, `twh_lmoveme` — read per-molecule storage
- `twh_rcut`, `twh_rcelect`, `twh_scp_num_foreign_lambda` — read scalars
- `twh_mimage` — pure geometric computation, no global state (still needs a SAVE/COMMON check — see Phase 3.0)

### GLB_SET audit

No `GLB_SET` calls are visible inside `do ichain = 1, nchain - 1` in the excerpt shown.
**This cannot be marked "confirmed" until the full 334-458 range is viewed directly**,
per the caveat in the Hot Loop Context section above.

## STEP-BY-STEP IMPLEMENTATION PLAN

### Phase 0: Source Verification (Estimated: 1 hour) — NEW, do this first

Before writing any directive, resolve the items in **OPEN QUESTIONS** at the end of this
document by reading the actual current source. This phase exists because several v1
audit claims (approximate line numbers, an incomplete excerpt, unverified variable
scoping) turned out to be assumptions rather than confirmed facts. Treat this phase as a
hard prerequisite, not optional due diligence.

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
Add OpenMP compiler **and linker** flags. The original spec only added the compiler
flags; if the build separates compiling from linking, omitting `-fopenmp` at link time
produces an "undefined reference" failure against `libgomp` even though compilation
succeeded cleanly:

```
if UseOpenMP
OMP_FFLAGS = -fopenmp
OMP_CPPFLAGS = -DUSEOPENMP
endif

AM_FFLAGS += $(OMP_FFLAGS)
AM_CPPFLAGS += $(OMP_CPPFLAGS)
AM_LDFLAGS += $(OMP_FFLAGS)
```

**Note:** `-fopenmp` is for gfortran. For ifort, use `-qopenmp`. Recommend a
configure-time compiler check that sets the right flag rather than hardcoding gfortran's.

### Phase 2: Insert OpenMP Directives in engtotal.F (Estimated: 2.5 hours)

#### Step 2.1 — Add the PARALLEL DO directive

Insert before the `do ichain` loop at line 335 in `twh_engtotal_twobody`. This corrects
three things from the original directive: `vftwo`/`vfelec` are no longer duplicated in
`SHARED`, `lovrlap` and `lfailure` are consistently `PRIVATE`, and an explicit
`SCHEDULE` clause is added because the workload is not evenly distributed across
`ichain` (see Guardrail 6):

```fortran
#ifdef USEOPENMP
!$OMP PARALLEL DO DEFAULT(NONE)                                         &
!$OMP& SCHEDULE(DYNAMIC, 8)                                              &
!$OMP& PRIVATE(ichain, jchain, testboxi, testboxj, imolty, jmolty,     &
!$OMP&         numiunit, numjunit, lcompute, lmovemei, lmovemej,        &
!$OMP&         xcmi, ycmi, zcmi, rcmi, xcmj, ycmj, zcmj, rcmj,         &
!$OMP&         rcm, rcmsq, rijsq, rij,                                  &
!$OMP&         rxuij, ryuij, rzuij,                                     &
!$OMP&         rxui, ryui, rzui, rxuj, ryuj, rzuj,                     &
!$OMP&         ii, jj, ntii, ntjj, luse_ci, luse_vi, luse_cj, luse_vj, &
!$OMP&         leval_v, leval_c, vcoulomb, dvdl, l,                     &
!$OMP&         vf, lovrlap, lfailure)                                    &
!$OMP& SHARED(nchain, ibox, crdflag, bxdflag, lmolecmove, lmovelist,    &
!$OMP&        cutstyle, lelect, lusecom, lmadelung,                      &
!$OMP&        scp_num_foreign_lambda, rcelect, lomp_failure)             &
!$OMP& REDUCTION(+:vtwointer, velectinter, dl_twointer, dl_elect)        &
!$OMP& REDUCTION(+:vftwo, vfelec)                                       &
!$OMP& REDUCTION(min:maddist)
#endif
      do ichain = 1, nchain - 1
```

`SCHEDULE(DYNAMIC, 8)` is a starting point, not a final answer — benchmark against
`GUIDED` and against a few chunk sizes once correctness is confirmed (Phase 4.3). The
important part is that it's *not* the default static chunking, which would badly
under-utilize threads handling high-`ichain` iterations.

#### Step 2.2 — Add OMP END PARALLEL DO

Insert after the `enddo` at line 458 (before the `do l=1,MAX_FOREIGN_LAMBDA` summation):

```fortran
#ifdef USEOPENMP
!$OMP END PARALLEL DO
#endif
```

#### Step 2.3 — Handle early exit (lfailure / lovrlap)

**Critical constraint:** `RETURN` and `STOP` cannot be executed from inside a parallel
region in standard OpenMP. The original spec's fix — set a shared flag and let every
thread finish its full remaining workload before checking it — is legal OpenMP, but it
silently changes program behavior: the original serial code aborts the *entire*
remaining O(N²) sweep the instant an overlap is found, specifically because continuing
to evaluate `twh_vtwobody`/`twh_vcoulomb` on near-zero-separation geometry risks
`NaN`/`Inf` propagation into the energy sums, or a `SIGFPE` crash on builds with
floating-point-exception trapping enabled (common in scientific Fortran debug builds).

The fix here narrows the blast radius instead of eliminating it entirely (a full
region-wide abort still isn't legal OpenMP): each thread stops working on its *current*
`ichain` as soon as it detects overlap, rather than working through the rest of that
`ichain`'s pairs on data already known to be bad.

1. Declare a shared flag variable `lomp_failure` (add to `twh_engtotal_twobody` locals):
   ```fortran
   logical lomp_failure
   ```

2. Initialize to `.false.` before the parallel region:
   ```fortran
   lomp_failure = .false.
   ```

3. Make the `jchain` loop a named construct (`douter`, as shown in the Hot Loop Context
   section above) so it can be exited from within the nested `ii`/`jj` loops without
   touching the OMP-controlled `ichain` loop itself. **Exiting the loop that's directly
   associated with `!$OMP PARALLEL DO` is illegal** — only the nested `douter` loop is
   safe to exit this way.

4. At each `if (lovrlap)` site, after `call twh_engtotal_overlap(...)`, check the
   resulting `lfailure` and, if set, atomically raise the shared flag and exit the
   named loop:
   ```fortran
   if ( lfailure ) then
      !$OMP ATOMIC WRITE
      lomp_failure = .true.
      exit douter
   endif
   ```
   The `ATOMIC WRITE` is technically cheap insurance here — every thread that reaches
   this point is writing the same value (`.true.`), so the practical risk of an
   unsynchronized write is low on real hardware, but there's no reason to rely on that
   when the atomic costs nothing measurable.

5. After the parallel region, check and propagate:
   ```fortran
   if ( lomp_failure ) then
      lfailure = .true.
      return
   endif
   ```

**On `lovrlap` specifically:** it does not need a `REDUCTION` clause. Based on the
visible code, it's produced by `twh_vtwobody`/`twh_vcoulomb` as a per-call output and
consumed in the same iteration — a `PRIVATE` variable, not a value that needs to be
combined across threads. Only its downstream effect (via `lfailure` → `lomp_failure`)
needs to cross the thread boundary, and that's handled by step 4 above. *Confirm this
matches the actual call signatures of `twh_vtwobody`/`twh_vcoulomb` before relying on
it — see Open Questions.*

### Phase 3: Variable Audit and globaldata.F Confirmation (Estimated: 2 hours)

#### Step 3.0 — Check for hidden shared state (NEW, do this before 3.1)

A `GLB_SET`-only search misses a classic Fortran thread-safety hazard: a subroutine can
avoid the `globaldata.F` global-data API entirely and still hold cross-call state via a
`SAVE` attribute on a local variable or a `COMMON` block — both invisible to a
GLB_GET/GLB_SET audit, and both would be silently shared (and raced on) across threads
under OpenMP. Before trusting any "read-only" classification:

```bash
grep -n "SAVE\|COMMON" Source/vtwobody.F Source/vcoulomb.F Source/mimage.F Source/engtotal.F
```

For every hit, confirm the variable is either (a) set once at initialization and never
written again during energy evaluation, or (b) not actually reused across calls (e.g., a
`SAVE`'d array that's fully overwritten before being read on every call, which is safe).
Anything that accumulates or persists partial state across calls is a race under
parallel `ichain` iteration and must be refactored to be call-local (or `THREADPRIVATE`)
before this spec can proceed.

#### Step 3.1 — Confirm all GLB_GET-only access

For each `twh_*` subroutine called from the loop body, verify the code path guarded by
`flag .eq. GLB_GET` only reads from SAVEd arrays. **The line numbers below are starting
points, not confirmed — re-verify each one against the current source:**

| Subroutine | globaldata.F line (unverified) | Expected behavior |
|-----------|-------------------|-------|
| `twh_nboxi` | ~8388 | GLB_GET reads `nboxi(array)` |
| `twh_moltyp` | ~11464 | GLB_GET reads `moltyp(array)` |
| `twh_nunit` | ~11629 | GLB_GET reads `nunit(array)` |
| `twh_ntype` | ~11568 | GLB_GET reads `ntype(array)` |
| `twh_luse_c` | ~10053 | GLB_GET reads `luse_c(array)` |
| `twh_luse_v` | ~10090 | GLB_GET reads `luse_v(array)` |
| `twh_coordinates` | ~4153 | GLB_GET reads coordinate arrays, read-only |
| `twh_com` | ~3756 | GLB_GET reads COM arrays, read-only |
| `twh_rcmu` | ~14847 | GLB_GET reads rcmu array, read-only |
| `twh_lmoveme` | ~9603 | GLB_GET reads moveme array, read-only |
| `twh_rcut` | ~14447 | GLB_GET reads scalar rcut |
| `twh_rcelect` | ~14407 | GLB_GET reads per-box scalar |
| `twh_scp_num_foreign_lambda` | ~14975 | GLB_GET reads scalar |

#### Step 3.2 — Verify called subroutines (vtwobody, vcoulomb, mimage)

- **`twh_vtwobody`** (`vtwobody.F`): expected to read nonbonded coefficients via
  `twh_nbcoeff(GLB_GET, ...)` only, returning `vf`, `dvdl`, and `lovrlap` as call-local
  outputs. Confirm no `SAVE`/`COMMON` state per Step 3.0.
- **`twh_vcoulomb`** (`vcoulomb.F`): expected to read charge/qqscale via GLB_GET only.
  Same confirmation required.
- **`twh_mimage`** (`mimage.F`): expected to be a pure geometric computation. Same
  confirmation required — "pure" claims in legacy F77 code are worth double-checking.
- **`twh_engtotal_overlap`** (`engtotal.F:1573`): expected to read `twh_rminsq(GLB_GET)`
  and write `lfailure` as a call-local output argument. Confirm it performs no I/O that
  would interleave garbled output across threads, and no `SAVE`/`COMMON` state.

**Do not mark this "all clear" until every item above has been checked against the
actual current source** — see Open Questions.

### Phase 4: Testing & Validation (Estimated: 3 hours)

#### Step 4.0 — Capture a pre-change baseline (NEW, do this first)

Before making any code changes, run the *unmodified* code (no OpenMP directives, no
flag-based early-exit refactor) on two inputs and save the output as ground truth:
1. The standard TraPPE-UA validation input (n-butane/n-hexane Gibbs ensemble).
2. A deliberately bad starting configuration that forces an overlap/`lfailure` path
   (e.g., two atoms placed well inside the hard-core cutoff). This exercises code that
   the standard benchmark likely never touches.

Both baselines are needed because Phase 4.2 below must confirm the refactored code
matches the *original* behavior, not just that it's self-consistent across thread
counts — self-consistency alone can't catch a bug the refactor itself introduces.

#### Step 4.1 — Build test
```bash
./configure --enable-openmp
make -j
```
Confirm the binary links (checks the Phase 1.2 `AM_LDFLAGS` fix) and `engtotal.o` was
compiled with `-fopenmp`.

#### Step 4.2 — Correctness test (comparison against baseline, not just cross-thread)
```bash
export OMP_NUM_THREADS=1
./towhee < test_input > output_1t.txt
export OMP_NUM_THREADS=4
./towhee < test_input > output_4t.txt
```
Compare both `output_1t.txt` and `output_4t.txt` against the Phase 4.0 baseline, and
against each other, to within 1×10⁻¹² relative error (floating-point reassociation
across threads can cause non-bitwise-identical results even for the same arithmetic).

#### Step 4.3 — Overlap / failure-path test (NEW)

Run the Phase 4.0 overlap-triggering configuration through the OpenMP build at
`OMP_NUM_THREADS=1` and `4`. Confirm:
- `lfailure` propagates correctly back to the caller in both cases.
- The result matches the pre-change baseline's behavior for the same input.
- No `NaN`/`Inf` appears in any energy accumulator, and no `SIGFPE` occurs, on a build
  with `-ffpe-trap=invalid,zero,overflow` enabled.

#### Step 4.4 — Performance test
Run a large system (10,000+ atoms, e.g., SPC/E water) and measure wall time at
`OMP_NUM_THREADS=1,2,4,8`. Expected scaling: near-linear up to available cores, with
diminishing returns when memory-bandwidth-bound — **conditional on the `SCHEDULE`
clause added in Step 2.1**; without it, expect the imbalance described in Guardrail 6 to
suppress scaling well below the target.

#### Step 4.5 — Race condition / thread-safety stress test
Run the correctness test 5 times with 4 threads — verify bitwise-identical output across
all 5 runs (within the 1e-12 tolerance). Non-deterministic results indicate a missing
`PRIVATE` clause or a data race.

#### Step 4.6 — valgrind DRD/Helgrind test
```bash
export OMP_NUM_THREADS=2
valgrind --tool=drd --check-stack-var=yes ./towhee < small_test_input
```
**Use a small, fast input here — not the 10,000+ atom system from Step 4.4.** DRD/Helgrind
commonly impose 10-100x slowdowns; running them against a production-scale system is
impractical. A small system exercises the same code paths for race detection purposes.

#### Step 4.7 — Relative hotspot check (optional but recommended)

`twh_engtotal_twobody` performs a full O(N²) recompute, which in a typical Monte Carlo
run is usually invoked far less often than the per-move delta-energy routines (the
`energy_change.F` work mentioned as a companion target). Before expecting this change
alone to deliver the overall simulation speedup implied by the Success Criteria below,
profile a representative production run to confirm what fraction of total wall time
`twh_engtotal_twobody` actually accounts for.

## CONSTRAINTS & GUARDRAILS

1. **No `GLB_SET` allowed inside the parallel region.** If future code changes introduce
   a `GLB_SET` call in the loop body, the parallelization must be re-audited.

2. **The `PRIVATE` list must be exhaustive.** Omitting a variable from `PRIVATE` that is
   written-to per-iteration causes undefined behavior (data race). `DEFAULT(NONE)` is
   used specifically to force the compiler to error on any unlisted variable — do not
   remove it.

3. **`lfailure` and `lovrlap` are both `PRIVATE`, not `REDUCTION`.** Only the derived
   shared flag `lomp_failure` crosses the thread boundary, via `ATOMIC WRITE`. This
   replaces the v1 guardrail that recommended `REDUCTION(.OR.:lovrlap)` — that
   recommendation directly contradicted the directive's own `PRIVATE` clause and, based
   on the visible code, wasn't needed in the first place.

4. **`STOP` and `RETURN` are forbidden in OpenMP parallel regions**, and so is exiting
   the loop directly associated with `!$OMP PARALLEL DO`. The named `douter` construct
   around `jchain` exists specifically so a thread can bail out of its *current*
   `ichain`'s remaining work without violating either rule.

5. **`maddist` min-reduction:** `maddist` is initialized to `1d10` before the parallel
   region. Note the actual mechanism: each thread's private copy is initialized by the
   compiler to the implementation's maximum representable value (not literally to
   `1d10`), accumulates its own minimum during the loop, and then **all thread results
   are combined with the original pre-region value of `maddist` using `min()`** after
   the region ends. Since `1d10` is safely larger than any real interatomic distance,
   the end result is correct — but don't reuse this "each thread starts at the initial
   value" mental model for a different reduction where the pre-loop value isn't a safe
   sentinel; it isn't literally what happens.

6. **Array reductions (`vftwo`, `vfelec`):** OpenMP 3.1+ supports array reductions.
   Fortran compilers (gfortran 5+, ifort 16+) support `REDUCTION(+:vftwo)` for
   explicit-shape arrays specifically — confirm `vftwo`/`vfelec` are declared with fixed
   dimensions (not allocatable) before relying on this, and confirm the target build's
   `gfortran --version`. For older compilers, manually sum into thread-local arrays and
   combine after the parallel region.

7. **`vforeign` recombination** happens AFTER the parallel region: `vforeign(l) =
   vfelec(l) + vftwo(l)`. Since `vftwo` and `vfelec` use `REDUCTION`, their values are
   correctly combined across threads before this summation. No change needed to this
   code, contingent on Guardrail 6 being satisfied.

8. **Coordinate data consistency:** the coordinates read via `twh_coordinates(GLB_GET,
   crdflag, ...)` are assumed to be set before `twh_engtotal_twobody` is called and not
   mutated during the loop. Confirm this assumption holds for every caller of
   `twh_engtotal_twobody`, not just the benchmark path being used for validation.

9. **OpenMP version:** target OpenMP 3.1 minimum (for array reductions). Fortran 2008 is
   not required.

10. **Load balancing (NEW):** the `ichain`/`jchain` workload is triangular (later
    `ichain` values have fewer `jchain` partners) and additionally filtered by
    `testboxi .eq. ibox` in multi-box GEMC runs, meaning many `ichain` iterations do
    zero inner work. Default static scheduling will assign a disproportionate share of
    real work to low-`ichain` threads. The `SCHEDULE(DYNAMIC, 8)` clause in Step 2.1 is
    a required starting point — benchmark it against `GUIDED` and alternate chunk sizes
    in Step 4.4, don't ship without an explicit `SCHEDULE`.

11. **Hidden shared state (NEW):** a `GLB_SET`-only audit does not catch `SAVE`
    variables or `COMMON` blocks used as scratch space inside called subroutines. Phase
    3, Step 3.0 must be completed — and must find nothing that persists mutable state
    across calls — before this spec's safety claims hold.

## SUCCESS CRITERIA

| Metric | Target | Measurement |
|--------|--------|-------------|
| Wall-clock speedup at 4 threads | 2.5-3.5x vs 1 thread | `time ./towhee` with `OMP_NUM_THREADS=1` vs `4`, contingent on Guardrail 10 (explicit `SCHEDULE`) |
| Per-step energy deviation vs. pre-change baseline | ≤1×10⁻¹² relative | Diff against Phase 4.0 baseline, not just against other thread counts |
| Overlap/failure-path behavior | Matches pre-change baseline | Phase 4.3 |
| Reproducibility (same thread count) | Bitwise identical | 5 runs at 4 threads, diff output |
| Data races | 0 | `valgrind --tool=drd` clean on a small test input |
| Compile with `--enable-openmp` | Clean build, links successfully | `make` with flag enabled |
| Compile without `--enable-openmp` | No regression | Default build unchanged |

## OPEN QUESTIONS — resolve before writing code

These require reading the actual current source, not this spec, to answer:

1. Does the full `engtotal.F:334-458` range match the excerpt above? What's in the gap
   between the visible `enddo` and line 458 (the `vforeign` summation is referenced at
   460-462, outside the stated range — confirm the boundaries)?
2. Is `cutstyle` actually used anywhere in this specific loop, or is it a holdover from
   a different code path? It doesn't appear in the excerpt.
3. Is `lfailure` a dummy argument of `twh_engtotal_twobody` (passed by the caller) or a
   purely local variable? This determines whether the `PRIVATE` classification is
   correct or whether special handling is needed for the pre/post-region values.
4. Do `twh_vtwobody` and `twh_vcoulomb` actually return `lovrlap` as a call-local output
   argument, as assumed? Confirm the real call signatures (the excerpt elides the
   arguments as `...`).
5. Are `vftwo` and `vfelec` declared as fixed-size (explicit-shape) arrays, or
   allocatable? This determines whether `REDUCTION(+:vftwo, vfelec)` is portable across
   the target compiler versions.
6. Do `vtwobody.F`, `vcoulomb.F`, or `mimage.F` contain any `SAVE` variables or `COMMON`
   blocks used as cross-call scratch state? (Phase 3, Step 3.0)
7. What fraction of total simulation wall time does `twh_engtotal_twobody` actually
   account for in a representative production run, relative to `energy_change.F`?
   (Phase 4.7 — informs whether this is the right first target for overall speedup.)
