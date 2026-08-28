# Towhee 8.2.3 Optimization Build Recipe

## Prerequisites
- gfortran (tested with GCC 11+)
- GNU make
- autotools (autoconf, automake) for full rebuild from Makefile.am

## Reproducible Build Commands

### Original (serial baseline)
```bash
cd src/towhee-8.2.3
./configure FFLAGS="-g -O2 -ffixed-line-length-none"
cd Source
make -j16 towhee
```

### OpenMP-only build (spec3: engtotal_twobody parallelization)
```bash
cd src/towhee-8.2.3
./configure --enable-openmp FFLAGS="-g -O2 -fopenmp -ffixed-line-length-none" LDFLAGS="-fopenmp"
cd Source
make -j16 towhee
```

### HOT_PATH_INLINE build (spec2 + spec3 combined)
```bash
cd src/towhee-8.2.3
./configure --enable-openmp FFLAGS="-g -O2 -fopenmp -ffixed-line-length-none -DHOT_PATH_INLINE" LDFLAGS="-fopenmp"
cd Source
make -j16 towhee
```

Note: HPI build requires `engtotal_globals.F` in the source list (add to Makefile.in):
```
energy_change.$(OBJEXT) engatom.$(OBJEXT) engtotal.$(OBJEXT) engtotal_globals.$(OBJEXT) \
...
engtotal.F engtotal_globals.F \
engtotal.o: preproc.h functions.h engtotal.F
engtotal_globals.o: preproc.h functions.h engtotal_globals.F hot_common.h
```

## Key Fixes Applied to Source

### engtotal.F (spec3 OpenMP region)
1. **lfailure**: Changed from `PRIVATE` → `FIRSTPRIVATE` (dummy arg was uninitialized)
2. **lovrlap**: Added explicit `lovrlap = .false.` before `twh_vtwobody` and `twh_vcoulomb` calls in non-HPI branch (function leaves it unset on beyond-cutoff return)
3. **lonefour_global**: Changed from `PRIVATE` → `FIRSTPRIVATE` (shared read-only flag)
4. Debug: Added `write(6,*) 'DEBUG twobody: lomp_failure=',lomp_failure` after parallel region (temporary)

### engtotal_globals.F (HPI module)
- `engtotal_globals_init()` must be called **after** `setmixterms`/`setclassical` complete — currently called at `readtowhee.F:2242` which is BEFORE mixing rules → nbcoeff copies are stale

## Benchmark Results Summary

| Build | 1T | 4T | 8T | 16T | Notes |
|-------|----|----|----|-----|-------|
| Original (serial) | 8.5s | 8.6s | 8.4s | 8.6s | Baseline |
| OMP-only (spec3) | 9.7s | 10.0s | 9.8s | 10.1s | **~17% SLOWER** - wrong loop parallelized |
| HPI+OMP | 9.7s | nondet. | nondet. | nondet. | nbcoeff snapshot timing bug |

## Root Cause: Profile Mismatch
gprof shows time is in accessor functions:
- `twh_hmatrix`: 123M calls (12%)
- `twh_mimage`: 40M calls
- `twh_coordinates`: 80M calls
- `twh_1dinteger/1ddouble` getters: ~12%

While `twh_engtotal_twobody` (the parallelized loop) is only ~3% of runtime.

## Recommendation
To get real speedup, implement **spec2 (HOT_PATH_INLINE)** properly:
1. Move `engtotal_globals_init()` call to **after** `twh_setclassical` in `readtowhee.F`
2. Verify nbcoeff copies stay in sync (no GLB_SET after init)
3. Test HPI determinism at 1T/16T before benchmarking