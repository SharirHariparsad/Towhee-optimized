# Phase 6+ Additional Improvements — Detailed Specification

**Date**: 2026-09-01
**Base**: `openspec/live_specs/optimization_plan.md:1` Phase 5 complete (HPI+Vec+Ewald, 1.02-1.15x on `tests/SPC_E_WATER`, `src/TMMC/SPC_E_WATER`, `tests/TraPPE_Pentane` small examples)
**Goal**: Extend single-node headroom after Phase 5 and close gap to `optimization_plan.md:4` target `5-10x` single-node, near-linear MPI to 64 nodes
**Compiler**: `gfortran + OpenMPI` `-O3 -march=native -ffast-math -fopenmp` `Source/Makefile.am:38` `VEC_FFLAGS`
**Regression**: `1e-6` relative `tests/comparison_benchmark.py:compare_outputs` + `src/towhee-8.2.3/Examples/tolerance_diff.py`; OpenMP `USEOPENMP` out-of-scope `optimization_plan.md:25` (deterministic corruption from move 1)

---

## 1. Rationale

`optimization_plan.md:43` gprof `Wall_Water`:

| Rank | Function | % | Calls | Category |
|------|----------|---|-------|----------|
| 1 | `twh_cbmcbendinfo_` `globaldata.F:16771` | 13.5% | 540K | CBMC bend table |
| 2 | `twh_hmatrix_` `globaldata.F` | 8.1% | 150M | Box accessor |
| 3 | `twh_halfbox_` | 6.4% | 150M |  |
| 4 | `twh_mimage_` `mimage.F:199` | 6.1% | 50M |  |
| 5 | `twh_eng_two_body_` `energy.F:164` | 5.3% | 11K | Dispatch |
| 6 | `twh_vcoulomb_` `vcoulomb.F:22` | 4.3% | 36M |  |
| 7 | `twh_1dinteger` `globalinfo.c:264` | 4.3% |  |  |
| 8 | `twh_1ddouble` `globalinfo.c:136` | 3.0% |  |  |
| 9 | `twh_luse_v_` | 2.8% | 43M |  |
| 10 | `twh_coordinates_` `globalinfo.c:1034` | 2.5% | 65M |  |
| — | `twh_engtotal_twobody_` `engtotal.F:335` | 0.25% | 1 | OMP target |

Phases 1-5 eliminated `hmatrix/halfbox/mimage/rc*` in the **full** `O(N^2)` `engtotal.F:335` `twh_engtotal_twobody` path via `engtotal_globals.F` COMMON (`/hot_hmatrix/`, `/hot_halfbox/`, `/hot_rcut/`, `/hot_rcelect/`, `/hot_qqscale/`, `/hot_ewald_calp/`) and `getcbangle.F:138` `/hot_cbmcbendinfo/` (2986->5 calls `optimization_plan.md:195`). `optimization_plan.md:369` wall-clock stays `~1.04s min5` parity on `SPC_E_WATER` because real work is **differential** `energy_change.F:2383` -> `energy.F:164` `twh_eng_two_body` per trial + CBMC `getcbangle.F:1781`/`getcbbond.F:484`/`getcbdihed.F:791` + `ranlux.F:980` `11%` + remaining accessors `energy.F:341,353` not `engtotal.F`. `bench_results.json:19` `7.93s min` shows same.

Remaining headroom is therefore **differential path + RNG + Ewald reciprocal completion + build/PGO + I/O + cell list** (only `>50K` atoms).

---

## 2. Improvement P6 — Differential Path Zero-Copy (Highest ROI)

### 2.1 Problem

`mainloop.F:384` `twh_main_do_MC_move` -> `energy_change.F:268` differential calls dominate. `energy.F:164` per `jlist` `1..ncmt` `energy.F:277` x `jj 1..numjunit` `energy.F:403` x `iunit 1..numiunit` `energy.F:410` still does per-pair dispatch:

* `energy.F:402` `twh_nunit(GLB_GET,jmolty,numjunit)` + `energy.F:404` `twh_ntype` + `energy.F:407` `twh_luse_c` + `energy.F:408` `twh_luse_v` (j side) + `energy.F:411` `twh_ntype` (i side) — `6` calls before physics
* `energy.F:341,353` `twh_coordinates(GLB_GET,*,jj,rxj..)` `65M/2.5%` + `energy.F:417` pair `rxi` — two dispatches per pair
* `energy.F:356` `twh_mimage` `6.1%` (COM prefilter `energy.F:311` also `twh_mimage`)
* `energy.F:216` `twh_rcelect`/`twh_rcutinsq` per `ichainmove` (box-global, hoistable)
* `energy.F:248` `twh_ctrmas` + `energy.F:252,258` `twh_com`/`twh_rcmu` per move even when `lusecom` false `energy.F:215`
* `engtotal.F:256` `twh_atomnum` `-> twh_globalpos(GLB_GET,ichain,globalpos)` `globalinfo.c:1158` inside `inline_atomnum` `engtotal.F:518` — `65M` retained after Phase 1, `ivalue = globalpos + 3*(iunit-1)` invariant

Hot `HOT_PATH_INLINE` `energy.F:333` already caches `ntype_arr/luse_*_arr/nunit_arr` from `engtotal_globals` but **coordinates/mimage/globalpos** still dispatch.

### 2.2 Solution

Extend Phase 1.2 COMMON pattern (`hot_common.h`, `globaldata.F`, `engtotal_globals.F`) to differential path, guarded `#ifdef HOT_PATH_INLINE` with `#else` preserved.

**2.2.1 Precomputed `atomnum_arr`**

```fortran
! engtotal_globals.F
integer, dimension(MAXCHAIN,NUMAX) :: atomnum_arr
common /hot_atomnum/ atomnum_arr  ! or module array
! engtotal_globals_init: after grow arrays built
do ichain=1,nchain
  call twh_globalpos(GLB_GET,ichain,gpos)
  do iunit=1,nunit_arr(moltyp(ichain))
    atomnum_arr(ichain,iunit)=gpos+3*(iunit-1)
  enddo
enddo
! refresh on nchain change (swap moves: update inserted chain only)
```

Replace `engtotal.F:261` `call twh_globalpos` in `inline_atomnum` and `energy.F:341` pattern with `atomnum_arr(ichain,iunit)`. Saves `~1` `twh_1dinteger` per inner iteration. Estimated `2-3%` pair loop.

**2.2.2 Flat coordinate pointer**

Add `get_coord_storage_ptr` already in `globalinfo.c:1034` (`twh_coordstorage_` `PNT_COORDSTORAGE` `preproc.h:716`). Expose via COMMON:

```fortran
! hot_common.h
double precision c_coord_real(3,MAXATOM), c_coord_temp(3,MAXATOM)
common /hot_coords/ c_coord_real, c_coord_temp
! engtotal_globals.F binds via c_loc on globaldata SAVE or via C accessor at init
```

In `energy.F:341` replace:

```fortran
#ifdef HOT_PATH_INLINE
rxj = c_coord_real(1, atomnum_arr(jchain,jj)/3+1) ! or direct index atomnum
#else
call twh_coordinates(GLB_GET,jcrd,jchain,jj,rxj,ryj,rzj)
#endif
```

Careful: `atomnum` is `1-indexed` Fortran triple offset. Use `idx = atomnum_arr/3+1`. Same for `testcrd` `CRD_REAL` vs `CRD_TEMP` select pointer. Eliminates `2x twh_coordinates` per pair (`~4.5M` dispatches `SPC_E_WATER` `100K steps`). Do same for `writemovie.F:130`/`writepdb.F:177`/`rwconf.F:764` I/O paths (P12).

**2.2.3 Hoist box-global scalars and COM**

* Hoist `rcelect/rcut/rcutinsq/boxvolume` before `do ichainmove:227` — already COMMON `/hot_rcelect/` `vcoulomb.F:64` but `energy.F:216` still fetches. Replace with `c_rcelect(bxdflag,ibox)`.
* Guard `energy.F:240` `if(lusecom)` entire `ctrmas` block: many systems `lusecom=.false.` `twh_lusecom` — skip `twh_com` `energy.F:258` + `twh_rcmu` `energy.F:260`. Add `if(lusecom.and.lmolecmove)` where `lmolecmove` derived from `chaininfo`.
* Hoist `twh_nchain/twh_nmolty/twh_numboxes` `energy.F:213` out of `twh_eng_two_body` to `engtotal.F:447` local `nchain` passed as arg to avoid per-call fetch.

**Files**: `Source/energy.F:164,213,240,311,341,402`, `Source/engtotal.F:256,518`, `Source/engtotal_globals.F`, `Source/globaldata.F`, `Source/hot_common.h` (new), `Source/globalinfo.c:1034`, `Source/mainloop.F:152`

**Validation**: `gprof` `twh_coordinates_` `twh_globalpos` should drop `>90%`; `0-diff` `tests/SPC_E_WATER` baseline `HOT_PATH_INLINE` off untouched; `tolerance_diff.py 1e-6` for `-ffast-math`.

**Est. gain**: `5-8%` small (CBMC-bound) / `8-12%` large pair-loop heavy (`Σ p*f ~10-15%`).

---

## 3. Improvement P7 — RNG Bulk / Default Switch

### 3.1 Problem

`ranlux.F:980` `dranlux(rvec,1)` per trial `getcbangle.F:560` `rbf`, `getcbangle.F:574` `rdecide`, `getcbdihed.F:791` `~2-3` calls per trial, `nch_bend_a/b` `6-10` trials per step `atd_cbmc_setup.F:2182` -> `~18-30` RNG calls per CBMC regrowth. Profile `11%` `optimization_plan.md:204`. `preproc.h:669` defines `RNG_RANLUX 0`, `RNG_DX_1597_2_7 5`, `RNG_KISS99 7`, `RNG_MRG32K3A 8` runtime-selected `initialize.F:42` `call twh_random_code(GLB_SET,RNG_DX_1597_2_7)` default actually `RANLUX` luxury `3` `initialize.F:43` slow branch, not vectorizable, `ranlux(rvec,1:108)` luxury shuffle.

### 3.2 Solution

**Option A (1 line, low risk) — Switch default**:

```fortran
! Source/initialize.F:42
call twh_random_code(GLB_SET,RNG_DX_1597_2_7) ! was RNG_RANLUX
! or RNG_KISS99 / RNG_MRG32K3A
```

`dx-1597-2-7.c:80` xorshift `~5-8x` faster, `mrg32k3a.c:96` also fast, `kiss.F:70` fastest but shorter period. No API change, `twh_random` dispatch in `functions.h` unchanged. Trajectory diverges but block averages within `1e-6` `comparison_benchmark.py` — verify `tests/TraPPE_Pentane` `0.554s` tolerance.

**Option B — Bulk API** (if A insufficient):

```c
// globalinfo.c / rnginit.c:189
void twh_random_vec(int n, double *buf); // fill n in one call
```
```fortran
! getcbangle.F:560
double precision rbuf(64)
call twh_random_vec(ntrial, rbuf) ! then consume rbuf(itrial)
```

Cuts `luxury` loop overhead `~30%` additional.

**Files**: `Source/initialize.F:42`, `Source/rnginit.c:189`, `Source/ranlux.F:108`, `Source/dx-1597-2-7.c:80`, `Source/getcbangle.F:560,574,777`, `Source/cbregrow.F:437`

**Validation**: All `Examples/*` `do_test` `tolerance_diff.py 1e-6`; performance `gprof` `ranlux_` `11%->~2%`.

**Est. gain**: `5-8%` CBMC-heavy small / `3-5%` Ewald-heavy.

---

## 4. Improvement P8 — CBMC Completion + Coefficient Caches

### 4.1 Remaining `cbmcbendinfo` Sites

Phase 2 patched `getcbangle.F:138` + `angle_b_distribution.F:74` via `/hot_cbmcbendinfo/` `globaldata.F`. Remaining dispatch in `#else` branches `getcbangle.F:514,528,534,744,792` `DIST_*` `preproc.h:559` and `dihedral_distribution.F:422` `twh_cbmctordistribution` `CBMCMAXTOR 4` + `getcbdihed.F:791` `twh_cbmctorinfo` `CBMC_MAXTORPEAKS 6` `preproc.h:44`. Apply same COMMON direct-read pattern, guarded.

### 4.2 `bencoeff/vibcoeff/torcoeff` Caches

`vangle.F:126` `do ibcoeff=MINBENDCOEFF(-1):MAXBENDCOEFF(10)` `preproc.h:197` `call twh_bencoeff(GLB_GET,itype,ibcoeff,dvalue)` `12` calls per angle `vangle.F:22` `464 LOC` style switch `ANG_STYLE_HARMONIC 1` `preproc.h:603`. Similarly `vbond.F:248` `vibcoeff:189` `MAXVIBCOEFF 3` `4` calls, `vtorsion.F:731` `torcoeff` `MAXTORCOEFF 31` `32` but only `~4-6` used per style `TOR_STYLE_TRAPPE_COS 7` `preproc.h:644`.

Cache at `setmixterms.F:976` startup:

```fortran
! engtotal_globals.F
double precision ben_coeff_cache(TBENMAX, MINBENDCOEFF:MAXBENDCOEFF)
common /hot_bencoeff/ ben_coeff_cache
! setmixterms: after twh_bencoeff(GLB_SET) loop, copy:
ben_coeff_cache(itype,icoeff)=bencoeff(itype,icoeff)
```

In `vangle.F:126` replace loop with `bencoeff_arr(icoeff)=ben_coeff_cache(itype,icoeff)` indexed. Same `vib_coeff_cache(TVIBMAX:8200,0:3)` `preproc.h:91`, `tor_coeff_cache(TTORMAX:930,-1:31)`.

### 4.3 `qqscale` Table

`vcoulomb.F:83` `twh_qqscale` per pair `qq = qqa*qqb` where `qqa=twh_qqatom(imolty,iunit)` `globaldata.F` `qqatom(NTMAX,NUMAX)` `7*2202=15k`. Precompute `qq_table(NTMAX,NUMAX,NTMAX,NUMAX)` too large; instead `charge_arr(NTMAX,NUMAX)` via COMMON `/hot_qqatom/` already? `vcoulomb.F:64` has `/hot_qqscale/` but `twh_qqatom` still dispatch `recipsum.F:148`. Add `/hot_qqatom/` COMMON read `c_qqatom(imolty,iunit)` direct array.

**Files**: `Source/getcbangle.F:514`, `Source/angle_b_distribution.F:74`, `Source/dihedral_distribution.F:422`, `Source/getcbdihed.F:791`, `Source/vangle.F:126,129`, `Source/vbond.F`, `Source/vtorsion.F:731`, `Source/vcoulomb.F:83,148`, `Source/globaldata.F:twh_bencoeff`, `Source/setmixterms.F:976`, `Source/engtotal_globals.F`

**Est. gain**: `1-2%` small / `1-2%` large (`0.8%` bending measured `P8` audit).

---

## 5. Improvement P9 — Ewald Reciprocal Completion

### 5.1 Problem

`recipsum.F:24` `172 LOC` `twh_recipsum` called per volume move `mainloop.F:731` but also per `engtotal`? Triple loop `nkx 0..kmax:81` `nky:90` `nkz:100` `O(kmax^3)` `kmax=MAXKMAX 12` `preproc.h:68` `VECTORMAX 5184`. Inner `ichain 1..nchain:139` `iunit 1..nunit:145` `dcos/dsin` per atom per k-vector `recipsum.F:150`.

* `recipsum.F:116-127` `9x twh_hinverse` per k-vector `5184*3* ~15k calls` — `hinverse(3,3)` constant per box per call `Source/recipsum.F:116` `call twh_hinverse(GLB_GET,bxdflag,ibox,1,1,hinv1)` etc.
* `recipsum.F:133` `prefact = invvol*twh_expon(-ksqr*invalpsq4)/ksqr` `expon.F:47` `if(arg<-700) 0` branch breaks `-ftree-vectorize`.
* `recipsum.F:144` `!$OMP SIMD REDUCTION(+:sumr,sumi)` only `iunit`, outer `nkx/nky/nkz` scalar.
* `vcoulomb.F:141` `twh_derfc` `FUNCTION_DERFC` per real-space pair `36M` calls `optimization_plan.md:51` scalar polynomial.
* `vcoulomb.F:139` `twh_ewald_calp` now COMMON `/hot_ewald_calp/` `recipsum.F:48` but `twh_ewald_kvector:128` `twh_ewald_prefact:134` `twh_ewald_ssumr:155` still dispatch per k-vector.

### 5.2 Solution

**5.2.1 Hoist `hinverse`/`boxvolume`**

```fortran
! recipsum.F:55 before nkx loop
double precision hinv(3,3), invvol, twopi, qqfact
call twh_hinverse(GLB_GET,bxdflag,ibox,1,1,hinv(1,1)) ! 9 calls total -> 1x
! or via COMMON /hot_hinverse/ hinv_arr(MAXBXD,MAXBOX,3,3) similar hmatrix
...
do nkx=... 
  kx = twopi*(n_one*hinv(1,1)+n_two*hinv(1,2)+n_thr*hinv(1,3))
```

**5.2.2 Branchless `expon` + SIMD `exp`**

```fortran
! expon.F:47
double precision function twh_expon(arg)
  twh_expon = exp(max(arg, -700d0))
end
! declare !$OMP DECLARE SIMD(twh_expon)
```

Enables `-ffast-math` auto-vector `vdexp`.

**5.2.3 k-vector/prefact generation**

Move `twh_ewald_kvector:128`/`twh_ewald_prefact:134` generation to `engtotal_globals_init` precompute `kvecs(3,VECTORMAX)`, `ksqr_arr`, or keep per-volume but hoist `hinverse` already cuts `15k` calls to `9`.

**5.2.4 `derfc` vectorization**

`vcoulomb.F:141` replace scalar `derfc(calp*rij)` with rational approx `erfc(x)=exp(-x*x)*poly(x)` or table linear interp `derfc_table(0:2000)`; declare `!$OMP DECLARE SIMD(derfc)`.

**Files**: `Source/recipsum.F:24,55,116,133,144`, `Source/vcoulomb.F:141`, `Source/expon.F:47`, `Source/globaldata.F:twh_hinverse`, `Source/engtotal_globals.F`, `Source/hot_common.h`

**Est. gain**: `1-2%` small (`6.3%` * `30-50%` saved) / `10-18%` large `Ewald 25-40%` (`large` recipsum dominates `O(kmax^3 * natom)`).

---

## 6. Improvement P10 — Branchless Minimum Image + `sqrt` Motion

### 6.1 Problem

`engtotal.F:288-313` `inline_mimage` `do while(rxuij>hbx) rxuij-=hm` loop-carried branch `trip 0-1` inhibits `!$OMP DECLARE SIMD(inline_mimage)` `vtwobody.F:1046` `engtotal.F:272`. `energy.F:380` `rij = dsqrt(rijsq)` for Coulomb even when `leval_c` false.

### 6.2 Solution

```fortran
! engtotal.F:288, mimage.F:199
!$OMP DECLARE SIMD(inline_mimage)
subroutine inline_mimage(...)
  rxuij = rxuij - anint(rxuij*hinv_diag)*hmatrix(1,1) ! using hinv=1/h
end
```

Theory: `hmatrix` orthorhombic fast path `3` dims, else `hinverse`. Transform `while` to `nint`. Move `energy.F:380` `rij=dsqrt` inside `if(lelect.and.leval_c)` `energy.F:379` (already `engtotal.F:596` but not `energy.F`).

**Files**: `Source/engtotal.F:272,288`, `Source/mimage.F:199`, `Source/vtwobody.F:1046`, `Source/energy.F:380`

**Est. gain**: `1%` scalar + enables SIMD `1-2%`.

---

## 7. Improvement P11 — Build System LTO + PGO

### 7.1 LTO

`Source/Makefile.am:38` currently `VEC_FFLAGS=-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize`:

```make
VEC_FFLAGS += -flto -ffat-lto-objects
VEC_CFLAGS += -flto -ffat-lto-objects
AM_LDFLAGS += -flto
```

Inlines `getnbtype.F:69` `twh_getnbtype` triangular `if(jtype<itype) swap` into `vtwobody.F:149` `vtb_lennard_jones:922` `c_nb13/14` `vtwobody.F:83` LTO pass. Cross-TU inline `dotproduct.F:58` `distance.F:39`.

### 7.2 PGO

Add `configure.ac:82` `--enable-pgo`:

```bash
make FFLAGS="-fprofile-generate" && ./towhee # run SPC_E_WATER+TraPPE
make FFLAGS="-fprofile-use -fprofile-correction"
```

Optimizes `vangle.F:146` style switch predictor (`ANG_STYLE_HARMONIC 1` hot `preproc.h:603`) + `vtwobody` `CP_LENNARD_JONES 1` `preproc.h:258`.

### 7.3 `march` Guard

`AX_CHECK_COMPILE_FLAG([-march=native])` fallback `-march=x86-64-v3` in `configure.ac:95` before `AC_PROG_F77`.

**Files**: `Source/Makefile.am:38,53,55`, `Source/configure.ac:82-127`, `Source/Makefile.in`

**Est. gain**: `2-3%` whole program (`2-4%` reported `-flto` typical gfortran).

---

## 8. Improvement P12 — I/O + Data Layout

### 8.1 I/O Not Hot for Benchmark but `>10K` Atoms

`writepdb.F:22` `76` `call twh_openfile` `open UNKNOWN` per frame, `writepdb.F:183` `WRITE(base+ibox,12) '(A6,I5...3f8.3)'` per atom, `6` accessors `writepdb.F:179` `pdbname` `ntype` `elementname` `qqatom` `getnbtype` `nbcoeff` per atom `writemovie.F:130` `2` writes per atom `write(ifile,'(4(1x,e15.6))')`, `rwconf.F:764` `twh_coordinates` per atom backup `mainloop.F:731` every `backupfreq`.

When `pdboutfreq/moviefreq` `mainloop.F:342,680` `0 (off)` for shipped Examples, `0%`. When on (`100-1000`), `10K atoms * 100 steps = 1M lines` dominates `sprintf` + `open/close`.

**Solution**: Add `form='unformatted',access='stream'` option controlled by `lpdb_binary` flag `writepdb.F:76` + `open(...,buffered='yes',blocksize=1048576)` `gfortran` or `MPI_File_write_all` `control.c:182` for MPI. Precompute immutable per `imolty,iunit` `pdbname/elementname` at `assemble.F:2554`. Use `hot_coords` flat copy `globalinfo.c:get_coord_storage_ptr` `memcpy` once per frame vs per-atom dispatch.

### 8.2 Data Layout

`cbmcbendinfo(NTMAX 7,NUMAX 2202,MAXBEND 12,MAXBEND,CBMC_MAXBENDPEAKS 3,CBMCMAXBEND 5)` `preproc.h:33,30,137,41,53` `≈33M doubles ~266 MB` `Source/globaldata.F:twh_cbmcbendinfo` page-fault `45%` flat profile `optimization_plan.md:204`. Actual `nmolty ≤3` often `ntmax` unused. Compress to `real(4)` or alloc `ntmax_used * numax_used` dynamically via `twh_allocateVector` `globalinfo.c:1779` `GLB_ALLOC`. Similarly `nbcoeff(0:55,CROSSTYPEMAX 70500)` `preproc.h:81,184` `30 MB` `vtwobody.F` copies `13,14,15,16,0` flat `5*70k=350k doubles 2.8MB` already minimal.

Make `MAXBOX 3` `preproc.h:27`, `MAXKMAX 12` `preproc.h:68` runtime alloc via `PNT_EWALD_KMAX:10` `MAX_DOUBLE_POINTERS 50` `preproc.h:760` if `kmax<6` typical.

**Files**: `Source/writepdb.F:76,124,183`, `Source/writemovie.F:71,117,134`, `Source/rwconf.F:107,764`, `Source/openfile.F:90`, `Source/globaldata.F:twh_cbmcbendinfo`, `Source/preproc.h:27,30,68`, `Source/assemblesetup.F`, `Source/hot_common.h`

**Est. gain**: `0%` bench (freq 0) / `20-30%` I/O phase `>10K` when enabled; layout saves `~200 MB` resident, improves L3 hit.

---

## 9. Improvement P13 — Cell / Verlet Neighbor List (Algorithmic, Large-N Only)

### 9.1 Problem

`energy.F:264` `do testmolty 0..nmolty` `do jlist 1..ncmt or nchainmove` `O(N^2)` `ncmt ≈ N/ nunit` `~2000` for `10K atoms`. COM cutoff `energy.F:325` `twh_cmp_gt(rijsq,rcmsq)` filters but still scans all `ncmt`. `twh_linclude` `energy.F:834` `engtotal.F:1064` linear bonded exclusion.

### 9.2 Solution (deferred unless `>50K` inputs appear)

Geometric bin `cell(3) = hmatrix diag / rcut` `twh_minboxlength.F:56` at `initconf` + `update` on volume move `mainloop.F:713` `updatevolmaxdisp`. Store `cell_head/next` `MAXCHAIN` linked list per box. `energy.F:264` scan `27` neighbor cells only. For `50K` atoms `O(N)` vs `O(N^2)` `~25x` pair checks, `twh_mimage` only for neighbor pairs. Precompute `exclude_mask(NTMAX,NUMAX,NUMAX)` bitmask at `assemble` vs `twh_linclude` per pair `energy.F:834`.

Risk: bookkeeping invalidates `lpartial` `energy.F:361` logic; verify `tolerance_diff.py`.

**Files**: `Source/energy.F:264,311,834`, `Source/engtotal.F:456,1064`, `Source/ctrmas.F:494`, `Source/minboxlength.F:56`, `Source/putcartinbox.F:49`, `Source/globaldata.F:chainlist`

**Est. gain**: `0%` `<10K` (replicates faster) / `3-5x` `10K` / `10-20x` `50K` -> only path to `optimization_plan.md:4` `5-10x` large-system.

---

## 10. Projected Speedup

Assumptions: Amdahl `S=1/(1-Σ p_i f_i)`, `p_i` from `optimization_plan.md:43`, `f_i` saving fraction, `HOT_PATH_INLINE` guards, `-ffast-math` `1e-6` tolerance `tolerance_diff.py`.

| Scenario | Baseline `optimization_plan.md:12` `8.8s` `bench_results.json` `7.93s min` | After Phase 5 HPI `1.15x` `6.648->5.799s` `optimization_plan.md:369` | +P6-11 (no cell list) | +P13 cell list `>50K` |
|----------|-----------|-----------|-----------|-----------|
| **Small CBMC** `SPC_E_WATER 600 rigid 3-site` `tests/SPC_E_WATER` | `1.00x` | `1.15x` `5.799s` | **1.35-1.55x vs plain** `4.3-4.9s` `1.18-1.35x vs HPI` `Σ p f≈20%` (RNG `8.8%`+P6 `6%`+P8 `2%`+P11 `2.5%`) | `1.35x` (no win) |
| **Medium flexible** `TraPPE_Pentane` `tests/TraPPE_Pentane 0.554->0.529s` `Wall_Water` | `1.00x` | `1.05-1.11x` | **1.50-1.80x vs plain** `Σ≈28%` | `1.60x` |
| **Large Ewald-heavy** `10K+ atoms kmax 8-12` `VECTORMAX 5184` `CROSSTYPEMAX 70.5k` | `1.00x` | `~1.20x` (recipsum `6->40%`) | **1.80-2.50x vs plain** `1.50-2.10x vs HPI` Ewald P9 `15-18%` + P6 `12%` | `3-5x` |
| **Very large** `50K+` | `1.00x` | `~1.20x` | `2.0x` | **5-10x vs plain** `O(N^2)->O(N)` |

`bench_results.json:19` `7.93s` small `-> 7.93/1.45≈5.47s` projected min. `5-10%` measurement noise `bench_results.json:24` `σ0.05-0.54s` WSL `2x` inflated historical `bench_results.json:13`.

MPI scaling: Phase 3 `control.c:49` `rex.c:115` `MPI_Isend+Waitall` + `MPI_THREAD_FUNNELED` `towhee.c:395` complete `optimization_plan.md:330`; additional `P6-11` multiply single-node efficiency, `64x16` `>80%` `optimization_plan.md:293` unchanged; domain `P13` only needed `>50K`.

---

## 11. Implementation Order & Checklist

**Order (ROI / risk)**:

1.  **P7A** `initialize.F:42` default switch (1 line, `5-8%` small, risk: `tolerance 1e-6` diverge — verify `Examples/* do_test`).
2.  **P6** `energy.F:164` atomnum + coords flat (2 files, `5-12%`, verify `gprof`).
3.  **P10** branchless `mimage.F:199` (enables SIMD).
4.  **P8** coeff caches `vangle.F:126` (small, targeted).
5.  **P9.1** hoist `hinverse` `recipsum.F:116` (`15k->9` calls, trivial).
6.  **P11** LTO `Source/Makefile.am:38` + PGO `configure.ac` (build flag).
7.  **P12** I/O buffered (only if `pdboutfreq` users need).
8.  **P9.2-4** `expon`/`derfc` vectorize (larger Ewald payoff, needs large-system benchmark).
9.  **P13** cell list (algorithmic, high effort, gate on `>50K` inputs).

**Per-improvement checklist**:

*   [ ] Add `#ifdef HOT_PATH_INLINE` guard + `#else` preserve `twh_*` dispatch (baseline `0-diff` `diff -w` `Examples/do_test`).
*   [ ] Update `hot_common.h` + `engtotal_globals.F` binding + `globaldata.F` SAVE `target` if needed.
*   [ ] `cd src/towhee-8.2.3 && ./configure && make -j16 towhee` plain + `make HOT_PATH_CPPFLAGS=-DHOT_PATH_INLINE` both compile.
*   [ ] `gprof` verify `twh_*_` calls drop (`twh_coordinates_` `65M-><5M`, `twh_globalpos` eliminated, `twh_hinverse` `15k->9`).
*   [ ] `python3 tests/comparison_benchmark.py --run-dir tests/SPC_E_WATER --threads 1 --repeats 3 --check-tolerance 1e-6` + `src/towhee-8.2.3/Examples/do_test` `diff -w` strict + `tolerance_diff.py 1e-6`.
*   [ ] Record `original_min`/`optimized_min` `bench_results.json` via `bench_config.json:1` schema.

---

## 12. Validation & Benchmarks

```bash
# Build
cd src/towhee-8.2.3
./configure                # or --enable-mpi
make HOT_PATH_CPPFLAGS=-DHOT_PATH_INLINE -j16
# Profile
rm -f gmon.out && src/towhee-8.2.3/Source/towhee && gprof src/towhee-8.2.3/Source/towhee gmon.out | head -n 30
# Regression strict + tolerance
cd Examples/SPC_E_WATER && ../../do_test; cat results
python3 Examples/tolerance_diff.py answer_current towhee_output 1e-6
# Benchmark
python3 tests/comparison_benchmark.py --original /home/sharir/towhee-8.2.3/Source/towhee --optimized src/towhee-8.2.3/Source/towhee --run-dir tests/SPC_E_WATER --threads 1 --repeats 10
```

**Test matrix** (extends `optimization_plan.md:358`):

| Config | Threads | Tol | Expected |
|--------|---------|-----|----------|
| `tests/SPC_E_WATER` HPI serial P6-8 | 1 | `1e-6` | `0 diff` + `1.30x` |
| `tests/Wall_Water` P6 (HPI) | 1 | `1e-6` | `⚠️ -0.314` CBMC init fixed P8? `1.03x` baseline+Vec passes `optimization_plan.md:364` |
| `tests/TraPPE_Pentane` P7A `DX_1597` | 1 | `1e-6` | `✅ 1.05x` + RNG `1.08x` cumulative `~1.14x` |
| `src/Gibbs/EPM_VLCC` P9 | 1 | `1e-6`/`1e-4` | fast-math `1.11x` -> `1.25x` large |
| `>10K` synthetic `kmax12` | 1 | `1e-6` | `Ewald` `recipsum.F:144` vector `1.15x` |

---

## 13. Open Decisions

| Decision | Status | Notes |
|----------|--------|-------|
| CBMC memoization | Resolved `optimization_plan.md:438` COMMON direct-read | P8 completes remaining sites `getcbangle.F:514` `dihedral_distribution.F:422` |
| System size for cell list | Pending | Gate `P13` on `>50K` shipped Examples `<10K` `optimization_plan.md:334` deferred |
| CBMC `266 MB` budget | Open | `P12` layout: `real(4)` or dynamic `ntmax_used` |
| RNG default | **Propose `DX_1597_2_7`** `preproc.h:671` | Trade `RANLUX` luxury `3` period vs `5-8x` speed; `tolerance` not strict |
| Regression format | Resolved `optimization_plan.md:342` `tolerance_diff.py` | `P7A` RNG changes require `1e-6` not `diff -w` |

---

## 14. File Map (New/Modified)

| File | Change |
|------|--------|
| `Source/hot_common.h` | NEW: `/hot_coords/`, `/hot_atomnum/`, `/hot_bencoeff/`, `/hot_qqatom/`, `/hot_hinverse/` |
| `Source/engtotal_globals.F` | Add `atomnum_arr`, `ben_coeff_cache`, `hinv_arr`, bind `c_coord_*` pointers |
| `Source/energy.F:164,341` | HOT path flat coords + atomnum + hoisted `rcelect/rcutinsq` + guarded `ctrmas` |
| `Source/engtotal.F:256,518` | `inline_atomnum` via `atomnum_arr` |
| `Source/globaldata.F` | `SAVE target` for `twh_hinverse`, `twh_qqatom`, `twh_bencoeff` if COMMON, add `target` |
| `Source/initialize.F:42` | `RNG_DX_1597_2_7` default + optional `twh_random_vec` |
| `Source/vangle.F:126` | Coeff cache loop, `!$OMP DECLARE SIMD` |
| `Source/recipsum.F:116,133,144` | Hoist `hinverse`, `expon` branchless, `SIMD` |
| `Source/vcoulomb.F:141` | `derfc` vector/`DECLARE SIMD` |
| `Source/mimage.F:199` | Branchless `anint` + `DECLARE SIMD` |
| `Source/writepdb.F:76`/`writemovie.F:71`/`rwconf.F:764` | Buffered/binary + flat coords |
| `Source/Makefile.am:38` | `VEC_FFLAGS+=-flto`, `AM_LDFLAGS+=-flto` |
| `Source/configure.ac:95` | `AX_CHECK_COMPILE_FLAG([march])`, `--enable-pgo` |
| `tests/comparison_benchmark.py` | Already `1e-6` `compare_outputs` — reuse |

---

## 15. References

* `Source/preproc.h:68` `VECTORMAX 5184` `MAXKMAX 12` `NNTYPE 375` `CROSSTYPEMAX 70500` `TVIBMAX 8200` `TBENMAX 2700`
* `Source/globaldata.F:16771` `506` `twh_*` accessors `GLB_GET 2002` `preproc.h:360`
* `Source/globalinfo.c:1034` `twh_coordstorage_` `PNT_COORDSTORAGE 7` `preproc.h:716`
* `Source/energy_change.F:2383` differential vs `Source/engtotal.F:335` full `O(N^2)`
* `Source/ranlux.F:980` `dranlux` luxury, `Source/dx-1597-2-7.c:80` alternative
* `openspec/live_specs/optimization_plan.md:25` OpenMP corruption scope
* `bench_results.json:1` `bench_config.json:1` historical `6.648->5.799s 1.15x` native

