#include "preproc.h"
C     engtotal_globals.F - Fortran module with pointer aliases
C     for bypassing the twh_* dispatch call chain in the hot path.
C
C     Guard: HOT_PATH_INLINE (compile-time CPP define)
C
C     ==================================================================
      module engtotal_globals
      implicit none
      save

C     --- COMMON block declarations shared with globaldata.F.
C         The TARGET attributes on each variable let the POINTER
C         aliases below bind to COMMON storage with no copies.  This
C         file must be compiled together with globaldata.F so that the
C         COMMON blocks are linked once.
#include "hot_common.h"

C     --- Atom type / flag / cutoff data: whole-array TARGET aliases ---
      integer, pointer, contiguous :: ntype_arr(:,:) => null()
      integer, pointer, contiguous :: nunit_arr(:) => null()
      logical, pointer, contiguous :: luse_c_arr(:,:) => null()
      logical, pointer, contiguous :: luse_v_arr(:,:) => null()
      real(8), pointer, contiguous :: qqscale_arr(:,:) => null()

C     --- Scalars / cutoffs ---
      real(8), pointer :: rcut_val => null()
      real(8), pointer :: rcutsq_val => null()
      logical, pointer :: lshift_val => null()
      logical, pointer :: lewald_val => null()
      logical, pointer :: lelect_val => null()
      real(8), pointer, contiguous :: rcelect_arr(:,:) => null()
      real(8), pointer :: rminsq_val => null()

C     --- Nonbonded coefficients: flat contiguous copies (see spec2) ---
      real(8), allocatable :: nbcoeff_0(:)
      real(8), allocatable :: nbcoeff_13(:)
      real(8), allocatable :: nbcoeff_14(:)
      real(8), allocatable :: nbcoeff_15(:)
      real(8), allocatable :: nbcoeff_16(:)

      contains

C     ==================================================================
C     engtotal_globals_init - one-time init of Fortran SAVE aliases
C     and nbcoeff copies. Called from readtowhee.F during startup.
C     ==================================================================
      subroutine engtotal_globals_init()
        integer :: ct

C       --- Fortran SAVE arrays with TARGET added (see globaldata.F).
        ntype_arr => ntype
        nunit_arr => nunit
        luse_c_arr => luse_c
        luse_v_arr => luse_v
        qqscale_arr => qqscale

C       --- Scalars ---
        rcut_val => rcut
        rcutsq_val => rcutsq
        lshift_val => lshift
        lewald_val => lewald
        lelect_val => lelect
        rcelect_arr => rcelect
        rminsq_val => rminsq

C       --- nbcoeff flat copies (cannot be CONTIGUOUS pointer targets
C           for fixed-first-index slices of column-major 2D arrays) ---
        allocate(nbcoeff_0(CROSSTYPEMAX))
        allocate(nbcoeff_13(CROSSTYPEMAX))
        allocate(nbcoeff_14(CROSSTYPEMAX))
        allocate(nbcoeff_15(CROSSTYPEMAX))
        allocate(nbcoeff_16(CROSSTYPEMAX))
        do ct = 1, CROSSTYPEMAX
          nbcoeff_0(ct)  = nbcoeff(0, ct)
          nbcoeff_13(ct) = nbcoeff(13, ct)
          nbcoeff_14(ct) = nbcoeff(14, ct)
          nbcoeff_15(ct) = nbcoeff(15, ct)
          nbcoeff_16(ct) = nbcoeff(16, ct)
        end do
      end subroutine engtotal_globals_init

      end module engtotal_globals
