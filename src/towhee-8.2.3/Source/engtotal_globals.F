#include "preproc.h"
C     engtotal_globals.F - Fortran module with pointer aliases
C     for bypassing the twh_* dispatch call chain in the hot path.
C
C     Guard: HOT_PATH_INLINE (compile-time CPP define)
C
C     The module declares COMMON blocks that share storage (by linker
C     association) with the SAVEd arrays inside the matching twh_*
C     accessor subroutines in globaldata.F.  Module POINTER aliases are
C     then bound to those COMMON arrays in engtotal_globals_init(),
C     giving the hot-path code direct, zero-copy array access without
C     any per-access subroutine dispatch overhead.
C
C     ==================================================================
      module engtotal_globals
      use iso_c_binding, only: c_ptr, c_f_pointer
      implicit none
      save

C     --- Parameter aliases for array bounds (kept in sync with
C         preproc.h via CPP; explicit integer,parameter avoids the
C         fixed-form module/COMMON macro-folding quirk) ---
      integer, parameter :: hp_NTMAX   = NTMAX
      integer, parameter :: hp_NUMAX   = NUMAX
      integer, parameter :: hp_CROSST  = CROSSTYPEMAX
      integer, parameter :: hp_MAXBXD  = MAXBXD
      integer, parameter :: hp_MAXBOX  = MAXBOX

C     --- COMMON storage shared with globaldata.F accessor subroutines ---
C         (linker-level COMMON association; TARGET so POINTERs can bind)
      integer, target :: ntype(hp_NTMAX, hp_NUMAX)
      common /hot_ntype/ ntype
      integer, target :: nunit(hp_NTMAX)
      common /hot_nunit/ nunit
      logical, target :: luse_c(hp_NTMAX, hp_NUMAX)
      common /hot_luse_c/ luse_c
      logical, target :: luse_v(hp_NTMAX, hp_NUMAX)
      common /hot_luse_v/ luse_v
      real(8), target :: qqscale(hp_NTMAX, hp_NUMAX)
      common /hot_qqscale/ qqscale
      real(8), target :: rcut
      common /hot_rcut/ rcut
      real(8), target :: rcutsq
      common /hot_rcutsq/ rcutsq
      real(8), target :: rminsq
      common /hot_rminsq/ rminsq
      real(8), target :: rcelect(hp_MAXBXD, hp_MAXBOX)
      common /hot_rcelect/ rcelect
      logical, target :: lshift
      common /hot_lshift/ lshift
      logical, target :: lewald
      common /hot_lewald/ lewald
      logical, target :: lelect
      common /hot_lelect/ lelect

C     --- hmatrix / halfbox (Phase 1.1 zero-copy mimage) ---
      real(8), target :: hmatrix(hp_MAXBXD, hp_MAXBOX, 3, 3)
      common /hot_hmatrix/ hmatrix
      real(8), target :: halfbox(hp_MAXBXD, hp_MAXBOX, 3)
      common /hot_halfbox/ halfbox

C     --- Ewald calp (Phase 5: precomputed, vectorized) ---
       real(8), target :: ewald_calp(hp_MAXBXD, 0:hp_MAXBOX)
       common /hot_ewald_calp/ ewald_calp

C     --- hinverse (Phase 6 P9.1) ---
       real(8), target :: hinverse(hp_MAXBXD, hp_MAXBOX, 3, 3)
       common /hot_hinverse/ hinverse

C     --- rcutinsq / lusecom (Phase 6 P6) ---
       real(8), target :: rcutinsq
       common /hot_rcutinsq/ rcutinsq
       logical, target :: lusecom
       common /hot_lusecom/ lusecom

C     --- POINTER aliases into COMMON storage (hot-path zero-copy) ---
      integer, pointer, contiguous :: ntype_arr(:,:) => null()
      integer, pointer, contiguous :: nunit_arr(:) => null()
      logical, pointer, contiguous :: luse_c_arr(:,:) => null()
      logical, pointer, contiguous :: luse_v_arr(:,:) => null()
      real(8), pointer, contiguous :: qqscale_arr(:,:) => null()
      real(8), pointer :: rcut_val => null()
      real(8), pointer :: rcutsq_val => null()
      real(8), pointer :: rminsq_val => null()
      real(8), pointer, contiguous :: rcelect_arr(:,:) => null()
      logical, pointer :: lshift_val => null()
      logical, pointer :: lewald_val => null()
      logical, pointer :: lelect_val => null()
       real(8), pointer, contiguous :: hmatrix_arr(:,:,:,:) => null()
       real(8), pointer, contiguous :: halfbox_arr(:,:,:) => null()
       real(8), pointer, contiguous :: ewald_calp_arr(:,:) => null()
       real(8), pointer, contiguous :: hinverse_arr(:,:,:,:) => null()
       real(8), pointer :: rcutinsq_val => null()
       logical, pointer :: lusecom_val => null()

C     --- Zero-copy coordinate access (C-allocated arrays) ---
      real(8), pointer, contiguous :: coord_real_arr(:,:) => null()
      real(8), pointer, contiguous :: coord_temp_arr(:,:) => null()

C     --- Nonbonded coefficients: flat contiguous copies (see spec2) ---
C         Stored in COMMON so the general (non-inline) two-body path in
C         vtwobody.F can read them with a matching extern COMMON
C         declaration (Phase 1.2 scalar getter batching).
      real(8), target :: nbcoeff_0(hp_CROSST)
      common /hot_nbcoeff0/ nbcoeff_0
      real(8), target :: nbcoeff_13(hp_CROSST)
      common /hot_nbcoeff13/ nbcoeff_13
      real(8), target :: nbcoeff_14(hp_CROSST)
      common /hot_nbcoeff14/ nbcoeff_14
      real(8), target :: nbcoeff_15(hp_CROSST)
      common /hot_nbcoeff15/ nbcoeff_15
      real(8), target :: nbcoeff_16(hp_CROSST)
      common /hot_nbcoeff16/ nbcoeff_16

      interface
        subroutine get_coord_storage_ptr(crdflag, cptr, natoms) bind(C, name='get_coord_storage_ptr')
          use iso_c_binding, only: c_ptr
          integer, value :: crdflag
          type(c_ptr), intent(out) :: cptr
          integer, intent(out) :: natoms
        end subroutine get_coord_storage_ptr
      end interface

      contains

C     ==================================================================
C     engtotal_globals_init - one-time init of Fortran SAVE aliases
C     and nbcoeff copies. Called from readtowhee.F during startup.
C     Uses twh_nbcoeff accessor to ensure post-mixing values.
C     Also initializes zero-copy coordinate pointers.
C     ==================================================================
      subroutine engtotal_globals_init()
        integer :: ct
        type(c_ptr) :: cptr
        integer :: natoms
        integer :: natoms_f, natoms_max

C       --- Bind POINTER aliases to the COMMON storage shared with
C           globaldata.F accessor subroutines (zero-copy) ---
        ntype_arr   => ntype
        nunit_arr   => nunit
        luse_c_arr  => luse_c
        luse_v_arr  => luse_v
        qqscale_arr => qqscale
        rcut_val    => rcut
        rcutsq_val  => rcutsq
        rminsq_val  => rminsq
        rcelect_arr => rcelect
        lshift_val  => lshift
        lewald_val  => lewald
        lelect_val  => lelect
         hmatrix_arr => hmatrix
         halfbox_arr => halfbox
         ewald_calp_arr => ewald_calp
         hinverse_arr => hinverse
         rcutinsq_val => rcutinsq
         lusecom_val => lusecom

C       --- nbcoeff flat copies via accessor (not direct COMMON read)
C           This guarantees post-mixing-rule values are captured.
        do ct = 1, CROSSTYPEMAX
          call twh_nbcoeff(GLB_GET, 0,  ct, nbcoeff_0(ct))
          call twh_nbcoeff(GLB_GET, 13, ct, nbcoeff_13(ct))
          call twh_nbcoeff(GLB_GET, 14, ct, nbcoeff_14(ct))
          call twh_nbcoeff(GLB_GET, 15, ct, nbcoeff_15(ct))
          call twh_nbcoeff(GLB_GET, 16, ct, nbcoeff_16(ct))
        end do
C       --- Runtime sanity check: nbcoeff should be non-zero after setup
        if (nbcoeff_13(1) == 0.0d0) then
           write(6,*) 'ERROR: engtotal_globals_init - nbcoeff appears uninitialized!'
           stop 'nbcoeff snapshot failed'
        endif

C       --- Initialize zero-copy coordinate pointers (CRD_REAL = 3001)
        call get_coord_storage_ptr(3001, cptr, natoms)
        natoms_f = natoms
        natoms_max = natoms
        call c_f_pointer(cptr, coord_real_arr, shape=[3, natoms_f])

C       --- Initialize zero-copy temporary coordinate pointers (CRD_TEMP = 3002)
        call get_coord_storage_ptr(3002, cptr, natoms)
        call c_f_pointer(cptr, coord_temp_arr, shape=[3, natoms_f])

      end subroutine engtotal_globals_init

C     ==================================================================
C     engtotal_globals_refresh - re-snapshot nbcoeff copies.
C     Call if nbcoeff is modified after init (future-proofing).
C     Coordinate pointers are live and don't need refresh.
C     ==================================================================
      subroutine engtotal_globals_refresh()
        integer :: ct
        do ct = 1, CROSSTYPEMAX
          call twh_nbcoeff(GLB_GET, 0,  ct, nbcoeff_0(ct))
          call twh_nbcoeff(GLB_GET, 13, ct, nbcoeff_13(ct))
          call twh_nbcoeff(GLB_GET, 14, ct, nbcoeff_14(ct))
          call twh_nbcoeff(GLB_GET, 15, ct, nbcoeff_15(ct))
          call twh_nbcoeff(GLB_GET, 16, ct, nbcoeff_16(ct))
        end do
      end subroutine engtotal_globals_refresh

      end module engtotal_globals
