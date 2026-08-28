#include "preproc.h"
C     engtotal_globals.F - Fortran module with pointer aliases
C     for bypassing the twh_* dispatch call chain in the hot path.
C
C     Guard: HOT_PATH_INLINE (compile-time CPP define)
C
C     ==================================================================
      module engtotal_globals
      use iso_c_binding, only: c_ptr, c_f_pointer
      implicit none
      save

C     --- Zero-copy coordinate access (C-allocated arrays) ---
      real(8), pointer, contiguous :: coord_real_arr(:,:) => null()
      real(8), pointer, contiguous :: coord_temp_arr(:,:) => null()

C     --- Nonbonded coefficients: flat contiguous copies (see spec2) ---
      real(8), allocatable :: nbcoeff_0(:)
      real(8), allocatable :: nbcoeff_13(:)
      real(8), allocatable :: nbcoeff_14(:)
      real(8), allocatable :: nbcoeff_15(:)
      real(8), allocatable :: nbcoeff_16(:)

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

C       --- nbcoeff flat copies via accessor (not direct COMMON read)
C           This guarantees post-mixing-rule values are captured.
        allocate(nbcoeff_0(CROSSTYPEMAX))
        allocate(nbcoeff_13(CROSSTYPEMAX))
        allocate(nbcoeff_14(CROSSTYPEMAX))
        allocate(nbcoeff_15(CROSSTYPEMAX))
        allocate(nbcoeff_16(CROSSTYPEMAX))
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
