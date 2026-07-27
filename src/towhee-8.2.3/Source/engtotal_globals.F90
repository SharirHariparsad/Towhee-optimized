!     engtotal_globals.F90 - Fortran module with pointer aliases
!     for bypassing the twh_* dispatch call chain in the hot path.
!
!     Guard: HOT_PATH_INLINE (compile-time CPP define)
!
!     ==================================================================
      module engtotal_globals
      use, intrinsic :: iso_c_binding
      implicit none
      save

      interface
        subroutine get_coord_storage_ptr(crdflag, cptr, natoms) bind(c)
          use, intrinsic :: iso_c_binding
          integer(c_int), value :: crdflag
          type(c_ptr), intent(out) :: cptr
          integer(c_int), intent(out) :: natoms
        end subroutine get_coord_storage_ptr
        subroutine get_integer_storage_ptr(ptr_index, cptr, length)
     &    bind(c)
          use, intrinsic :: iso_c_binding
          integer(c_int), value :: ptr_index
          type(c_ptr), intent(out) :: cptr
          integer(c_int), intent(out) :: length
        end subroutine get_integer_storage_ptr
      end interface

!     --- Coordinates: zero-copy views into C-owned storage ---
      real(8), pointer :: crd2d_real(:,:) => null()
      real(8), pointer :: crd2d_temp(:,:) => null()
      real(8), pointer :: crd_x(:) => null()
      real(8), pointer :: crd_y(:) => null()
      real(8), pointer :: crd_z(:) => null()
      real(8), pointer :: crd_temp_x(:) => null()
      real(8), pointer :: crd_temp_y(:) => null()
      real(8), pointer :: crd_temp_z(:) => null()

!     --- Nonbonded coefficients: flat contiguous copies (see spec2) ---
      real(8), allocatable :: nbcoeff_0(:)
      real(8), allocatable :: nbcoeff_13(:)
      real(8), allocatable :: nbcoeff_14(:)
      real(8), allocatable :: nbcoeff_15(:)
      real(8), allocatable :: nbcoeff_16(:)

!     --- Atom type / flag data: whole-array TARGET aliases ---
      integer, pointer, contiguous :: ntype_arr(:,:) => null()
      integer, pointer, contiguous :: nunit_arr(:) => null()
      integer, pointer, contiguous :: moltyp_arr(:) => null()
      logical, pointer, contiguous :: luse_c_arr(:,:) => null()
      logical, pointer, contiguous :: luse_v_arr(:,:) => null()
      real(8), pointer, contiguous :: qqscale_arr(:,:) => null()

!     --- C-based integer arrays (bound via c_f_pointer) ---
      integer, pointer, contiguous :: globalpos_arr(:) => null()
      integer, pointer, contiguous :: nboxi_arr(:) => null()

!     --- Scalars / cutoffs ---
      real(8), pointer :: rcut_val => null()
      real(8), pointer :: rcutsq_val => null()
      logical, pointer :: lshift_val => null()
      logical, pointer :: lewald_val => null()
      logical, pointer :: lelect_val => null()
      real(8), pointer :: rcelect_val => null()
      real(8), pointer :: rminsq_val => null()

      integer :: natoms_cached = 0
      integer :: classical_pot_num_val = 0
      logical :: module_initialized = .false.
      logical :: int_ptrs_initialized = .false.

      contains

!     ==================================================================
!     engtotal_globals_init - one-time init of Fortran SAVE aliases
!     and nbcoeff copies. Called from allocate.F after all allocations.
!     ==================================================================
      subroutine engtotal_globals_init()
        integer :: ct

!       --- Fortran SAVE arrays with TARGET added ---
!       ntype, nunit, luse_c, luse_v, qqscale are declared in globaldata.F
!       with TARGET attribute. These are genuinely contiguous.
        ntype_arr => ntype
        nunit_arr => nunit
        luse_c_arr => luse_c
        luse_v_arr => luse_v
        qqscale_arr => qqscale

!       --- Scalars ---
        rcut_val => rcut
        rcutsq_val => rcutsq
        lshift_val => lshift
        lewald_val => lewald
        lelect_val => lelect
        rcelect_val => rcelect
        rminsq_val => rminsq

!       --- Cache classical potential number for inline getnbtype ---
        call twh_cp_number(GLB_GET, classical_pot_num_val)

!       --- nbcoeff flat copies (cannot be CONTIGUOUS pointer targets
!           for fixed-first-index slices of column-major 2D arrays) ---
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

!     ==================================================================
!     engtotal_sync_coord_ptrs - O(1) rebind of coordinate pointers.
!     Called at the top of engtotal.F and energy.F hot-path entries.
!     Safe/idempotent to call every invocation.
!     ==================================================================
      subroutine engtotal_sync_coord_ptrs()
        type(c_ptr) :: cptr_real, cptr_temp
        integer(c_int) :: n

        call get_coord_storage_ptr(0, cptr_real, n)
        call c_f_pointer(cptr_real, crd2d_real, [3, int(n)])
        crd_x => crd2d_real(1, :)
        crd_y => crd2d_real(2, :)
        crd_z => crd2d_real(3, :)

        call get_coord_storage_ptr(1, cptr_temp, n)
        call c_f_pointer(cptr_temp, crd2d_temp, [3, int(n)])
        crd_temp_x => crd2d_temp(1, :)
        crd_temp_y => crd2d_temp(2, :)
        crd_temp_z => crd2d_temp(3, :)

        natoms_cached = int(n)
      end subroutine engtotal_sync_coord_ptrs

!     ==================================================================
!     engtotal_sync_int_ptrs - bind C-based integer array pointers
!     (globalpos, nboxi, moltyp). Called once during init after
!     nchain is known.
!     ==================================================================
      subroutine engtotal_sync_int_ptrs(nchain)
        use, intrinsic :: iso_c_binding
        integer, intent(in) :: nchain
        type(c_ptr) :: cptr
        integer(c_int) :: length

        call get_integer_storage_ptr(5, cptr, length)
        call c_f_pointer(cptr, globalpos_arr, [int(nchain)])

        call get_integer_storage_ptr(14, cptr, length)
        call c_f_pointer(cptr, nboxi_arr, [int(nchain)])

        call get_integer_storage_ptr(13, cptr, length)
        call c_f_pointer(cptr, moltyp_arr, [int(nchain)])
      end subroutine engtotal_sync_int_ptrs

      end module engtotal_globals
