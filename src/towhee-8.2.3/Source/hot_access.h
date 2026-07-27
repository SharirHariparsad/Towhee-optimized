C     hot_access.h - local COMMON accessors for hot-path inline code.
C     Uses _hot suffix to avoid naming conflicts with routine-local
C     variables in engtotal.F and energy.F.
C     Included only under #ifdef HOT_PATH_INLINE.
C
C     These COMMON blocks map to the same storage as hot_common.h
C     (which is included in engtotal_globals.F90).  The blocks share
C     storage by linker-level COMMON association.

C     --- NBCOEFF: nonbonded coefficients ---
      double precision nbcoeff_hot
      common /hot_nbcoeff/ nbcoeff_hot(MINNBCOEFF:MAXNBCOEFF
     &     ,CROSSTYPEMAX)

C     --- NTYPE: atom type per molecule type and unit ---
      integer ntype_hot
      common /hot_ntype/ ntype_hot(NTMAX,NUMAX)

C     --- NUNIT: number of units per molecule type ---
      integer nunit_hot
      common /hot_nunit/ nunit_hot(NTMAX)

C     --- LUSE_C: Coulomb flag per molecule type and unit ---
      logical luse_c_hot
      common /hot_luse_c/ luse_c_hot(NTMAX,NUMAX)

C     --- LUSE_V: VDW flag per molecule type and unit ---
      logical luse_v_hot
      common /hot_luse_v/ luse_v_hot(NTMAX,NUMAX)

C     --- QQSCALE: charge scaling per molecule type and unit ---
      double precision qqscale_hot
      common /hot_qqscale/ qqscale_hot(NTMAX,NUMAX)

C     --- LSHIFT: shift potential flag ---
      logical lshift_hot
      common /hot_lshift/ lshift_hot

C     --- LEWALD: Ewald summation flag ---
      logical lewald_hot
      common /hot_lewald/ lewald_hot

C     --- LELECT: electrostatics flag ---
      logical lelect_hot
      common /hot_lelect/ lelect_hot

C     --- RCUTSQ: squared VDW cutoff ---
      double precision rcutsq_hot
      common /hot_rcutsq/ rcutsq_hot

C     --- RMINSQ: squared hard-core inner radius ---
      double precision rminsq_hot
      common /hot_rminsq/ rminsq_hot
