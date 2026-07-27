C     hot_common.h - COMMON block declarations for hot-path variable
C     sharing between globaldata.F and engtotal_globals.F90.
C     Included in both files to guarantee identical layout.
C
C     The globaldata.F subroutines place their original variable name
C     in these COMMON blocks (replacing their SAVE).  This means NO
C     changes are needed to the body of each subroutine.

C     --- NBCOEFF: nonbonded coefficients (MINNBCOEFF:MAXNBCOEFF, CROSSTYPEMAX) ---
      double precision nbcoeff
      common /hot_nbcoeff/ nbcoeff(MINNBCOEFF:MAXNBCOEFF
     &     ,CROSSTYPEMAX)

C     --- NTYPE: atom type per molecule type and unit ---
      integer ntype
      common /hot_ntype/ ntype(NTMAX,NUMAX)

C     --- NUNIT: number of units per molecule type ---
      integer nunit
      common /hot_nunit/ nunit(NTMAX)

C     --- LUSE_C: Coulomb flag per molecule type and unit ---
      logical luse_c
      common /hot_luse_c/ luse_c(NTMAX,NUMAX)

C     --- LUSE_V: VDW flag per molecule type and unit ---
      logical luse_v
      common /hot_luse_v/ luse_v(NTMAX,NUMAX)

C     --- QQSCALE: charge scaling per molecule type and unit ---
      double precision qqscale
      common /hot_qqscale/ qqscale(NTMAX,NUMAX)

C     --- RCUT: VDW cutoff ---
      double precision rcut
      common /hot_rcut/ rcut

C     --- RCUTSQ: squared VDW cutoff ---
      double precision rcutsq
      common /hot_rcutsq/ rcutsq

C     --- LSHIFT: shift potential flag ---
      logical lshift
      common /hot_lshift/ lshift

C     --- LEWALD: Ewald summation flag ---
      logical lewald
      common /hot_lewald/ lewald

C     --- LELECT: electrostatics flag ---
      logical lelect
      common /hot_lelect/ lelect

C     --- RCELECT: electrostatic cutoff ---
      double precision rcelect
      common /hot_rcelect/ rcelect(MAXBXD,MAXBOX)

C     --- RMINSQ: squared hard-core inner radius ---
      double precision rminsq
      common /hot_rminsq/ rminsq
