C     hot_common.h - COMMON block declarations for hot-path variable
C     sharing between globaldata.F and engtotal_globals.F.
C     Included in both files to guarantee identical layout.
C     TARGET attributes added so module POINTER aliases can bind to
C     COMMON storage.
C
C     The globaldata.F subroutines place their original variable name
C     in these COMMON blocks (with TARGET).  This means NO changes are
C     needed to the body of each subroutine beyond the declaration area.

C     --- NBCOEFF: nonbonded coefficients (MINNBCOEFF:MAXNBCOEFF, CROSSTYPEMAX) ---
      double precision nbcoeff
      common /hot_nbcoeff/ nbcoeff(MINNBCOEFF:MAXNBCOEFF
     &     ,CROSSTYPEMAX)
      target nbcoeff

C     --- NTYPE: atom type per molecule type and unit ---
      integer ntype
      common /hot_ntype/ ntype(NTMAX,NUMAX)
      target ntype

C     --- NUNIT: number of units per molecule type ---
      integer nunit
      common /hot_nunit/ nunit(NTMAX)
      target nunit

C     --- LUSE_C: Coulomb flag per molecule type and unit ---
      logical luse_c
      common /hot_luse_c/ luse_c(NTMAX,NUMAX)
      target luse_c

C     --- LUSE_V: VDW flag per molecule type and unit ---
      logical luse_v
      common /hot_luse_v/ luse_v(NTMAX,NUMAX)
      target luse_v

C     --- QQSCALE: charge scaling per molecule type and unit ---
      double precision qqscale
      common /hot_qqscale/ qqscale(NTMAX,NUMAX)
      target qqscale

C     --- RCUT: VDW cutoff ---
      double precision rcut
      common /hot_rcut/ rcut
      target rcut

C     --- RCUTSQ: squared VDW cutoff ---
      double precision rcutsq
      common /hot_rcutsq/ rcutsq
      target rcutsq

C     --- LSHIFT: shift potential flag ---
      logical lshift
      common /hot_lshift/ lshift
      target lshift

C     --- LEWALD: Ewald summation flag ---
      logical lewald
      common /hot_lewald/ lewald
      target lewald

C     --- LELECT: electrostatics flag ---
      logical lelect
      common /hot_lelect/ lelect
      target lelect

C     --- RCELECT: electrostatic cutoff ---
      double precision rcelect
      common /hot_rcelect/ rcelect(MAXBXD,MAXBOX)
      target rcelect

C     --- RMINSQ: squared hard-core inner radius ---
      double precision rminsq
      common /hot_rminsq/ rminsq
      target rminsq

C     --- EWALD_CALP: Ewald real-space cutoff parameter calp=kalp/boxlen ---
       double precision ewald_calp
       common /hot_ewald_calp/ ewald_calp(MAXBXD,0:MAXBOX)
       target ewald_calp

C     --- HINVERSE: inverse box matrix (Phase 6 P9.1) ---
       double precision hinverse
       common /hot_hinverse/ hinverse(MAXBXD,MAXBOX,3,3)
       target hinverse

C     --- RCUTINSQ: inner cutoff squared (Phase 6 P6) ---
       double precision rcutinsq
       common /hot_rcutinsq/ rcutinsq
       target rcutinsq

C     --- LUSECOM: COM cutoff flag (Phase 6 P6) ---
       logical lusecom
       common /hot_lusecom/ lusecom
       target lusecom
