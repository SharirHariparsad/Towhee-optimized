#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2020 Marcus G. Martin                       *
c     * see the file license.gpl for the full license information      *
c     *                                                                *
c     * This program is free software; you can redistribute it and/or  *
c     * modify it under the terms of the GNU General Public License    *
c     * as published by the Free Software Foundation; either version 2 *
c     * of the License, or (at your option) any later version.         *
c     *                                                                *
c     * This program is distributed in the hope that it will be useful,*
c     * but WITHOUT ANY WARRANTY; without even the implied warranty of *
c     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
c     * GNU General Public License for more details.                   *
c     *                                                                *
c     * You should have received a copy of the GNU General Public      *
c     * License along with this program; if not, write to the Free     *
c     * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
c     * MA  02111-1307, USA.                                           *
c     *                                                                *
c     * See the file towhee.F for more information about the code      *
c     ******************************************************************
      subroutine twh_ffsum2003(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the Sum2003      *
c     * forcefield that relies heavily on the NERD version 2 forcefield*
c     * NOTE: some of the new style torsions were set up improperly    *
c     * prior to 09-06-2018 (all of the TOR_STYLE_SUM_2003)            *
c     * substantial rewrite of this routine for version 8.1.0          *
c     *                                                                *
c     * originally written 08-20-2004 by M.G. Martin                   *
c     * last modified 03-02-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer ii,type
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
      double precision onepi

      ffname = 'NERD'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      type = 0

c     ******************************************************************
c     * Sum2003                                                        *
c     ******************************************************************
c     This force field supplements the other NERD parameters with some
c     additional parameters for triglycerols
c     It is not really designed to stand on its own, but will usually
c     need to be combined with one of the versions of NERD.
c     -Sum2003------------------------------------------------
c     jpcb2003b: A.K. Sum; M.J. Biddy; J.J. de Pablo "Predictive 
c     Molecular Model for the Thermodynamic and Transport Properties of 
c     Triacylglycerols", J. Phys. Chem. B, 107, 14443-14451 (2003)
c     Several parameters are not included in that paper, yet appear in
c     this file.  I believe they are from personal communication with
c     Mary Biddy, but the documentation of that has been lost to time
c     M.G. Martin 09-10-2018
      
c     ** Carbon (plus bonded Hydrogen) **

c     methyl group in methylacetate (CH3) united-atom bonded to 
c     oxygen
c     Methyl Acetate CH3- [jcp2003b, Table 2]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.830d0
      nbclocal(2) = (0.6347d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3sp3oac'
      cbnlocal = 'Csp3oac'
      canlocal = 'CH3sp3oac'
      ctnlocal = 'CH3sp3oac'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methyl group in methylacetate (CH3) united-atom bonded to 
c     carbon
c     Methyl Acetate -CH3 [jcp2003b Table 2]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.835d0
      nbclocal(2) = (0.6253d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3sp3cac'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     acetal carbon in acetate
c     Methyl Acetate -C- [jpc2003b, Table 2]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.312d0
      nbclocal(2) = (0.7185d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Csp2ac'
      cbnlocal = 'Csp2ac'
      canlocal = 'Csp2ac'
      ctnlocal = 'Csp2ac'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methylene group (CH2) united atom bonded to oxygen in a
c     triglyceride [jpcb2003]
c     This parameter does not appear in the jpcb2003 paper and I no
c     longer recall its source [M.G. Martin, 09-07-2018]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.830d0
      nbclocal(2) = (0.6347d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 14.0268d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2sp3oac'
      cbnlocal = 'Csp3oac'
      canlocal = 'CH2sp3oac'
      ctnlocal = 'CH2sp3oac'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methylene group (CH2) united atom bonded to acetal carbon in
c     triglyceride [jpcb2003]
c     This parameter is not listed the jpcb2003 paper, but its use is
c     implied.
c     I no longer recall the source [M.G. Martin, 09-07-2018]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.835d0
      nbclocal(2) = (0.6253d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 14.0268d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2sp3cac'
      cbnlocal = 'Csp3'
      canlocal = 'CH2sp3'
      ctnlocal = 'CH2sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methine group (CH) united-atom bonded to an oxygen in a 
c     triglyceride
c     This is implied in [jpcb2003], but the parameters are not listed
c     No longer have a documented reference for this [MG Martin 2018]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.830d0
      nbclocal(2) = (0.6347d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0189d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHsp3oac'
      cbnlocal = 'Csp3oac'
      canlocal = 'CHsp3oac'
      ctnlocal = 'CHsp3oac'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Oxygen **

c     oxygen with 2 single bonds in acetate [jpcb2003b, Table 2]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.472d0
      nbclocal(2) = (0.6736d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Osp3'
      cbnlocal = 'Osp3'
      canlocal = 'Osp3'
      ctnlocal = 'Osp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     oxygen with 1 double bond in acetate [jpcb2003b, Table 2]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.840d0
      nbclocal(2) = (0.7922d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Osp2'
      cbnlocal = 'Osp2'
      canlocal = 'Osp2'
      ctnlocal = 'Osp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * NERD bond lengths - simple harmonic                            *
c     ******************************************************************

c     CH3 - O in methyl acetate
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.48d0
      vclocal(1) = (0.5d0)*(2200d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Csp3oac - Osp3 [jpcb2003b, Table 9, C1-O2]
      vnlocal(1,1) = 'Csp3oac'
      vnlocal(1,2) = 'Osp3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     O - C single bond in acetate
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (0.5d0)*(2700d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Osp3 - Csp2ac [jpcb2003b, Table 9, O2-C3]
      vnlocal(1,1) = 'Osp3'
      vnlocal(1,2) = 'Csp2ac'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     O = C double bond in acetate
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.26d0
      vclocal(1) = (0.5d0)*(3000d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Osp2 = Csp2ac [jpcb2003b, Table 9, C3-O4]
      vnlocal(1,1) = 'Osp2'
      vnlocal(1,2) = 'Csp2ac'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     C - CH3 in acetate or C - CH2 in glyceride
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.53d0
      vclocal(1) = (0.5d0)*(2000d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Csp2ac - Csp3 [jpcb2003b, Table 9, C3-C5]
      vnlocal(1,1) = 'Csp2ac'
      vnlocal(1,2) = 'Csp3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     C-C single bond between carbons
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.54d0
      vclocal(1) = (0.5d0)*(802.3d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C(sp3) - C(sp3) [jpcb2003b, Table 9, C5-C6]
      ii = 1
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'Csp3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     Coac - Coac in triacylglycerol [jpcb2003b]
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.60d0
      vclocal(1) = (0.5d0)*(2000d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Csp3oac - Csp3oac [jpcb2003b, Table 9, C1-C1]
      vnlocal(1,1) = 'Csp3oac'
      vnlocal(1,2) = 'Csp3oac'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * NERD angles - simple harmonic                                  *
c     ******************************************************************

c     CH3 - O - C acetal [jpcb2003b]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.0d0
      bclocal(1) = (0.5d0)*(923d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     C1-O2-C3 glycerols [jpcb2003, Table 9, C1-O2-C3]
c     CH2sp3oac - Osp3 - Csp2ac
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3oac'
      bnlocal(ii,2) = 'Osp3'
      bnlocal(ii,3) = 'Csp2ac'
c     CHsp3oac - Osp3 - Csp2ac
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp3oac'
      bnlocal(ii,2) = 'Osp3'
      bnlocal(ii,3) = 'Csp2ac'
c     Csp3oac - Osp3 - Csp2ac
      ii = ii + 1
      bnlocal(ii,1) = 'Csp3oac'
      bnlocal(ii,2) = 'Osp3'
      bnlocal(ii,3) = 'Csp2ac'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     O - C - O acetal [jpcb2003b]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.0d0
      bclocal(1) = (0.5d0)*(1520d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     Osp3 - Csp2ac - Osp2 [jpcb2003b, Table 9, O2-C3-O4]
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'Csp2ac'
      bnlocal(ii,3) = 'Osp2'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     O - C - CH3 acetal [jpcb2003b]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.0d0
      bclocal(1) = (0.5d0)*(1200d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     Osp3 - Csp2ac - CH2sp3 triglycerols [jpcb2003, Table 9, O2-C3-C5]
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'Csp2ac'
      bnlocal(ii,3) = 'CH2sp3'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     O = C - CH3 acetal [jpcb2003b]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.0d0
      bclocal(1) = (0.5d0)*(1520d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     Osp2 - Csp2ac - CH2sp3 gycerols [jpcb2003, Table 9, O4-C3-C5]
      ii = ii + 1
      bnlocal(ii,1) = 'Osp2'
      bnlocal(ii,2) = 'Csp2ac'
      bnlocal(ii,3) = 'CH2sp3'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.0d0
      bclocal(1) = (0.5d0)*(864d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     C3-C5-C6 in glycerols [jpcb2003]
c     Csp2ac - CH2sp3 - CH3sp3
      ii = ii + 1
      bnlocal(ii,1) = 'Csp2ac'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'CH3sp3'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.0d0
      bclocal(1) = (0.5d0)*(1000d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     O2-C1-C1 in glycerols [jpcb2003]
c     Osp3 - CH2sp3oac - CH2sp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'CH2sp3oac'
      bnlocal(ii,3) = 'CH2sp3oac'
c     Osp3 - CHsp3oac - CH2sp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'CHsp3oac'
      bnlocal(ii,3) = 'CH2sp3oac'
c     Osp3 - Csp3oac - CH2sp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'Csp3oac'
      bnlocal(ii,3) = 'CH2sp3oac'
c     Osp3 - CH2sp3oac - CHsp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'CH2sp3oac'
      bnlocal(ii,3) = 'CHsp3oac'
c     Osp3 - CHsp3oac - CHsp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'CHsp3oac'
      bnlocal(ii,3) = 'CHsp3oac'
c     Osp3 - Csp3oac - CHsp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'Csp3oac'
      bnlocal(ii,3) = 'CHsp3oac'
c     Osp3 - CH2sp3oac - Csp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'CH2sp3oac'
      bnlocal(ii,3) = 'Csp3oac'
c     Osp3 - CHsp3oac - Csp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'CHsp3oac'
      bnlocal(ii,3) = 'Csp3oac'
c     Osp3 - Csp3oac - Csp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'Osp3'
      bnlocal(ii,2) = 'Csp3oac'
      bnlocal(ii,3) = 'Csp3oac'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.0d0
      bclocal(1) = (0.5d0)*(1200d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     C1-C1-C1 in glycerols [jpcb2003]
c     CH2sp3oac - CHsp3oac - CH2sp3oac
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3oac'
      bnlocal(ii,2) = 'CHsp3oac'
      bnlocal(ii,3) = 'CH2sp3oac'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * NERD torsions using standard OPLS cosine series                *
c     * k0 + k1*(1 + Cos(phi)) + k2*(1 - Cos(2*phi))                   *
c     * + k3*(1 + Cos(3*phi))                                          *
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)
      
c     --- using the new torsion potential introduced in Sum et al 2003
      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_SUM_COS)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(0) = 0.0d0*(onepi/180.0d0)
      tclocal(1) = 0.5d0*(0.0d0)*CON_KJTOK
      tclocal(2) = 0.5d0*(52.650d0)*CON_KJTOK
      tclocal(3) = 0.5d0*(30.262d0)*CON_KJTOK
      tclocal(4) = 0.5d0*(7.061d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     CH3 - O - C = O [jpcb2003b, Table 1]
c     CH3sp3oac - Osp3 - Csp2ac - Osp2 (explicitly listed)
      ii = ii + 1
      tnlocal(ii,1) = 'CH3sp3oac'
      tnlocal(ii,2) = 'Osp3'
      tnlocal(ii,3) = 'Csp2ac'
      tnlocal(ii,4) = 'Osp2'
c     C1 - O2 - C3 = O4 [jpcb2003b, Table 9]
c     CH2sp3oac - Osp3 - Csp2ac - Osp2 (assumed C1 is any CH**oac)
      ii = ii + 1
      tnlocal(ii,1) = 'CH2sp3oac'
      tnlocal(ii,2) = 'Osp3'
      tnlocal(ii,3) = 'Csp2ac'
      tnlocal(ii,4) = 'Osp2'
c     CH3 - O - C = O [jpcb2003b]
c     CHsp3oac - Osp3 - Csp2ac - Osp2 (assumed C1 is any CH**oac)
      ii = ii + 1
      tnlocal(ii,1) = 'CHsp3oac'
      tnlocal(ii,2) = 'Osp3'
      tnlocal(ii,3) = 'Csp2ac'
      tnlocal(ii,4) = 'Osp2'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_SUM_COS)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(0) = 180.0d0*(onepi/180.0d0)
      tclocal(1) = 0.5d0*(0.0d0)*CON_KJTOK
      tclocal(2) = 0.5d0*(52.074d0)*CON_KJTOK
      tclocal(3) = 0.5d0*(33.679d0)*CON_KJTOK
      tclocal(4) = 0.5d0*(6.493d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     CH3 - O - C - CH3 [jpcb2003b, Table 2]
c     CH3sp3oac - Osp3 - Csp2ac - Cwild (listed, wild implied)
      ii = ii + 1
      tnlocal(ii,1) = 'CH3sp3oac'
      tnlocal(ii,2) = 'Osp3'
      tnlocal(ii,3) = 'Csp2ac'
      tnlocal(ii,4) = 'Cwild'
c     C1-O2-C3-C5 glycerols [jpcb2003, Table 9]
c     CH2sp3oac - Osp3 - Csp2ac - Cwild (assumed)
      ii = ii + 1
      tnlocal(ii,1) = 'CH2sp3oac'
      tnlocal(ii,2) = 'Osp3'
      tnlocal(ii,3) = 'Csp2ac'
      tnlocal(ii,4) = 'Cwild'
c     C1-O2-C3-C5 glycerols [jpcb2003, Table 9]
c     CHsp3oac - Osp3 - Csp2ac - Cwild (assumed)
      ii = ii + 1
      tnlocal(ii,1) = 'CHsp3oac'
      tnlocal(ii,2) = 'Osp3'
      tnlocal(ii,3) = 'Csp2ac'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_SUM_COS)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(0) = 180.0d0*(onepi/180.0d0)
      tclocal(1) = 0.5d0*(0.0d0)*CON_KJTOK
      tclocal(2) = 0.5d0*(3.076d0)*CON_KJTOK
      tclocal(3) = 0.5d0*(1.459d0)*CON_KJTOK
      tclocal(4) = 0.5d0*(2.928d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     O2-C3-C5-C6 [jpcb2003, Table 9]
c     Osp3 - Csp2ac - CH2sp3 - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Osp3'
      tnlocal(ii,2) = 'Csp2ac'
      tnlocal(ii,3) = 'CH2sp3'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_SUM_COS)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(0) = 0.0d0*(onepi/180.0d0)
      tclocal(1) = 0.5d0*(0.0d0)*CON_KJTOK
      tclocal(2) = 0.5d0*(3.076d0)*CON_KJTOK
      tclocal(3) = 0.5d0*(1.459d0)*CON_KJTOK
      tclocal(4) = 0.5d0*(2.928d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     O4-C2-C5-C6 [jpcb2003, Table 9]
c     Osp2 - Csp2ac - CH2sp3 - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Osp2'
      tnlocal(ii,2) = 'Csp2ac'
      tnlocal(ii,3) = 'CH2sp3'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_SUM_COS)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(0) = 180.0d0*(onepi/180.0d0)
      tclocal(1) = 0.5d0*(0.0d0)*CON_KJTOK
      tclocal(2) = 0.5d0*(51.596d0)*CON_KJTOK
      tclocal(3) = 0.5d0*(-44.053d0)*CON_KJTOK
      tclocal(4) = 0.5d0*(32.759d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     C3-C5-C6-C7 [jpcb2003, Table 9]
c     Csp2ac - CH2sp3 - CH2sp3 - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Csp2ac'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'CH2sp3'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_SUM_COS)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(0) = 180.0d0*(onepi/180.0d0)
      tclocal(1) = 0.5d0*(0.0d0)*CON_KJTOK
      tclocal(2) = 0.5d0*(32.302d0)*CON_KJTOK
      tclocal(3) = 0.5d0*(-15.707d0)*CON_KJTOK
      tclocal(4) = 0.5d0*(26.392d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     C1-C1-O2-C3 [jpcb2003, Table 9]
c     CH2sp3oac - CHsp3oac - Osp3 - Csp2ac
      ii = ii + 1
      tnlocal(ii,1) = 'CH2sp3oac'
      tnlocal(ii,2) = 'CHsp3oac'
      tnlocal(ii,3) = 'Osp3'
      tnlocal(ii,4) = 'Csp2ac'
c     CHsp3oac - CH2sp3oac - Osp3 - Csp2ac
      ii = ii + 1
      tnlocal(ii,1) = 'CHsp3oac'
      tnlocal(ii,2) = 'CH2sp3oac'
      tnlocal(ii,3) = 'Osp3'
      tnlocal(ii,4) = 'Csp2ac'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_SUM_COS)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(0) = 180.0d0*(onepi/180.0d0)
      tclocal(1) = 0.5d0*(0.0d0)*CON_KJTOK
      tclocal(2) = 0.5d0*(8.608d0)*CON_KJTOK
      tclocal(3) = 0.5d0*(42.768d0)*CON_KJTOK
      tclocal(4) = 0.5d0*(20.228d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     C1-C1-C1-O2 [jpcb2003, Table 9]
c     CH2sp3oac - CHsp3oac - CH2sp3oac - Osp3
      ii = ii + 1
      tnlocal(ii,1) = 'CH2sp3oac'
      tnlocal(ii,2) = 'CHsp3oac'
      tnlocal(ii,3) = 'CH2sp3oac'
      tnlocal(ii,4) = 'Osp3'
c     O2-C1-C1-O2 [jpcb2003, Table 9]
c     Osp3 - CH2sp3oac - CHsp3oac - Osp3
      ii = ii + 1
      tnlocal(ii,1) = 'Osp3'
      tnlocal(ii,2) = 'CH2sp3oac'
      tnlocal(ii,3) = 'CHsp3oac'
      tnlocal(ii,4) = 'Osp3'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_SUM_COS)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(0) = 180.0d0*(onepi/180.0d0)
      tclocal(1) = 0.5d0*(0.0d0)*CON_KJTOK
      tclocal(2) = 0.5d0*(0.130d0)*CON_KJTOK
      tclocal(3) = 0.5d0*(-1.837d0)*CON_KJTOK
      tclocal(4) = 0.5d0*(-8.687d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     C=C-C-C [jpcb2003, Table 9]
c     CH2(sp2) - CH(sp2) CH2(sp3) - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'CH2sp2'
      tnlocal(ii,2) = 'CHsp2'
      tnlocal(ii,3) = 'CH2sp3'
      tnlocal(ii,4) = 'Cwild'
c     C=C-C-C [jpcb2003, Table 9]
c     CH(sp2) - CH(sp2) CH2(sp3) - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'CHsp2'
      tnlocal(ii,2) = 'CHsp2'
      tnlocal(ii,3) = 'CH2sp3'
      tnlocal(ii,4) = 'Cwild'
c     C=C-C-C [jpcb2003, Table 9]
c     CH2(sp2) - CH(sp2) CH(sp3) - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'CH2sp2'
      tnlocal(ii,2) = 'CHsp2'
      tnlocal(ii,3) = 'CHsp3'
      tnlocal(ii,4) = 'Cwild'
c     C=C-C-C [jpcb2003, Table 9]
c     CH(sp2) - CH(sp2) CH(sp3) - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'CHsp2'
      tnlocal(ii,2) = 'CHsp2'
      tnlocal(ii,3) = 'CHsp3'
      tnlocal(ii,4) = 'Cwild'
c     C=C-C-C [jpcb2003, Table 9]
c     CH2(sp2) - CH(sp2) C(sp3) - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'CH2sp2'
      tnlocal(ii,2) = 'CHsp2'
      tnlocal(ii,3) = 'Csp3'
      tnlocal(ii,4) = 'Cwild'
c     C=C-C-C [jpcb2003, Table 9]
c     CH(sp2) - CH(sp2) C(sp3) - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'CHsp2'
      tnlocal(ii,2) = 'CHsp2'
      tnlocal(ii,3) = 'Csp3'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,TOR_STYLE_SUM_COS)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,4)
      tclocal(0) = 180.0d0*(onepi/180.0d0)
      tclocal(1) = 0.5d0*(0.0d0)*CON_KJTOK
      tclocal(2) = 0.5d0*(47.853d0)*CON_KJTOK
      tclocal(3) = 0.5d0*(264.781d0)*CON_KJTOK
      tclocal(4) = 0.5d0*(15.129d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     C-C=C-C [jpcb2003, Table 9]
c     Cwild - CH(sp2) CH(sp2) - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'CHsp2'
      tnlocal(ii,3) = 'CHsp2'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
