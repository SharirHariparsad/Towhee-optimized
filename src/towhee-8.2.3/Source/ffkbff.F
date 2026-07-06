#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2021-2021 Marcus G. Martin                       *
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
c     ******************************************************************
      subroutine twh_ffkbff(lfailure,ntypes,nbondtypes,nangletypes
     &     ,nimprotypes,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the KBFF         *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 07-19-2021 by M.G. Martin                   *
c     * last modified 07-22-2021 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,nimprotypes,nbitypes
c     --- local variables
      integer ii,type,ivalue
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal,testmass
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
      character*(FFNAMELEN) ninlocal 
      integer isame,nislocal
      dimension ninlocal(MAXSAMEIMP,MAXIMPNAME)
      double precision onepi
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)
      double precision qbivlocal
      character*(FFNAMELEN) qbinlocal 
      dimension qbinlocal(2)
      
      ffname = 'KBFF'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Geometric')

c     ******************************************************************
c     * KBFF: designed as an acetone forcefield to be used with SPC-E  *
c     * water.                                                         *
c     ******************************************************************
c     --- S. Weerasinge, P. E. Smith
c     --- "Kirkwood-Buff derived force field for mixtures of acetone
c     --- and water",
c     --- J. Chem. Phys. 118, 10663-10670 (2003).
c     --- native units are sigma in nm and epsilon in kJ/mol
      
      type = 0
c     (C) carbonyl carbon in acetone
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.336d0)*10.0d0
      nbclocal(2) = (0.330)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      elelocal = ' C'
      call twh_getelement(lfailure,'getmass',testmass,elelocal,ivalue)
      masslocal = testmass
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (O) carbonyl oxygen in acetone
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.310d0)*10.0d0
      nbclocal(2) = (0.560)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      elelocal = ' O'
      call twh_getelement(lfailure,'getmass',testmass,elelocal,ivalue)
      if ( lfailure ) return
      masslocal = testmass
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CH3) methyl group in acetone
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.3748d0)*10.0d0
      nbclocal(2) = (0.8672)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 0.0d0
c     --- carbon and 3 hydrogen
      elelocal = ' C'
      call twh_getelement(lfailure,'getmass',testmass,elelocal,ivalue)
      if ( lfailure ) return
      masslocal = masslocal + testmass
      elelocal = ' H'
      call twh_getelement(lfailure,'getmass',testmass,elelocal,ivalue)
      if ( lfailure ) return
      masslocal = masslocal + 3.0d0*testmass
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
c     --- use the carbon for the methyl group
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Fixed bond length                                              *
c     ******************************************************************

      type =  type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = (0.1222d0)*10.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'O'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type =  type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = (0.1507d0)*10.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CH3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * harmonic angles                                                *
c     * paper uses V = 1/2 k (theta - theta_0)^2                       *
c     * paper units kJ/mol/rad^2 and degrees
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.44d0
      bclocal(1) = (0.5d0)*(730.0d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CH3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.12d0
      bclocal(1) = (0.5d0)*(670.0d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CH3'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CH3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     --- not completely clear which improper torsion version is
c     --- used in the paper.  Since it is enforcing planar, I
c     --- believe it is the charmm version.
      
      call twh_constant_pi(GLB_GET,onepi)
c     --- initialize type count
      type = 0

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.5d0)*(167.36d0)*(CON_KJTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- C-O-CH3-CH3
      isame = isame + 1
      ninlocal(isame,1) = 'C'
      ninlocal(isame,2) = 'O'
      ninlocal(isame,3) = 'CH3'
      ninlocal(isame,4) = 'CH3'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return

c     --- bond increments inferred from the paper
      
      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.565d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     O - C(=O)  acetone
      qbinlocal(1) = 'O'
      qbinlocal(2) = 'C'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C - C(=O)  acetone
      qbinlocal(1) = 'CH3'
      qbinlocal(2) = 'C'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return
      
      return
      end
