#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000-2020 Marcus G. Martin                       *
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
      function twh_vbond(lfailure,lreport,bxdflag,simbox,vibtype,length
     &     ,imolty,iunit,junit,lovrlap)
c     ******************************************************************
c     * Computes the bond energy for a given type and length           *
c     * lreport: if true we are storing the componenent energies for   *
c     * a full report (generally used at the start and end)            *
c     * simbox: the simulation box number (required for some vibtypes) *
c     * vibtype: the forcefield number of the bond vibration           *
c     * length: the distance between the two atoms                     *
c     * imolty: the molecule type - used for nonbonded vib style       *
c     * iunit: the first unit number - used for nonbonded vib style    *
c     * junit: the second unit number - used for nonboned vib style    *
c     * lovrlap: returned as true if there is a hard overlap in the    *
c     *          bond interaction - this occurs with nonbonded vib     *
c     *          style, or if a fixed bond length is > 1% from its mean*
c     *                                                                *
c     * originally written 1-24-2000 M.G. Martin                       *
c     * last modified 03-30-2020 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_EXPON
#define FUNCTION_VBOND
#include "functions.h"
c     --- variables passed to/from function
      logical lfailure,lreport,lovrlap
      integer bxdflag,simbox,vibtype,imolty,iunit,junit
      double precision length
c     --- local variables
      logical lnonbond
      integer style,icoeff
c     --- double precision scalars
      double precision dist,distsq,max_bond_length,onethird
      double precision vvdw,vcharge
c     --- double precisionarrays
      double precision vibcoeff
      dimension vibcoeff(MINVIBCOEFF:MAXVIBCOEFF)
      
c     --- retrieve constants
      call twh_vibstyle(GLB_GET,vibtype,style)
      call twh_max_bond_length(GLB_GET,max_bond_length)
c     --- initialize nonbonded logical
      lnonbond = .false.
c     --- for all bond lengths, check against maximum bond length
      if ( length .gt. max_bond_length ) then
         twh_vbond=0.0d0
         lovrlap = .true.
         return
      endif
c     --- retrieve all vibcoeff values for this type
      do icoeff = MINVIBCOEFF,MAXVIBCOEFF
         call twh_vibcoeff(GLB_GET,vibtype,icoeff,vibcoeff(icoeff))
      enddo
c     --- for valid bond lengths, check the style to determine the
c     --- bond energy
      if ( style .eq. BOND_STYLE_FIXED ) then
c        --- fixed bond length
c        --- zero bond stretch energy unless outside tolerance
         dist = dabs( (length-vibcoeff(0))
     &        /vibcoeff(0)) 
         if ( twh_cmp_gt(dist, 0.01d0) ) then
c           --- this fixed bond length differs too much from its 
c           --- target value
            lovrlap = .true.
         endif
         twh_vbond = 0.0d0
      elseif (style .eq. BOND_STYLE_HARMONIC ) then
c        --- standard harmonic
         twh_vbond = vibcoeff(1) 
     &        * (length - vibcoeff(0))**2
      elseif ( style .eq. BOND_STYLE_GROMOS ) then
c        --- GROMOS bond stretch potential
         twh_vbond = vibcoeff(1)*(length*length 
     &        - vibcoeff(0)
     &        *vibcoeff(0) )**2

      elseif ( style .eq. BOND_STYLE_NONLINEAR ) then
c        --- nonlinear (same as LAMMPS 2001 implementation)
c        --- LAMMPS references Rector, van Swol, Henderson
c        --- Molecular Physics 82 1009 (1994)
c        --- U = vibcoeff(1) * (r-r0)^2 / ( vibcoeff(2) - (r-r0)^2 )
         dist = (length - vibcoeff(0))**2
         twh_vbond = vibcoeff(1)*dist 
     &        / (vibcoeff(2)-dist)

      elseif ( style .eq. BOND_STYLE_MMTWO ) then
c        --- MM2 style
c        ---- this bond potential is not well formed and has a 
c        --- sensible minimum when the length = vibcoeff(0), but then
c        --- has a nonsense maximum when the length = vibcoeff + 1/3
c        --- and then decreases in energy towards -infinity afterwards
c        --- In order to make this sensible enforce a +infinity energy
c        --- beyond that distance
         dist = length - vibcoeff(0)
         call twh_constant_onethird(GLB_GET,onethird)
         if ( dist .gt. onethird ) then
c           --- infinite energy
            twh_vbond=0.0d0
            lovrlap = .true.
         else
            twh_vbond = vibcoeff(1)*dist*dist
     &           *(1.0d0 - 2.0d0*dist)
         endif
      elseif ( style .eq. BOND_STYLE_COMPASS ) then
c        --- Compass quartic potential
c        --- also used by MMFF94
         dist = length - vibcoeff(0)
         distsq = dist*dist
         twh_vbond = vibcoeff(1)*(distsq)
     &        + vibcoeff(2)*(dist*distsq)
     &        + vibcoeff(3)*(distsq*distsq)
      elseif ( style .eq. BOND_STYLE_NONBOND ) then
c        --- nonbonded interactions used as a bonding potential
c        --- need to initialize bond energy and set logical
         twh_vbond = 0.0d0
         lnonbond = .true.
      elseif ( style .eq. BOND_STYLE_NONE ) then
c         --- no interaction
         twh_vbond = 0.0d0
      elseif ( style .eq. BOND_STYLE_MORSE ) then
c         --- Morse potential
         dist = vibcoeff(2)
     &        *(length - vibcoeff(0))
         twh_vbond = vibcoeff(1)
     &        *((twh_expon(dist) - 1.0d0)**2)
      elseif ( style .eq. BOND_STYLE_SQUARE_WELL ) then
c         --- infinite square well
         if ( twh_cmp_lt(length, vibcoeff(1)) .or.
     &        twh_cmp_gt(length, vibcoeff(2)) ) then
c            --- outside of the square well so overlap
            twh_vbond = 0.0d0
            lovrlap = .true.
         else
c           ---- inside the square well so assign energy
            twh_vbond = vibcoeff(3)
         endif
      elseif ( style .eq. BOND_STYLE_HARMONIC_NONBOND ) then
c        --- standard harmonic
         twh_vbond = vibcoeff(1) 
     &        * (length - vibcoeff(0))**2
c        --- harmonic plus nonbond
         lnonbond = .true.
      elseif ( style .eq. BOND_STYLE_FENE ) then
c        --- Frank Heilmann, BASF AG
c        --- FENE bond potential for CG models
c        --- similar to to P. Virnau's thesis:
c        --- V_fene = - 33.75 * eps * ln( 1 - ( r/(1.5*sig) )**2 )
         if ( length .le. vibcoeff(1) ) then
c           --- distance based FENE term
c           --- vibcoeff(3) = -0.5 * vibcoeff(2) * vibcoeff(1)^2
            twh_vbond = vibcoeff(3)
     &           *dlog(1.0d0 -(length/(vibcoeff(1)))**2)
c           --- nonbond portion
            lnonbond = .true.
         else
c           --- infinite energy for overextended bond
            twh_vbond = 0.0d0
            lovrlap = .true.
         endif

      else
         write(6,*) 'vibration style ',style,' is unknown.',vibtype
         twh_vbond = 0.0d0
         lfailure = .true.
         return
      endif

      if ( lreport ) then
c        --- store the regular bond energy
         call twh_report(GLB_INCR,REP_BON_REG,twh_vbond)
      endif


      if ( lnonbond ) then
         call twh_vnonbondinbond(lfailure,lovrlap,bxdflag,simbox
     &        ,imolty,iunit,junit,length,vvdw,vcharge)
         twh_vbond = twh_vbond + vvdw + vcharge
         if ( lfailure .or. lovrlap ) return
         if ( lreport ) then
            call twh_report(GLB_INCR,REP_BON_VDW,vvdw)
            call twh_report(GLB_INCR,REP_BON_CHG,vcharge)
         endif
      endif


c      write(6,*) 'iunit,junit,length,vbond',iunit,junit
c     &     ,length,twh_vbond/CON_KCALTOK

      return
      end

      subroutine twh_vnonbondinbond(lfailure,lovrlap,bxdflag,simbox
     &     ,imolty,iunit,junit,length,vvdw,vcharge)
c     ******************************************************************
c     * computes nonbond terms in the bond energy                      *
c     *                                                                *
c     * split from vbond 07-18-2014 by M.G. Martin                     *
c     * last modified 04-02-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lovrlap
      integer bxdflag,simbox,imolty,iunit,junit
      double precision length,vvdw,vcharge
c     --- local variables
      integer ntii,ntjj
      double precision distsq,dv
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c     --- compute the nonbonded energy between these atoms
      call twh_ntype(GLB_GET,imolty,iunit,ntii)
      call twh_ntype(GLB_GET,imolty,junit,ntjj)
      distsq = length*length
c     --- nonbond pair potential
      vvdw = twh_vtwobody(lfailure,distsq,ntii,ntjj,0,lovrlap,vf,dv)
      if ( lfailure .or. lovrlap ) return

      call twh_vcoulomb(0,bxdflag,simbox,length,imolty,iunit,imolty
     &     ,junit,vf,dv,lovrlap,vcharge)

      return
      end
