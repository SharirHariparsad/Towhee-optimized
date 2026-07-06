#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2020 Marcus G. Martin                       *
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
      subroutine twh_writeimproper(lfailure,index)
c     ******************************************************************
c     * writes all of the improper torsion parameters used in this     *
c     * simulation to the file index                                   *
c     * index: file index for output                                   *
c     *                                                                *
c     * rewritten from writetorsion 03-21-2002 by M.G. Martin          *
c     * last modified 03-30-2020 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LIMPHERE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer index
c     --- local variables
      logical lfound
      integer type,itemp,istyle,icoeff
      double precision impcoeff
      dimension impcoeff(MINIMPCOEFF:MAXIMPCOEFF)

      write(index,'(a22)') 'Improper Torsion Types'
      lfound = .false.
      do type = 1,TIMPMAX
         do icoeff = MINIMPCOEFF,MAXIMPCOEFF
            call twh_impcoeff(GLB_GET,type,icoeff,impcoeff(icoeff))
         enddo
         if ( twh_get_limphere(type) ) then
            lfound = .true.
            call twh_impstyle(GLB_GET,type,istyle)
            if ( istyle .eq. 1 ) then
               write(index,10) type,impcoeff(1),impcoeff(0)
            elseif ( istyle .eq. 2 ) then
               write(index,20) type
     &              ,(impcoeff(itemp),itemp=1,3)
            elseif ( istyle .eq. 3 ) then
               write(index,30) type
            elseif ( istyle .eq. 4 ) then
               write(index,40) type,impcoeff(1),impcoeff(2),impcoeff(3)
            elseif ( istyle .eq. 5 ) then
               write(index,50) type,impcoeff(1),impcoeff(0)
            elseif ( istyle .eq. 6 ) then
               write(index,60) type,(impcoeff(itemp),itemp=0,3)
            elseif ( istyle .eq. 7 ) then
               write(index,70) type,impcoeff(1)
            elseif ( istyle .eq. 8 ) then
               write(index,80) type,impcoeff(1)
            else
               write(index,*) 'nondefined improper torsion type of '
     &              ,istyle
     &              ,'for torsion number ',type
               lfailure = .true.
               return
            endif
         endif
      enddo
      if ( .not. lfound ) write(index,*) '  No Improper Types'
 10   format('Type: ',i3,1x,'Style: Harmonic  Angle:'
     &     ,f6.4,1x,'k ',f10.1)
 20   format('Type: ',i3
     &     ,1x,'Style: OPLS Cosine Series',/
     &     ,10x,'k1:',f10.1,1x,'k2:',f10.1,1x,'k3:'
     &     ,f10.1)
 30   format('Type: ',i3,1x,'Style: MGM Stereochemistry Enforcer')
 40   format('Type: ',i3,1x,'Style: Amber Improper Torsion  k:'
     &     ,f10.1,' m:',f3.1,' d:',f6.3)
 50   format('Type: ',i3,1x,'Style: Wilson Harmonic  Angle:'
     &     ,f6.4,1x,'k ',f10.1)
 60   format('Type: ',i3,1x,'Style: UFF OOP  C0:',f10.4,1x,'C1:'
     &     ,f10.4,1x,'C2:',f10.4,1x,'K:',f10.4)
 70   format('Type: ',i3,1x,'Style: k[1 - Cos(phi)]   k:',f10.4)
 80   format('Type: ',i3,1x,'Style: Wilson triple Harmonic  k:',f10.4)
      return
      end
    
