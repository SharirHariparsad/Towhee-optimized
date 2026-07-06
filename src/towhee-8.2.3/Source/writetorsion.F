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
      subroutine twh_writetorsion(index)
c     ******************************************************************
c     * writes all of the torsion parameters used in this simulation to*
c     * the file index                                                 *
c     * index: file index for output                                   *
c     *                                                                *
c     * originally written 02-22-2002 by M.G. Martin                   *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LOFTOR
#define FUNCTION_GET_LTORHERE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer index
c     --- local variables
      logical lfound
      character*(FFNAMELEN) torff
      character*(FFNAMELEN) torstring
      dimension torstring(4)
      integer type,i,itemp,tstyle,ntorloop,icoeff
      double precision radtodeg
      double precision torcoeff
      dimension torcoeff(MINTORCOEFF:MAXTORCOEFF)
      
c     --- retrieve constants
      call twh_constant_radtodeg(GLB_GET,radtodeg)
c     --- output torsion information
      write(index,'(a13)') 'Torsion Types'
      lfound = .false.
      do type = 1,TTORMAX
         if ( twh_get_ltorhere(type) ) then
            do icoeff = MINTORCOEFF,MAXTORCOEFF
               call twh_torcoeff(GLB_GET,type,icoeff,torcoeff(icoeff))
            enddo
            lfound = .true.
            call twh_torstyle(GLB_GET,type,tstyle)
            if ( tstyle .eq. TOR_STYLE_SIM_HARM ) then
               write(index,10) type,torcoeff(1)
     &              ,torcoeff(0)
 10            format('Type: ',i3,1x,'Style: Harmonic  Angle:'
     &              ,f6.4,1x,'k ',f10.1)
            elseif ( tstyle .eq. TOR_STYLE_OPLS_COS ) then
               write(index,20) type,(torcoeff(i),i=1,3)
 20            format('Type: ',i3
     &              ,1x,'Style: OPLS Cosine Series',/
     &              ,10x,'k1:',f10.1,1x,'k2:',f10.1,1x,'k3:'
     &              ,f10.1)
            elseif ( tstyle .eq. TOR_STYLE_GCA_COS ) then
               call twh_ntorloop(GLB_GET,type,ntorloop)
               write(index,30) type,ntorloop
 30            format('Type: ',i3,1x
     &              ,'Style: Amber/Charmm/Gromos Cosine Series'
     &              ,' n: ',i3)
               do i = 1,ntorloop
                  itemp = 3*(i-1)
                  write(index,33) i,torcoeff(itemp+1)
     &                 ,torcoeff(itemp+2)
     &                 ,torcoeff(itemp+3)
               enddo
 33            format(10x,'Index',i3,3x,' k:',f10.1,' m:',f3.1
     &              ,' d:',f10.7)
            elseif ( tstyle .eq. TOR_STYLE_GCA_CUB ) then
               call twh_ntorloop(GLB_GET,type,ntorloop)
               write(index,40) type,ntorloop
 40            format('Type: ',i3,1x
     &              ,'Style: Amber/Charmm/Gromos Cosine Series'
     &              ,' n: ',i3)
               do i = 1,ntorloop
                  itemp = 3*(i-1)
                  write(index,33) i,torcoeff(itemp+1)
     &                 ,torcoeff(itemp+2)
     &                 ,torcoeff(itemp+3)
               enddo
               itemp = 3*ntorloop
               write(index,41) torcoeff(itemp+2)*radtodeg
     &              ,torcoeff(itemp+1)
 41            format('With additional Harmonic   Angle:'
     &              ,f8.4,1x,'k ',f10.1)
            elseif ( tstyle .eq. TOR_STYLE_COMP_X_AD ) then
               do i = 1,4
                  call twh_torstring(GLB_GET,type,i,torstring(i))
               enddo
               write(index,62) type,(torstring(i),i=1,4)
     &              ,(torcoeff(i),i=0,19)
            elseif ( tstyle .eq. TOR_STYLE_COMP_R_AD ) then
               write(index,63) type,(torcoeff(i),i=0,2)
               
            elseif ( tstyle .eq. TOR_STYLE_TRAPPE_COS ) then
               write(index,70) type,(torcoeff(i),i=0,1)
 70            format('Type: ',i3,1x,'Style: TraPPE simple cosine '
     &              ,/,10x,'e0 ',f10.2,'e1 ',f10.2)
            elseif ( tstyle .eq. TOR_STYLE_NB ) then
               call twh_torff(GLB_GET,type,torff)
               write(index,73) type,torff
               
            elseif ( tstyle .eq. TOR_STYLE_COMP_X_EX ) then
c               --- class 2 with explicit terms
               write(index,91) type,(torcoeff(i),i=0,31)
               
            elseif ( tstyle .eq. TOR_STYLE_COS_POW ) then
               call twh_ntorloop(GLB_GET,type,ntorloop)
               write(index,100)type,ntorloop
 100           format('Type: ',i3,1x,'Style: Cosine Power Series '
     &              ,'Sum x=1,Ntor a(x) Cos(phi)^x'
     &              ,/,10x,'Ntor:',1x,i3)
               do i=0,ntorloop
                  write(index,101) i,torcoeff(i)
               enddo
 101           format (10x,'a(',i1,'):',1x,g12.6)
            elseif ( tstyle .eq. TOR_STYLE_O_OPLS_C ) then
               write(index,50) type,(torcoeff(i),i=0,3)
            elseif ( tstyle .eq. TOR_STYLE_SUM_COS ) then
               call twh_ntorloop(GLB_GET,type,ntorloop)
               write(index,'(a5,1x,i3,1x,a28,1x,a2,1x,i3)') 
     &              'Type:',type
     &              ,'Style: Sum2003 Cosine Series','n:',ntorloop
               write(index,'(a34,1x,f6.2)') 
     &              'Torsion shift angle [torcoeff(0)]:'
     &              ,torcoeff(0)*radtodeg
               do i = 0,ntorloop-1
                  itemp = i+1
                  write(index,'(a2,i3,1x,a9,f15.5)')
     &                 'n:',i,'Constant:',torcoeff(itemp)
               enddo
            elseif ( tstyle .eq. TOR_STYLE_O_OPLS_T ) then
               write(index,'(a5,1x,i3,1x,a27)') 'Type:',type
     &              ,'Style: 2 term Cosine Series'
               do i = 1,2
                  write(index,'(a9,i1,a2,f15.5)') 
     &                 'torcoeff(',i,'):',torcoeff(i)
               enddo
            elseif ( tstyle .eq. TOR_STYLE_UFF ) then
               write(index,140) 'Type:',type
     &              ,'Style: UFF K[1 - Cn Cos(n phi)]'
     &              ,'K:',torcoeff(1),'Cn:'
     &              ,torcoeff(2)
     &              ,'n:',torcoeff(3)
 140           format(a5,1x,i4,1x,a31,/,11x,a2,1x,f10.3,1x,a3,f10.3,1x
     &              ,a2,1x,f5.2)
            elseif ( tstyle .eq. TOR_STYLE_DREIDING ) then
               write(index,150) 'Type:',type
     &              ,'Style: K[1 - Cos(n (phi - phi_0))]'
     &              ,'K:',torcoeff(1),'n:'
     &              ,torcoeff(2)
     &              ,'phi_0:',torcoeff(3)
 150           format(a5,1x,i4,1x,a34,/,11x,a2,1x,f10.3,1x,a2,f10.3,1x
     &              ,a6,1x,f10.3)
            elseif ( tstyle .eq. TOR_STYLE_TF_COS ) then
               write(index,160) 'Type:',type
     &              ,'Style: K[1 - Cos(2 phi)]','K:'
     &              ,torcoeff(1)
 160           format(a5,1x,i4,1x,a24,/,10x,a2,1x,f10.3)
            elseif ( tstyle .eq. TOR_STYLE_TRAPPE_PLAN ) then
               write(index,170) type,'c1:',torcoeff(1)
     &              ,'c2:',torcoeff(2),'f:'
     &              ,torcoeff(3)
 170           format('Type:',i4,1x
     &              ,'Style: c1[1 + Cos(phi + f)] + c2[1 - Cos^2(phi)]'
     &              ,/,10x,2(a3,1x,f10.3,1x),a2,1x
     &              ,f10.5)
            elseif ( tstyle .eq. TOR_STYLE_SQUARE_WELL ) then
               write(index,180) type,torcoeff(1)
     &              ,torcoeff(2)
     &              ,torcoeff(3)
 180           format('Type:',i4,1x,'Style: Square Well',/,10x
     &              ,'Low Wall:',1x,f10.3,1x,'High Wall:',1x,f10.3,1x
     &              ,'Well Depth:',f10.3)
            elseif ( tstyle .eq. TOR_STYLE_AMBER ) then
               call twh_ntorloop(GLB_GET,type,ntorloop)
               write(index,190) type,ntorloop
 190           format('Type: ',i3,1x
     &              ,'Style: Amber Cosine Series / Total Torsions'
     &              ,' n: ',i3)
               do i = 1,ntorloop
                  itemp = 3*(i-1)
                  write(index,33) i,torcoeff(itemp+1)
     &                 ,torcoeff(itemp+2)
     &                 ,torcoeff(itemp+3)
               enddo
            elseif ( tstyle .eq. TOR_STYLE_OPLS_FL ) then
               write(index,200) type,(torcoeff(i),i=1,4)
 200           format('Type: ',i3
     &              ,1x,'Style: OPLS 4 Parameter Cosine Series',/
     &              ,10x,'k1:',f10.1,1x,'k2:',f10.1,1x,'k3:'
     &              ,f10.1,1x,'k4:',f10.1)
            elseif ( tstyle .eq. TOR_STYLE_MRD ) then
               call twh_ntorloop(GLB_GET,type,ntorloop)
               write(index,210) type,ntorloop
 210           format('Type:',1x,i3,1x,'Style: Multiple Rigid Torsions'
     &              ,/,10x,'Number of allowed Torsions:',1x,i3)
               do i = 1,ntorloop
                  write(index,211) i,torcoeff(i)
               enddo
 211           format(10x,'Torsion(',i1,'):',1x,f10.5)
            elseif ( tstyle .eq. TOR_STYLE_FL ) then
               write(index,220) type,(torcoeff(i),i=0,5)
 220           format('Type:',1x,i3,1x,'Style: Siepmann Fluoroalkane'
     &              ,/,'a0 + a1(1-cos(phi)) + a2(1-Cos(3 phi))',1x
     &              ,'+ a3(1-Cos(phi))^5 + b0 exp[-b1 phi^2]'
     &              ,/,10x,'a0:',f10.1,1x,'a1:',f10.1,1x,'a2:',f10.1
     &              ,1x,'a3:',f10.1
     &              ,/,10x,'b0:',f10.1,1x,'b1:',f10.1)
            elseif ( tstyle .eq. TOR_STYLE_DRUKKER ) then
               write(index,230) type,(torcoeff(i),i=0,1)
 230           format('Type:',1x,i3,1x,'Style: k1 (1 - Cos(phi - k0))'
     &              ,/,'k0:',1x,g12.6
     &              ,/,'k1:',1x,g12.6)
            else
               write(6,*) 'nondefined torsion style of '
     &              ,tstyle
     &              ,'for torsion number ',type
               stop
            endif
            if ( twh_get_loftor(type) ) then
               write(6,11) torcoeff(-1)
            endif
 11         format(10x,'with 1-4 vdw and scaled (',f5.3,') 1-4 coulomb')
         endif
      enddo
      if ( .not. lfound ) write(index,*) '  No Torsion Types'

c     --- format statements
      
 50   format('Type: ',i3
     &     ,1x,'Style: Old UA OPLS Cosine Series',/
     &     ,10x,'k0:',f10.1,1x,'k1:',f10.1,1x,'k2:'
     &     ,f10.1,1x,'k3:',f10.1)
 62   format('Type: ',i3,1x,'Style: Compass with Class 2 cross terms'
     &     ,/,10x,4(a,1x)
     &     ,/,10x,'Normal Torsion Constants',3f10.2
     &     ,/,10x,'Bond(ij) - Torsion Constants',3f10.2
     &     ,/,10x,'Bond(jk) - Torsion Constants',3f10.2
     &     ,/,10x,'Bond(kl) - Torsion Constants',3f10.2
     &     ,/,10x,'Angle(ijk) - Torsion Constants',3f10.2
     &     ,/,10x,'Angle(jkl) - Torsion Constants',3f10.2
     &     ,/,10x,'Angle(ijk) - Angle(jkl) - Torsion Constant',f10.2
     &     ,/,10x,'Bond(ij) - Bond(kl) Constant',f10.2
     &     )
 63   format('Type: ',i3,1x,'Style: Compass without cross terms'
     &     ,/,10x,'Normal Torsion Constants',3f10.2)
 73   format('Type: ',i3,1x,'Style: ',a10,' Nonbonded terms only')

 91   format('Type: ',i3,1x,'Style: Class 2 with explicit cross terms'
     &     ,/,10x,'Normal Torsion Constants',6f10.2
     &     ,/,10x,'Bond(jk) - Torsion Constants',4f10.2
     &     ,/,10x,'Bond(ij) - Torsion Constants',4f10.2
     &     ,/,10x,'Bond(kl) - Torsion Constants',4f10.2
     &     ,/,10x,'Angle(ijk) - Torsion Constants',4f10.2
     &     ,/,10x,'Angle(jkl) - Torsion Constants',4f10.2
     &     ,/,10x,'Angle(ijk) - Angle(jkl) - Torsion Constants',3f10.2
     &     ,/,10x,'Bond(ij) - Bond(kl) Constants',3f10.2
     &     )

      return
      end
    
