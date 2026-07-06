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
      subroutine twh_writeangle(index)
c     ******************************************************************
c     * writes all of the angle parameters used in this simulation to  *
c     * the file index                                                 *
c     * index: file index for output                                   *
c     *                                                                *
c     * originally written 02-22-2002 by M.G. Martin                   *
c     * last modified 03-30-2020 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LBENDHERE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer index
c     --- local variables
      logical lfound,lbonang
c     --- integer scalars
      integer type,style,icoeff
c     --- character string scalars
      character*(FFNAMELEN) namea,nameb
      character*15 order
c     --- double precision scalars
      double precision radtodeg
c     --- double precision arrays
      double precision bencoeff
      dimension bencoeff(MINBENDCOEFF:MAXBENDCOEFF)
c     --- retrieve constant
      call twh_constant_radtodeg(GLB_GET,radtodeg)

      write(index,'(a11)') 'Angle Types'
      lfound = .false.
      do type = 1,TBENMAX
         if ( twh_get_lbendhere(type) )  then
c           --- retrieve all bencoeffs
            do icoeff = MINBENDCOEFF,MAXBENDCOEFF
               call twh_bencoeff(GLB_GET,type,icoeff,bencoeff(icoeff))
            enddo
            lfound = .true.
            call twh_anglestyle(GLB_GET,type,style)
            if ( style .eq. ANG_STYLE_RIGID ) then
               write(index,10) type,radtodeg*bencoeff(0)
 10            format('Type: ',i3,1x,'Style: Fixed Angle  Angle: '
     &              ,f10.3)
            elseif ( style .eq. ANG_STYLE_HARMONIC ) then
               write(index,11) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1)
 11            format('Type: ',i3,1x,'Style: Standard Harmonic'
     &              ,'  Angle: ',f10.3,1x,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_DREIDING ) then
               write(index,21) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1)
 21            format('Type: ',i3,1x,'Style: 1 - Cosine'
     &              ,'  Angle: ',f10.3,1x,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_HARM_COS ) then
               write(index,31) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1)
 31            format('Type: ',i3,1x,'Style: Harmonic Cosine Angle: '
     &              ,f10.3,1x,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_COMPASS_AUTO ) then
               call twh_bendnames(GLB_GET,type,1,1,namea)
               call twh_bendnames(GLB_GET,type,1,3,nameb)
               write(index,41) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1),bencoeff(2)
     &              ,bencoeff(3),namea
     &              ,bencoeff(4),nameb
     &              ,bencoeff(5),bencoeff(6)
 41            format('Type: ',i3,1x,'Style: Compass  Angle: ',f10.3
     &              ,/,10x,'Constants: ',3f10.1,/,11x,a,f10.3,2x,a,f10.3
     &              ,/,10x,'Bond-Bond Constant ',f10.3)
            elseif ( style .eq. ANG_STYLE_CHARMM_UB ) then
               write(index,51) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1),bencoeff(2)
     &              ,bencoeff(3)
 51            format('Type: ',i3,1x,'Style: Charmm   Angle: ',f10.3
     &              ,1x,'Constant: ',f10.1,/,10x,'UB Length: ',f6.4
     &              ,1x,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_NONBOND ) then
               write(index,61) type
 61            format('Type: ',i3,1x,'Style: Nonbonded Interaction')
            elseif ( style .eq. ANG_STYLE_NB_HARM ) then
               write(index,71) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1)
 71            format('Type: ',i3,1x,'Style: Harmonic plus Nonbonded '
     &              ,'Interaction'
     &              ,/,10x,'Angle: ',f10.3,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_COMPASS_EX ) then
               write(index,81) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1),bencoeff(2)
     &              ,bencoeff(3),bencoeff(4)
     &              ,bencoeff(5),bencoeff(6)
     &              ,bencoeff(7),bencoeff(8)
     &              ,bencoeff(9),bencoeff(10)
 81            format('Type: ',i3,1x,'Style: Compass  Angle: ',f10.3
     &              ,/,10x,'Constants: ',3f10.1
     &              ,/,10x,'Bond-Angle Constants '
     &              ,'k1 :',f10.3,' d1:',f10.3,' k2 :',f10.3
     &              ,' d2:',f10.3,/
     &              ,10x,'Bond-Bond Constants k: ',f10.3,' d1:',f10.3
     &              ,' d2:',f10.3)
            elseif ( style .eq. ANG_STYLE_FOURIER ) then
               write(index,91)'Type:',type
     &              ,'Style: UFF One Term Fourier Expansion'
     &              ,'K (1 + m Cos(n theta))'
     &              ,'Angle:'
     &              ,radtodeg*bencoeff(0)
     &              ,'n:',bencoeff(1)
     &              ,'m:',bencoeff(2)
 91            format(a5,1x,i4,1x,a37,1x,a22,/
     &              ,11x,a6,1x,f10.3,1x,a2,1x,f5.2,1x,a2,1x,f4.1)
               call twh_writeanglek(type)
            elseif ( style .eq. ANG_STYLE_TTFOURIER ) then
               write(index,101) 'Type:',type
     &              ,'Style: UFF Three Term Fourier Expansion','Angle:'
     &              ,radtodeg*bencoeff(0)
     &              ,'C0:',bencoeff(1)
     &              ,'C1:',bencoeff(2)
     &              ,'C2:',bencoeff(3)
 101           format(a5,1x,i4,1x,a39,/,11x,a6,1x,f7.3
     &              ,3(1x,a3,1x,f10.3))
               call twh_writeanglek(type)
            elseif ( style .eq. ANG_STYLE_NONE ) then
               write(index,111) type
 111            format('Type: ',i3,1x,'Style: No Interaction')
            elseif ( style .eq. ANG_STYLE_SIXTH_AUTO ) then
               call twh_bendnames(GLB_GET,type,1,1,namea)
               call twh_bendnames(GLB_GET,type,1,3,nameb)
               call twh_lbonang(GLB_GET,type,lbonang)
               write(index,121) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1),bencoeff(2)
     &              ,lbonang
     &              ,namea,nameb
     &              ,bencoeff(3)
 121           format('Type: ',i3,1x,'Style: MM2  Angle: ',f10.3,/
     &              ,10x,'Constants: ',2f10.3,/
     &              ,10x,'Angle-Bond Logical',1x,l1,1x
     &              ,'Atom A:',1x,a10,1x,'Atom B:',1x,a10,1x,/
     &              ,10x,'Angle-Bond Constant:',f10.3)
            elseif ( style .eq. ANG_STYLE_SQUARE_WELL ) then
               write(index,131) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1),bencoeff(2)
     &              ,bencoeff(3)
 131           format('Type:',1x,i3,1x,'Style: Infinite Square Well'
     &              ,1x,'Angle:',1x,f10.3,/
     &              ,10x,'Minimum Distance:',1x,f6.3
     &              ,1x,'Maximum Distance:',1x,f6.3,1x,'Well Depth:'
     &              ,1x,f10.5)
            elseif ( style .eq. ANG_STYLE_MULTI_RIGID ) then
               write(index,141) type,radtodeg*bencoeff(0)
     &              ,radtodeg*bencoeff(1)
     &              ,bencoeff(2)
 141           format('Type:',1x,i3,1x,'Style: Multiple Rigid',/
     &              ,10x,'Angle(1):',1x,f10.3,/
     &              ,10x,'Angle(2):',1x,f10.3,/
     &              ,10x,'Tolerance:',1x,e10.5)
            elseif ( style .eq. ANG_STYLE_MMFF ) then
               write(index,151) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1),bencoeff(2)
 151           format('Type:',1x,i4,1x,'Style: MMFF Cubic',/
     &              ,10x,'Angle:',1x,f10.3,1x,/
     &              ,10x,'Constants:',2(1x,f10.3))
               call twh_lbonang(GLB_GET,type,lbonang)
               write(index,152) lbonang
 152           format(10x,'Bond-Angle Logical',1x,l1)
               if ( lbonang ) then
                  write(index,153) bencoeff(3)
     &                 ,bencoeff(4)
 153              format(10x,'Bond-Angle Constants:',2(1x,f10.2))
               endif
               call twh_angleorder(GLB_GET,type,order)
               write(index,154) order
 154           format(10x,'Angle Order:',1x,a15)
            elseif ( style .eq. ANG_STYLE_HARM_COS_OT ) then
               write(index,161) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1)
 161           format('Type: ',i3,1x
     &              ,'Style: Harmonic Cosine Angle plus 1-3 nonbond: '
     &              ,f10.3,1x,'Constant: ',f10.1)
            elseif ( style .eq. ANG_STYLE_MMFFC_ANG ) then
               write(index,251) type,radtodeg*bencoeff(0)
     &              ,bencoeff(1),bencoeff(2)
 251           format('Type:',1x,i4,1x,'Style: MMFF Cubic Combo Angle'
     &              ,/,10x,'Angle:',1x,f10.3,1x,/
     &              ,10x,'Constants:',2(1x,f10.3))
               call twh_angleorder(GLB_GET,type,order)
               write(index,254) order
 254           format(10x,'Angle Order:',1x,a15)

            elseif ( style .eq. ANG_STYLE_MMFFC_BA ) then
               write(index,351) type
 351           format('Type:',1x,i4,1x
     &              ,'Style: MMFF Cubic Combo Bond-Angle')
               write(index,353) bencoeff(1)
     &              ,bencoeff(2)
 353           format(10x,'Bond-Angle Constants:',2(1x,f10.2))
               call twh_angleorder(GLB_GET,type,order)
               write(index,354) order
 354           format(10x,'Angle Order:',1x,a15)
            else
               write(index,*) 'WRITEANGLE: angle type ',type,
     &              'has an invalid anglestyle',style
               stop
            endif
         endif
      enddo

      if ( .not. lfound ) write(index,*) '  No Angle Types'

      return
      end
    
      subroutine twh_writeanglek(type)
c     *****************************************************************
c     * writes out the energetic parameters for the angles            *
c     *                                                               *
c     * last modified 03-27-2020 by M.G. Martin                       *
c     *****************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer type
c     --- local variables 
      integer imolty,iunit,ibend,nmolty,nunit,itest,ijben2,ijben3
      integer inben
      double precision benprefact
c     --- the prefactor for this depend upon the actual atoms
c     --- run through all of the molecules and print out the
c     --- constants for each instance
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_inben(GLB_GET,imolty,iunit,inben)
            do ibend = 1,inben
               call twh_itben(GLB_GET,imolty,iunit,ibend,itest)
               if ( itest .eq. type ) then
c                 --- output the atoms and prefactor
                  call twh_benprefact(GLB_GET,imolty,iunit,ibend
     &                 ,benprefact)
                  call twh_ijben2(GLB_GET,imolty,iunit,ibend,ijben2)
                  call twh_ijben3(GLB_GET,imolty,iunit,ibend,ijben3)
                  write(6,102) 'Molecule:',imolty,'Atoms:',
     &                 iunit,ijben2,ijben3,'K:',benprefact
 102              format(11x,a9,1x,i3,1x,a6,3(1x,i4),1x,a2,1x
     &                 ,f10.3)
               endif
            enddo
         enddo
      enddo
      return
      end
