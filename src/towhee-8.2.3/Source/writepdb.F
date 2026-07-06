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
      subroutine twh_writepdb(lfailure,step)
c     ******************************************************************
c     * makes a pdb file of the system                                 *
c     * the pdb file is written to boxone.pdb and boxtwo.pdb  etc      *
c     *                                                                *
c     * subroutine originally written in 2000 by M.G. Martin           *
c     * last modified 03-31-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_DISTANCE
#define FUNCTION_ARCCOS
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer step
c     --- local variables
c     --- logical scalars
      logical lrectbox
c     --- character string scalars
      character*(ELENAMLEN) elename
      character*2 numbers
      character*3 shortname
      character*4 pdbname
      character*14 stepchar
      character*30 filename
c     --- integer scalars
      integer imolty,base,ibox,ichain,iunit,ipeptide,type,ii
      integer nchain,nunit,numboxes,cpnumber,ntype
c     --- integer arrays
      integer counter
      dimension counter(MAXBOX)
c     --- double precision scalars
      double precision xmin,ymin,zmin,sigma,xtemp,ytemp,ztemp
      double precision radtodeg,hm1,hm2,hm3,qqatom
c     --- double precision arrays
      double precision vector_len,vector_ang
      dimension  vector_len(3), vector_ang(3)

c      write(6,*) 'starting writepdb'

c     --- retrieve constants
      call twh_constant_radtodeg(GLB_GET,radtodeg)
      call twh_cp_number(GLB_GET,cpnumber)
      call twh_numboxes(GLB_GET,numboxes)

      base = 80
c     --- turn step into a character string
      write(stepchar,'(i14.14)') step
      filename = 'box_xx_step_'//stepchar//'.pdb'
      do ibox = 1,numboxes
c        --- open the output files
         write(numbers,'(i2.2)') ibox
         filename(5:6) = numbers
         call twh_openfile(lfailure,80+ibox,filename,'DEFAULT','UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) return
         counter(ibox) = 0
      enddo

c     --- output CRYST1 records to pdb files for non-rectangular boxes
      do ibox = 1,numboxes
         call twh_lrectbox(GLB_GET,ibox,lrectbox)
        if ( .not. lrectbox ) then
          do ii = 1,3
             call twh_hmatrix(GLB_GET,BXD_REAL,ibox,ii,1,hm1)
             call twh_hmatrix(GLB_GET,BXD_REAL,ibox,ii,2,hm2)
             call twh_hmatrix(GLB_GET,BXD_REAL,ibox,ii,3,hm3)
            vector_len(ii) = twh_distance(hm1,hm2,hm3)
            vector_ang(ii) = 0.0d0
          enddo
          do ii = 1,3
             call twh_hmatrix(GLB_GET,BXD_REAL,ibox,1,ii,hm1)
             call twh_hmatrix(GLB_GET,BXD_REAL,ibox,2,ii,hm2)
             call twh_hmatrix(GLB_GET,BXD_REAL,ibox,3,ii,hm3)
            vector_ang(1) = vector_ang(1) + hm2*hm3
            vector_ang(2) = vector_ang(2) + hm1*hm3
            vector_ang(3) = vector_ang(3) + hm1*hm2
          enddo
          vector_ang(1) = vector_ang(1) / (vector_len(2)*vector_len(3))
          vector_ang(2) = vector_ang(2) / (vector_len(1)*vector_len(3))
          vector_ang(3) = vector_ang(3) / (vector_len(1)*vector_len(2))
          do ii = 1,3
            vector_ang(ii) = twh_arccos(vector_ang(ii)) * radtodeg
          enddo
          write(base+ibox,'(A6,3f9.3,3f7.2,A11,i4)')
     &              'CRYST1'
     &              ,vector_len(1), vector_len(2), vector_len(3)
     &              ,vector_ang(1), vector_ang(2), vector_ang(3)
     &              ,'P 1', 1
        endif
      enddo

c     --- if we are just playing around with one protein then move
c     --- it so it is in the bottom of the box
      call twh_nchain(GLB_GET,nchain)
      if ( nchain .eq. 1 ) then
         ichain = 1
         call twh_coordinates(GLB_GET,CRD_REAL,ichain,1,xmin,ymin,zmin)
         call twh_moltyp(GLB_GET,ichain,imolty)
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 2,nunit
            call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &           ,xtemp,ytemp,ztemp)
            if ( xtemp .lt. xmin ) xmin = xtemp
            if ( ytemp .lt. ymin ) ymin = ytemp
            if ( ztemp .lt. zmin ) zmin = ztemp
         enddo
      else
c        --- do not adjust the coordinates when outputting
         xmin = 0.0d0
         ymin = 0.0d0 
         zmin = 0.0d0
      endif

      do ichain = 1,nchain
         call twh_nboxi(GLB_GET,ichain,ibox)
         if ( ibox .ne. 0 ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1,nunit
c              --- increment counter
               counter(ibox) = counter(ibox) + 1
c              --- find peptide name
               call twh_aminonum(GLB_GET,imolty,iunit,ipeptide)
               if ( ipeptide .eq. 0 ) then
                  shortname = '   '
               else
                  call twh_aminoshort(GLB_GET,imolty,iunit,shortname)
               endif
               if ( cpnumber .eq. CP_LENNARD_JONES .or.
     &              cpnumber .eq. CP_SCALED_LENNARD_JONES .or.
     &              cpnumber .eq. CP_UFF_12_6 ) then
c                 --- output the sigma
                  call twh_ntype(GLB_GET,imolty,iunit,type)
                  type = twh_getnbtype(type,type)
                  call twh_nbcoeff(GLB_GET,1,type,sigma)
               elseif ( cpnumber .eq. CP_HARD_SPHERE .or. 
     &                 cpnumber .eq. CP_REPULSIVE_SPHERE ) then 
c                 --- output hard sphere radius
                  call twh_ntype(GLB_GET,imolty,iunit,type)
                  type = twh_getnbtype(type,type)
                  call twh_nbcoeff(GLB_GET,1,type,sigma)
               elseif ( cpnumber .eq. CP_SQUARE_WELL .or.
     &                 cpnumber .eq. CP_REPULSIVE_WELL ) then
c                 --- output hard sphere radius
                  call twh_ntype(GLB_GET,imolty,iunit,type)
                  type = twh_getnbtype(type,type)
                  call twh_nbcoeff(GLB_GET,1,type,sigma)
                  sigma = dsqrt(sigma)
               else
c                 --- just put in a dummy factor
                  sigma = 1.0d0
               endif
c              --- write to the appropriate file (base + box)
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               call twh_pdbname(GLB_GET,imolty,iunit,pdbname)
               call twh_ntype(GLB_GET,imolty,iunit,ntype)
               call twh_elementname(GLB_GET,ntype,elename)
               call twh_qqatom(GLB_GET,imolty,iunit,qqatom)
               write(base+ibox,12)
     &              'ATOM  ',counter(ibox),pdbname,shortname,ipeptide
     &              ,xtemp-xmin,ytemp-ymin,ztemp-zmin
     &              ,sigma,1.0,elename,qqatom
            enddo
         endif
      enddo
 12   format(A6,I5,1x,A4,1x,A3,2x,i4,4x,3f8.3,2f6.2,4x,5x,A2,1x,f7.4)
      do ibox = 1,numboxes
         write(base+ibox,*) 'TER'
         write(base+ibox,*) 'END'
         close(base+ibox)
         rewind(base+ibox)
      enddo
 
c      write(6,*) 'done with writepdb'

      return
      end
    
