#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2020 Marcus G. Martin                       *
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
      subroutine twh_updatevolmaxdisp()
c     ******************************************************************
c     * updated the maximum volume and cell displacements              *
c     *                                                                *
c     * split from old pupdate 10-13-2005 by M.G. Martin               *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- local variables
      character*20 routput
      integer idim,jdim,imove,movevector,atemp,ntemp,ivalue,jvalue
      integer nvmmax,ensemble,ipair,numboxes
      double precision ratvol,ratcell,dtemp,pmcell,tacell,tavol
      double precision pmvol,rmvol

c     --- retrieve constants
      call twh_runoutput(GLB_GET,routput)
      call twh_nvmmax(GLB_GET,nvmmax)
      call twh_ensemble(GLB_GET,ensemble)
      call twh_numboxes(GLB_GET,numboxes)

c     --- perform periodic operations 
      call twh_pmvol(GLB_GET,pmvol)
      if ( pmvol .gt. 0.0d0 ) then
c        --- adjust 3D maximum volume displacements 
         if ( routput .eq. 'full' .or. routput .eq. 'updates' 
     &        ) then
            write(6,*) 'Updating 3D volume maximum displacements'
         endif
c        --- get the target acceptance rate
         do imove = 1,nvmmax
            call twh_bavol(GLB_GET,imove,atemp)
            call twh_bnvol(GLB_GET,imove,ntemp)
            if ( ntemp .ne. 0 ) then
               ratvol = dble(atemp) / dble(ntemp)
               if ( ratvol .eq. 0.0d0 ) then
                  call twh_rmvol(GLB_SCALE,imove,0.1d0)
               else
                  call twh_tavol(GLB_GET,tavol)
                  dtemp = ratvol/tavol
                  call twh_rmvol(GLB_SCALE,imove,dtemp)
               endif
            endif
            if ( routput .eq. 'full' .or. routput .eq. 'updates'
     &           ) then
               call twh_rmvol(GLB_GET,imove,rmvol)
               if ( imove .le. numboxes ) then
                  write(6,60) imove,ntemp,atemp,rmvol
               else
                  ipair = imove - numboxes
                  call twh_pairbox(GLB_GET,ipair,1,ivalue)
                  call twh_pairbox(GLB_GET,ipair,2,jvalue)
                  write(6,61) ivalue,jvalue,ntemp,atemp,rmvol
               endif
            endif
            call twh_acnvol(GLB_INCR,imove,ntemp)
            call twh_acsvol(GLB_INCR,imove,atemp)
         enddo
c        --- reset volume accumulators
         call twh_bnvol(GLB_INIT,nvmmax,0)
         call twh_bavol(GLB_INIT,nvmmax,0)

 60      format(' Box ',i2,' Tries: ',i8,' Accepted: ',i8
     &        ,' Max Disp.: ',e9.3)
 61      format(' Boxes ',i2,' and ',i2,' Tries: ',i8,
     &        ' Accepted: ',i8,' Max Disp.: ',e9.3)
      endif
      call twh_pmcell(GLB_GET,pmcell)
      if ( pmcell .gt. 0.0d0 ) then
c        --- adjust unit cell maximum displacements
         if ( routput .eq. 'full' .or. routput .eq. 'updates' )
     &        write(6,*) 'Updating unit cell maximum displacements'
         do imove = 1,nvmmax
            do idim = 1,3
               do jdim = 1,3
                  movevector = 3*(idim-1)+jdim
                  call twh_bncell(GLB_GET,imove,movevector,ntemp)
                  call twh_bacell(GLB_GET,imove,movevector,atemp)
                  if ( ntemp .ne. 0 ) then
                     ratcell = dble(atemp) / dble(ntemp)
                     if ( ratcell .eq. 0.0d0 ) then
                        call twh_rmcell(GLB_SCALE,imove,idim,jdim,0.1d0)
                     else
                        call twh_tacell(GLB_GET,tacell)
                        dtemp = ratcell/tacell
                        call twh_rmcell(GLB_SCALE,imove,idim,jdim,dtemp)
                     endif
                  endif
                  if ( routput .eq. 'full' 
     &                 .or. routput .eq. 'updates' ) then
c                    --- output the new displacements
                     call twh_rmcell(GLB_GET,imove,idim,jdim,dtemp)
                     if ( imove .le. numboxes ) then
                        write(6,70) imove,idim,jdim,ntemp,atemp,dtemp
                     else
                        ipair = imove - numboxes
                        call twh_pairbox(GLB_GET,ipair,1,ivalue)
                        call twh_pairbox(GLB_GET,ipair,2,jvalue)
                        write(6,71) ivalue,jvalue,idim,jdim,ntemp
     &                       ,atemp,dtemp
                     endif
                  endif
                  call twh_acncell(GLB_INCR,imove,movevector,ntemp)
                  call twh_acscell(GLB_INCR,imove,movevector,atemp)
               enddo
            enddo
         enddo
c        --- reinitialize cell accumulators
         call twh_bncell(GLB_INIT,nvmmax,9,0)
         call twh_bacell(GLB_INIT,nvmmax,9,0)
 70      format(' Box ',i2,' H-matrix(',i1,',',i1,') Tries: ',i12
     &        ,' Accepted: ',i8,' Max Disp.: ',e9.3)
 71      format(' Boxes ',i2,' and ',i2,' H-matrix(',i1,','
     &           ,i1,') Tries: ',i12
     &        ,' Accepted: ',i12,' Max Disp.: ',e9.3)
      endif

      return
      end
