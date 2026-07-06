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
      subroutine twh_comparestruc(imolty,jmolty,lsame)
c     ******************************************************************
c     * this subroutine checks to see if two molecules have identical  *
c     * molecule templates                                             *
c     *                                                                *
c     * originally written 07-03-2002 by M.G. Martin                   *
c     * last modified 03-28-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lsame
      integer imolty,jmolty,nmaxi,nmaxj
c     --- local variables
      integer unit,ivib,iben,itor,iaa,iunit,junit
      integer itbentest,jtbentest
      integer ijaa0i,ijaa0j,ijaa1i,ijaa1j,ijaa2i,ijaa2j
      integer ijben2i,ijben2j,ijben3i,ijben3j
      integer ijimp2i,ijimp2j,ijimp3i,ijimp3j,ijimp4i,ijimp4j
      integer ijtor2i,ijtor2j,ijtor3i,ijtor3j,ijtor4i,ijtor4j
      integer ijvibi,ijvibj,inaai,inaaj,inbendi,inbendj
      integer inimproperi,inimproperj,intorsioni,intorsionj
      integer invibi,invibj,itaai,itaaj,itimpropi,itimpropj
      integer ittori,ittorj,itvibi,itvibj
c     --- initialize lsame
      lsame = .true.

      call twh_nunit(GLB_GET,imolty,iunit)
      call twh_nunit(GLB_GET,jmolty,junit)
      if ( iunit .ne. junit ) then
         lsame = .false.
         return
      endif

      call twh_nmaxcbmc(GLB_GET,imolty,nmaxi)
      call twh_nmaxcbmc(GLB_GET,jmolty,nmaxj)
      if ( nmaxi .ne. nmaxj ) then
         lsame = .false.
         return
      endif

      do unit = 1,iunit
c        --- vibrations
         call twh_invib(GLB_GET,imolty,unit,invibi)
         call twh_invib(GLB_GET,jmolty,unit,invibj)
         if ( invibi .ne. invibj ) then
            lsame = .false.
            return
         endif
         do ivib = 1,invibi
            call twh_ijvib(GLB_GET,imolty,unit,ivib,ijvibi)
            call twh_ijvib(GLB_GET,jmolty,unit,ivib,ijvibj)
            call twh_itvib(GLB_GET,imolty,unit,ivib,itvibi)
            call twh_itvib(GLB_GET,jmolty,unit,ivib,itvibj)
            if ( ijvibi .ne. ijvibj .or. itvibi .ne. itvibj ) then
               lsame = .false.
               return
            endif
         enddo
c        --- angles
         call twh_inben(GLB_GET,imolty,unit,inbendi)
         call twh_inben(GLB_GET,jmolty,unit,inbendj)
         if ( inbendi .ne. inbendj ) then
            lsame = .false.
            return
         endif
         do iben = 1,inbendi
            call twh_itben(GLB_GET,imolty,unit,iben,itbentest)
            call twh_itben(GLB_GET,jmolty,unit,iben,jtbentest)
            call twh_ijben2(GLB_GET,imolty,unit,iben,ijben2i)
            call twh_ijben2(GLB_GET,jmolty,unit,iben,ijben2j)
            call twh_ijben3(GLB_GET,imolty,unit,iben,ijben3i)
            call twh_ijben3(GLB_GET,jmolty,unit,iben,ijben3j)
            if ( (ijben2i .ne. ijben2j) .or. ( ijben3i .ne. ijben3j )
     &           .or. ( itbentest .ne. jtbentest ) ) then
               lsame = .false.
               return
            endif
         enddo
c        --- torsions
         call twh_intor(GLB_GET,imolty,unit,intorsioni)
         call twh_intor(GLB_GET,jmolty,unit,intorsionj)
         if ( intorsioni .ne. intorsionj ) then
            lsame = .false.
            return
         endif
         do itor = 1,intorsioni
            call twh_ijtor2(GLB_GET,imolty,unit,itor,ijtor2i)
            call twh_ijtor2(GLB_GET,jmolty,unit,itor,ijtor2j)
            call twh_ijtor3(GLB_GET,imolty,unit,itor,ijtor3i)
            call twh_ijtor3(GLB_GET,jmolty,unit,itor,ijtor3j)
            call twh_ijtor4(GLB_GET,imolty,unit,itor,ijtor4i)
            call twh_ijtor4(GLB_GET,jmolty,unit,itor,ijtor4j)
            call twh_ittor(GLB_GET,imolty,unit,itor,ittori)
            call twh_ittor(GLB_GET,jmolty,unit,itor,ittorj)
            if ( (ijtor2i .ne. ijtor2j) .or. ( ijtor3i .ne. ijtor3j )
     &           .or. ( ijtor4i .ne. ijtor4j )
     &           .or. ( ittori .ne. ittorj ) ) then
               lsame = .false.
               return
            endif
         enddo
c        --- improper torsions
         call twh_inimprop(GLB_GET,imolty,unit,inimproperi)
         call twh_inimprop(GLB_GET,jmolty,unit,inimproperj)
         if ( inimproperi .ne. inimproperj ) then
            lsame = .false.
            return
         endif
         do itor = 1,inimproperi
            call twh_ijimprop2(GLB_GET,imolty,unit,itor,ijimp2i)
            call twh_ijimprop2(GLB_GET,jmolty,unit,itor,ijimp2j)
            call twh_ijimprop3(GLB_GET,imolty,unit,itor,ijimp3i)
            call twh_ijimprop3(GLB_GET,jmolty,unit,itor,ijimp3j)
            call twh_ijimprop4(GLB_GET,imolty,unit,itor,ijimp4i)
            call twh_ijimprop4(GLB_GET,jmolty,unit,itor,ijimp4j)
            call twh_itimprop(GLB_GET,imolty,unit,itor,itimpropi)
            call twh_itimprop(GLB_GET,jmolty,unit,itor,itimpropj)
            if ( (ijimp2i .ne. ijimp2j) .or. (ijimp3i .ne. ijimp3j)
     &           .or. (ijimp4i .ne. ijimp4j)
     &           .or. (itimpropi .ne. itimpropj) ) then
               lsame = .false.
               return
            endif
         enddo
c        --- angle-angle
         call twh_inaa(GLB_GET,imolty,unit,inaai)
         call twh_inaa(GLB_GET,jmolty,unit,inaaj)
         if ( inaai .ne. inaaj ) then
            lsame = .false.
            return
         endif
         do iaa = 1,inaai
            call twh_ijaa0(GLB_GET,imolty,unit,iaa,ijaa0i)
            call twh_ijaa0(GLB_GET,jmolty,unit,iaa,ijaa0j)
            call twh_ijaa1(GLB_GET,imolty,unit,iaa,ijaa1i)
            call twh_ijaa1(GLB_GET,jmolty,unit,iaa,ijaa1j)
            call twh_ijaa2(GLB_GET,imolty,unit,iaa,ijaa2i)
            call twh_ijaa2(GLB_GET,jmolty,unit,iaa,ijaa2j)
            call twh_itaa(GLB_GET,imolty,unit,iaa,itaai)
            call twh_itaa(GLB_GET,jmolty,unit,iaa,itaaj)
            if ( ijaa0i .ne. ijaa0j .or. ijaa1i .ne. ijaa1j 
     &           .or. ijaa2i .ne. ijaa2j .or. itaai .ne. itaaj ) then
               lsame = .false.
               return
            endif
         enddo
      enddo

      return
      end

