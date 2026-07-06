#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2020-2020 Marcus G. Martin                       *
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
      subroutine twh_writechempotdata(lfailure,datacode,ibox,imolty
     &     ,ncmolty,boxvolume,nvtweight)
c     ******************************************************************
c     * this routine outputs a single data point of insertion chemical *
c     * potential information to the towhee_chempotdata file           *
c     * for post-processing analysis                                   *
c     *                                                                *
c     * originally written 03-11-2020 by M.G. Martin                   *
c     * last modified 03-11-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer datacode,ibox,imolty,ncmolty
      double precision boxvolume,nvtweight
c     --- local variables
      character*30 filename
      logical loutchempotdata
      integer ifile,nmolty,iversion,itemp
      double precision debroglie
c     --- set the movie file number
      ifile = FILE_CHEMPOTDATA
c     --- chemical potential version
      iversion = 1

      call twh_loutchempotdata(GLB_GET,loutchempotdata)
      if ( .not. loutchempotdata ) return

      if ( datacode .eq. GLB_OUTPUT ) then
c        --- regular output of chemical potential data
         write(ifile,'(i13,1x,i11,1x,i12,1x,f20.10,1x,e20.10)')
     &        ibox,imolty,ncmolty,boxvolume,nvtweight
      elseif ( datacode .eq. GLB_INIT ) then
c        --- open the chempotdata file
         filename = 'towhee_chempotdata'
         call twh_openfile(lfailure,ifile,filename,'DEFAULT'
     &        ,'UNKNOWN','DEFAULT')
         if ( lfailure ) return
c        --- output version
         write(ifile,'(a32,1x,i10)') 'Chemical Potential Data Version:'
     &        ,iversion
         call twh_nmolty(GLB_GET,nmolty)
         write(ifile,'(a15,1x,i10)') 'Molecule Types:',nmolty
         write(ifile,'(a38)') 'Molecule Type     Debroglie Wavelength'
         do itemp = 1,nmolty
            call twh_debroglie(GLB_GET,itemp,debroglie)
            write(ifile,*) itemp,debroglie
         enddo
         write(ifile,'(a13,1x,a11,1x,a12,11x,a10,11x,a10)')
     &        'Insertion Box','Molec. Type'
     &        ,'Molec. Count','Box Volume','Rosenbluth' 
      elseif ( datacode .eq. GLB_CLOSE ) then
         close(ifile)
      else
         write(6,*) 'WRITECHEMPOTDATA: unknown datacode:',datacode
         lfailure = .true.
      endif

      return
      end
