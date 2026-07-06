#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2018 Marcus G. Martin                       *
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
      subroutine twh_resetcoulomb(lfailure,bxdflag,ibox)
c     ******************************************************************
c     * resets all of the coulomb parameters in the specified box      *
c     * complicated dependencies on the bxdflag.  We are always        *
c     * fetching from BXD_REAL, but set to the bxdflag                 *
c     *                                                                *
c     * originally resetewald 08-13-2004 by M.G. Martin                *
c     * last modified 04-03-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_MAXBOXLENGTH
#define FUNCTION_MINBOXLENGTH
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer bxdflag,ibox
c     --- local variables
      logical lelect,lewald
      character*20 coulombstyle
      integer kmax
      double precision minlength,onepi,ewald_prec,kalp,calp,rcelect
c     --- check if electrostatics is used
      call twh_lelect(GLB_GET,lelect)
      if ( .not. lelect ) return
c     --- get the coulombstyle
      call twh_coulombstyle(GLB_GET,coulombstyle)
      if ( coulombstyle .eq. 'ewald_fixed_kmax' ) then
         if ( bxdflag .eq. BXD_T_TO_R ) then
c           --- set real to the value of temp
            call twh_ewald_calp(GLB_STORE,bxdflag,ibox,calp)
            call twh_rcelect(GLB_STORE,bxdflag,ibox,rcelect)
            return
         else
c           --- set kmax to the real value in box 1
            call twh_ewald_kmax(GLB_GET,BXD_REAL,1,kmax)
            call twh_ewald_kmax(GLB_SET,bxdflag,ibox,kmax)
c           --- set calp based upon the minimum box lengths
c           --- determine the minimum box length
            minlength = twh_minboxlength(bxdflag,ibox)
            call twh_ewald_kalp(GLB_GET,kalp)
            calp = kalp/minlength
            call twh_ewald_calp(GLB_SET,bxdflag,ibox,calp)
            rcelect = 0.5d0*minlength
            call twh_rcelect(GLB_SET,bxdflag,ibox,rcelect)
         endif
      elseif ( coulombstyle .eq. 'ewald_fixed_cutoff' ) then
         if ( bxdflag .eq. BXD_T_TO_R ) then
c           --- set real to the value stored in temp
            call twh_ewald_calp(GLB_STORE,bxdflag,ibox,calp)
            call twh_ewald_kmax(GLB_GET,BXD_TEMP,ibox,kmax)
            call twh_ewald_kmax(GLB_SET,BXD_REAL,ibox,kmax)
            return
         else
            call twh_constant_pi(GLB_GET,onepi)
            call twh_ewald_prec(GLB_GET,ewald_prec)
            call twh_rcelect(GLB_GET,BXD_REAL,ibox,rcelect)
            calp = (1.35D0 - 0.15D0*log(ewald_prec))/rcelect
            call twh_ewald_calp(GLB_SET,bxdflag,ibox,calp)
            kmax = int( (calp*twh_maxboxlength(bxdflag,ibox)/onepi)
     &           *dsqrt(-dlog(ewald_prec)))
            call twh_ewald_kmax(GLB_SET,bxdflag,ibox,kmax)
         endif
      elseif ( coulombstyle .eq. 'minimum image' ) then
         if ( bxdflag .eq. BXD_T_TO_R ) then
            call twh_rcelect(GLB_STORE,bxdflag,ibox,rcelect)
            return
         else
c           --- we interact with the closest atom regardless of cutoff
            rcelect = 10.0d0*twh_maxboxlength(bxdflag,ibox)
            call twh_rcelect(GLB_SET,bxdflag,ibox,rcelect)
c           --- disable the center of mass cutoff
            call twh_lusecom(GLB_SET,.false.)
         endif
      else
         write(6,*) 'RESETCOULOMB: unknown coulombstyle'
         lfailure = .true.
         return
      endif

      call twh_lewald(GLB_GET,lewald)
      if ( lewald ) then
c        --- safety check kmax
         call twh_ewald_kmax(GLB_GET,bxdflag,ibox,kmax)
         if ( kmax .lt. 0 ) then
            write(6,*) 'RESETCOULOMB: negative kmax ',kmax
     &           ,' in box ',ibox
            lfailure = .true.
            return
         elseif ( kmax .gt. MAXKMAX ) then
            write(6,*) 'RESETCOULOMB: new kmax exceeds array dimensions'
            write(6,*) 'Please set MAXKMAX to at least ',kmax
            lfailure = .true.
            return
         endif
      endif

      return
      end
