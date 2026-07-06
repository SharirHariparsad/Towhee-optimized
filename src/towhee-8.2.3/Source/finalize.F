#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2008-2018 Marcus G. Martin                       *
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
      subroutine twh_finalize( )
c     ******************************************************************
c     * frees all of the allocated memory                              *
c     *                                                                *
c     * originally written 12-03-2008 M.G. Martin                      *
c     * last modified 08-11-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- local variables
      logical lvalue
      integer itemp,jtemp,ktemp,ivalue,nboxpair,nvmmax,nvmsim
      double precision dvalue
c     --- constants
      call twh_nboxpair(GLB_GET,nboxpair)
      call twh_nvmmax(GLB_GET,nvmmax)
      call twh_nvmsim(GLB_GET,nvmsim)
c     --- set temp variables to null values
      lvalue = .false.
      itemp = 0
      jtemp = 0
      ktemp = 0
      ivalue = 0
      dvalue = 0.0d0
c     --- free arrays based only on maxmolty
      call twh_growvalidcount(GLB_FREE,itemp,jtemp)
c     --- free arrays based only on maxunit
      call twh_glist(GLB_FREE,itemp,jtemp)
      call twh_growfrom(GLB_FREE,itemp,jtemp)
      call twh_grownum(GLB_FREE,itemp,jtemp)
      call twh_growprev(GLB_FREE,itemp,jtemp)
      call twh_lexsched(GLB_FREE,itemp,lvalue)
      call twh_lexist(GLB_FREE,itemp,lvalue)
c     --- free arrays based only on maxchain
      call twh_eam_rho_real(GLB_FREE,itemp,dvalue)
      call twh_eam_rho_temp(GLB_FREE,itemp,dvalue)
      call twh_gyration(GLB_FREE,itemp,dvalue)
      call twh_tmmc_weight(GLB_FREE,itemp,dvalue)
      call twh_c_matrix(GLB_FREE,itemp,dvalue)
      call twh_v_semigrand(GLB_FREE,itemp,dvalue)
      call twh_rcmu(GLB_FREE,itemp,dvalue)
      call twh_globalpos(GLB_FREE,itemp,jtemp)
      call twh_logical_moveme(GLB_FREE,itemp,lvalue)
      call twh_logical_periodic(GLB_FREE,itemp,lvalue)
      call twh_moltyp(GLB_FREE,itemp,jtemp)
      call twh_nboxi(GLB_FREE,itemp,jtemp)
c     --- free arrays based only on nboxpair
      if ( nboxpair .ne. 0 ) then
         call twh_pairbox(GLB_FREE,itemp,jtemp,ivalue)
         call twh_pm2cbswpr(GLB_FREE,itemp,dvalue)
         call twh_pm2rbswpr(GLB_FREE,itemp,dvalue)
      endif
c     --- free arrays based only on nvmmax
      if ( nvmmax .ne. 0 ) then
         call twh_acnvol(GLB_FREE,itemp,jtemp)
         call twh_acsvol(GLB_FREE,itemp,jtemp)
         call twh_rmvol(GLB_FREE,itemp,dvalue)
         call twh_bavol(GLB_FREE,itemp,ivalue)
         call twh_bnvol(GLB_FREE,itemp,ivalue)
         call twh_acncell(GLB_FREE,itemp,9,ivalue)
         call twh_acscell(GLB_FREE,itemp,9,ivalue)
         call twh_bacell(GLB_FREE,itemp,9,ivalue)
         call twh_bncell(GLB_FREE,itemp,9,ivalue)
         call twh_rmcell(GLB_FREE,itemp,3,3,dvalue)
      endif
c     --- free arrays based only on nvmsim
      if ( nvmsim .ne. 0 ) then
         call twh_pmcellpr(GLB_FREE,itemp,dvalue)
         call twh_pmcellpt(GLB_FREE,itemp,dvalue)
         call twh_pmvlpr(GLB_FREE,itemp,dvalue)
      endif
c     --- free arrays based on cubes
      call twh_cubeletweight(GLB_FREE,itemp,dvalue)
c     --- free arrays for coordinate storage
      call twh_coordfield(GLB_FREE,itemp,dvalue)
      call twh_coordstorage(GLB_FREE,itemp,dvalue)
      call twh_coordtemp(GLB_FREE,itemp,dvalue)
      call twh_arbcomfield(GLB_FREE,itemp,dvalue)
      call twh_comfield(GLB_FREE,itemp,dvalue)
      call twh_comtempfield(GLB_FREE,itemp,dvalue)
c     --- free arrays based on numboxes and nmolty
      call twh_acsrot(GLB_FREE,itemp,jtemp,dvalue)
      call twh_acnrot(GLB_FREE,itemp,jtemp,dvalue)
      call twh_barot(GLB_FREE,itemp,jtemp,dvalue)
      call twh_bnrot(GLB_FREE,itemp,jtemp,dvalue)
      call twh_rmrot(GLB_FREE,itemp,jtemp,dvalue)
      call twh_acncomp(GLB_FREE,itemp,jtemp,dvalue)
      call twh_rmcomrot(GLB_FREE,itemp,jtemp,dvalue)
      call twh_rmcomtra(GLB_FREE,itemp,jtemp,dvalue)
      call twh_acnswitch(GLB_FREE,itemp,jtemp,dvalue)
      call twh_acsswitch(GLB_FREE,itemp,jtemp,dvalue)
      call twh_bncom2bswitch(GLB_FREE,itemp,jtemp,ktemp,dvalue)
      call twh_acntraa(GLB_FREE,itemp,jtemp,dvalue)
      call twh_acstraa(GLB_FREE,itemp,jtemp,dvalue)
      call twh_batraa(GLB_FREE,itemp,jtemp,dvalue)
      call twh_bntraa(GLB_FREE,itemp,jtemp,dvalue)
      call twh_acntrac(GLB_FREE,itemp,jtemp,dvalue)
      call twh_acstrac(GLB_FREE,itemp,jtemp,dvalue)
      call twh_batrac(GLB_FREE,itemp,jtemp,dvalue)
      call twh_bntrac(GLB_FREE,itemp,jtemp,dvalue)
      call twh_rmtraa(GLB_FREE,itemp,jtemp,dvalue)
      call twh_rmtrac(GLB_FREE,itemp,jtemp,dvalue)
c     --- free based on molty and chain
      call twh_parall(GLB_FREE,itemp,jtemp,ivalue)
      call twh_chainlist(GLB_FREE,itemp,jtemp,ktemp,ivalue)
c     --- free based on molty, unit and torsion
      call twh_torofcode(GLB_FREE,itemp,jtemp,ktemp,ivalue)
      return
      end
