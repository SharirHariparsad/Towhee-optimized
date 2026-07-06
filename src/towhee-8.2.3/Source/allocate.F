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
      subroutine twh_allocate_blocks(maxprop,maxblock)
c     ******************************************************************
c     * allocates the array for block averages                         *
c     *                                                                *
c     * originally written 02-03-2011 by M.G. Martin                   *
c     * last modified 01-16-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxprop,maxblock
c     --- local variables
      double precision dvalue

      if ( maxprop*maxblock .ne. 0 ) then
         dvalue = 0.0d0
         call twh_blockvalue(GLB_ALLOC_INIT,maxprop,maxblock,dvalue)
      endif
      return
      end

      subroutine twh_allocate_coords(natoms,maxunit)
c     ******************************************************************
c     * allocates arrays for coordinate storage                        *
c     *                                                                *
c     * originally written 12-05-2008 M.G. Martin                      *
c     * last modified 01-16-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer natoms,maxunit
c     --- local variables
      integer total
      double precision dvalue
c     --- set to null values
      dvalue = 0.0d0

c     --- total size for coordstorage is 3*natoms
      total = 3*natoms
      call twh_coordstorage(GLB_ALLOC_INIT,total,dvalue)
c     --- total size for coordtemp also must include an additional 
c     --- maxunit allocation for swapmoves
      total = 3*natoms + 3*maxunit
      call twh_coordtemp(GLB_ALLOC_INIT,total,dvalue)

      return
      end

      subroutine twh_allocate_eam(maxchain)
c     ******************************************************************
c     * allocates arrays used with the EAM potentials                  *
c     *                                                                *
c     * originally written 12-05-2008 M.G. Martin                      *
c     * last modified 01-16-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxchain
c     --- local variables
      double precision dvalue

c     --- set to null values
      dvalue = 0.0d0
c     --- array allocations only required for EAM potential
c     --- allocate real and temp
      call twh_eam_rho_real(GLB_ALLOC_INIT,maxchain,dvalue)
      call twh_eam_rho_temp(GLB_ALLOC_INIT,maxchain,dvalue)
c     --- initialize real
      call twh_eam_rho(GLB_INIT,CRD_REAL,maxchain,dvalue)

      return
      end

      subroutine twh_allocate_maxchain(maxchain)
c     ******************************************************************
c     * allocates arrays based soley upon maxchain                     *
c     *                                                                *
c     * originally written 12-03-2008 M.G. Martin                      *
c     * last modified 01-16-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxchain
c     --- local variables
      logical lvalue
      integer ivalue,total
      double precision dvalue
      double precision vector
      dimension vector(3)

c     --- set ivalue to null values
      lvalue = .false.
      ivalue = 0
      dvalue = 0.0d0
c     --- arrays depending only on maxchain
      call twh_gyration(GLB_ALLOC,maxchain,dvalue)
      call twh_rcmu(GLB_ALLOC_INIT,maxchain,dvalue)
      call twh_globalpos(GLB_ALLOC,maxchain,ivalue)
      call twh_lmoveme(GLB_ALLOC,maxchain,lvalue)
      call twh_lperiodic(GLB_ALLOC_INIT,maxchain,lvalue)
      call twh_moltyp(GLB_ALLOC_INIT,maxchain,ivalue)
      call twh_nboxi(GLB_ALLOC,maxchain,ivalue)
c     --- allocate COM storage arrays for volume moves
      total = 3*maxchain-1
      call twh_arbcomfield(GLB_ALLOC,total,vector)
      call twh_comfield(GLB_ALLOC,total,vector)
      call twh_comtempfield(GLB_ALLOC,total,vector)

      return
      end

      subroutine twh_allocate_maxchain_maxmolty(maxchain,maxmolty)
c     ******************************************************************
c     * allocate arrays that depend on the maximum number of chains and*
c     * the maximum number of molecule types                           *
c     *                                                                *
c     * originally written 06-24-2009 by M.G. Martin                   *
c     * last modified 01-16-2018 M.G. Martin                           *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxchain,maxmolty
c     --- local variables
      integer ivalue

c     --- set to null values
      ivalue = 0

c     --- allocate arrays
      call twh_parall(GLB_ALLOC_INIT,maxmolty,maxchain,ivalue)

      return
      end

      subroutine twh_allocate_maxmolty(maxmolty)
c     ******************************************************************
c     * allocates arrays based soley upon maxmolty                     *
c     *                                                                *
c     * originally written 12-03-2008 M.G. Martin                      *
c     * last modified 01-16-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxmolty
c     --- local variables
      integer ivalue

c     --- set ivalue to null values
      ivalue = 0
         
c     --- arrays depending only on maxmolty
      call twh_growvalidcount(GLB_ALLOC_INIT,maxmolty,ivalue)

      return
      end

      subroutine twh_allocate_maxmolty_maxunit_maxnumtor(maxmolty
     &     ,maxunit,maxnumtor)
c     ******************************************************************
c     * allocate arrays that depend on the maximum number of molecule  *
c     * types, maximum units per molecule and maximum number of        *
c     * torsions on any unit                                           *
c     *                                                                *
c     * originally written 01-30-2011 by M.G. Martin                   *
c     * last modified 03-23-2018 M.G. Martin                           *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxmolty,maxunit,maxnumtor
c     --- local variables
      integer ivalue

      if ( maxmolty*maxunit*maxnumtor .eq. 0 ) return

c     --- default to CODE_ONE
      ivalue = CODE_ONE

c     --- allocate arrays
      call twh_torofcode(GLB_ALLOC_INIT,maxmolty,maxunit,maxnumtor
     &     ,ivalue)

      return
      end

      subroutine twh_allocate_maxmolty_numboxes(maxmolty,numboxes)
c     ******************************************************************
c     * allocate arrays based on maximum number of molecule types and  *
c     * the maximum number of simulation boxes                         *
c     *                                                                *
c     * originally written 06-17-2009 by M.G. Martin                   *
c     * last modified 08-11-2018 M.G. Martin                           *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxmolty,numboxes
c     --- local variables
      integer itemp,npair
      double precision dvaluezero,dvaluehalf
c     --- set to null values
      itemp = 0
      dvaluezero = 0.0d0
      dvaluehalf = 0.5d0
c     --- allocate and initialize acnrot
      call twh_acnrot(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
      call twh_acsrot(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
      call twh_rmrot(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluehalf)
      call twh_barot(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
      call twh_bnrot(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
c     --- allocate and initialize composite move
      call twh_acncomp(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
      call twh_acscomp(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
      call twh_rmcomrot(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluehalf)
      call twh_rmcomtra(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluehalf)
c     --- accumulators for the single-atom translation move
c     --- allocate and initialize acstraa
      call twh_acstraa(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
c     --- allocate and initialize acntraa
      call twh_acntraa(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
c     --- allocate and initialize batraa
      call twh_batraa(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
c     --- allocate and initialize bntraa
      call twh_bntraa(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
      call twh_rmtraa(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluehalf)
c     --- accumulators for the COM translation move
c     --- allocate and initialize acntrac
      call twh_acntrac(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
c     --- allocate and initialize acstrac
      call twh_acstrac(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
c     --- allocate and initialize batrac
      call twh_batrac(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
c     --- allocate and initialize bntrac
      call twh_bntrac(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluezero)
      call twh_rmtrac(GLB_ALLOC_INIT,maxmolty,numboxes,dvaluehalf)
      npair = ((maxmolty-1)*(maxmolty))/2
      if ( npair .gt. 0 ) then
c        --- allocate and initialize acnswitch
         call twh_acnswitch(GLB_ALLOC_INIT,npair,numboxes,dvaluezero)
c        --- allocate and initialize acsswitch
         call twh_acsswitch(GLB_ALLOC_INIT,npair,numboxes,dvaluezero)
c        --- allocate and initialize bncom2bswitch
         call twh_bncom2bswitch(GLB_ALLOC_INIT,npair,numboxes,numboxes
     &        ,dvaluezero)
      endif
      return
      end

      subroutine twh_allocate_maxmolty_numboxes_nchain(maxmolty,maxboxes
     &     ,maxchain)
c     ******************************************************************
c     * allocate arrays that depend on the maximum number of molecule  *
c     * types, boxes, and chains                                       *
c     *                                                                *
c     * originally written 06-25-2009 by M.G. Martin                   *
c     * last modified 03-23-2018 M.G. Martin                           *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxmolty,maxboxes,maxchain
c     --- local variables
      integer ivalue

c     --- set to null values
      ivalue = 0
c     --- allocate arrays
      call twh_chainlist(GLB_ALLOC_INIT,maxchain,maxboxes,maxmolty
     &     ,ivalue)

      return
      end

      subroutine twh_allocate_maxunit(maxunit)
c     ******************************************************************
c     * allocates arrays based soley upon maxunit                      *
c     *                                                                *
c     * originally written 12-03-2008 M.G. Martin                      *
c     * last modified 01-16-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxunit
c     --- local variables
      logical lvalue
      integer ivalue

c     --- set ivalue to null values
      ivalue = 0
      lvalue = .false.

c     --- arrays depending only on maxunit
      call twh_glist(GLB_ALLOC,maxunit,ivalue)
      call twh_growfrom(GLB_ALLOC,maxunit,ivalue)
      call twh_grownum(GLB_ALLOC,maxunit,ivalue)
      call twh_growprev(GLB_ALLOC,maxunit,ivalue)
      call twh_lexsched(GLB_ALLOC,maxunit,lvalue)
      call twh_lexist(GLB_ALLOC,maxunit,lvalue)

      return
      end

      subroutine twh_allocate_nboxpair(nboxpair)
c     ******************************************************************
c     * allocate arrays that depend on the number of box pairs         *
c     *                                                                *
c     * originally written 01-14-2018 by M.G. Martin                   *
c     * last modified 04-12-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nboxpair
      double precision dvalue

      if ( nboxpair .eq. 0 ) return
      call twh_pairbox(GLB_ALLOC,nboxpair,2,0)
      dvalue = 1.0d0
      call twh_pm2cbswpr(GLB_ALLOC_INIT,nboxpair,dvalue)
      call twh_pm2rbswpr(GLB_ALLOC_INIT,nboxpair,dvalue)
      call twh_pm2comswboxpair(GLB_ALLOC_INIT,nboxpair,dvalue)

      return
      end

      subroutine twh_allocate_numboxes(numboxes)
c     ******************************************************************
c     * allocate arrays based on numboxes                              *
c     *                                                                *
c     * originally written 12-04-2008 by M.G. Martin                   *
c     * last modified 04-03-2018 M.G. Martin                           *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer numboxes
c     --- local variables
      integer itemp

      itemp = 5
c     --- ewald kmax
      call twh_ewald_kmax(GLB_ALLOC_INIT,MAXBXD,numboxes,itemp)

      return
      end

      subroutine twh_allocate_nvmmax(nvmmax)
c     ******************************************************************
c     * allocate arrays based on nvmmax                                *
c     *                                                                *
c     * originally written 12-04-2008 by M.G. Martin                   *
c     * last modified 01-17-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nvmmax
c     --- local variables
      integer ivalue
      double precision dvalue

c     --- no need to allocate if the value is 0
      if ( nvmmax .eq. 0 ) return

      ivalue = 0
c     --- isotropic volume move
      call twh_acnvol(GLB_ALLOC_INIT,nvmmax,ivalue)
      call twh_acsvol(GLB_ALLOC_INIT,nvmmax,ivalue)
      dvalue = 0.1d0
      call twh_rmvol(GLB_ALLOC_INIT,nvmmax,dvalue)
      call twh_bavol(GLB_ALLOC_INIT,nvmmax,ivalue)
      call twh_bnvol(GLB_ALLOC_INIT,nvmmax,ivalue)
c     --- cell volume move
      dvalue = 0.0d0
      call twh_acncell(GLB_ALLOC_INIT,nvmmax,9,ivalue)
      call twh_acscell(GLB_ALLOC_INIT,nvmmax,9,ivalue)
      call twh_bacell(GLB_ALLOC_INIT,nvmmax,9,ivalue)
      call twh_bncell(GLB_ALLOC_INIT,nvmmax,9,ivalue)
      dvalue = 1.0d0
      call twh_rmcell(GLB_ALLOC_INIT,nvmmax,3,3,dvalue)

      return
      end

      subroutine twh_allocate_nvmsim(nvmsim)
c     ******************************************************************
c     * allocate arrays based on nvmsim                                *
c     *                                                                *
c     * originally written 01-15-2018 by M.G. Martin                   *
c     * last modified 01-16-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nvmsim
c     --- local variables
      double precision dvalue
      if ( nvmsim .eq. 0 ) return

      dvalue = 1.0d0
      call twh_pmcellpr(GLB_ALLOC_INIT,nvmsim,dvalue)
      call twh_pmvlpr(GLB_ALLOC_INIT,nvmsim,dvalue)
      dvalue = 0.5d0
      call twh_pmcellpt(GLB_ALLOC_INIT,nvmsim,dvalue)

      return
      end

      subroutine twh_allocate_tmmc(maxchain)
c     ******************************************************************
c     * allocate arrays for TMMC                                       *
c     * most of these arrays run from 0 to nchain, where maxchain is   *
c     * nchain+1                                                       *
c     *                                                                *
c     * originally written 12-05-2008 by M.G. Martin                   *
c     * last modified 06-23-2009 M.G. Martin                           *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer maxchain
c     --- local variables
      double precision dvalue

c     --- set to null values
      dvalue = 0.0d0

c     --- allocate arrays
      call twh_tmmc_weight(GLB_ALLOC,maxchain,dvalue)
      call twh_v_semigrand(GLB_ALLOC,maxchain,dvalue)
      call twh_c_matrix(GLB_ALLOC,maxchain,1,dvalue)

      return
      end

