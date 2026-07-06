#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
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
      subroutine twh_engtotal(lfailure,lovrlap
     &     ,vclassical,vinter,vthrinter,vtail,vintra,vvib,vbend,vtg
     &     ,vext,velect,vsolvation,vforeign,dl_inter, dl_elect
     &     ,crdflag,bxdflag,ibox
     &     ,lmolecmove,lmovelist,lreport,lbonded)
c     ******************************************************************
c     * calculates the total potential energy for a configuration.     *
c     *                                                                *
c     * The following are set by the subroutine                        *
c     * lovrlap is true if this subroutine finds substantial overlap   *
c     * vclassical is the total classical energy                       *
c     * vinter is the total intermolecular two-body, three-body and    *
c     *        tail correction non-coulombic energy                    *
c     * vtail: just the tail correction terms of the intermolecular    *
c     *        energy                                                  *
c     * vthrinter: just the 3-body terms of the intermolecular energy  *
c     * vintra is the intramolecular (non-coulombic) non-bonded energy *
c     * vvib is the bond vibration energy                              *
c     * vbend is the bond bending energy                               *
c     * vtg is the dihedral angle energy - including improper torsions *
c     * vext is the energy of external fields/walls                    *
c     * velect is the total coulombic energy - execpt 1-4 coulombics   *
c     *        which are included in the torsional potential           *
c     * vsolvation is the solvation energy                             *
c     * vforeign is total classical energy evaluated at different      *
c     *     values of lambda for 'Scaled Lennard-Jones' potential.     *
c     *     See scaled_potential.F                                     *
c     * dl_inter -- derivative of two body potential with respect      *
c     *     to lambda_lj, calculated if scp_eval_derivative true       *
c     *     (else 0)                                                   *
c     * dl_elect -- derivative of inter-molecule coulomb potential     *
c     *     with respect to lambda_c                                   *
c     *                                                                *
c     * The following are set by the caller                            *
c     * crdflag is the coordinate set (real or temp) to use here       *
c     * bxdflag is the box dimensions (real or temp) to use here       *
c     * ibox is the box we are computing the energy in                 *
c     * lmolecmove is true if this call is from a Monte Carlo move,    *
c     *            otherwise this is the call at the start or end of   *
c     *            the simulation.  If this is a move then we use the  *
c     *            additional COM cutoff.                              *
c     * lmovelist is true if we are performing a move and only need the*
c     *            interactions for a subset of the molecules          *
c     * lreport is true if we want to print out all of the energies    *
c     * lbonded is true if we are computing bonded interactions        *
c     *            (vibration, bending, torsion etc)                   *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May exit early if overlap found (error may be printed)       *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modified 03-29-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lovrlap, lmolecmove,lmovelist,lreport,lbonded
      integer crdflag,bxdflag,ibox
      double precision vclassical,vinter,vintra,vthrinter,vtail,vvib
      double precision vbend,vtg,velect,vsolvation,vext
      double precision vforeign,dl_inter, dl_elect
      dimension vforeign(MAX_FOREIGN_LAMBDA)
c     --- local variables
c     --- logical scalars
      logical lmadelung,ltailc,lelect,lewald
c     --- integer scalars
      integer ichain, nmcount, cutstyle,lambda,testboxi
      integer classical_potnum,nchain,nchbox
      integer scp_num_foreign_lambda
c     --- double precision scalars
      double precision maddist,sself,correct,vtmp
      double precision velectinter,velectintra,vtwointer,vangleangle
      double precision vimpr,vrealelect,vrecipsum,vbondbond

c      write(6,*) 'start ENGTOTAL'

c     --- retrieve constants
      call twh_ltailc(GLB_GET,ltailc)
      call twh_lelect(GLB_GET,lelect)
      call twh_lewald(GLB_GET,lewald)

      if ( lreport ) then
c        --- initialize report variables
         call twh_report(GLB_INIT,REP_MAX,0.0d0)
      endif

      cutstyle = 0
      lovrlap = .false.

      vclassical = 0.0d0
      vinter = 0.0d0
      vtwointer = 0.0d0
      vthrinter = 0.0d0
      vintra = 0.0d0
      vtail = 0.0d0
      vtg = 0.0d0
      vimpr = 0.0d0
      vbend = 0.0d0
      vangleangle = 0.0d0
      vvib = 0.0d0
      vbondbond = 0.0d0 
      vext = 0.0d0
      velect = 0.0d0
      sself = 0.0d0
      correct = 0.0d0
      vrealelect = 0.0d0
      vrecipsum = 0.0d0
      vsolvation = 0.0d0
      velectinter = 0.0d0
      velectintra = 0.0d0
      dl_inter = 0.0d0
      dl_elect = 0.0d0
      do lambda=1,MAX_FOREIGN_LAMBDA
         vforeign(lambda) = 0.0d0
      enddo
c     --- determine whether we want to compute the madelung constant
      if ( lreport ) then
         call twh_lcompute_madelung(GLB_GET,lmadelung)
      else
         lmadelung = .false.
      endif
c     ---  check the molecule count
      call twh_nchain(GLB_GET,nchain)
      nmcount = 0
      do ichain = 1, nchain
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox ) then
            nmcount = nmcount + 1
         endif
      enddo

      call twh_nchbox(GLB_GET,ibox,nchbox)
      if ( nmcount .ne. nchbox ) then
         write(6,*) 'ENGTOTAL: nmcount ne nchbox'
         write(6,*) 'ibox nmcount ',ibox,nmcount
         write(6,*) 'nchbox ',nchbox
         write(6,*) 'nchain ',nchain
         lfailure = .true.
         return
      endif
 
      call twh_engtotal_twobody(lfailure,crdflag
     &     ,velectinter,vtwointer
     &     ,vforeign,maddist,dl_inter,dl_elect,lovrlap,lmadelung
     &     ,lmolecmove,lmovelist,cutstyle,bxdflag,ibox )
      if ( lfailure ) return
      if ( lovrlap ) then
         return
      endif
      if ( ltailc ) then
         call twh_engtotal_tail(vtail,bxdflag,ibox)
      endif
      call twh_cp_number(GLB_GET,classical_potnum)
      if ( classical_potnum .eq. CP_STILLINGER_WEBER ) then
         call  twh_engtotal_tb_sw(lfailure,crdflag,vthrinter
     &        ,bxdflag,ibox)
         if ( lfailure ) return
      elseif ( classical_potnum .eq. CP_EMBEDDED_ATOM_METHOD ) then
         call twh_engtotal_tb_eam(lfailure,crdflag
     &        ,vthrinter, lovrlap,bxdflag,ibox)
         if ( lfailure ) return
      endif
      if ( lewald ) then
c        --- TODO: add foreign energy terms here
         call twh_engtotal_ewald(lfailure,crdflag,bxdflag,ibox
     &        ,vrecipsum,sself,correct)
         if ( lfailure ) return
      endif
      call twh_engtotal_intrachain(lfailure,crdflag
     &     ,vintra, velectintra
     &     , lovrlap,cutstyle,bxdflag,ibox, lmolecmove, lmovelist)
      if ( lfailure ) return
      if ( lovrlap ) then
         return
      endif   
      if ( lbonded ) then
         call twh_engtotal_vibbendtors(lfailure,crdflag
     &        ,vvib,vbend,vtg,vangleangle,vimpr,vbondbond
     &        ,lovrlap,lreport,bxdflag,ibox)
         if ( lfailure ) return
         if ( lovrlap ) then
            return
         endif
      endif 
      call twh_engtotal_ext(lfailure,bxdflag,crdflag
     &     ,vext,lovrlap,ibox,lmolecmove
     &     ,lmovelist )
      if ( lfailure ) return
      if ( lovrlap ) return
      call twh_engtotal_solv(lfailure,crdflag,bxdflag,ibox
     &     ,vsolvation)
      if ( lfailure ) return
      if ( lelect ) then
         vrealelect = velectintra + velectinter
         velect = vrealelect + sself + correct + vrecipsum
      endif
      vinter = vtwointer + vthrinter + vtail
      vtmp = vintra + vvib + vbondbond + vbend 
     &     + vangleangle + vtg + vimpr + vext + vsolvation
      vclassical = vinter + velect + vtmp
c     --- foreign energies add terms which do not change with lambda
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      do lambda = 1, scp_num_foreign_lambda
         vforeign(lambda) = vforeign(lambda) + vtmp
      enddo

      if ( lreport ) then
         call twh_eng_print(nmcount, vvib, vbondbond, vbend, 
     &       vangleangle, vtg, vimpr, vintra, vinter, vtwointer, 
     &       vthrinter, velect, vrealelect, sself, correct, vrecipsum, 
     &       vext, vsolvation, vclassical, vtail, velectinter, 
     &       velectintra,vforeign,dl_inter,dl_elect, lmadelung,
     &       maddist, ibox)
      endif

c     --- add the bond-bond terms into the regular bond
      vvib = vvib + vbondbond
c     --- add the angleangle terms into the regular angle
      vbend = vbend + vangleangle
c     --- add the improper terms into the regular torsions
      vtg = vtg + vimpr
c      write(6,*) 'end ENGTOTAL'
      return
      end

      subroutine twh_engtotal_twobody(lfailure,crdflag
     &     ,velectinter,vtwointer
     &     ,vforeign,maddist,dl_twointer,dl_elect,lovrlap,lmadelung
     &     ,lmolecmove,lmovelist,cutstyle,bxdflag,ibox)
c     ******************************************************************
c     * Evaluate two body non-bonded interactions                      *
c     * Arguments set by subroutine:                                   *
c     *  velectinter -- the intermolecular component of real space     *
c     *     coulombic potential                                        *
c     *  vtwointer -- the two-body non-bond potential                  *
c     *  vforeign -- coulomb and twobody energies at different lambda  *
c     *  maddist -- Distance for Madelung constant                     *
c     *  dl_twointer -- derivative of two body potential with respect  *
c     *     to lambda_lj, calculated if scp_eval_derivative true       *
c     *  dl_elect -- derivative of inter-molecule coulomb potential    *
c     *     with respect to lambda_lj                                  *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  lmadelung -- compute madelung constant?                       *
c     *  lmolecmove -- is true if this call is from a Monte Carlo move *
c     *  lmovelist -- true if we are performing a move and only need   *
c     *            the interactions for a subset of the molecules      *
c     *  cutstyle -- cutoff style.  See twh_vtwobody for definition    *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * All energies and derivatives in units of [K], with the         *
c     * exception of velectinter, which is in units of [EE]            *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May halt program                                             *
c     *   May exit early (sets lovrlap=true)                           *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 03-29-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision velectinter, vtwointer, maddist
      double precision dl_twointer, dl_elect
      double precision vforeign(MAX_FOREIGN_LAMBDA)
      logical lfailure,lmadelung, lmolecmove, lmovelist, lovrlap
      integer crdflag,cutstyle,bxdflag,ibox
c     --- local variables
c     --- logical scalars
      logical lcompute,leval_c,leval_v,lusecom,lmovemei,lmovemej
      logical luse_ci,luse_cj,luse_vi,luse_vj,lelect
c     --- integer scalars
      integer ichain,imolty,jchain,jmolty,ii,ntii,jj,ntjj,l,nchain
      integer numiunit,numjunit,scp_num_foreign_lambda
      integer testboxi,testboxj
c     --- double precision scalars
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,rijsq,rij,dvdl
      double precision xcmj,ycmj,zcmj,rcmj,rcelect,vcoulomb
c     --- vf holds incremental foreign energies from vtwobody, vcoulomb
c     --- calls
c     --- vfelec and vftwo keep running totals, their sum becomes
c     --- vforeign
      double precision vf, vfelec, vftwo
      dimension vf(MAX_FOREIGN_LAMBDA), vfelec(MAX_FOREIGN_LAMBDA), 
     &      vftwo(MAX_FOREIGN_LAMBDA)

c     --- retrieve constants
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      call twh_lusecom(GLB_GET,lusecom)
      call twh_rcelect(GLB_GET,bxdflag,ibox,rcelect)
      call twh_lelect(GLB_GET,lelect)

      vtwointer = 0.0d0
      velectinter = 0.0d0
      do l=1,MAX_FOREIGN_LAMBDA
         vfelec(l) = 0.0d0
         vftwo(l) = 0.0d0
      enddo
      maddist = 1d10
      dl_twointer = 0.0d0
      dl_elect = 0.0d0
      rcmi = 0.0d0

c     --- two-body nonbonded interactions
c     --- loop over all chains i 
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1, nchain - 1
c        --- check if i is in relevant box
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,numiunit)
            if ( lmolecmove ) then
               call twh_com(GLB_GET,crdflag,ichain,xcmi,ycmi,zcmi)
               call twh_rcmu(GLB_GET,ichain,rcmi)
            endif
c           --- loop over all chains j with j>i 
            do jchain = ichain + 1, nchain
c              --- check for simulation box
               call twh_nboxi(GLB_GET,jchain,testboxj)
               if ( testboxj .eq. ibox ) then
                  call twh_moltyp(GLB_GET,jchain,jmolty)
                  call twh_nunit(GLB_GET,jmolty,numjunit)
                  lcompute = .true.
                  if ( lmovelist ) then
c                    --- see if at least one of these molecules is on
c                    --- the movelist, otherwise skip
                     call twh_lmoveme(GLB_GET,ichain,lmovemei)
                     call twh_lmoveme(GLB_GET,jchain,lmovemej)
                     if ( .not. (lmovemei .or. lmovemej) ) then
c                       --- neither of these molecules moved, do not
c                       --- compute energy
                        lcompute = .false.
                     endif
                  endif
                  if ( lmolecmove .and. lusecom ) then
c                    --- check if ctrmas within rcmsq
c                    --- minimum image the ctrmas pair separations
                     call twh_com(GLB_GET,crdflag,jchain
     &                    ,xcmj,ycmj,zcmj)
                     call twh_mimage(bxdflag,ibox,xcmi,ycmi,zcmi,xcmj
     &                    ,ycmj,zcmj,rxuij,ryuij,rzuij)
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     call twh_rcmu(GLB_GET,jchain,rcmj)
c                    --- start with the vdw cutoff
                     call twh_rcut(GLB_GET,rcm)
                     if ( lelect ) then
c                       --- check the coulombic cutoff
                        if ( rcelect .gt. rcm ) rcm = rcelect
                     endif
c                    --- add in the com radii
                     rcm = rcm + rcmi + rcmj
                     rcmsq = rcm*rcm
                     if ( rijsq .gt. rcmsq ) lcompute = .false.
                  endif
                  if ( lcompute ) then
                     do ii = 1,numiunit
                        call twh_ntype(GLB_GET,imolty,ii,ntii)
                        call twh_coordinates(GLB_GET,crdflag,ichain,ii
     &                       ,rxui,ryui,rzui)
                        call twh_luse_c(GLB_GET,imolty,ii,luse_ci)
                        call twh_luse_v(GLB_GET,imolty,ii,luse_vi)
c                       --- loop over all beads jj of chain j 
                        do jj = 1, numjunit
                           call twh_ntype(GLB_GET,jmolty,jj,ntjj)
                           call twh_coordinates(GLB_GET,crdflag,jchain
     &                          ,jj,rxuj,ryuj,rzuj)
                           call twh_luse_c(GLB_GET,jmolty,jj,luse_cj)
                           call twh_luse_v(GLB_GET,jmolty,jj,luse_vj)
c                          --- check if nonzero interaction sites
                           leval_v = luse_vi .and. luse_vj
                           leval_c = luse_ci .and. luse_cj
                           if ( leval_v .or. leval_c ) then
c                             --- minimum image the pair separation
                              call twh_mimage(bxdflag,ibox
     &                             ,rxui,ryui,rzui,rxuj,ryuj,rzuj
     &                             ,rxuij,ryuij,rzuij )
                              rijsq = (rxuij*rxuij)+(ryuij*ryuij)
     &                          + (rzuij*rzuij)
                              if (leval_v) then
                                 vtwointer = vtwointer 
     &                              + twh_vtwobody(lfailure,rijsq,ntii
     &                                ,ntjj,cutstyle,lovrlap,vf,dvdl)
                                 if ( lfailure ) return
                                 dl_twointer = dl_twointer + dvdl
                                 do l = 1, scp_num_foreign_lambda
                                    vftwo(l)=vftwo(l) + vf(l)
                                 enddo

                                 if ( lovrlap ) then
                                    call twh_engtotal_overlap(lfailure
     &                                   ,lmolecmove,ichain,ii,jchain,jj
     &                                   ,rijsq,rxui,ryui,rzui,rxuj,ryuj
     &                                   ,rzuj)
                                    if ( lfailure ) return
                                 endif
                              endif
c                             --- electrostatics.
                              if ( lelect .and. leval_c) then
                                 rij = dsqrt(rijsq)
                                 call twh_vcoulomb(-1,bxdflag,ibox,rij
     &                                ,imolty,ii
     &                                ,jmolty,jj,vf,dvdl,lovrlap
     &                                ,vcoulomb)
                                 velectinter = velectinter + vcoulomb
                                 if ( lovrlap ) then
                                    call twh_engtotal_overlap(lfailure
     &                                   ,lmolecmove,ichain,ii,jchain,jj
     &                                   ,rijsq,rxui,ryui,rzui,rxuj,ryuj
     &                                   ,rzuj)
                                    if ( lfailure ) return
                                 endif
                                 dl_elect = dl_elect + dvdl
                                 do l = 1, scp_num_foreign_lambda
                                    vfelec(l) = vfelec(l) + vf(l)
                                 enddo
                                 if ( lmadelung ) then
                                   if ( rij .lt. maddist ) maddist = rij
                                 endif
                              endif
                           endif
c                          --- loop over all beads jj of chain j 
                        enddo
                     enddo
c                    --- lcompute
                  endif
               endif
            enddo
         endif
      enddo

      do l=1,MAX_FOREIGN_LAMBDA
         vforeign(l) = vfelec(l) + vftwo(l)
      enddo

      return
      end


      subroutine twh_engtotal_tail(vtail,bxdflag,ibox)
c     ******************************************************************
c     * Calculate tail corrections beyond the cutoff                   *
c     * Arguments set by subroutine:                                   *
c     *  vtail -- just the tail correction terms of the intermolecular *
c     *        energy                                                  *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * Error conditions:                                              *
c     *   None                                                         *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 08-04-2016 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision vtail
      integer bxdflag,ibox
c     --- local variables
      integer imolty,jmolty,nmolty,ncmt
      double precision rho,coru,boxvolume

c     --- retrieve constants
      call twh_boxvolume(GLB_GET,bxdflag,ibox,boxvolume)
      vtail = 0.0d0
c     --- add tail corrections beyond the cutoff
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1, nmolty
         do jmolty = 1, nmolty
            call twh_ncmt(GLB_GET,ibox,jmolty,ncmt)
            rho = dble(ncmt) / boxvolume
            call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
            call twh_coru(GLB_GET,imolty,jmolty,coru)
            vtail = vtail + dble(ncmt)*rho*coru
         enddo
      enddo
      return
      end

      subroutine twh_engtotal_tb_sw(lfailure,crdflag,vthrinter
     &     ,bxdflag,ibox)
c     ******************************************************************
c     * calculates the three-body energy contribution                  *
c     * Arguments set by subroutine:                                   *
c     *  vthrinter -- just the 3-body terms of the intermolecular      *
c     *     energy                                                     *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * Error conditions:                                              *
c     *  None                                                          *
c     *                                                                *
c     * split from original engtotal 08-12-2005 by MAW                 *
c     * last modified 03-29-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_VTHREEBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      double precision vthrinter
      integer crdflag,bxdflag,ibox
c     --- local variables
      integer ichain,imolty,ntii,jchain,jmolty,ntjj,kchain,kmolty,ntkk
      integer nchain,testboxi,testboxj,testboxk
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rxuk,ryuk,rzuk,rij,bij,duma,thetac
      double precision rxujk,ryujk,rzujk,bjk,rjk
      logical lzero

      vthrinter = 0.0d0

c     --- stillinger-weber three body potential
c     --- this only works for chains that are monatomic 
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1,nchain
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_ntype(GLB_GET,imolty,1,ntii)
            call twh_coordinates(GLB_GET,crdflag,ichain,1
     &           ,rxui,ryui,rzui)
            do jchain = 1,nchain
               call twh_nboxi(GLB_GET,jchain,testboxj)
               if ( testboxj .eq. ibox .and. ichain .ne. jchain) then
                  call twh_moltyp(GLB_GET,jchain,jmolty)
                  call twh_ntype(GLB_GET,jmolty,1,ntjj)
                  call twh_coordinates(GLB_GET,crdflag,jchain,1
     &                 ,rxuj,ryuj,rzuj)
c                    --- compute i-j bond order
                  call twh_mimage(bxdflag,ibox,rxui,ryui,rzui
     &                 ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
                  rij = twh_distance(rxuij,ryuij,rzuij)
                  call twh_bondorder(lfailure,ntii,ntjj,rij,bij,duma
     &                 ,lzero)
                  if ( lfailure ) return
                  if ( .not. lzero ) then
                     do kchain = 1,nchain
                        call twh_nboxi(GLB_GET,kchain,testboxk)
                        if ( testboxk .eq. ibox 
     &                       .and. ichain .lt. kchain 
     &                       .and. jchain .ne. kchain ) then
                           call twh_moltyp(GLB_GET,kchain,kmolty)
                           call twh_ntype(GLB_GET,kmolty,1,ntkk)
                           call twh_coordinates(GLB_GET,crdflag
     &                          ,kchain,1,rxuk,ryuk,rzuk)
c                             --- compute j-k bondorder
                           call twh_mimage(bxdflag,ibox,rxuk,ryuk,rzuk
     &                          ,rxuj,ryuj,rzuj,rxujk,ryujk,rzujk)
                           rjk = twh_distance(rxujk,ryujk,rzujk)
                           call twh_bondorder(lfailure,ntjj,ntkk,rjk,bjk
     &                          ,duma,lzero)
                           if ( lfailure ) return
                           if ( .not. lzero ) then
c                                --- compute angle
                              thetac = (rxuij*rxujk + ryuij*ryujk 
     &                             + rzuij*rzujk)/(rij*rjk)
c                                --- compute 3-body energy
                              vthrinter = vthrinter + 
     &                             twh_vthreebody(lfailure,rij,rjk
     &                             ,thetac,ntii,ntjj,ntkk)
                              if ( lfailure ) return
                           endif
                        endif
                     enddo
                  endif
               endif
            enddo
         endif
      enddo
      return
      end

      subroutine twh_engtotal_tb_eam(lfailure,crdflag
     &     ,vthrinter,lovrlap,bxdflag, ibox)
c     ******************************************************************
c     * calculates the three-body energy contribution for Embedded     *
c     * Atom Method                                                    *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  vthrinter -- just the 3-body terms of the intermolecular      *
c     *     energy                                                     *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May halt program                                             *
c     *                                                                *
c     * split from original engtotal 08-12-2005 by MAW                 *
c     * last modified 03-29-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_VEMBED
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision vthrinter
      logical lfailure,lovrlap
      integer crdflag,bxdflag,ibox
c     --- local variables
      logical lzero
      integer ichain,imolty,ntii,jchain,jmolty,ntjj,nchain
      integer testboxi,testboxj
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rij,bij,bji,density

      vthrinter = 0.0d0

c     --- Embedded atom potential
c     --- reinitialize all the local densities to zero
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1,nchain
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox ) then
            call twh_eam_rho(GLB_SET,crdflag,ichain,0.0d0)
         endif
      enddo

c        --- need to compute local densities for each atom and then
c        --- use that information to compute the embedding energy
      do ichain = 1,nchain-1
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_ntype(GLB_GET,imolty,1,ntii)
            call twh_coordinates(GLB_GET,crdflag,ichain,1
     &           ,rxui,ryui,rzui)
            do jchain = ichain+1,nchain
               call twh_nboxi(GLB_GET,jchain,testboxj)
               if ( testboxj .eq. ibox .and. ichain .ne. jchain) then
                  call twh_moltyp(GLB_GET,jchain,jmolty)
                  call twh_ntype(GLB_GET,jmolty,1,ntjj)
                  call twh_coordinates(GLB_GET,crdflag,jchain,1
     &                 ,rxuj,ryuj,rzuj)
c                    --- compute i-j bond density
                  call twh_mimage(bxdflag,ibox,rxui,ryui,rzui
     &                 ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
                  rij = twh_distance(rxuij,ryuij,rzuij)
                  call twh_bondorder(lfailure,ntii,ntjj,rij,bij,bji
     &                 ,lzero)
                  if ( lfailure ) return
                  call twh_eam_rho(GLB_INCR,crdflag,ichain,bij)
                  call twh_eam_rho(GLB_INCR,crdflag,jchain,bji)
                  if ( lovrlap ) then
                     write(6,*) 'ENGTOTAL: undershoot in density'
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
         endif
      enddo
c        --- add up the embedding energy
      do ichain = 1,nchain
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_ntype(GLB_GET,imolty,1,ntii)
            call twh_eam_rho(GLB_GET,crdflag,ichain,density)
            vthrinter = vthrinter + twh_vembed(lfailure,ntii,density)
            if ( lfailure ) return
            if ( lovrlap ) then
               write(6,*) 'ENGTOTAL: undershoot in embedding'
               write(6,*) 'density',density
               lfailure = .true.
               return
            endif
         endif
      enddo
      return
      end


      subroutine twh_engtotal_ewald(lfailure,crdflag,bxdflag,ibox
     &     ,vrecipsum,sself,correct)
c     ******************************************************************
c     * calculates the self and correction ewald terms.                *
c     * Note that coulomb potential scaling (scp_scale_coulomb) is not *
c     * yet implemented here                                           *
c     * Ewald corrections do not distinguish between inter- and intra- *
c     * molecular terms.  'Correct' term explicitly excludes real-space*
c     * terms                                                          *
c     *                                                                *
c     * Arguments set by subroutine:                                   *
c     *  vrecipsum -- reciprocal space ewald sum in units [K]          *
c     *  sself -- Coulombic 'self' in units [EE]                       *
c     *  correct -- Coulombic 'correction' in units [EE]               *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     * bxdflag: which box dimensions to use                           *
c     *                                                                *
c     * Error conditions:                                              *
c     *  none                                                          *
c     *                                                                *
c     * split from original engtotal 08-12-2005 by MAW                 *
c     * last modified 03-29-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_EWALD_CORRECT
#define FUNCTION_EWALD_SELF
#define FUNCTION_LINCLUDE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer crdflag,bxdflag,ibox
      double precision vrecipsum, sself, correct
c     --- local variables
      integer ichain,imolty,ii,jj,nchain,numiunit,testboxi
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rij

      sself = 0.0d0
      correct = 0.0d0

      call twh_recipsum(lfailure,bxdflag,ibox,crdflag,vrecipsum)
      if ( lfailure ) return
c        --- update self terms and correction terms
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1,nchain
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,numiunit)
            do ii = 1,numiunit
               sself = sself + twh_ewald_self(bxdflag,ibox,imolty,ii)
               call twh_coordinates(GLB_GET,crdflag,ichain,ii
     &              ,rxui,ryui,rzui)
               do jj = ii+1,numiunit
                  if ( .not. twh_linclude(imolty,ii,jj) ) then
c                       --- compute the correction term assuming
c                       --- these atoms have no coulombic interations
c                       --- with each other.  For those with 1-4 
c                       --- interactions this assumption is fixed in 
c                       --- vtorsion
                     call twh_coordinates(GLB_GET,crdflag,ichain,jj
     &                    ,rxuj,ryuj,rzuj)
                     call twh_mimage(bxdflag,ibox,rxui,ryui,rzui
     &                    ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
                     rij = twh_distance(rxuij,ryuij,rzuij)
                     correct = correct
     &                    + twh_ewald_correct(bxdflag,ibox,imolty,ii,jj
     &                    ,rij)
                  endif
               enddo
            enddo
         endif
      enddo
      
      return
      end


      subroutine twh_engtotal_intrachain(lfailure,crdflag
     &     ,vintra, velectintra
     &     ,lovrlap,cutstyle,bxdflag,ibox, lmolecmove, lmovelist)
c     ******************************************************************
c     * calculates the intrachain interactions.                        *
c     * Arguments set by subroutine:                                   *
c     *  vintra -- the intramolecular (non-coulombic) non-bonded energy*
c     *  velectintra -- the intramolecular component of real space     *
c     *     coulombic potential                                        *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  cutstyle -- cutoff style.  See twh_vtwobody for definition    *
c     *  ibox -- the box we are computing the energy in                *
c     *  lmolecmove -- is true if this call is from a Monte Carlo move *
c     *  lmovelist -- true if we are performing a move and only need   *
c     *            the interactions for a subset of the molecules      *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May exit early (sets lovrlap=true)                           *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_LINCLUDE
#define FUNCTION_VONEFIVE
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision vintra, velectintra
      logical lfailure,lovrlap,lmolecmove,lmovelist
      integer crdflag,cutstyle,bxdflag,ibox
c     --- local variables
      logical lcompute,lmoveme,lelect
      integer ichain,imolty,ii,ntii,jj,ntjj,ioftype,nchain,nunit
      integer testboxi
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision rijsq,rij,dvdl,rminsq,vcoulomb
c     --- foreign energies ignored for intramolecular evaluations
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)

c     --- retrieve constants
      call twh_lelect(GLB_GET,lelect)

      vintra = 0.0d0 
      velectintra = 0.0d0 

c     --- loop over all chains
      call twh_nchain(GLB_GET,nchain)
      do ichain = 1, nchain
c        --- check if i is in relevant box
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox ) then
            if ( lmovelist ) then
               call twh_lmoveme(GLB_GET,ichain,lmoveme)
               if ( lmoveme ) then
                  lcompute = .true.
               else
                  lcompute = .false.
               endif
            else
               lcompute = .true.
            endif
            
            if ( lcompute ) then 
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               do ii = 1, nunit-1
                  call twh_ntype(GLB_GET,imolty,ii,ntii)
                  call twh_coordinates(GLB_GET,crdflag,ichain,ii
     &                 ,rxui,ryui,rzui)
                  do jj = ii+1, nunit
                     if ( twh_linclude(imolty,ii,jj) ) then
                        call twh_ntype(GLB_GET,imolty,jj,ntjj)
                        call twh_coordinates(GLB_GET,crdflag,ichain,jj
     &                       ,rxuj,ryuj,rzuj)
                        call twh_mimage(bxdflag,ibox,rxui,ryui,rzui
     &                       ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij)
                        rijsq = rxuij*rxuij + ryuij*ryuij 
     &                       + rzuij*rzuij
c                       --- get one-five type
                        call twh_onefivetype(imolty,ii,jj,ioftype)
                        if ( ioftype .eq. 0 ) then
c                          --- compute nonbond intramolecular energy
                           vintra = vintra + twh_vtwobody(lfailure,rijsq
     &                          ,ntii,ntjj,cutstyle,lovrlap,vf,dvdl)
                           if ( lfailure ) return
                        else
c                          --- use special one-five terms for energy
                           vintra = vintra + twh_vonefive(lfailure,rijsq
     &                          ,ioftype,cutstyle,lovrlap)
                           if ( lfailure ) return
                        endif
                        if ( lovrlap ) then
                           if ( .not. lmolecmove ) then
                              write(6,*) 'overlap intra'
                              call twh_rminsq(GLB_GET,rminsq)
                              write(6,*) 'rijsq rminsq', rijsq,rminsq
                              write(6,*) 'chain,unit', ichain, ii
                              write(6,*) 'i-pos', rxui,ryui,rzui
                              write(6,*) 'jj', jj
                              write(6,*) 'j-pos',rxuj,ryuj,rzuj
                           endif
                           return
                        endif
                        if ( lelect ) then
c                          --- for volume moves we only need the intra
c                          --- coulomb
                           rij = dsqrt(rijsq)
                           call twh_vcoulomb(-1,bxdflag,ibox
     &                          ,rij,imolty,ii,imolty,jj,vf,dvdl
     &                          ,lovrlap,vcoulomb)
                           velectintra = velectintra + vcoulomb
                           if ( lovrlap ) then
                              if ( .not. lmolecmove ) then
                                 write(6,*) 'overlap intra coulomb'
                                 call twh_rminsq(GLB_GET,rminsq)
                                 write(6,*) 'rijsq rminsq',rijsq,rminsq
                                 write(6,*) 'chain,unit', ichain, ii
                                 write(6,*) 'i-pos', rxui,ryui,rzui
                                 write(6,*) 'jj', jj
                                 write(6,*) 'j-pos',rxuj,ryuj,rzuj
                              endif
                              return
                           endif
                        endif
                     endif
                  enddo
               enddo
            endif
         endif
      enddo
      return
      end


      subroutine twh_engtotal_vibbendtors(lfailure,crdflag
     &     ,vvib,vbend,vtg,vangleangle,vimpr,vbondbond
     &     ,lovrlap,lreport,bxdflag,ibox)
c     ******************************************************************
c     * Calculate vibrational, bending and torsional energies          *
c     * Arguments set by subroutine:                                   *
c     *  vvib -- the bond vibration energy                             *
c     *  vbend -- the bond bending energy                              *
c     *  vtg -- the dihedral angle energy - including improper torsions*
c     *  vangleangle -- angle-angle potential                          *
c     *  vimpr -- the improper torsional energy                        *
c     *  vbondbond -- the IJ-KL bond-bond cross term                   *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *  lreport -- print out all of the energies?                      *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May exit early (sets lovrlap=true)                           *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 03-29-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_VANGANG
#define FUNCTION_VANGLE
#define FUNCTION_VBONBON
#define FUNCTION_VBOND
#define FUNCTION_VIMPROPER
#define FUNCTION_VTORSION
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision vvib, vbend, vtg, vangleangle, vimpr, vbondbond
      logical lfailure,lovrlap, lreport
      integer crdflag,bxdflag,ibox
c     --- local variables
      logical lbadstereo
c     --- integer scalars
      integer ichain,imolty,ii,jjvib,ip1,ip2,it,jjben,ibead,iitor
      integer nchain,nunit,inaa,inbend,inbondbond,inimproper,intorsion
      integer invibration,testboxi
c     --- integer arrays
      integer ttype
      dimension ttype(4)
c     --- double precision scalars
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij
      double precision dist,duma,dumb,dumc
c     --- double precision arrays
      double precision xtmp,ytmp,ztmp,length,angle
      dimension xtmp(4),ytmp(4),ztmp(4),length(4),angle(2)

      vtg = 0.0d0
      vimpr = 0.0d0
      vbend = 0.0d0
      vangleangle = 0.0d0
      vvib = 0.0d0

      call twh_nchain(GLB_GET,nchain)
      do ichain = 1, nchain
         call twh_moltyp(GLB_GET,ichain,imolty)
         call twh_nunit(GLB_GET,imolty,nunit)
c        --- check if i is in relevant box
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ibox ) then
c              --- calculate all bonds vectors and lengths
c              --- calculate all stretching, bending, and torsional 
c              --- potentials that have an end-bead with an index less
c              --- than the current bead
c              --- bond vibrations
            do ii = 2, nunit
               call twh_coordinates(GLB_GET,crdflag,ichain,ii
     &              ,rxui,ryui,rzui)
               call twh_invib(GLB_GET,imolty,ii,invibration)
               do jjvib = 1,invibration
                  call twh_ijvib(GLB_GET,imolty,ii,jjvib,ip1)
                  call twh_itvib(GLB_GET,imolty,ii,jjvib,it)
c                    --- increment vibrational energy
                  if ( ip1 .lt. ii ) then
                     call twh_coordinates(GLB_GET,crdflag,ichain,ip1
     &                    ,rxuj,ryuj,rzuj)
c                       --- compute distance
c                       --- minimum image the bond vectors. this is 
c                       --- required in case the molecule network wraps
c                       --- completely through the periodic box
                     call twh_mimage(bxdflag,ibox,rxui,ryui,rzui
     &                    ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij )
                     dist = twh_distance(rxuij,ryuij,rzuij)
                     vvib = vvib + twh_vbond(lfailure,lreport,bxdflag
     &                    ,ibox,it,dist,imolty,ii,ip1,lovrlap)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        write(6,*) 'overlap vibration'
                        write(6,*) 'chain,molecule,atom A,atom B'
     &                          ,ichain,imolty,ii,ip1
                        write(6,*) 'bond length:',dist
                        write(6,*) 'bond type:',it
                        return
                     endif
                  endif
               enddo
            enddo

c           --- molecule with bending angles
            do ii = 2, nunit
               call twh_coordinates(GLB_GET,crdflag,ichain,ii
     &                 ,xtmp(1),ytmp(1),ztmp(1))
               call twh_inben(GLB_GET,imolty,ii,inbend)
               do jjben = 1, inbend
                  call twh_ijben3(GLB_GET,imolty,ii,jjben,ip2)
                  if ( ip2 .lt. ii ) then
                     call twh_ijben2(GLB_GET,imolty,ii,jjben,ip1)
c                       --- store the coordinates
                     call twh_coordinates(GLB_GET,crdflag,ichain,ip1
     &                       ,xtmp(2),ytmp(2),ztmp(2))
                     call twh_coordinates(GLB_GET,crdflag,ichain,ip2
     &                       ,xtmp(3),ytmp(3),ztmp(3))
                     vbend = vbend + twh_vangle(lfailure,lreport
     &                    ,lovrlap,bxdflag,ibox
     &                    ,imolty,ii,jjben,1,duma,dumb,dumc,xtmp,ytmp
     &                    ,ztmp)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        write(6,*) 'ENGTOTAL: overlap bending'
                        write(6,*) 'imolty,ii,ip1,ip2',imolty,ii
     &                          ,ip1,ip2
                        return
                     endif
                  endif
               enddo
            enddo

c           --- bond bond cross terms
            do ii = 1,nunit
               call twh_inbnbn(GLB_GET,imolty,ii,inbondbond)
               if ( inbondbond .ne. 0 ) then
                  call twh_coordinates(GLB_GET,crdflag,ichain,ii
     &                 ,rxui,ryui,rzui)
                  do jjben = 1,inbondbond
                     call twh_ijbnbn1(GLB_GET,imolty,ii,jjben,ip1)
                     call twh_ijbnbn2(GLB_GET,imolty,ii,jjben,ip2)
c                    --- compute and store first bond length
                     call twh_coordinates(GLB_GET,crdflag,ichain,ip1
     &                    ,rxuj,ryuj,rzuj)
                     call twh_mimage(bxdflag,ibox,rxui,ryui,rzui
     &                    ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij )
                     length(1) = twh_distance(rxuij,ryuij,rzuij)
c                    --- compute and store second bond length
                     call twh_coordinates(GLB_GET,crdflag,ichain,ip2
     &                    ,rxuj,ryuj,rzuj)
                     call twh_mimage(bxdflag,ibox,rxui,ryui,rzui
     &                    ,rxuj,ryuj,rzuj,rxuij,ryuij,rzuij )
                     length(2) = twh_distance(rxuij,ryuij,rzuij)
c                    --- compute energy as a bond vibration term
                     vbondbond = vbondbond+twh_vbonbon(imolty,ii
     &                       ,jjben,length)
                  enddo
               endif
            enddo

c           --- regular torsions 
            do ii = 1, nunit-1
c              --- store the charge on atom ii
               ttype(1) = ii
               call twh_intor(GLB_GET,imolty,ii,intorsion)
               do iitor = 1,intorsion
                  call twh_ijtor4(GLB_GET,imolty,ii,iitor,ttype(4))
                  if ( ttype(4) .gt. ii ) then
c                    --- store unit numbers
                     call twh_ijtor2(GLB_GET,imolty,ii,iitor,ttype(2))
                     call twh_ijtor3(GLB_GET,imolty,ii,iitor,ttype(3))
c                    --- store atom positions
                     do ibead = 1,4
                        call twh_coordinates(GLB_GET,crdflag,ichain
     &                       ,ttype(ibead)
     &                       ,xtmp(ibead),ytmp(ibead),ztmp(ibead))
                     enddo
c                    --- add in the torsional energy
                     vtg = vtg + twh_vtorsion(lfailure,lreport
     &                    ,bxdflag,ibox,1,imolty,ii,iitor
     &                    ,xtmp,ytmp,ztmp,duma,lovrlap)
                     if ( lfailure ) return
                     if ( lovrlap ) then
                        write(6,*) 'overlap torsion'
                        write(6,*) 'imolty,ttype',imolty,ttype(2)
                        return
                     endif
                  endif
               enddo
            enddo

c           --- angle angle cross terms
            do ii = 1,nunit
               ttype(1) = ii
               call twh_inaa(GLB_GET,imolty,ii,inaa)
               do iitor = 1,inaa
                  call twh_ijaa0(GLB_GET,imolty,ii,iitor,ttype(2))
                  call twh_ijaa1(GLB_GET,imolty,ii,iitor,ttype(3))
                  call twh_ijaa2(GLB_GET,imolty,ii,iitor,ttype(4))
                  do ibead = 1,4
                     call twh_coordinates(GLB_GET,crdflag,ichain
     &                    ,ttype(ibead)
     &                    ,xtmp(ibead),ytmp(ibead),ztmp(ibead))
                  enddo
c                    --- compute the angle-angle energy
                  vangleangle = vangleangle 
     &                 + twh_vangang(AA_COORD,bxdflag,ibox,imolty
     &                 ,ii,iitor,angle,xtmp,ytmp,ztmp)
               enddo
            enddo                   

c           --- improper torsions 
c           --- these are only stored on a single atom as they are 
c           --- not symmetrically defined
            do ii = 1,nunit
               ttype(1) = ii
               call twh_inimprop(GLB_GET,imolty,ii,inimproper)
               do iitor = 1,inimproper
                  call twh_ijimprop2(GLB_GET,imolty,ii,iitor,ttype(2))
                  call twh_ijimprop3(GLB_GET,imolty,ii,iitor,ttype(3))
                  call twh_ijimprop4(GLB_GET,imolty,ii,iitor,ttype(4))
c                    --- store the coordinates
                  do ibead = 1,4
                     call twh_coordinates(GLB_GET,crdflag,ichain
     &                    ,ttype(ibead)
     &                       ,xtmp(ibead),ytmp(ibead),ztmp(ibead))
                  enddo
                  vimpr = vimpr + twh_vimproper(bxdflag,ibox,imolty
     &                 ,ii,iitor,xtmp,ytmp,ztmp,duma,lbadstereo)
                  if ( lbadstereo .and. lreport ) then
                     write(6,*) 'ENGTOTAL: Bad stereochemistry '
                     write(6,*) 'moltyp chain unit imptor'
     &                       ,imolty,ichain,ii,iitor
                  endif
               enddo
            enddo
         endif
      enddo
      return
      end

      subroutine twh_engtotal_ext(lfailure,bxdflag,crdflag
     &     ,vext,lovrlap,ibox,lmolecmove
     &     ,lmovelist)
c     ******************************************************************
c     * calculates the interaction energy with external surface        *
c     * Arguments set by subroutine:                                   *
c     *  vext -- the energy of external fields/walls                   *
c     *                                                                *
c     * Arguments modified (but not initialized) by subroutine:        *
c     *  lovrlap -- Atomic overlap found?                              *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *  lmolecmove -- is true if this call is from a Monte Carlo move *
c     *  lmovelist -- true if we are performing a move and only need   *
c     *            the interactions for a subset of the molecules      *
c     *                                                                *
c     * Error conditions:                                              *
c     *   May exit early (sets lovrlap=true)                           *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 03-29-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_VFIELD
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision vext
      logical lfailure,lovrlap, lmolecmove, lmovelist
      integer bxdflag,crdflag,ibox
c     --- local variables
      logical lcompute,lmoveme
      integer ichain,imolty,ii,ntii,nchain,nunit,nfield,testboxi
      double precision cord
      dimension cord(3)

      vext = 0.0d0

      call twh_nfield(GLB_GET,nfield)
      if ( nfield .ne. 0 ) then
         call twh_nchain(GLB_GET,nchain)
         do ichain = 1, nchain
            call twh_nboxi(GLB_GET,ichain,testboxi)
            if ( testboxi .eq. ibox ) then
               if ( lmovelist ) then
                  call twh_lmoveme(GLB_GET,ichain,lmoveme)
                  if ( lmoveme ) then
                     lcompute = .true.
                  else
                     lcompute = .false.
                  endif
               else
                  lcompute = .true.
               endif
               
               if ( lcompute ) then
                  call twh_moltyp(GLB_GET,ichain,imolty)
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do ii = 1, nunit
                     call twh_coordinates(GLB_GET,crdflag,ichain,ii
     &                    ,cord(1),cord(2),cord(3))
                     call twh_ntype(GLB_GET,imolty,ii,ntii)
                     vext = vext + twh_vfield(lfailure,lovrlap
     &                    ,bxdflag,ibox,ntii,imolty,ichain,ii,cord)
                     if ( lfailure ) return
                     if (lovrlap) then
                        if ( .not. lmolecmove ) then
                           write(6,*) 'overlap with field'
                           write(6,*) 'chain,unit',ichain,ii
                           write(6,*) 'x,y,z',cord(1),cord(2),cord(3)
                        endif
                        return
                     endif
                  enddo
               endif
            endif
         enddo
      endif
      return
      end


      subroutine twh_engtotal_solv(lfailure,crdflag,bxdflag,ibox
     &     ,vsolvation)
c     ******************************************************************
c     * calculates the solvation energy for a configuration.           *
c     * Arguments set by subroutine:                                   *
c     *  vsolvation -- solvation energy term                           *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * Error conditions:                                              *
c     *   none                                                         *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 03-29-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_VEEFONE
#define FUNCTION_VSASA
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure 
      double precision vsolvation
      integer bxdflag,ibox
c     --- local variables
      integer ichain,imolty,totunit,isolvtype,nchain,testboxi
      double precision vnt
      integer crdflag,chaininfo
      dimension chaininfo(1,3)

      vsolvation = 0.0d0

      call twh_isolvtype(GLB_GET,isolvtype)
      if ( isolvtype .eq. SOLV_SASA ) then
         call twh_nchain(GLB_GET,nchain)
         do ichain = 1, nchain  
            call twh_nboxi(GLB_GET,ichain,testboxi)
            if ( testboxi .eq. ibox ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,totunit)
               vsolvation = vsolvation + twh_vsasa(ichain,bxdflag,ibox
     &              ,crdflag)
            endif
         enddo
      elseif( isolvtype .eq. SOLV_EEF1 ) then
         call twh_nchain(GLB_GET,nchain)
         do ichain = 1, nchain
            call twh_nboxi(GLB_GET,ichain,testboxi)
            if ( testboxi .eq. ibox ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,totunit)
               vsolvation = vsolvation + twh_veefone(ichain,bxdflag,ibox
     &              ,crdflag)
            endif
         enddo
      elseif ( isolvtype .eq. SOLV_TRAMONTO ) then
         chaininfo(1,1)=-1
         chaininfo(1,2)=-1
         chaininfo(1,3)=-1

         call twh_fielddft(lfailure,0,ibox,1,chaininfo,vnt)
c         call twh_fielddft(lfailure,0,0,chaininfo,ibox,vnt)
         if ( lfailure ) return
         vsolvation = vsolvation + vnt
      endif
      return
      end
     

      subroutine twh_eng_print(nmcount, vvib, vbondbond, vbend, 
     &     vangleangle, vtg, vimpr, vintra, vinter, vtwointer, 
     &     vthrinter, velect, vrealelect, sself, correct, vrecipsum, 
     &     vext, vsolvation, vclassical, vtail, velectinter, 
     &     velectintra, vforeign, dl_inter, dl_elect, 
     &     lmadelung, maddist, ibox)
c     ******************************************************************
c     * print the output                                               *
c     *                                                                *
c     * Arguments set by caller:                                       *
c     *  nmcount -- number molecules in this box                       *
c     *  vvib -- the bond vibration energy                             *
c     *  vbondbond -- the IJ-KL bond-bond cross term                   *
c     *  vbend -- the bond bending energy                              *
c     *  vangleangle -- angle-angle potential                          *
c     *  vtg -- the dihedral angle energy - including improper torsions*
c     *  vimpr -- the improper torsional energy                        *
c     *  vintra -- the intramolecular (non-coulombic) non-bonded energy*
c     *  vinter -- the total intermolecular two-body, three-body and   *
c     *        tail correction non-coulombic energy                    *
c     *  vtwointer -- the two-body non-bond potential                  *
c     *  vthrinter -- just the 3-body terms of the intermolecular      *
c     *        energy                                                  *
c     *  velect -- the total coulombic energy - except 1-4 coulombics  *
c     *        which are included in the torsional potential           *
c     *  vrealelect -- real space component of coulombic potential     *
c     *  sself -- Coulombic 'self'                                     *
c     *  correct -- Coulombic 'correction'                             *
c     *  vrecipsum -- reciprocal space ewald sum                       *
c     *  vext -- the energy of external fields/walls                   *
c     *  vsolvation -- solvation energy term                           *
c     *  vclassical -- the total classical energy                      *
c     *  vtail -- just the tail correction terms of the intermolecular *
c     *        energy                                                  *
c     *  velectinter -- the intermolecular component of real space     *
c     *     coulombic potential                                        *
c     *  velectintra -- the intramolecular component of real space     *
c     *     coulombic potential                                        *
c     *  vforeign -- classical potential at different values of lambda *
c     *     for Scaled Lennard-Jones potential only; see               *
c     *     scaled_potential.F                                         *
c     *  dl_inter -- derivative of two body potential with respect     *
c     *     to lambda_lj, calculated if scp_eval_derivative true       *
c     *  dl_elect -- derivative of inter-molecule coulomb potential    *
c     *     with respect to lambda_lj                                  *
c     *  lmadelung -- compute madelung constant?                       *
c     *  maddist -- distance for Madelung constant                     *
c     *  ibox -- the box we are computing the energy in                *
c     *                                                                *
c     * All energies and derivatives are in units of [K], with the     *
c     * exception of: vrealelect,velectintra,velectinter,sself,correct *
c     * which are in units of [EE].                                    *
c     *                                                                *
c     * Error conditions:                                              *
c     *  None                                                          *
c     *                                                                *
c     * split from original engtotal 08-11-2005 by MAW                 *
c     * last modified 03-29-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_FOREIGN_LAMBDA
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision vvib, vbondbond, vbend, vangleangle, vtg, vimpr  
      double precision vintra, vinter, vtwointer, vthrinter, velect  
      double precision vrealelect, sself, correct, vrecipsum, vext 
      double precision vsolvation, vclassical, maddist, vtail
     &     ,velectinter,velectintra,dl_inter,dl_elect
      double precision vforeign(MAX_FOREIGN_LAMBDA)
      logical lmadelung
      integer nmcount, ibox, l
c     --- local variables
c     --- character string
      character*54 outtext
c     --- logical scalars
      logical scp_eval_derivative,ltailc,lelect,lewald
c     --- integer scalars
      integer imolty,classical_potnum,nmolty,scp_num_foreign_lambda
      integer ncmt
c     --- double precision scalars
      double precision qqfact,vreport,sumvdw,sumchg

c     --- retrieve constants
      call twh_nmolty(GLB_GET,nmolty)
      call twh_ltailc(GLB_GET,ltailc)
      call twh_lelect(GLB_GET,lelect)
      call twh_lewald(GLB_GET,lewald)

      sumvdw = 0.0d0
      sumchg = 0.0d0

      write(6,'(a27,1x,i10)') 'Total molecules in this box', nmcount
      do imolty = 1, nmolty
         call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
         write(6,8) 'Molecules of type',imolty,': ',ncmt
 8       format(a17,1x,i3,1x,a2,i10)
      enddo

c     --- global output format
 111  format(a54,f15.3,' [K] ',f15.5,' [kcal/mol]')

c     --- output classical energies
      outtext = ' Total vibration (atoms ij bonded together)'
      write(6,111) outtext
     &     ,vvib+vbondbond,(vvib+vbondbond)/CON_KCALTOK
      call twh_report(GLB_GET,REP_BON_REG,vreport)
      outtext = '  bond (ij) regular'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_BON_VDW,vreport)
      sumvdw = sumvdw + vreport
      outtext = '  bond (ij) van der waals'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_BON_CHG,vreport)
      sumchg = sumchg + vreport
      outtext = '  bond (ij) Electric'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      outtext = '  bond-bond (ij-jk)'
      write(6,111) outtext,vbondbond,vbondbond/CON_KCALTOK
      outtext = ' Total angle (atoms ijk bonded in sequence)'
      write(6,111) outtext,vbend+vangleangle
     &     ,(vbend+vangleangle)/CON_KCALTOK
      call twh_report(GLB_GET,REP_ANG_REG,vreport)
      outtext = '  angle (ijk) regular'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_ANG_BA,vreport)
      outtext = '  bond-angle (ij-ijk)'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_ANG_VDW,vreport)
      sumvdw = sumvdw + vreport
      outtext = '  angle (ik) van der waals'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_ANG_CHG,vreport)
      sumchg = sumchg + vreport
      outtext = '  angle (ik) Electric'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      outtext = '  angle-angle (ijk-ijm)'
      write(6,111) outtext,vangleangle,vangleangle/CON_KCALTOK
      outtext = ' Total torsion (ijkl bonded in sequence)'
      write(6,111) outtext,vtg,vtg/CON_KCALTOK
      call twh_report(GLB_GET,REP_TOR_REG,vreport)
      outtext = '  torsion (ijkl) regular'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_TOR_BB,vreport)
      outtext = '  torsion (ij-kl) bond-bond'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_TOR_BT,vreport)
      outtext = '  torsion (ij-ijkl) bond-torsion'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_TOR_AT,vreport)
      outtext = '  torsion (ijk-ijkl) angle-torsion'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_TOR_AAT,vreport)
      outtext = '  torsion (ijk-jkl-ijkl) angle-angle-torsion'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_TOR_VDW,vreport)
      sumvdw = sumvdw + vreport
      outtext = '  torsion (i-l) one-four van der waals'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      call twh_report(GLB_GET,REP_TOR_CHG,vreport)
      sumchg = sumchg + vreport
      outtext = '  torsion (i-l) one-four Electric'
      write(6,111) outtext,vreport,vreport/CON_KCALTOK
      outtext = ' Total improper (ijkm, i,k,m bonded to j)'
      write(6,111) outtext,vimpr,vimpr/CON_KCALTOK
      outtext = ' Total nonbond (more than 3 bonds apart) vdw'
      sumvdw = sumvdw + vintra + vinter
      write(6,111) outtext,vintra+vinter,(vintra+vinter)/CON_KCALTOK
      outtext = '  intramolecular'
      write(6,111) outtext,vintra,vintra/CON_KCALTOK
      outtext = '  2-body vdw intermolecular'
      write(6,111) outtext,vtwointer,vtwointer/CON_KCALTOK
      outtext = '  3-body vdw intermolecular'
      write(6,111) outtext,vthrinter,vthrinter/CON_KCALTOK
      if ( ltailc ) then
         outtext = '  tail correction'
         write(6,111) outtext,vtail,vtail/CON_KCALTOK
      endif
      if ( lelect ) then
         outtext =
     &        ' Total nonbond (more than 3 bonds apart) Electric'
         sumchg = sumchg + velect
         write(6,111) outtext,velect,velect/CON_KCALTOK
         outtext = '  total real space Electric'
         write(6,111) outtext,vrealelect,vrealelect/CON_KCALTOK
         outtext = '   intramolecular Electric'
         write(6,111) outtext,velectintra,velectintra/CON_KCALTOK
         outtext = '   intermolecular Electric'
         write(6,111) outtext,velectinter,velectinter/CON_KCALTOK
         if ( lewald ) then
            outtext = '  self'
            write(6,111) outtext,sself,sself/CON_KCALTOK
            outtext = '  correction'
            write(6,111) outtext,correct,correct/CON_KCALTOK
            outtext = '  reciprocal sum'
            write(6,111) outtext,vrecipsum,vrecipsum/CON_KCALTOK
         endif
      endif
      outtext = ' Total external field'
      write(6,111) outtext,vext,vext/CON_KCALTOK
      outtext = ' Total solvation'
      write(6,111) outtext,vsolvation,vsolvation/CON_KCALTOK
      outtext = 'Total classical'
      write(6,111) outtext,vclassical,vclassical/CON_KCALTOK
      outtext = ' Aggregate van der waals (bond and nonbond)'
      write(6,111) outtext,sumvdw,sumvdw/CON_KCALTOK
      outtext = ' Aggregate Electric (bond and nonbond)'
      write(6,111) outtext,sumchg,sumchg/CON_KCALTOK
      call twh_scp_eval_derivative(GLB_GET,scp_eval_derivative)
      if (scp_eval_derivative) then
         write(6,91) dl_inter,dl_inter/CON_KCALTOK
         write(6,92) dl_elect,dl_elect/CON_KCALTOK
 91      format(' dU/d[lambda_lj]',f17.3,' [K] ',f15.5,
     &        ' [kcal/mol]')
 92      format(' dU/d[lambda_c] ',f17.3,' [K] ',f15.5,
     &        ' [kcal/mol]')
      endif
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      if (scp_num_foreign_lambda .ne. 0) then
         write(6,'(1x,a,i1,a)')
     &        "Foreign Lambda energy calculations (box ",ibox,")"
         write(6,901) "index","lambda_lj","lambda_c",
     &        "Total Energy [K]", "Total Energy [kcal/mol]"
         do l = 1,scp_num_foreign_lambda
            write(6,902) l,scp_get_foreign_lambda(LAMBDA_LJ,l),
     &           scp_get_foreign_lambda(LAMBDA_C,l)
     &           ,vforeign(l),vforeign(l)/CON_KCALTOK
         enddo
 901     format(a6,a11,a11,a20,a25)
 902     format(i6,f11.3,f11.3,f20.4,f20.4)
      endif
      
c        --- hack to print out the madelung constant
c        --- only want to do this for a system of sodium chloride
      if ( lmadelung ) then
         call twh_constant_qqfact(GLB_GET,qqfact)
         call twh_ncmt(GLB_GET,ibox,1,ncmt)
         write(6,'(1x,a,f15.10)') 'Madelung Constant [unitless]: '
     &        ,(-velect*maddist)/(qqfact*dble(ncmt))
      endif
      call twh_cp_number(GLB_GET,classical_potnum)
      if ( classical_potnum .eq. CP_EMBEDDED_ATOM_METHOD ) then
c           --- print total energy in eV for EAM potential
         write(6,'(1x,a,1x,g12.6)') 'Total energy [eV]'
     &        ,vclassical/CON_EVTOK
      endif
      return
      end

      subroutine twh_engtotal_overlap(lfailure,lmolecmove,ichain,ii
     &     ,jchain,jj,rijsq,rxui,ryui,rzui,rxuj,ryuj,rzuj)
c     ******************************************************************
c     * handles the overlap error conditions and assigns the proper    *
c     * degree of problem for an overlap that occurs in engtotal       *
c     *                                                                *
c     * split from eng_twobody 07-17-2006 by M.G. Martin               *
c     * last modified 03-20-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lmolecmove
      integer ichain,ii,jchain,jj,ivalue,jvalue
      double precision rijsq,rxui,ryui,rzui,rxuj,ryuj,rzuj,rminsq

      if ( .not. lmolecmove ) then
         write(6,*) 'intermolecular overlap'
         call twh_rminsq(GLB_GET,rminsq)
         write(6,*) 'rijsq rminsq',rijsq,rminsq
         write(6,*) 'rij:',dsqrt(rijsq),' rmin:',dsqrt(rminsq)
         write(6,*) 'chain, unit',ichain,ii
         write(6,*) 'i-pos', rxui,ryui,rzui
         write(6,*) 'jchain jj',jchain,jj
         write(6,*) 'j-pos',rxuj,ryuj,rzuj
         call twh_atomnum(GLB_GET,ichain,ii,ivalue)
         call twh_atomnum(GLB_GET,jchain,jj,jvalue)
         write(6,*) 'getatompos',ivalue,jvalue
         lfailure = .true.
         return
      endif
      return
      end
