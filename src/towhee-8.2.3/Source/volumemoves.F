#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2018 Marcus G. Martin                       *
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
      subroutine twh_volumemoves(lfailure,vmcode)
c     ******************************************************************
c     * control subroutine for all of the volume change moves.         *
c     * vmcode is defined in preproc.h and controls move details       *
c     *                                                                *
c     * originally written 07-14-2016 by M.G. Martin                   *
c     * last modified 04-05-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#include "functions.h"
c     --- debug toggle
#define DEBUG_VOLMOVE 0
c     --- variables passed to/from the subroutine
      logical lfailure
      integer vmcode
c     --- local logicals
      logical lfound
c     --- local integers
      integer boxcode,coordtype,ifind,nboxpair,movevch
      integer numboxes,ibox
c     --- local double precision
      double precision pcheck,rantest
#if DEBUG_VOLMOVE
      write(6,*) 'VOLUMEMOVE start'
#endif
c     --- retrieve constants
      call twh_numboxes(GLB_GET,numboxes)

c     --- determine move type
      if ( vmcode .eq. VMC_NVT_ISO ) then
         boxcode = BCODE_EXCHANGE
         coordtype = VMC_THREED
      elseif ( vmcode .eq. VMC_NVT_ONED ) then
         boxcode = BCODE_EXCHANGE
         coordtype = VMC_ONED
      elseif ( vmcode .eq. VMC_NPT_ISO ) then
         boxcode = BCODE_SINGLE
         coordtype = VMC_THREED
      elseif ( vmcode .eq. VMC_NPT_ONED ) then
         boxcode = BCODE_SINGLE
         coordtype = VMC_ONED
      elseif ( vmcode .eq. VMC_TPRESS ) then
         boxcode = BCODE_GHOST
         coordtype = VMC_THREED
      else
         write(6,*) 'VOLUMEMOVES unknown vmcode:',vmcode
         lfailure = .true.
         return
      endif

c     --- choose boxes
      if ( boxcode .eq. BCODE_SINGLE ) then
c        --- Select a box at  random to change the volume of box
         lfound = .false.
         rantest = twh_random()
         do ifind = 1,numboxes
            if ( coordtype .eq. VMC_THREED ) then
c              --- isotropic volume move, use pmvlpr
               call twh_pmvlpr(GLB_GET,ifind,pcheck)
            elseif ( coordtype .eq. VMC_ONED ) then
c              --- unit cell perturbation, use pmcellpr
               call twh_pmcellpr(GLB_GET,ifind,pcheck)
            else
               write(6,*) 'VOLUMEMOVES: invalid coordtype ',coordtype
               lfailure = .true.
               return
            endif
            if ( .not. lfound .and. rantest .lt. pcheck ) then
               lfound = .true.
               movevch = ifind
            endif
         enddo
      elseif ( boxcode .eq. BCODE_EXCHANGE ) then
         call twh_nboxpair(GLB_GET,nboxpair)
         if ( nboxpair .eq. 0 ) then
            write(6,*) 'VOLUMEMOVE: not enough simulation boxes for NVT'
            lfailure = .true.
            return
         endif
c        --- determine which boxes to perform the volume exchange
         rantest = twh_random()
         lfound = .false.
         do ifind = 1,nboxpair
            if ( coordtype .eq. VMC_THREED ) then
               call twh_pmvlpr(GLB_GET,ifind,pcheck)
            elseif ( coordtype .eq. VMC_ONED ) then
               call twh_pmcellpr(GLB_GET,ifind,pcheck)
            else
               write(6,*) 'VOLUMEMOVES: invalid coordtype ',coordtype
               lfailure = .true.
               return
            endif
            if ( .not. lfound .and. rantest .lt. pcheck ) then
               lfound = .true.
               movevch = ifind + numboxes
            endif
         enddo
      endif
      
      if ( boxcode .eq. BCODE_SINGLE 
     &     .or. boxcode .eq. BCODE_EXCHANGE )then
         call twh_volume_control(lfailure,coordtype
     &        ,movevch,boxcode)
         if ( lfailure ) return
      elseif ( boxcode .eq. BCODE_GHOST ) then
         do ibox = 1, numboxes
            movevch = ibox
            call twh_volume_control(lfailure,coordtype
     &           ,movevch,boxcode)
            if ( lfailure ) return
         enddo
      else
         write(6,*) 'VOLUMEMOVES unknown boxcode:',boxcode
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_volume_control(lfailure,coordtype
     &     ,movevch,boxcode)
c     ******************************************************************
c     * controls the volume move subroutines                           *
c     *                                                                *
c     * originally written 01-17-2018 by M.G. Martin                   *
c     * last modified 04-05-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer coordtype,movevch,boxcode
c     --- local variables
      integer movevector,ivector,jvector

      if ( coordtype .eq. VMC_ONED ) then
c        --- determine which of the 9 quantities to alter
         ivector = idint(twh_random()*3.0d0) + 1
         jvector = idint(twh_random()*3.0d0) + 1
         movevector = 3*(ivector-1) + jvector
      endif
      
      if ( boxcode .ne. BCODE_GHOST ) then
         if ( coordtype .eq. VMC_THREED ) then
c           --- increment the number of attempted 3d volume moves
            call twh_bnvol(GLB_INCR,movevch,1)
         elseif ( coordtype .eq. VMC_ONED ) then
c           --- increment the number of attemped unit cell changes
            call twh_bncell(GLB_INCR,movevch,movevector,1)
         endif
      endif

c     --- call moves
      if ( boxcode .eq. BCODE_SINGLE .or.
     &     boxcode .eq. BCODE_GHOST ) then
c        --- single box
         call twh_volnpt(lfailure,boxcode,coordtype,movevch
     &        ,ivector,jvector,movevector)
         if ( lfailure ) return
      elseif ( boxcode .eq. BCODE_EXCHANGE ) then
c        --- pair of boxes
         call twh_volnvt(lfailure,coordtype,movevch
     &        ,ivector,jvector,movevector)
         if ( lfailure ) return
      else
         write(6,*) 'VOLUMECONTROL: unknown boxcode:',boxcode
         lfailure = .true.
         return
      endif

      return
      end

      subroutine twh_volnvt(lfailure,coordtype
     &     ,movevch,ivector,jvector,movevector)
c     ******************************************************************
c     * coordtype VMC_THREED: makes an isotropic volume change over a  *
c     *   range of ln(vol).  The maximum change is controled by rmvol  *
c     *                                                                *
c     * coordtype VMC_ONED: makes an adjustment of a single hmatrix    *
c     *   element in one of the boxes, and a corresponding isotropic   *
c     *   volume move in the other box.  The maximum change is         *
c     *   controled by rmcell.                                         *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables that are passed to/from the subroutines
      logical lfailure
      integer coordtype,movevch,ivector,jvector,movevector
c     --- local variables
c     --- logical scalars
      logical lovrlap,ladjust,laccept,lreject,lbad,lperiodic
c     --- logical arrays
      logical lbonded
      dimension lbonded(2)
c     --- integer scalars
      integer ibox,imolty,ilabel,ichain,iunit,keycode,cpnumber
      integer idim,jdim,pairvch,numboxes
      integer boxiso,boxhmat,labiso,labhmat,nchain
      integer ivalue,jvalue,nunit,nchbox,testboxi
c     --- double precision scalars
      double precision volt,expdv,dx,dy,dz,vtotal,vvibb
      double precision dele,vinter,vtail,vext,velect,vintra,vvib,vbend
      double precision vbendb,vtorb,boxvolume,rmvol
      double precision vtor,delecurrent
      double precision vsolv,vduma
      double precision xtemp,ytemp,ztemp,value,dvl,dvc
      double precision beta,temperature,xcmi,ycmi,zcmi,pmtest,rmcell
      double precision arbacom,arbbcom,arbccom,dvalue,density
c     --- double precision arrays
      double precision vf(MAX_FOREIGN_LAMBDA)
      double precision hmold
      dimension hmold(2,3,3)
      double precision volo,voln,sfact
      dimension volo(2),voln(2),sfact(2)
      double precision vboxn,vintern,vtailn,vextn,velectn,vintran
      double precision vsolvn
      dimension vboxn(2),vintern(2),vtailn(2),vextn(2),velectn(2)
      dimension vintran(2),vsolvn(2)

c     --- hack
      integer bxdflag
      bxdflag = BXD_REAL
      call twh_numboxes(GLB_GET,numboxes)
      pairvch = movevch - numboxes
      lbonded(1) = .false.
      lbonded(2) = .false.

c     --- retrieve constants
      call twh_nchain(GLB_GET,nchain)
      call twh_cp_number(GLB_GET,cpnumber)

c     --- store old box dimensions, energy, configuration etc
      do ilabel = 1,2
         call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
         do idim = 1,3
            do jdim = 1,3
               call twh_hmatrix(GLB_GET,bxdflag,ibox,idim,jdim
     &              ,hmold(ilabel,idim,jdim))
            enddo
         enddo
         call twh_boxvolume(GLB_GET,bxdflag,ibox,boxvolume)
         volo(ilabel) = boxvolume
c        --- store old k vectors and reciprocal sum
         call twh_recip(lfailure,ibox,RCP_STORE)
         if ( lfailure ) return
      enddo

      do ichain = 1, nchain
         call twh_pairbox(GLB_GET,pairvch,1,ivalue)
         call twh_pairbox(GLB_GET,pairvch,2,jvalue)
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. ivalue .or. testboxi .eq. jvalue ) then
c           --- make sure this molecule is not bonded through the
c           --- periodic boundary
            call twh_lperiodic(GLB_GET,ichain,lperiodic)
            if ( lperiodic ) then
               write(6,*) 'VOLNVT: attempted move in a box that has'
     &              ,' a molecule bonded through the periodic boundary'
               write(6,*) 'Molecule type:',imolty
               write(6,*) 'Chain number:',ichain
               write(6,*) 'Box number:',testboxi
c              --- automatic rejection
               return
            endif
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
c           --- store the old COM in cartesian coords
            call twh_com(GLB_STORE,CRD_R_TO_T,ichain,xcmi,ycmi,zcmi)
c           --- store the old COM in arb coords
            call twh_com(GLB_GET,CRD_REAL,ichain,xtemp,ytemp,ztemp)
            call twh_nboxi(GLB_GET,ichain,ibox)
            call twh_carttoarb(bxdflag,ibox,xtemp,ytemp,ztemp,arbacom
     &           ,arbbcom,arbccom)
            call twh_arbcom(GLB_SET,ichain,arbacom,arbbcom,arbccom)
c           --- store the old cartesian coords for each unit
            do iunit = 1, nunit
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
            enddo
            if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c              --- store the old embedding density
               call twh_eam_rho(GLB_STORE,CRD_R_TO_T,ichain,density)
            endif
         endif
      enddo

c     --- calculate total volume (using label indexes)
      volt = volo(1) + volo(2)

c     --- initialize lreject, this is used to let us know if one of the
c     --- box dimensions, or box volumes has gone somewhere we will 
c     --- automatically reject.  If so then we skip all the energy 
c     --- computataions and just reset the box dimensions
      lreject = .false.
      if ( coordtype .eq. VMC_THREED ) then
c        --- calculate new volume
         call twh_rmvol(GLB_GET,movevch,rmvol)
         expdv = twh_expon(dlog(volo(1)/volo(2))
     &        + rmvol*(2.0d0*twh_random()-1.0d0))
         voln(1)= expdv*volt/(1+expdv)
         voln(2)= volt-voln(1)
         if ( voln(1) .lt. 0.0d0 .or. voln(2) .lt. 0.0d0 ) then
            write(6,*) 'VOLNVT: a trial volume went negative'
            write(6,*) 'should not happen in the NVT version of volume'
            lfailure = .true.
            return
         endif
         sfact(1) = (voln(1)/volo(1))**(1.0d0/3.0d0)
         sfact(2) = (voln(2)/volo(2))**(1.0d0/3.0d0)

c        --- adjust the box dimensions of the two boxes
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_SCALE,bxdflag,ibox,idim,jdim
     &                 ,sfact(ilabel))
               enddo
            enddo
         enddo
c        --- set the new hinverse and boxvolume for the two boxes
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
            call twh_inverthmatrix(bxdflag,ibox)
         enddo
      elseif ( coordtype .eq. VMC_ONED ) then
c        --- we simply adjust a single box dimension for one of the 
c        --- boxes, and cause a corresponding volume change 
c        --- isotropically in the other box
c        --- randomly select which box to change first
         call twh_pmcellpt(GLB_GET,pairvch,pmtest)
         if ( twh_random() .lt. pmtest ) then
            labhmat = 1
            labiso = 2
         else
            labhmat = 2
            labiso = 1
         endif
c        --- get the box numbers
         call twh_pairbox(GLB_GET,pairvch,labhmat,boxhmat)
         call twh_pairbox(GLB_GET,pairvch,labiso,boxiso)
c        --- modify a single vector of the hmatrix in boxhmat
         call twh_rmcell(GLB_GET,boxhmat,ivector,jvector,rmcell)
         call twh_hmatrix(GLB_GET,bxdflag,boxhmat,ivector,jvector
     &        ,dvalue)
         dvalue = dvalue + (2.0d0*twh_random()-1.0d0)*rmcell
         call twh_hmatrix(GLB_SET,bxdflag,boxhmat,ivector,jvector
     &        ,dvalue)
c        --- compute the new volume for this box
         call twh_inverthmatrix(bxdflag,boxhmat)
c        --- make sure the new volume is still positive
         call twh_boxvolume(GLB_GET,bxdflag,boxhmat,boxvolume)
         if ( boxvolume .le. 0.0d0 ) then
            lreject = .true.
         else
c           --- compute the volume change and transfer it to the other
c           --- box
            voln(labhmat) = boxvolume
            voln(labiso) = volo(labiso) - (voln(labhmat)-volo(labhmat))
c           --- make sure this new volume is still positive
            if ( voln(labiso) .lt. 0.0d0 ) then
               lreject = .true.
            else
c              --- compute the isotropic scaling factor for hmatrix
               sfact(labiso) =(voln(labiso)/volo(labiso))**(1.0d0/3.0d0)
c              --- scale the hmatrix
               do idim = 1,3
                  do jdim = 1,3
                     call twh_hmatrix(GLB_SCALE,bxdflag,boxiso,idim,jdim
     &                    ,sfact(labiso)) 
                  enddo
               enddo
c              --- compute the new hinverse
               call twh_inverthmatrix(bxdflag,boxiso)
            endif
         endif
      endif

      do ilabel = 1,2
         call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
c        --- make sure the box dimensions are still valid for a right
c        --- handed coordinate system, and for rcut
         call twh_checkhmatrix(bxdflag,ibox,.false.,lbad)
         if ( lbad ) lreject = .true.
      enddo

      if ( lreject ) then
         if ( coordtype .eq. VMC_THREED ) then 
            write(6,*) 'VOLNVT: move rejected, trial volume too low'
         else
            write(6,*) 'VOLNVT: cell move rejected, '
     &           ,'invalid trial hmatrix or volume'
         endif
c        --- restore old hmatrix
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_SET,bxdflag,ibox,idim,jdim
     &                 ,hmold(ilabel,idim,jdim))
               enddo
            enddo
c           --- reset the hinverse and boxvolume
            call twh_inverthmatrix(bxdflag,ibox)
         enddo
         return
      endif

c     --- determine new positions of the molecules
c     --- calculate centre of mass and its displacement
      do ichain = 1, nchain
         call twh_nboxi(GLB_GET,ichain,ibox)
         ladjust = .false.
         call twh_pairbox(GLB_GET,pairvch,1,ivalue)
         call twh_pairbox(GLB_GET,pairvch,2,jvalue)
         if ( ibox .eq. ivalue ) then
            ilabel = 1
            ladjust = .true.
         elseif (ibox .eq. jvalue ) then
            ilabel = 2
            ladjust = .true.
         endif
         if ( ladjust ) then
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
c           --- transform the old COM in arbitrary coordinates into
c           --- the new COM in cartesian coordinates
            call twh_arbcom(GLB_GET,ichain,arbacom,arbbcom,arbccom)
            call twh_arbtocart(bxdflag,ibox,arbacom,arbbcom,arbccom
     &           ,xtemp,ytemp,ztemp)
c           --- compute the COM displacement
            call twh_com(GLB_GET,CRD_TEMP,ichain,xcmi,ycmi,zcmi)
            dx = xtemp - xcmi
            dy = ytemp - ycmi
            dz = ztemp - zcmi
c           --- increment the difference for the new COM 
            call twh_com(GLB_INCR,CRD_REAL,ichain,dx,dy,dz)
c           --- adjust all of the unit coordinates
            do iunit = 1, nunit
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               xtemp = xtemp + dx
               ytemp = ytemp + dy
               ztemp = ztemp + dz
               call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
            enddo
         endif
      enddo

c     --- compute the energy change
      dele = 0.0d0
      do ilabel = 1,2
         call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
c        --- set coulombic parameters for this new box length
         call twh_resetcoulomb(lfailure,bxdflag,ibox)
         if ( lfailure ) return
c        --- calculate the new intermolecular energies
c        --- note that the intramolecular nonbonded term could change
c        --- if the molecule is large and changes from interacting 
c        --- through real space to through an image
         call twh_engtotal(lfailure,lovrlap,vtotal,vinter,vduma,vtail
     &        ,vintra,vvib,vbend,vtor,vext,velect,vsolv,vf,dvl
     &        ,dvc
     &        ,CRD_REAL,bxdflag,ibox
     &        ,.true.,.false.,.false.,lbonded(ilabel))
         if ( lfailure ) return
         if ( lovrlap ) then
            laccept = .false.
            goto 500
         endif
         vintern(ilabel) = vinter
         vtailn(ilabel)  = vtail
         vextn(ilabel)   = vext  
         velectn(ilabel) = velect
         vintran(ilabel) = vintra
         vsolvn(ilabel)  = vsolv
c        --- did not recompute bonded energies
         call twh_boxvbond(GLB_GET,ibox,vvibb)
         call twh_boxvangle(GLB_GET,ibox,vbendb)
         call twh_boxvtorsion(GLB_GET,ibox,vtorb)
         vboxn(ilabel)   = vtotal + vtorb + vbendb + vvibb
         call scp_set_derivative(ilabel, LAMBDA_C, dvc)
         call scp_set_derivative(ilabel, LAMBDA_LJ, dvl)
         call scp_set_foreign_energies(ilabel, vf)
c        --- compute current energy change for this box
         call twh_boxvclassic(GLB_GET,ibox,vtotal)
         delecurrent = vboxn(ilabel) - vtotal
c        --- increment total dele
         dele = dele + delecurrent

         if ( coordtype .eq. VMC_THREED ) then
c           --- store dU/dV in the averages
            keycode = AC_DU_DV
            value = delecurrent /(voln(ilabel)-volo(ilabel))
            call twh_averages(lfailure,AVG_INCREMENT,keycode,ibox,0
     &           ,value,1)
            if ( lfailure ) return
         endif
      enddo

      call twh_temperature(GLB_GET,temperature)
      if ( coordtype .eq. VMC_THREED ) then
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
            call twh_nchbox(GLB_GET,ibox,nchbox)
            dele = dele - (dble(nchbox+1)*dlog(voln(ilabel)
     &           /volo(ilabel))*temperature)
         enddo
      elseif ( coordtype .eq. VMC_ONED ) then
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
            call twh_nchbox(GLB_GET,ibox,nchbox)
            dele = dele - (dble(nchbox)*dlog(voln(ilabel)/volo(ilabel))
     &           *temperature )
         enddo
      endif

c     --- acceptance test
      laccept = .false.
      if ( dele .lt. SMALLEST ) then
c        --- automatic acceptance
         laccept = .true.
      else
         call twh_beta(GLB_GET,beta)
         if ( twh_random() .lt. twh_expon(-beta*dele) ) then
            laccept = .true.
         endif
      endif

 500  if ( laccept ) then
#if DEBUG_VOLMOVE
         write(6,*) 'volume move accepted',dele
#endif
c        --- move accepted
         if ( coordtype .eq. VMC_THREED ) then
            call twh_bavol(GLB_INCR,movevch,1)
         elseif ( coordtype .eq. VMC_ONED ) then
            call twh_bacell(GLB_INCR,movevch,movevector,1)
         endif
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
c           --- update the energies
            call twh_boxvclassic(GLB_SET,ibox,vboxn(ilabel))
            call twh_boxvinter(GLB_SET,ibox,vintern(ilabel))
            call twh_boxvintra(GLB_SET,ibox,vintran(ilabel))
            call twh_boxvtail(GLB_SET,ibox,vtailn(ilabel))
            call twh_boxvfield(GLB_SET,ibox,vextn(ilabel))
            call twh_boxvelect(GLB_SET,ibox,velectn(ilabel))
            call twh_boxvsolvation(GLB_SET,ibox,vsolvn(ilabel))
         enddo
      else
c        --- volume move was rejected: restore old box lengths
         do ilabel = 1,2
            call twh_pairbox(GLB_GET,pairvch,ilabel,ibox)
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_SET,bxdflag,ibox,idim,jdim
     &                 ,hmold(ilabel,idim,jdim))
               enddo
            enddo
c           --- reset old volume and hinverse
            call twh_inverthmatrix(bxdflag,ibox)

            call twh_resetcoulomb(lfailure,bxdflag,ibox)
            if ( lfailure ) return
            call twh_recip(lfailure,ibox,RCP_RESTORE)
            if ( lfailure ) return
         enddo

c        --- restore old coordinates
         do ichain = 1, nchain
            call twh_pairbox(GLB_GET,pairvch,1,ivalue)
            call twh_pairbox(GLB_GET,pairvch,2,jvalue)
            call twh_nboxi(GLB_GET,ichain,testboxi)
            if ( testboxi .eq. ivalue .or. testboxi .eq. jvalue ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               call twh_com(GLB_STORE,CRD_T_TO_R,ichain,xcmi,ycmi,zcmi)
               do iunit = 1, nunit
                  call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
                  call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
               enddo
               if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c                 --- restore the old embedding density
                  call twh_eam_rho(GLB_STORE,CRD_T_TO_R,ichain,density)
               endif
            endif
         enddo
      endif

      return
      end

      subroutine twh_volnpt(lfailure,boxcode,coordtype,boxvch
     &        ,ivector,jvector,movevector)
c     ******************************************************************
c     * coordtype: VMC_THREED                                          *
c     *    makes an isotropic volume change under const. pressure      *
c     *    the maximum change is controlled by rmvol and the           *
c     *    number of successful trial moves is stored in bsvol         *
c     * coordtype: VMC_ONED                                            *
c     *    changes a single element of the hmatrix which controls the  *
c     *    unit cell size and shape.  Move is performed on a single    *
c     *    simulation box under constant pressure                      *
c     *    this move was first described by Yashonath and Rao in       *
c     *    Mol. Phys. 54, 1, 245-251 (1985).                           *
c     *                                                                *
c     * rewritten from old volnpt 01-22-2018 by M.G. Martin            *
c     * last modified 03-30-2020 by M.G.Martin                         *
c     ******************************************************************
      implicit none
c     --- functions
#define FUNCTION_EXPON
#define FUNCTION_RANDOM
#include "functions.h"
c     --- constants
#define NEWVOL_RMVOL 0
#define NEWVOL_RANGE 1
c     --- variables passed to/from the subroutine
      logical lfailure
      integer boxcode,coordtype,boxvch,ivector,jvector,movevector
c     --- local variables
c     --- logical scalars
      logical lovrlap,lvol,laccept,lreject,lbad
c     --- integer scalars
      integer ichain,iunit,imolty,cpnumber,newcode
      integer idim,jdim,nchain
      integer nunit,nchbox,keycode
      integer nfield,ifield,fieldbox,testboxi
c     --- double precision scalars
      double precision volo,deltav,sfact,vinterb,vintrab,rmvol
      double precision voln,vintern,vtailn,vextn,velectn,vintran
      double precision vsolvn,vextb,velectb,vsolvb
      double precision dx,dy,dz,v,dele,vduma
      double precision xtemp,ytemp,ztemp,external_pressure
      double precision vvibn,vbenddum,vtorn
      double precision dvl,dvc,beta,temperature,xcmi,ycmi,zcmi,rmcell
      double precision arbacom,arbbcom,arbccom,dvalue,density
      double precision rangelow,rangehigh,deltavolume
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)
      character*20 pressure_thermo_style
      
c     --- retrieve constants
      call twh_nchain(GLB_GET,nchain)
      call twh_cp_number(GLB_GET,cpnumber)
c     --- safety check for external fields in this box
      call twh_nfield(GLB_GET,nfield)
      do ifield = 1,nfield
         call twh_fieldint(GLB_GET,ifield,1,fieldbox)
         if ( fieldbox .eq. boxvch ) then
c           --- automatic rejection, print warning
            write(6,*) 'Cannot perform a volume change move in a box'
     &           ,' that has an external field'
            write(6,*) 'Automatic rejection in box:',boxvch
            return
         endif
      enddo
c     --- get the old volume
      call twh_boxvolume(GLB_GET,BXD_REAL,boxvch,volo)
      lreject = .false.
      if ( coordtype .eq. VMC_THREED ) then
         if ( boxcode .eq. BCODE_GHOST ) then
c           --- check the pressure_thermo_style
            call twh_pressure_thermo_style(GLB_GET
     &           ,pressure_thermo_style)
            if ( pressure_thermo_style .eq. 'move data' ) then
               newcode = NEWVOL_RMVOL
            elseif ( pressure_thermo_style .eq. 'range' ) then
               newcode = NEWVOL_RANGE
            else
               write(6,*) 'VOLNPT: unknown pressure_thermo_style'
               write(6,*) pressure_thermo_style
               lfailure = .true.
               return
            endif
         else
c           --- regular moves use rmvol
            newcode = NEWVOL_RMVOL
         endif
         if ( newcode .eq. NEWVOL_RMVOL ) then
c        --- calculate new volume         
            call twh_rmvol(GLB_GET,boxvch,rmvol)
            voln = volo + rmvol*( 2.0d0*twh_random() - 1.0d0 )
         elseif ( newcode .eq.NEWVOL_RANGE ) then
            call twh_pressure_thermo_range(GLB_GET,1,rangelow)
            call twh_pressure_thermo_range(GLB_GET,2,rangehigh)
            voln = volo + volo*(
     &           (rangehigh-rangelow)*twh_random() + rangelow )
         else
            write(6,*) 'VOLNPT: unknown newcode:',newcode
            lfailure = .true.
            return
         endif
c        --- make sure the new volume is positive
         if ( voln .lt. 0.0d0 ) then
            lreject = .true.
         else
c           --- isotropic volume change
            sfact = (voln/volo)**(1.0d0/3.0d0)
c           --- set new box dimensions
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_SCALE,BXD_R_TO_T,boxvch,idim,jdim
     &                 ,sfact)
               enddo
            enddo
         endif
      elseif ( coordtype .eq. VMC_ONED ) then
c        --- store the hmatrix in temp
         call twh_hmatrix(GLB_STORE_ALL,BXD_R_TO_T,boxvch,3,3,dvalue)
c        --- modify a single vector of the hmatrix
         call twh_rmcell(GLB_GET,boxvch,ivector,jvector,rmcell)
         call twh_hmatrix(GLB_GET,BXD_TEMP,boxvch,ivector,jvector
     &        ,dvalue)
         dvalue = dvalue + (2.0d0*twh_random()-1.0d0)*rmcell
         call twh_hmatrix(GLB_SET,BXD_TEMP,boxvch,ivector,jvector
     &        ,dvalue)
      endif
c     --- check the new box dimensions
      call twh_checkhmatrix(BXD_TEMP,boxvch,.false.,lbad)
      if ( lbad ) lreject = .true.

      if ( .not. lreject ) then
c        --- set the new boxvolume and hinverse
         call twh_inverthmatrix(BXD_TEMP,boxvch)
c        --- make sure the volume is still positive
         call twh_boxvolume(GLB_GET,BXD_TEMP,boxvch,voln)
         if ( voln .lt. 0.0d0 ) then
            lreject = .true.
         endif
      endif
#if DEBUG_VOLMOVE
      write(6,*) 'voln,volo:',voln,volo
#endif
      if ( lreject ) then
         write(6,*) 'VOLNPT: move rejected, negative volume or bad '
     &        ,'trial hmatrix'
         return
      endif
c     --- determine new positions of the molecules
c     --- calculate centre of mass and its displacement
      do ichain = 1, nchain
         call twh_moltyp(GLB_GET,ichain,imolty)
         call twh_nunit(GLB_GET,imolty,nunit)
         call twh_nboxi(GLB_GET,ichain,testboxi)
         if ( testboxi .eq. boxvch) then
c           --- compute the shift from the old to the new coords
            call twh_com(GLB_GET,CRD_REAL,ichain,xcmi,ycmi,zcmi)
c           --- convert this to arbitrary coordinates
            call twh_carttoarb(BXD_REAL,boxvch,xcmi,ycmi,zcmi
     &           ,arbacom,arbbcom,arbccom)
c           --- convert back to cart using the new box settings
            call twh_arbtocart(BXD_TEMP,boxvch,arbacom,arbbcom,arbccom
     &           ,xtemp,ytemp,ztemp)
c           --- store the new COM
            call twh_com(GLB_SET,CRD_TEMP,ichain,xtemp,ytemp,ztemp)
c           --- compute shift for this molecule
            dx = xtemp - xcmi
            dy = ytemp - ycmi
            dz = ztemp - zcmi
c           --- adjust all of the unit coordinates
            do iunit = 1, nunit
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
               xtemp = xtemp + dx
               ytemp = ytemp + dy
               ztemp = ztemp + dz
               call twh_coordinates(GLB_SET,CRD_TEMP,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
            enddo
         endif
      enddo
      lvol = .true.
c     --- set new coulombic parameters for the new box length
      call twh_resetcoulomb(lfailure,BXD_TEMP,boxvch)
      if ( lfailure ) return
      call twh_engtotal(lfailure,lovrlap,v,vintern,vduma,vtailn,vintran
     &     ,vvibn,vbenddum,vtorn,vextn,velectn,vsolvn,vf,dvl,dvc
     &     ,CRD_TEMP,BXD_TEMP,boxvch,lvol,.false.,.false.,.false.)
      if ( lfailure ) return
#if DEBUG_VOLMOVE
      write(6,*) 'v:',v
      write(6,*) 'vintern:',vintern
      write(6,*) 'vtailn:',vtailn
      write(6,*) 'vintran:',vintran
      write(6,*) 'vvibn:',vvibn
      write(6,*) 'vtorn:',vtorn
      write(6,*) 'vextn:',vextn
      write(6,*) 'velectn:',velectn
      write(6,*) 'vsolvn:',vsolvn
#endif
      if ( lovrlap ) then
#if DEBUG_VOLMOVE
         write(6,*) 'move rejected due to overlap in VOLNPT'
#endif
c        --- hard overlap: automatic reject
         return
      endif
c     --- initialize laccept
      laccept = .false.
c     -- retrieve the old total box energy from accumulators
      call twh_boxvinter(GLB_GET,boxvch,vinterb)
      call twh_boxvintra(GLB_GET,boxvch,vintrab)
      call twh_boxvfield(GLB_GET,boxvch,vextb)
      call twh_boxvelect(GLB_GET,boxvch,velectb)
      call twh_boxvsolvation(GLB_GET,boxvch,vsolvb)
c     --- compute energy change for move
      deltav = (vintern - vinterb) + (vextn  - vextb)
     &     + (velectn - velectb) + (vintran - vintrab)
     &     + (vsolvn - vsolvb)
      if ( coordtype .eq. VMC_THREED ) then
c        --- store dU/dV in the averages
         keycode = AC_DU_DV
         deltavolume = voln - volo
         dvalue = deltav/deltavolume
         call twh_averages(lfailure,AVG_INCREMENT,keycode,boxvch,0
     &        ,dvalue,1)
         if ( lfailure ) return
c        --- split the averaging depending on deltavolume
         if ( deltavolume .lt. 0.0d0 ) then
            keycode = AC_DU_DV_CONTRACT
         else
            keycode = AC_DU_DV_EXPAND
         endif
         call twh_averages(lfailure,AVG_INCREMENT,keycode,boxvch,0
     &        ,dvalue,1)
      endif
      if ( boxcode .eq. BCODE_GHOST ) then
c        --- ghost volume change, du/dv is computed, we are done
         return
      elseif ( boxcode .eq. BCODE_SINGLE ) then
c        --- isobaric-isothermal ensemble move
c        --- add in the number density and pressure terms
         call twh_express(GLB_GET,external_pressure)
         call twh_nchbox(GLB_GET,boxvch,nchbox)
         call twh_temperature(GLB_GET,temperature)
         dele = deltav + external_pressure*(voln - volo) 
     &        - (dble(nchbox)*dlog(voln/volo)*temperature )
c        --- acceptance test 
         if ( dele .lt. SMALLEST ) then
c           --- accept move
            laccept = .true.
         else
            call twh_beta(GLB_GET,beta)
            if ( twh_random() .lt. twh_expon(-(beta*dele)) ) then
               laccept = .true.
            endif
         endif
      endif
#if DEBUG_VOLMOVE
         write(6,*) 'dele,laccept:',dele,laccept
         write(6,*) 'deltav:',deltav
#endif
      if ( laccept ) then
         if ( coordtype .eq. VMC_THREED ) then
c           --- increment the number of accepted volume moves
            call twh_bavol(GLB_INCR,boxvch,1)
         elseif ( coordtype .eq. VMC_ONED ) then
c           --- increment the number of accepted unit cell moves
            call twh_bacell(GLB_INCR,boxvch,movevector,1)
         endif
c        --- store the new energies
         call twh_boxvclassic(GLB_INCR,boxvch,deltav)
         call twh_boxvinter(GLB_SET,boxvch,vintern)
         call twh_boxvintra(GLB_SET,boxvch,vintran)
         call twh_boxvtail(GLB_SET,boxvch,vtailn)
         call twh_boxvfield(GLB_SET,boxvch,vextn)
         call twh_boxvelect(GLB_SET,boxvch,velectn)
         call twh_boxvsolvation(GLB_SET,boxvch,vsolvn)
         call scp_set_derivative(boxvch, LAMBDA_C, dvc)
         call scp_set_derivative(boxvch, LAMBDA_LJ, dvl)
         call scp_set_foreign_energies(boxvch, vf)
c        --- store new box dimensions
         call twh_hmatrix(GLB_STORE_ALL,BXD_T_TO_R,boxvch,3,3,dvalue)
c        --- store the inverse hmatrix information
         call twh_inverthmatrix(BXD_T_TO_R,boxvch)
c        --- store new electrostatic information
         call twh_resetcoulomb(lfailure,BXD_T_TO_R,boxvch)
         if ( lfailure ) return
         call twh_recip(lfailure,boxvch,RCP_ACCEPT_NEW_VOL)
         if ( lfailure ) return
c        --- store temp coordinates as real
         do ichain = 1, nchain
            call twh_nboxi(GLB_GET,ichain,testboxi)
            if ( testboxi .eq. boxvch) then
               call twh_com(GLB_STORE,CRD_T_TO_R,ichain,xcmi,ycmi,zcmi)
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               do iunit = 1, nunit
                  call twh_coordinates(GLB_GET,CRD_TEMP,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
                  call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
               enddo
               if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c                 --- store the current eam_rho
                  call twh_eam_rho(GLB_STORE,CRD_T_TO_R,ichain,density)
               endif
            endif
         enddo
      endif

      return
      end
