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
      subroutine twh_setmixterms(lfailure,ntypes)
c     ******************************************************************
c     * sets up all of the nonbonded force field cross terms with the  *
c     * appropriate constants                                          *
c     * also converts degrees into radians for the bending terms       *
c     *                                                                *
c     * originally written 03-01-2002 by M.G. Martin                   *
c     * last modified 03-31-2020 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_GETNBTYPE
#define FUNCTION_VTWOBODY
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
c     --- logical scalars
      logical lovrlap,lproblem,lhbond,lfound,lshift,loutput
c     --- character string scalars
      character*1 dai,daj
      character*(ELENAMLEN) twoelei,twoelej
      character*(FFNAMELEN) testone,testtwo,elei,elej,key
      character*20 interpolatestyle,eampairstyle,eamembedstyle
      character*20 eamdensstyle
      character*(CLAPOTLEN) localcp
      character*(MIXRULLEN) localcm
      character*(FFNAMELEN) nbnamei,nbnamej
c     --- character string arrays
      character*(FFNAMELEN) hbondnames
      dimension hbondnames(2)
c     --- integer scalars
      integer idist,imax,itype,jtype,ii,jj,ij,itrial,jlo,jhi
      integer hbtype,id,nadjust,iadjust,icoeff,itest,style
      integer anglestyle,cpnumber
c     --- double precision scalars
      double precision maxval,test,dist,rpos,dv,dtemp,coefftemp
      double precision esigsix,esigtwelve,polari,polarj,denom
      double precision valuei,valuej,twopi
      double precision ai,aj,alphai,alphaj,gi,gj,nni,nnj,power,bbb
      double precision epsscale,radscale,gamma,nbbeta,rii,rjj,rij,rcut
      double precision dvibone,dvibtwo,dvibthree,rijfinal,epsfinal
      double precision nbc_1_ii,nbc_1_jj,nbc_2_ii,nbc_2_jj
      double precision nbc_1id_ii,nbc_2id_ii,nbc_1id_jj,nbc_2id_jj
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)
      double precision nbcoeff
      dimension nbcoeff(MINNBCOEFF:MAXNBCOEFF)
      
c     --- get the classical potential
      call twh_classical_potential(GLB_GET,localcp)
      call twh_cp_number(GLB_GET,cpnumber)
      call twh_classical_mixrule(GLB_GET,localcm)
c     --- retrieve constants
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_lshift(GLB_GET,lshift)
      call twh_rcut(GLB_GET,rcut)

c     --- initialize parameters
      lovrlap = .false.
      lproblem = .false.
      lhbond = .false.
 111  format(1x,a)
c     --- output method in which cross terms are computed
      if ( localcm .eq. 'Lorentz-Berthelot' ) then
         write(6,111) 'Lorentz-Berthelot Mixing rules'
         write(6,111) 'Arithmetic mean of sigma terms'
         write(6,111) 'Geometric mean of epsilon term'
      elseif ( localcm .eq. 'Geometric' ) then
         write(6,111) 'Geometric Mixing rules'
         write(6,111) 'Geometric mean of all nonbond coefficients'
      elseif ( localcm .eq. 'Sixth Power' ) then
         write(6,111) 'Sixth Power Mixing rules'
         write(6,111) 'Sixth order combination of sigma'
         write(6,111) 'Sixth order combination of epsilon'
      elseif ( localcm .eq. 'Explicit' ) then
         call twh_loutput(GLB_GET,loutput)
         if ( loutput ) then
            write(6,111) 'Explicit declaration of Cross Terms'
            write(6,111) 'Please check the forcefield file carefully'
         endif
      elseif ( localcm .eq. 'Arithmetic' ) then
         write(6,111) 'Arithmetic mean of all nonbond parameters'
      elseif ( localcm .eq. 'Shukla' ) then
         write(6,111) 'Mixing rules from Equations 19 and 20 of'
         write(6,111) 'K.P. Shukla Mol. Phys. 62, 5, 1143-1163 (1987)'
         write(6,111) 'with cross terms adjusted by polarizabilities'
      elseif ( localcm .eq. 'LB plus manual' ) then
         write(6,'(1x,a,1x,a)')
     &        'Lorentz-Berthelot Mixing rules plus manual setting'
     &        ,'of certain specified cross terms'
         write(6,111) 'Arithmetic mean of the other sigma terms'
         write(6,111) 'Geometric mean of the other epsilon term'
      elseif ( localcm .eq. 'MMFF' ) then
         write(6,'(1x,a,1x,a)')
     &        'MMFF mixing rules from Equations 2-5 of'
     &        ,'T.A. Halgren J. Comp. Chem. 17, 5&6, 52-552 (1996)'
         write(6,111) 'Augmented Arithmetic mean of distances'
         write(6,111) 'Slater-Kirkwood combination of epsilons'
         write(6,111) 'using B=0.2, beta=12 in most cases'
         write(6,111) 'terms involving polar hydrogens use B=0.0'
         write(6,'(1x,a,1x,a,1x,a)')
     &        'terms involving a donor/acceptor pair subsequently'
     &        ,'scale the radius by DARAD=0.8'
     &        ,'and epsilon by DAEPS=0.5'
      else
         write(6,111) 'SETMIXTERMS: undefined mixing rule style'
         write(6,111) localcm
         lfailure = .true.
         return
      endif

c     --- classical_potential specific section
      if ( cpnumber .eq. CP_EXPONENTIAL_6 ) then
c        --- need to avoid some divide by zero problems.  If term 3
c        --- is zero we have to reset it to be 1.0d0
c        --- term 3 needs to be reinitialized to 1.0 to aviod 
c        --- division by zero in the tail corrections
         do itype = 1,ntypes
            do jtype = 1,ntypes
               ij = twh_getnbtype(itype,jtype)
               call twh_nbcoeff(GLB_GET,3,ij,nbcoeff(3))
               if ( nbcoeff(3) .eq. 0.0d0 ) then
                  call twh_nbcoeff(GLB_SET,3,ij,1.0d0)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_EXPONENTIAL_12_6 ) then
c        --- if term 4 is zero then it needs to be reinitialized to 1.0 
c        --- to avoid problem in the exp-6 force fields
         do itype = 1,ntypes
            do jtype = 1,ntypes
               ij = twh_getnbtype(itype,jtype)
               call twh_nbcoeff(GLB_GET,4,ij,nbcoeff(4))
               if ( nbcoeff(4) .eq. 0.0d0 ) then
                  call twh_nbcoeff(GLB_SET,4,ij,1.0d0)
               endif
            enddo
         enddo
      elseif ( cpnumber .eq. CP_12_6_PLUS_12_10_H_BOND ) then
c        --- need to do hbonds
         lhbond = .true.
      endif

      if ( lhbond ) then
         do itype = 1,ntypes
            do jtype = itype,ntypes
               ij = twh_getnbtype(itype,jtype)
c              --- search through the hbond data structures to try and
c              --- find a match
               lfound = .false.
               hbtype = 0
               do while ( hbtype  .lt. THBONDMAX .and. .not. lfound )
                  hbtype = hbtype + 1
                  call twh_hbondnames(GLB_GET,hbtype,1,hbondnames(1))
                  call twh_hbondnames(GLB_GET,hbtype,2,hbondnames(2))
                  call twh_nbname(GLB_GET,itype,nbnamei)
                  call twh_nbname(GLB_GET,jtype,nbnamej)
                  if ( (nbnamei .eq. hbondnames(1)
     &                 .and. nbnamej .eq. hbondnames(2))
     &                 .or. (nbnamei .eq. hbondnames(2)
     &                 .and. nbnamej .eq. hbondnames(1) ) ) then
c                    --- match found, set nb_to_hb array
                     call twh_nb_to_hbtype(GLB_SET,ij,hbtype)
                     lfound = .true.
                  endif
               enddo
            enddo
         enddo
      endif

      if ( cpnumber .eq. CP_SCALED_LENNARD_JONES .or.
     &     cpnumber .eq. CP_WEEKS_CHANDLER_ANDERSON ) then
c        --- we wish to mark which atoms will be rescaled.  Rescaling of
c        --- LJ terms takes place during energy evaluation since lambdas
c        --- may vary at runtime.  
c        --- for the rest of setmixterms, scaled LJ is treated like LJ
c        --- for WCA this keeps track of repulsive-only atomtypes
         call scp_initialize_atom_scaling(ntypes)
      endif
            
c     --- set all of the nbcoeffs according to mixing rules and cp
      do itrial = 1,2
         do itype = 1,ntypes
            if ( itrial .eq. 1 ) then
c               --- we are setting up all of the dissimilar cross terms
               jlo = itype + 1
               jhi = ntypes
            else
c              --- we are setting up all of the self terms
               jlo = itype
               jhi = itype
            endif
            do jtype = jlo,jhi
               ii = twh_getnbtype(itype,itype)
               jj = twh_getnbtype(jtype,jtype)
               ij = twh_getnbtype(itype,jtype)
               if ( localcm .eq. 'Lorentz-Berthelot' ) then
                  call twh_set_lorentz_berthelot(lproblem,cpnumber
     &                 ,localcm,ii,jj,ij)
               elseif ( localcm .eq. 'LB plus manual' ) then
c                 --- check to see if this is one that is manually 
c                 --- adjusted, or if we should just use LB
                  call twh_mixrule_adjust_total(GLB_GET,nadjust)
                  iadjust = 0
                  lfound = .false.
                  do while ( iadjust .lt. nadjust .and. .not. lfound )
                     iadjust = iadjust + 1
                     call twh_mixrule_adjust_key(GLB_GET,iadjust,key)
                     if ( key .eq. 'element' ) then
                        call twh_mixrule_adjust_keynames(GLB_GET,iadjust
     &                       ,1,testone)
                        call twh_mixrule_adjust_keynames(GLB_GET,iadjust
     &                       ,2,testtwo)
                        call twh_elementname(GLB_GET,itype,elei)
                        call twh_elementname(GLB_GET,jtype,elej)
                        if ( (testone .eq. elei .and. testtwo .eq. elej)
     &                       .or.
     &                       (testone .eq. elej .and. testtwo .eq. elei)
     &                       ) then
                           lfound = .true.
                        endif
                     elseif ( key .eq. 'nbname' ) then
                        call twh_mixrule_adjust_keynames(GLB_GET,iadjust
     &                       ,1,testone)
                        call twh_mixrule_adjust_keynames(GLB_GET,iadjust
     &                       ,2,testtwo)
                        call twh_nbname(GLB_GET,itype,nbnamei)
                        elei = nbnamei(1:2)
                        call twh_nbname(GLB_GET,jtype,nbnamej)
                        elej = nbnamej(1:2)
                        if ( (testone .eq. elei .and. testtwo .eq. elej)
     &                       .or.
     &                       (testone .eq. elej .and. testtwo .eq. elei)
     &                       ) then
                           lfound = .true.
                        endif
                     else
                        write(6,*) 'SETMIXTERMS: unknown '
     &                       ,'mixrule_adjust_key for iadjust:',iadjust
                        write(6,'(a)') key
                        lfailure = .true.
                        return
                     endif
                     if ( lfound ) then
                        if ( cpnumber .eq. CP_LENNARD_JONES ) then
                           icoeff = 1
                           call twh_mixrule_adjustments(GLB_GET
     &                          ,iadjust,icoeff,dtemp)
c                          --- store sigma^6
                           dtemp = dtemp**6
                           call twh_nbcoeff(GLB_SET,icoeff+12,ij
     &                          ,dtemp)
                           icoeff = 2
                           call twh_mixrule_adjustments(GLB_GET
     &                          ,iadjust,icoeff,dtemp)
c                          --- store 4*epslion
                           dtemp = 4.0d0*dtemp
                           call twh_nbcoeff(GLB_SET,icoeff+12,ij
     &                          ,dtemp)
                           icoeff = 3
                           call twh_mixrule_adjustments(GLB_GET
     &                          ,iadjust,icoeff,dtemp)
c                          --- store sigma^6 (1-4)
                           dtemp = dtemp**6
                           call twh_nbcoeff(GLB_SET,icoeff+12,ij
     &                          ,dtemp)
                           icoeff = 4
                           call twh_mixrule_adjustments(GLB_GET
     &                          ,iadjust,icoeff,dtemp)
c                          --- store 4*epslion (1-4)
                           dtemp = 4.0d0*dtemp
                           call twh_nbcoeff(GLB_SET,icoeff+12,ij
     &                          ,dtemp)
                        else
                           write(6,*) 'SETMIXTERMS: classical'
     &                          ,' potential of :',localcp
                           write(6,*) 'not currently supported with '
     &                          ,localcm
                           lfailure = .true.
                           return
                        endif
                     endif
                  enddo
                  if ( .not. lfound ) then
c                    --- did not find a manual match, use LB
                     call twh_set_lorentz_berthelot(lproblem,cpnumber
     &                    ,localcm,ii,jj,ij)
                  endif
               elseif ( localcm .eq. 'Geometric' ) then
                  if ( cpnumber .eq. CP_LENNARD_JONES .or.
     &                 cpnumber .eq. CP_SCALED_LENNARD_JONES .or.
     &                 cpnumber .eq. CP_UFF_12_6 .or.
     &                 cpnumber .eq. CP_WEEKS_CHANDLER_ANDERSON ) then
c                    --- geometric mean of sigma values
c                    --- store sigma**6
                     call twh_nbcoeff(GLB_GET,1,ii,valuei)
                     call twh_nbcoeff(GLB_GET,1,jj,valuej)
                     call twh_geometricmean(lfailure,valuei,valuej
     &                    ,dtemp)
                     if ( lfailure ) return
                     dtemp = dtemp**6
                     call twh_nbcoeff(GLB_SET,13,ij,dtemp)
c                    --- store 1-4 sigma**6
                     call twh_nbcoeff(GLB_GET,3,ii,valuei)
                     call twh_nbcoeff(GLB_GET,3,jj,valuej)
                     call twh_geometricmean(lfailure,valuei,valuej
     &                    ,dtemp)
                     if ( lfailure ) return
                     dtemp = dtemp**6
                     call twh_nbcoeff(GLB_SET,15,ij,dtemp)
c                    --- geometric mean of epsilon
c                    --- store 4*epsilon
                     call twh_nbcoeff(GLB_GET,2,ii,valuei)
                     call twh_nbcoeff(GLB_GET,2,jj,valuej)
                     call twh_geometricmean(lfailure,valuei,valuej
     &                    ,dtemp)
                     if ( lfailure ) return
                     dtemp = 4.0d0*dtemp
                     call twh_nbcoeff(GLB_SET,14,ij,dtemp)
c                    --- store 4*epsion 1-4
                     call twh_nbcoeff(GLB_GET,4,ii,valuei)
                     call twh_nbcoeff(GLB_GET,4,jj,valuej)
                     call twh_geometricmean(lfailure,valuei,valuej
     &                    ,dtemp)
                     if ( lfailure ) return
                     dtemp = 4.0d0*dtemp
                     call twh_nbcoeff(GLB_SET,16,ij,dtemp)
                  elseif ( cpnumber .eq. CP_12_9_6 ) then
c                    --- geometric mean of the three nbcoeffs
c                    --- geometric mean of the "12" values
                     call twh_nbcoeff(GLB_GET,1,ii,valuei)
                     call twh_nbcoeff(GLB_GET,1,jj,valuej)
                     call twh_geometricmean(lfailure,valuei,valuej
     &                    ,dtemp)
                     if ( lfailure ) return
                     call twh_nbcoeff(GLB_SET,1,ij,dtemp)
c                    --- geometric mean of the "9" values
                     call twh_nbcoeff(GLB_GET,2,ii,valuei)
                     call twh_nbcoeff(GLB_GET,2,jj,valuej)
                     call twh_geometricmean(lfailure,valuei,valuej
     &                    ,dtemp)
                     if ( lfailure ) return
                     call twh_nbcoeff(GLB_SET,2,ij,dtemp)
c                    --- geometric mean of the "6" parameters
                     call twh_nbcoeff(GLB_GET,3,ii,valuei)
                     call twh_nbcoeff(GLB_GET,3,jj,valuej)
                     call twh_geometricmean(lfailure,valuei,valuej
     &                    ,dtemp)
                     if ( lfailure ) return
                     call twh_nbcoeff(GLB_SET,3,ij,dtemp)
                  else
                     lproblem = .true.
                  endif
               elseif ( localcm .eq. 'Sixth Power' ) then
c                 --- sixth order combination of epsilon and sigma
c                 --- note that nbcoeff(1) contains r^3
                  if ( cpnumber .eq. CP_9_6 ) then
                     call twh_nbcoeff(GLB_GET,1,ii,nbc_1_ii)
                     call twh_nbcoeff(GLB_GET,1,jj,nbc_1_jj)
                     if ( nbc_1_ii*nbc_1_jj .gt. 1.0d-5 ) then
                        dtemp = dsqrt( (nbc_1_ii**6+nbc_1_jj**6)
     &                       /(2.0d0))
                        call twh_nbcoeff(GLB_SET,1,ij,dtemp)
                        call twh_nbcoeff(GLB_GET,2,ii,nbc_2_ii)
                        call twh_nbcoeff(GLB_GET,2,jj,nbc_2_jj)
                        dtemp = 2.0d0*dsqrt(nbc_2_ii*nbc_2_jj)*
     &                       (  (nbc_1_ii**3)*(nbc_1_jj**3) )
     &                       /(nbc_1_ii**6 + nbc_1_jj**6)
                        call twh_nbcoeff(GLB_SET,2,ij,dtemp)
                     else
                        call twh_nbcoeff(GLB_SET,1,ij,0.0d0)
                        call twh_nbcoeff(GLB_SET,2,ij,0.0d0)
                     endif
                  else
                     lproblem = .true.
                  endif
               elseif ( localcm .eq. 'Explicit' ) then
                  if ( cpnumber .eq. CP_LENNARD_JONES .or.
     &                 cpnumber .eq. CP_SCALED_LENNARD_JONES .or.
     &                 cpnumber .eq. CP_12_6_PLUS_SOLVATION .or.
     &                 cpnumber .eq. CP_WEEKS_CHANDLER_ANDERSON
     &                 ) then
c                    --- explicit terms for LJ 12-6 
c                    --- need to convert
c                    --- store sixth power of sigma
                     call twh_nbcoeff(GLB_GET,1,ij,dtemp)
                     dtemp = (dtemp)**6
                     call twh_nbcoeff(GLB_SET,13,ij,dtemp)
c                    --- store 4*epsilon
                     call twh_nbcoeff(GLB_GET,2,ij,dtemp)
                     dtemp = 4.0d0*dtemp
                     call twh_nbcoeff(GLB_SET,14,ij,dtemp)
c                    --- sixth power 1-4 sigma
                     call twh_nbcoeff(GLB_GET,3,ij,dtemp)
                     dtemp = (dtemp)**6
                     call twh_nbcoeff(GLB_SET,15,ij,dtemp)
c                    --- 4*epsilon 1-4
                     call twh_nbcoeff(GLB_GET,4,ij,dtemp)
                     dtemp = 4.0d0*dtemp
                     call twh_nbcoeff(GLB_SET,16,ij,dtemp)
                  elseif ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD .or.
     &                    cpnumber .eq. CP_EAM_PAIR_ONLY ) then
                     call twh_interpolatestyle(GLB_GET,interpolatestyle)
                     if ( interpolatestyle .eq. 'cubicspline' ) then
c                       --- need to set up the spline routines for EAM
                        if ( itype .eq. jtype ) then
                           call twh_eam_embed_style(GLB_GET,itype
     &                          ,eamembedstyle)
                           if ( eamembedstyle .eq. 'table' ) then
c                             --- set up eam_embed splin knots
                              call twh_splineset(2,itype,jtype)
                           endif
                        endif
                        call twh_eam_dens_style(GLB_GET,itype,jtype
     &                       ,eamdensstyle)
                        if ( eamdensstyle .eq. 'table' ) then
c                          --- set up eam_dens spline knots
                           call twh_splineset(1,itype,jtype)
                        endif
                        if ( itype .ne. jtype .and. 
     &                       eamdensstyle .eq. 'table' ) then
c                          --- set up eam_dens spline knots
                           call twh_splineset(1,jtype,itype)
                        endif
                        call twh_eam_pair_style(GLB_GET,ij,eampairstyle)
                        if ( eampairstyle .eq. 'table' ) then
c                          --- set up table_pair
                           call twh_splineset(3,itype,jtype)
                        endif
                     endif
                  elseif ( cpnumber .eq. CP_TABULATED_PAIR ) then
                     call twh_interpolatestyle(GLB_GET,interpolatestyle)
                     if ( interpolatestyle .eq. 'cubicspline' ) then
c                       --- need to set up the spline routines
c                       --- set up table_pair
                        call twh_splineset(3,itype,jtype)
                     endif
                  else
c                    --- for explicit mixrules then for any other 
c                    --- potential there is just nothing to be done
                  endif
               elseif ( localcm .eq. 'Arithmetic' ) then
c                 --- arithmetic mean
                  if ( cpnumber .eq. CP_HARD_SPHERE .or.
     &                 cpnumber .eq. CP_REPULSIVE_SPHERE ) then
c                    --- store the actual distance in the 1 spot
                     call twh_nbcoeff(GLB_GET,1,ii,nbc_1_ii)
                     call twh_nbcoeff(GLB_GET,1,jj,nbc_1_jj)
                     dtemp = 0.5d0*(nbc_1_ii + nbc_1_jj)
                     call twh_nbcoeff(GLB_SET,1,ij,dtemp)
c                    --- store the square of the distance in the 2 spot
                     dtemp = (0.5d0*(nbc_1_ii + nbc_1_jj) )**2
                     call twh_nbcoeff(GLB_SET,2,ij,dtemp)
                  else
                     lproblem = .true.
                  endif
               elseif ( localcm .eq. 'Shukla' ) then
                  if ( cpnumber .eq. CP_LENNARD_JONES .or.
     &                 cpnumber .eq. CP_SCALED_LENNARD_JONES .or.
     &                 cpnumber .eq. CP_WEEKS_CHANDLER_ANDERSON .or.
     &                 cpnumber .eq. CP_12_6_PLUS_12_10_H_BOND ) then
c                    --- store nbcoeff(13) and nb(15) with sigma^6 
c                    --- store nbcoeff(14) and nb(16) with 4*epsilon
                     call twh_polarizability(GLB_GET,itype,polari)
                     call twh_polarizability(GLB_GET,jtype,polarj)
                     do id = 0,2,2
c                       --- first compute the epsilon*sigma^6 term
                        call twh_nbcoeff(GLB_GET,2+id,ii,nbc_2id_ii)
                        call twh_nbcoeff(GLB_GET,1+id,ii,nbc_1id_ii)
                        call twh_nbcoeff(GLB_GET,2+id,jj,nbc_2id_jj)
                        call twh_nbcoeff(GLB_GET,1+id,jj,nbc_1id_jj)
                        denom = ( nbc_2id_ii*(nbc_1id_ii**6)*(polarj**2)
     &                       + nbc_2id_jj*(nbc_1id_jj**6)*(polari**2))
                        if ( denom .ne. 0.0d0 ) then
                           esigsix = 2.0d0*(nbc_2id_ii*(nbc_1id_ii**6)
     &                          *nbc_2id_jj*(nbc_1id_jj**6)
     &                          *polari*polarj) / denom
                           esigtwelve = (((nbc_2id_ii
     &                          *(nbc_1id_ii**12))**(1.0d0/13.0d0)
     &                          + (nbc_2id_jj*(nbc_1id_jj**12)
     &                          )**(1.0d0/13.0d0)) /2.0d0 )**13
c                          --- store sigma^6
                           dtemp = esigtwelve/esigsix
                           call twh_nbcoeff(GLB_SET,13+id,ij,dtemp)
c                          --- store 4*epsilon
                           dtemp = 4.0d0*(esigsix/dtemp)
                           call twh_nbcoeff(GLB_SET,14+id,ij,dtemp)
                        else
c                          --- set to zero
                           dtemp = 0.0d0
                           call twh_nbcoeff(GLB_SET,13+id,ij,dtemp)
                           call twh_nbcoeff(GLB_SET,14+id,ij,dtemp)
                        endif
                     enddo
                  else
                     lproblem = .true.
                  endif
               elseif ( localcm .eq. 'MMFF' ) then
                  if ( cpnumber .eq. CP_BUFFERED_14_7 ) then
c                    --- get the donor/acceptor flag for both atoms
                     call twh_donor_acceptor(GLB_GET,itype,dai)
                     call twh_donor_acceptor(GLB_GET,jtype,daj)
c                    --- get the element for both atoms
                     call twh_elementname(GLB_GET,itype,twoelei)
                     call twh_elementname(GLB_GET,jtype,twoelej)
c                    --- check to see if there is a polar hydrogen
                     if ( ( dai .ne. '-' .and. twoelei .eq. ' H' ) .or.
     &                    ( daj .ne. '-' .and. twoelej .eq. ' H' )) then
c                       --- polar hydrogen so cross-term B is 0.0
                        bbb = 0.0d0
                     else
c                       --- use the normal cross-term B of 0.2
                        bbb = 0.2d0
                     endif
c                    --- check for a donor/acceptor pairing
                     if ( ( dai .eq. 'A' .and. daj .eq. 'D' ) .or.
     &                    ( dai .eq. 'D' .and. daj .eq. 'A' ) ) then
c                       --- use the DARAD and DAEPS values
                        epsscale = 0.5d0
                        radscale = 0.8d0
                     else
c                       --- no scaling
                        epsscale = 1.0d0
                        radscale = 1.0d0
                     endif
                     power = 0.25d0
                     nbbeta = 12.0d0
c                    --- determine the self parameters
                     call twh_nbcoeff(GLB_GET,3,ii,ai)
                     call twh_nbcoeff(GLB_GET,3,jj,aj)
                     call twh_nbcoeff(GLB_GET,1,ii,alphai)
                     call twh_nbcoeff(GLB_GET,1,jj,alphaj)
                     rii = ai*alphai**power
                     rjj = aj*alphaj**power
                     if ( itype .eq. jtype ) then
c                       --- set rij to either value
                        rij = rii
                     else
c                       --- augmented arithmetic for R terms
                        gamma = (rii - rjj)/(rii + rjj)
                        rij = 0.5d0*(rii + rjj) * (1.0d0 
     &                       + bbb*(1.0d0 
     &                       - twh_expon(-nbbeta*gamma*gamma)))
                     endif
                     rijfinal = radscale*rij
c                    --- slater-kirkwood for epsilon
                     call twh_nbcoeff(GLB_GET,4,ii,gi)
                     call twh_nbcoeff(GLB_GET,4,jj,gj)
                     call twh_nbcoeff(GLB_GET,2,ii,nni)
                     call twh_nbcoeff(GLB_GET,2,jj,nnj)
c                    --- use the unscaled rij to compute epsilon
                     epsfinal = epsscale
     &                    *((181.16d0*CON_KCALTOK*gi*gj*alphai*alphaj)
     &                    /( dsqrt(alphai/nni) + dsqrt(alphaj/nnj) ))
     &                    *(1.0d0/(rij**6))
c                    --- store the scaled rij in coeff 5
                     dtemp = rijfinal
                     call twh_nbcoeff(GLB_SET,5,ij,dtemp)
c                    --- store 1.07 rij in coeff 7
                     dtemp = 1.07d0*rijfinal
                     call twh_nbcoeff(GLB_SET,7,ij,dtemp)
c                    --- store 0.07 rij in coeff 8
                     dtemp = 0.07d0*rijfinal
                     call twh_nbcoeff(GLB_SET,8,ij,dtemp)
c                    --- store 1.12*rij^7 in coeff 9
                     dtemp = 1.12d0*(rijfinal**7)
                     call twh_nbcoeff(GLB_SET,9,ij,dtemp)
c                    --- store 0.12*rij^7 in coeff 10
                     dtemp = 0.12d0*(rijfinal**7)
                     call twh_nbcoeff(GLB_SET,10,ij,dtemp)
c                    --- store epsilon in coeff 6
                     call twh_nbcoeff(GLB_SET,6,ij,epsfinal)
                  else
                     lproblem = .true.
                  endif
               else
c                 --- non-defined mixing rule
                  write(6,*)'classical mixrule of ',localcm
     &                 ,' is not defined'
                  lfailure = .true.
                  return
               endif
               if ( lproblem ) then
                  write(6,*) 'SETMIXTERMS: invalid combination of'
     &                 ,' classical potential and classical mixrule'
                  write(6,*) 'potential:',localcp
                  write(6,*) 'mixrule:',localcm
                  lfailure = .true.
                  return
               endif
c              --- set the shift potential (nbcoeff(0,ij))
               if ( lshift )  then
                  dtemp = twh_vtwobody(lfailure,rcut*rcut,itype,jtype,-1
     &                 ,lovrlap,vf,dv)
                  if ( lfailure ) return
                  call twh_nbcoeff(GLB_SET,0,ij,dtemp)
               endif
               if ( lovrlap ) then
                  write(6,*) 'strange trouble in ffnonbond'
                  lfailure = .true.
                  return
               endif
            enddo
         enddo
      enddo

      if ( cpnumber .eq. CP_EXPONENTIAL_6 ) then
         do itype = 1,ntypes
            do jtype = itype,ntypes
               ij = twh_getnbtype(itype,jtype)
c              --- have to set up the nbcoeff(5,ij) for the exp-6
c              --- this is needed to avoid those annoying -infinite
c              --- energies that occur as r get small.  
               imax = 1000
               maxval = -10.0d0
               rpos = 0.0d0
               do idist = 1,imax
                  dist = dble(idist)*rcut/dble(imax)
                  test = twh_vtwobody(lfailure,dist*dist,itype,jtype,0
     &                 ,lovrlap,vf,dv)
                  if ( lfailure ) return
                  if ( lovrlap ) then
                     write(6,*) 'SETMIXTERMS: trouble with expminsq'
                     write(6,*) 'itype,jtype,ij',itype,jtype,ij
                     write(6,*) 'dist ',dist
                     lfailure = .true.
                     return
                  endif
                  if ( test .gt. maxval ) then
                     maxval = test
                     rpos = dist
                  endif
               enddo
c              --- set exp-min-sq 
               dtemp = rpos*rpos
               call twh_nbcoeff(GLB_SET,5,ij,dtemp)
            enddo
         enddo
      endif

c     --- setup any additional required vibrations
      do itest = 1,TVIBMAX
         call twh_vibstyle(GLB_GET,itest,style)
         if ( style .eq. BOND_STYLE_FENE ) then
c           --- setup vib(3)
            call twh_vibcoeff(GLB_GET,itest,1,dvibone)
            call twh_vibcoeff(GLB_GET,itest,2,dvibtwo)
            dvibthree = -0.5d0 * dvibtwo * dvibone*dvibone
            call twh_vibcoeff(GLB_SET,itest,3,dvibthree)
         endif
      enddo
      
c     --- change degrees into radians and assign cosine of angle
      do itest = 1, TBENMAX
         call twh_bencoeff(GLB_GET,itest,0,coefftemp)
         dtemp = coefftemp*twopi/360.0d0
         call twh_anglestyle(GLB_GET,itest,anglestyle)
         if ( anglestyle .eq. ANG_STYLE_NONBOND ) then
c           --- there is no bencoeff 0 value so use a default of Pi/2
            dtemp = twopi/4.0d0
         endif
         call twh_bencoeff(GLB_SET,itest,0,dtemp)
         call twh_bencoeff(GLB_GET,itest,0,coefftemp)
         dtemp = dcos(coefftemp)
         call twh_bencoeff(GLB_SET,itest,-1,dtemp)
         if ( anglestyle .eq. ANG_STYLE_MULTI_RIGID ) then
c           --- also need to convert term 1 for multiple angle
            call twh_bencoeff(GLB_GET,itest,1,coefftemp)
            dtemp = coefftemp*twopi/360.0d0
            call twh_bencoeff(GLB_SET,itest,1,dtemp)
         endif
      enddo

      return
      end

      function twh_getsign(test)
c     ******************************************************************
c     * returns the sign of a double precision argument                *
c     *                                                                *
c     * originally written 04-13-2005 by M.G. Martin                   *
c     * last modified 03-03-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETSIGN
#include "functions.h"
c     --- variables passed to/from the function
      double precision test
      if ( test .lt. 0.0d0 ) then
         twh_getsign = -1.0d0
      else
         twh_getsign = 1.0d0
      endif
      return
      end

      subroutine twh_getcrosssign(lfailure,testa,testb,goodsign)
c     ******************************************************************
c     * take two double precision numbers and returns the sign of both *
c     * of them if the signs are the same.  If the signs are different *
c     * this this returns an error condition                           *
c     *                                                                *
c     * originally written 12-02-2005 by M.G. Martin                   *
c     * last modified 11-14-2006 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETSIGN
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      double precision testa,testb,goodsign
c     --- local variables
      double precision sign_i,sign_j

c     --- combine the signs
      sign_i = twh_getsign(testa)
      sign_j = twh_getsign(testb)
      if ( sign_i .ne. sign_j ) then
c        --- incompatable signs
         write(6,*) 'GETCROSSSIGN: incompatiable signs'
         write(6,*) 'sign_i,sign_j',sign_i,sign_j
         write(6,*) 'testa,testb',testa,testb
         lfailure = .true.
      else
         goodsign = sign_i
      endif
      return
      end

      subroutine twh_set_lorentz_berthelot(lfailure,cpnumber,localcm
     &     ,ii,jj,ij)
c     ******************************************************************
c     * sets the nbcoeffs according to the Lorentz-Berthelot mixing    *
c     * rules                                                          *
c     *                                                                *
c     * split from setmixterms 04-26-2006 by M.G. Martin               *
c     * last modified 03-31-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETSIGN
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      character*(*) localcm
      integer cpnumber,ii,jj,ij
c     --- local variables
      integer icoeff,tabnpairi,tabnpairj
      double precision dtemp,gorfact,sigma,soverr,rminimum
      double precision valuei,valuej,tpairii,tpairij,tpairjj
      double precision powrn
      double precision nbc_1_ii,nbc_1_jj,nbc_2_ii,nbc_2_jj
      double precision nbc_3_ii,nbc_3_jj
      double precision nbc_1_ij,nbc_2_ij,nbc_4_ij

c     --- Lorentz-Berthelot rules 
      if ( cpnumber .eq. CP_LENNARD_JONES
     &     .or. cpnumber .eq. CP_SCALED_LENNARD_JONES
     &     .or. cpnumber .eq. CP_WEEKS_CHANDLER_ANDERSON
     &     .or. cpnumber .eq. CP_12_6_PLUS_12_10_H_BOND
     &     .or. cpnumber .eq. CP_DRUKKER_H_BOND ) then
c        --- store nbcoeff(13) and nb(15) with sigma^6 
c        --- store nbcoeff(14) and nb(16) with 4*epsilon
c        --- arimethic mean of sigma (including 1-4)
         call twh_nbcoeff(GLB_GET,1,ii,nbc_1_ii)
         call twh_nbcoeff(GLB_GET,1,jj,nbc_1_jj)
         dtemp = (0.5d0*(nbc_1_ii + nbc_1_jj))**6
         call twh_nbcoeff(GLB_SET,13,ij,dtemp)
         call twh_nbcoeff(GLB_GET,3,ii,nbc_3_ii)
         call twh_nbcoeff(GLB_GET,3,jj,nbc_3_jj)
         dtemp = (0.5d0*(nbc_3_ii + nbc_3_jj))**6
         call twh_nbcoeff(GLB_SET,15,ij,dtemp)
c        --- geometric mean of epsilon (including 1-4)
c        --- store 4 times the geometric mean of epsilon
         call twh_nbcoeff(GLB_GET,2,ii,valuei)
         call twh_nbcoeff(GLB_GET,2,jj,valuej)
         call twh_geometricmean(lfailure,valuei,valuej,dtemp)
         if ( lfailure ) return
         dtemp = 4.0d0*dtemp
         call twh_nbcoeff(GLB_SET,14,ij,dtemp)
c        --- store 4 times the geometric mean of 1-4 epsilon
         call twh_nbcoeff(GLB_GET,4,ii,valuei)
         call twh_nbcoeff(GLB_GET,4,jj,valuej)
         call twh_geometricmean(lfailure,valuei,valuej,dtemp)
         if ( lfailure ) return
         dtemp = 4.0d0*dtemp
         call twh_nbcoeff(GLB_SET,16,ij,dtemp)
      elseif ( cpnumber .eq. CP_SQUARE_WELL .or.
     &        cpnumber .eq. CP_REPULSIVE_WELL ) then
c        --- arithmetic mean of hard sphere and attractive 
c        --- well diameters
c        --- store the square of the distance
         call twh_nbcoeff(GLB_GET,1,ii,nbc_1_ii)
         call twh_nbcoeff(GLB_GET,1,jj,nbc_1_jj)
         dtemp = (0.5d0*(nbc_1_ii + nbc_1_jj) )**2
         call twh_nbcoeff(GLB_SET,1,ij,dtemp)
         call twh_nbcoeff(GLB_GET,2,ii,nbc_2_ii)
         call twh_nbcoeff(GLB_GET,2,jj,nbc_2_jj)
         dtemp = (0.5d0*(nbc_2_ii + nbc_2_jj) )**2
         call twh_nbcoeff(GLB_SET,2,ij,dtemp)
c        --- geometric mean of the attractive well depth
         call twh_nbcoeff(GLB_GET,3,ii,valuei)
         call twh_nbcoeff(GLB_GET,3,jj,valuej)
         call twh_geometricmean(lfailure,valuei,valuej,dtemp)
         if ( lfailure ) return
         call twh_nbcoeff(GLB_SET,3,ij,dtemp)
      elseif ( cpnumber .eq. CP_GORDON_N_6 ) then
c        --- sigma with arithmetic mean
         call twh_nbcoeff(GLB_GET,1,ii,nbc_1_ii)
         call twh_nbcoeff(GLB_GET,1,jj,nbc_1_jj)
         dtemp = 0.5d0*(nbc_1_ii + nbc_1_jj)
         call twh_nbcoeff(GLB_SET,1,ij,dtemp)
c        --- epsilon with geometric mean
         call twh_nbcoeff(GLB_GET,2,ii,valuei)
         call twh_nbcoeff(GLB_GET,2,jj,valuej)
         call twh_geometricmean(lfailure,valuei,valuej,dtemp)
         if ( lfailure ) return
         call twh_nbcoeff(GLB_SET,2,ij,dtemp)
c        --- n with geometric mean
         call twh_nbcoeff(GLB_GET,3,ii,valuei)
         call twh_nbcoeff(GLB_GET,3,jj,valuej)
         call twh_geometricmean(lfailure,valuei,valuej,powrn)
         if ( lfailure ) return
         call twh_nbcoeff(GLB_SET,3,ij,powrn)
c        --- compute the gordon (n/6)^(1/(n-6)) factor
         gorfact = (powrn/6.0d0)**(1.0d0/(powrn-6.0d0))
c        --- store a(n)
         call twh_nbcoeff(GLB_GET,1,ij,nbc_1_ij)
         dtemp = nbc_1_ij*(2.0d0**(1.0d0/6.0d0) - gorfact)
         call twh_nbcoeff(GLB_SET,4,ij,dtemp)
c        --- compute rminimum
         call twh_nbcoeff(GLB_GET,4,ij,nbc_4_ij)
         rminimum = nbc_4_ij + gorfact*nbc_1_ij
         soverr = nbc_1_ij/(rminimum-nbc_4_ij)
c        --- store 4 epsilon c(n)
         call twh_nbcoeff(GLB_GET,2,ij,nbc_2_ij)
         dtemp = -nbc_2_ij*(1.0d0/(soverr**powrn - soverr**6))
         call twh_nbcoeff(GLB_SET,5,ij,dtemp)
      elseif ( cpnumber .eq. CP_MULTIWELL .or.
     &        cpnumber .eq. CP_REPULSIVE_MULTIWELL
     &        ) then
c        --- run through the table_npair and combine the
c        --- distances with arithmetic means and the well 
c        --- depths with a geometric mean
c        --- this only works if both potentials have the 
c        --- same number of table_npair
         call twh_table_npair(GLB_GET,ii,tabnpairi)
         call twh_table_npair(GLB_GET,jj,tabnpairj)
         if ( tabnpairi .ne. tabnpairj ) then
c           --- problem
            write(6,*) 'SETMIXTERMS: cannot use '
     &           ,localcm,' with classical potential number '
     &           ,cpnumber,' unless the '
     &           ,'table_npair are the same'
            write(6,*) 'ii,table_npair',ii,tabnpairi
            write(6,*) 'jj,table_npair',jj,tabnpairj
            lfailure = .true.
            return
         endif
c        --- set the cross term table_npair
         call twh_table_npair(GLB_SET,ij,tabnpairi)
         do icoeff = 1,tabnpairi
c           --- arithmetic mean of the first parameter
            call twh_table_pair(GLB_GET,1,ii,icoeff,tpairii)
            call twh_table_pair(GLB_GET,1,jj,icoeff,tpairjj)
            tpairij = 0.5d0*(tpairii + tpairjj)
            call twh_table_pair(GLB_SET,1,ij,icoeff,tpairij)
c           --- geometric mean of the second parameter
            call twh_table_pair(GLB_GET,2,ii,icoeff,tpairii)
            call twh_table_pair(GLB_GET,2,jj,icoeff,tpairjj)
            valuei = twh_getsign(tpairii)
            valuej = twh_getsign(tpairjj)
            call twh_geometricmean(lfailure,valuei,valuej,dtemp)
            if ( lfailure ) return
            call twh_table_pair(GLB_SET,2,ij,icoeff,dtemp)
         enddo
      elseif ( cpnumber .eq. CP_HARD_2580_MULTISTEP
     &        .or. cpnumber .eq. CP_REPULSIVE_2580_MULTISTEP ) then
c        --- take the arithmetic average of the sigmas
         call twh_nbcoeff(GLB_GET,1,ii,nbc_1_ii)
         call twh_nbcoeff(GLB_GET,1,jj,nbc_1_jj)
         sigma = 0.5d0*(nbc_1_ii + nbc_1_jj)
         call twh_nbcoeff(GLB_SET,1,ij,sigma)
c        --- store the square of the distances
c        --- hard core distance
         dtemp = sigma**2
         call twh_nbcoeff(GLB_SET,7,ij,dtemp)
c        --- the 1.2 distance
         dtemp = (1.2d0*sigma)**2
         call twh_nbcoeff(GLB_SET,8,ij,dtemp)
c        --- the 1.5 distance
         dtemp = (1.5d0*sigma)**2
         call twh_nbcoeff(GLB_SET,9,ij,dtemp)
c        --- the 1.8 distance
         dtemp = (1.8d0*sigma)**2
         call twh_nbcoeff(GLB_SET,10,ij,dtemp)
c        --- the 2.0 distance
         dtemp = (2.0d0*sigma)**2
         call twh_nbcoeff(GLB_SET,11,ij,dtemp)
c        --- geometric mean of the well depths
c        --- the 1.2 well
         call twh_nbcoeff(GLB_GET,2,ii,valuei)
         call twh_nbcoeff(GLB_GET,2,jj,valuej)
         call twh_geometricmean(lfailure,valuei,valuej,dtemp)
         if ( lfailure ) return
         call twh_nbcoeff(GLB_SET,2,ij,dtemp)
c        --- the 1.5 well
         call twh_nbcoeff(GLB_GET,3,ii,valuei)
         call twh_nbcoeff(GLB_GET,3,jj,valuej)
         call twh_geometricmean(lfailure,valuei,valuej,dtemp)
         if ( lfailure ) return
         call twh_nbcoeff(GLB_SET,3,ij,dtemp)
c        --- the 1.8 well
         call twh_nbcoeff(GLB_GET,4,ii,valuei)
         call twh_nbcoeff(GLB_GET,4,jj,valuej)
         call twh_geometricmean(lfailure,valuei,valuej,dtemp)
         if ( lfailure ) return
         call twh_nbcoeff(GLB_SET,4,ij,dtemp)
c        --- the 2.0 well
         call twh_nbcoeff(GLB_GET,5,ii,valuei)
         call twh_nbcoeff(GLB_GET,5,jj,valuej)
         call twh_geometricmean(lfailure,valuei,valuej,dtemp)
         if ( lfailure ) return
         call twh_nbcoeff(GLB_SET,5,ij,dtemp)
      else
         lfailure = .true.
      endif

      return
      end

      subroutine twh_geometricmean(lfailure,valuei,valuej,mean)
c     ******************************************************************
c     * returns the geometric mean of two variables preserving their   *
c     * sign                                                           *
c     *                                                                *
c     * originally written 11-14-2006 by M.G. Martin                   *
c     * last modified 11-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      double precision valuei,valuej,mean
c     --- local variables
      double precision thesign

      call twh_getcrosssign(lfailure,valuei,valuej,thesign)
      if ( lfailure ) return

      mean =  thesign*dsqrt(valuei*valuej)

      return
      end
