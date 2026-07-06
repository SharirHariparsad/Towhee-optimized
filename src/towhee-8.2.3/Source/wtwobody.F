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
      function twh_wtwobody(lfailure,rijsq,ntii,ntjj,lovrlap )
c     ******************************************************************
c     * Computes the intermolecular pair virial function divided by r^2*
c     * as defined by equation 2.60 in Allen and Tildesley             *
c     * w(r) = r dv(r)/dr                                              *
c     * w(r)/r^2 = (dv(r)/dr) /r                                       *
c     *                                                                *
c     * originally written (as ftwobody) 07-03-2002 by M.G. Martin     *
c     * last update 03-31-2020 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_EXPON
#define FUNCTION_EXTRACTPAIR
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_NATIVE_LAMBDA
#define FUNCTION_SCALE_ATOMS
#define FUNCTION_WTWOBODY
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure,lovrlap
      integer ntii,ntjj
      double precision rijsq
c     --- local variables
c     --- logical scalars
      logical lcontinue
c     --- integer scalars
      integer ntij,icoeff,ipower,hbtype,cpnumber,aindex,eampairnum
c     --- double precision scalars
      double precision srsix,srthree,rij,lambda,rdiff,rtest,rvalue
      double precision rinvfour,rinveight,soverr
      double precision rinvthree,rinvsix,rinvnine,rinvtwelve,wca_cutsq
      double precision x,acklandphi,acklandphiprime,xprime,ainv,onethird
      double precision cbrt_two,tpairone,tpairtwo,rcutsq,rminsq
      double precision avalue,powrn
c     --- double precision arrays
      double precision hbondcoeff
      dimension hbondcoeff(2)
      double precision nbcoeff
      dimension nbcoeff(MINNBCOEFF:MAXNBCOEFF)
c     --- functions
      double precision wtb_lennard_jones

c      write(6,*) 'starting WTWOBODY'
c     --- retrieve constants
      call twh_cp_number(GLB_GET,cpnumber)
      call twh_rcutsq(GLB_GET,rcutsq)
      call twh_rminsq(GLB_GET,rminsq)

c     --- initialize wtwobody
      twh_wtwobody = 0.0d0
      lovrlap = .false.
c     --- all use the same logic for finding pairs
      ntij = twh_getnbtype(ntii,ntjj)

c     --- compute two-body forces according to the classical_potential
      if ( cpnumber .eq. CP_LENNARD_JONES .or.
     &     cpnumber .eq. CP_12_6_PLUS_SOLVATION .or.
     &     cpnumber .eq. CP_UFF_12_6 ) then
c        --- LJ 12-6 forces
         if ( twh_cmp_lt(rijsq, rminsq) ) then
            lovrlap = .true.
         elseif ( twh_cmp_lt(rijsq, rcutsq) ) then
            twh_wtwobody =  wtb_lennard_jones(ntij, rijsq)
         endif
      elseif ( cpnumber .eq. CP_SCALED_LENNARD_JONES ) then
c        --- soft-core lennard-jones
         if ( rijsq .lt. rminsq ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rcutsq ) then
c        --- v = 4 lambda**4 epsilon (ainv**2 - ainv)
c        --- srsix here is inverse of what it is for above case
c        --- for sigma (nbcoeff(13)) = 0, wtwobody = 0
            do icoeff = 13,14
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            if ( nbcoeff(13) .eq. 0.0d0) then
               twh_wtwobody = 0.0d0
            else
               lambda = scp_get_native_lambda(LAMBDA_LJ)
               if (scp_scale_atoms(ntii,ntjj)) then 
                  srsix = (rijsq**3) / nbcoeff(13) 
                  ainv = 1.0d0 / (scp_get_native_lambda(ALPHA_LJ) * 
     &               (1 - lambda)**2 + srsix )
                  twh_wtwobody = 6.0d0 * lambda**4 * 
     &               nbcoeff(14)*(1.0d0-2.0d0*ainv)*ainv**2*srsix/rijsq
               else
                  twh_wtwobody =  wtb_lennard_jones(ntij, rijsq)
               endif 
            endif
         endif
      elseif ( cpnumber .eq. CP_WEEKS_CHANDLER_ANDERSON ) then
c        --- 1/r dv(r)/dr is same as LJ 12-6 except for range
         if ( rijsq .lt. rminsq ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rcutsq ) then
            do icoeff = 13,14
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            srsix = nbcoeff(13) / (rijsq**3)
            if (scp_scale_atoms(ntii, ntjj)) then
               call twh_constant_onethird(GLB_GET,onethird)
               call twh_constant_cuberoottwo(GLB_GET,cbrt_two)
               wca_cutsq = cbrt_two * (nbcoeff(13)**onethird)
               if (rijsq .lt. wca_cutsq) then
                  twh_wtwobody = nbcoeff(14)*12.0d0*srsix
     &                 *(0.5d0-srsix)/rijsq
               else 
                  twh_wtwobody = 0.0d0
c              --- 0 if rijsq .ge. wca_cutsq
               endif
            else
c           --- if not lscale then normal LJ interaction
            endif
            twh_wtwobody = nbcoeff(14)*12.0d0*srsix
     &           *(0.5d0-srsix)/ rijsq
         endif
      elseif ( cpnumber .eq. CP_12_6_PLUS_12_10_H_BOND ) then
c        --- LJ 12-6 forces as above, plus a 12-10 repulsion
         if ( rijsq .lt. rminsq ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rcutsq ) then
            do icoeff = 13,14
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            srsix = nbcoeff(13) / (rijsq**3)
            twh_wtwobody = nbcoeff(14)*12.0d0*srsix
     &           *(0.5d0-srsix)/ rijsq
            call twh_nb_to_hbtype(GLB_GET,ntij,hbtype)
            if ( hbtype .ne. 0 ) then
               call twh_hbondcoeff(GLB_GET,hbtype,1,hbondcoeff(1))
               call twh_hbondcoeff(GLB_GET,hbtype,2,hbondcoeff(2))
               twh_wtwobody = twh_wtwobody 
     &              - 12.0d0*hbondcoeff(1)/rijsq**7
     &              + 10.0d0*hbondcoeff(2)/rijsq**6
            endif
         endif
      elseif ( cpnumber .eq. CP_9_6 ) then
c        --- LJ 9-6 forces
         if ( rijsq .lt. rminsq ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rcutsq ) then
            rij = dsqrt(rijsq)
            do icoeff = 1,2
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            srthree = nbcoeff(1) / (rijsq*rij)
            twh_wtwobody = (18.0d0*nbcoeff(2)*
     &           srthree*srthree*(1.0d0-srthree)
     &           / rijsq)
         endif
      elseif ( cpnumber .eq. CP_12_9_6 ) then
c        --- 12-9-6 potential
         if ( rijsq .lt. rminsq ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rcutsq ) then
            rij = dsqrt(rijsq)
            rinvthree = 1.0d0/(rijsq*rij)
            rinvsix = rinvthree*rinvthree
            rinvnine = rinvsix*rinvthree
            rinvtwelve = rinvsix*rinvsix
            do icoeff = 1,3
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            twh_wtwobody = -12.0d0*nbcoeff(1)*rinvtwelve 
     &           - 9.0d0*nbcoeff(2)*rinvnine
     &           - 6.0d0*nbcoeff(3)*rinvsix
         endif
      elseif ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD .or.
     &        cpnumber .eq. CP_EAM_PAIR_ONLY ) then
c        --- the pair potential part of Embedded Atom
         call twh_eam_pair_number(GLB_GET,ntij,eampairnum)
         if ( eampairnum .eq. EAM_PAIR_ACKLANDTP ) then
            rij = dsqrt(rijsq)
            do icoeff = 1,9
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            if ( rij .lt. nbcoeff(1) ) then
c                 --- coulombic type repulsion
               x = rij/nbcoeff(4)
               acklandphi = 0.1818d0*twh_expon(-3.2d0*x) 
     &              + 0.5099d0*twh_expon(-0.9423d0*x)
     &              + 0.2802d0*twh_expon(-0.4029d0*x)
     &              + 0.02817d0*twh_expon(-0.2016d0*x)
               acklandphiprime = 
     &              0.1818d0*(-3.2d0)*twh_expon(-3.2d0*x)/nbcoeff(4)
     &              + 0.5099d0*(-0.9423d0)*twh_expon(-0.9423d0*x)
     &              /nbcoeff(4)
     &              + 0.2802d0*(-0.4029d0)*twh_expon(-0.4029d0*x)
     &              /nbcoeff(4)
     &              + 0.02817d0*(-0.2016d0)*twh_expon(-0.2016d0*x)
     &              /nbcoeff(4)
               twh_wtwobody = nbcoeff(3)
     &              *( acklandphiprime/rij - acklandphi/rijsq ) /rij
            elseif ( rij .lt. nbcoeff(2) ) then
c                 --- exponential repulsion
               x = nbcoeff(5) + nbcoeff(6)*rij + nbcoeff(7)*rijsq
     &              + nbcoeff(8)*rij*rijsq
               xprime = nbcoeff(6) + 2.0d0*nbcoeff(7)*rij
     &              + 3.0d0*nbcoeff(8)*rijsq
               twh_wtwobody = nbcoeff(9)*xprime*twh_expon(x)/rij
            else
c                 --- cubic sum with heaviside (using a while loop to 
c                 --- recreate the heaviside and save some expense)
               call twh_table_npair(GLB_GET,ntij,icoeff)
               if ( icoeff .gt. 0 ) then
                  lcontinue = .true.
               else
                  lcontinue = .false.
               endif
               do while ( lcontinue )
                  call twh_table_pair(GLB_GET,1,ntij,icoeff,tpairone)
                  x = tpairone - rij
                  call twh_table_pair(GLB_GET,2,ntij,icoeff,tpairtwo)
                  twh_wtwobody = twh_wtwobody 
     &                 + tpairtwo*(-3.0d0)*x*x/rij
                  icoeff = icoeff - 1
                  if ( icoeff .eq. 0 ) then
                     lcontinue = .false.
                  else
                     call twh_table_pair(GLB_GET,1,ntij,icoeff,tpairone)
                     if ( tpairone .le. rij ) lcontinue = .false.
                  endif
               enddo
            endif
         elseif ( eampairnum .eq. EAM_PAIR_ACKLANDPOW ) then
c           --- complicated pair potential that has 3 different forms
c           --- depending upon the distances
            rij = dsqrt(rijsq)
            do icoeff = 1,16
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            if ( rij .lt. nbcoeff(1) ) then
c              --- coulombic type repulsion
               acklandphi 
     &              = nbcoeff(4)
     &              *twh_expon(nbcoeff(5)*rij) 
     &              + nbcoeff(6)
     &              *twh_expon(nbcoeff(7)*rij)
     &              + nbcoeff(8)
     &              *twh_expon(nbcoeff(9)*rij)
     &              + nbcoeff(10)
     &              *twh_expon(nbcoeff(11)*rij)
               acklandphiprime
     &              = nbcoeff(4)*nbcoeff(5)
     &              * twh_expon(nbcoeff(5)*rij)
     &              + nbcoeff(6)*nbcoeff(7)
     &              * twh_expon(nbcoeff(7)*rij)
     &              + nbcoeff(8)*nbcoeff(9)
     &              * twh_expon(nbcoeff(9)*rij)
     &              + nbcoeff(10)*nbcoeff(11)
     &              * twh_expon(nbcoeff(11)*rij)
               twh_wtwobody = nbcoeff(3)
     &              *(acklandphiprime - acklandphi/rij )/rijsq
            elseif ( rij .lt. nbcoeff(2) ) then
c              --- exponential repulsion
               x = nbcoeff(12) 
     &              + nbcoeff(13)*rij
     &              + nbcoeff(14)*rijsq
     &              + nbcoeff(15)*rij*rijsq
               xprime = nbcoeff(13) 
     &              + nbcoeff(14)*(2.0d0)*rij
     &              + nbcoeff(15)*(3.0d0)*rijsq
               twh_wtwobody = xprime*nbcoeff(16)
     &              *twh_expon(x)/rij
            else
c              --- power series with heaviside (using a while loop to
c              --- recreate the heaviside and save some expense)
               call twh_table_npair(GLB_GET,ntij,icoeff)
               if ( icoeff .gt. 0 ) then
                  lcontinue = .true.
               else
                  lcontinue = .false.
               endif
               do while ( lcontinue ) 
                  call twh_table_pair(GLB_GET,1,ntij,icoeff,tpairone)
                  x = tpairone - rij
                  ipower = icoeff + 3
                  call twh_table_pair(GLB_GET,2,ntij,icoeff,tpairtwo)
                  twh_wtwobody = twh_wtwobody 
     &                 - tpairtwo*dble(ipower)*(x**(ipower-1))
                  icoeff = icoeff - 1
                  if ( icoeff .eq. 0 ) then
                     lcontinue = .false.
                  else
                     call twh_table_pair(GLB_GET,1,ntij,icoeff,tpairone)
                     if ( tpairone .le. rij ) lcontinue = .true.
                  endif
               enddo
            endif
         elseif ( eampairnum .eq. EAM_PAIR_BELASHFE ) then
c           --- Belashchenko Heaviside power series
c           --- doing this in declining distance in an attempt to only
c           --- fetch the minimum required parameters
            rij = dsqrt(rijsq)
            do icoeff = 1,10
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            if ( twh_cmp_gt(rij,nbcoeff(6)) ) then
c              --- beyond the cutoff distance, zero interaction
               twh_wtwobody = 0.0d0
            elseif ( twh_cmp_lt(rij,nbcoeff(1)) ) then
c              --- below r1, use a different potential
               rtest = nbcoeff(1)
               twh_wtwobody = (- nbcoeff(8) - nbcoeff(9)*nbcoeff(10)
     &              *twh_expon(nbcoeff(10)*(rtest-rij)))/rij
            else
c              --- between r1 and r6, power series
               if ( twh_cmp_lt(rij,nbcoeff(2)) ) then
c                 --- between r1 and r2
                  rvalue = nbcoeff(2)
                  aindex = 11
               elseif ( twh_cmp_lt(rij, nbcoeff(3)) ) then
c                 --- between r2 and r3
                  rvalue = nbcoeff(3)
                  aindex = 20
               elseif ( twh_cmp_lt(rij, nbcoeff(4)) ) then
c                 --- between r3 and r4
                  rvalue = nbcoeff(4)
                  aindex = 29
               elseif ( twh_cmp_lt(rij, nbcoeff(5)) ) then
c                 --- between r4 and r5
                  rvalue = nbcoeff(5)
                  aindex = 38
               else
c                 --- between r5 and r6
                  rvalue = nbcoeff(6)
                  aindex = 47
               endif
c              --- compute the power series
               rdiff = rij - rvalue
c              --- can skip power 0 as it has a 0 deriviative
c              --- power 1 rdiff is 1.0 (rdiff**(1-1))
               ipower = 1
               call twh_nbcoeff(GLB_GET,aindex+ipower,ntij,avalue)
               twh_wtwobody = twh_wtwobody + avalue/rij
               do ipower = 2,8
                  call twh_nbcoeff(GLB_GET,aindex+ipower,ntij,avalue)
                  twh_wtwobody = twh_wtwobody
     &                 + dble(ipower)*avalue*(rdiff**(ipower-1))/rij
               enddo
            endif   
         elseif ( eampairnum .eq. EAM_PAIR_TABLE ) then
c          --- compute the derivatives from the tabulated potential
            rij = dsqrt(rijsq)
            twh_wtwobody = twh_extractpair(lfailure,INTERP_DERIV
     &           ,ntij,rij,lovrlap)/rij
            if ( lfailure ) return
         elseif ( eampairnum .eq. EAM_PAIR_NONE ) then
c           --- no energy, no forces
            twh_wtwobody = 0.0d0
         else
            write(6,*) 'WTWOBODY: unknown eam_pair_num'
            write(6,*) eampairnum
            lfailure = .true.
            return
         endif
      elseif ( cpnumber .eq. CP_EXPONENTIAL_6 ) then
c        --- exponential-6 forces
         do icoeff = 1,5
            call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
         enddo
         if ( rijsq .lt. nbcoeff(5) ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rcutsq ) then
            rij = dsqrt(rijsq)
            rinvfour = 1.0d0/(rijsq*rijsq)
            rinveight = rinvfour*rinvfour
            twh_wtwobody =  -6.0d0*nbcoeff(1)*rinveight
     &           + nbcoeff(2)*nbcoeff(3)*twh_expon(nbcoeff(3)*rij)/rij
         endif
      elseif ( cpnumber .eq. CP_EXPONENTIAL_12_6 ) then
c        --- exp-6 with LJ 12
         if ( rijsq .lt. rminsq ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rcutsq ) then
            rij = dsqrt(rijsq)
            rinvfour = 1.0d0/(rijsq*rijsq)
            rinveight = rinvfour*rinvfour
            do icoeff = 1,4
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            twh_wtwobody = -6.0d0*nbcoeff(1)*rinveight
     &           + nbcoeff(3)*nbcoeff(4)*twh_expon(nbcoeff(4)*rij)/rij
     &           - 12.0*nbcoeff(2)*rinveight*rinveight*rijsq
         endif
      elseif ( cpnumber .eq. CP_GORDON_N_6 ) then
         if ( rijsq .lt. rminsq ) then
            lovrlap = .true.
         elseif ( rijsq .lt. rcutsq ) then
            do icoeff = 1,5
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            powrn = nbcoeff(3)
            rij = dsqrt(rijsq)
            soverr = nbcoeff(1)/(rij - nbcoeff(4))
            twh_wtwobody = (nbcoeff(5)
     &           /(nbcoeff(1)*rij))*(6.0d0*(soverr**7)
     &           -powrn*(soverr**(powrn+1.0d0)))
         endif
      else
         write(6,*) 'WTWOBODY: undefined classical potential'
     &        ,cpnumber
         lfailure = .true.
         return
      endif

      return
      end

      double precision function wtb_lennard_jones(ntij, rijsq)
c     ******************************************************************
c     * Computes the Lennard Jones intermolecular pair virial function *
c     * divided by r^2                                                 *
c     *                                                                *
c     * originally split out of wtwobody 03-01-2006 by MAW             *
c     * last modified 03-31-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ntij
      double precision rijsq
c     --- local variables
      double precision srsix,nbcoeff_13,nbcoeff_14

      call twh_nbcoeff(GLB_GET,13,ntij,nbcoeff_13)
      call twh_nbcoeff(GLB_GET,14,ntij,nbcoeff_14)
      srsix = nbcoeff_13/(rijsq**3)
      wtb_lennard_jones = nbcoeff_14*12.0d0*srsix*(0.5d0-srsix)/rijsq

      return
      end
