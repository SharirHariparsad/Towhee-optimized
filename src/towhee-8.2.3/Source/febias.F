#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2020 Marcus G. Martin                       *
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
      function twh_febias(lfailure,bxdflag,ibox,biastype,imolty,iunit
     &     ,index,iufrom,xcord,ycord,zcord)
c     ******************************************************************
c     * computes a biasing probability for use in fixed-endpoint cbmc  *
c     * ibox: the box the molecule is in - used for angles             *
c     * biastype is the type of biasing to compute (angle,torsion,etc) *
c     *          1: 1-3 interaction bias (angle)                       *
c     *          2: 1-4 interaction bias (dihedral)                    *
c     * imolty is the molecule type                                    *
c     * iunit is the unit number of the first atom                     *
c     * index is the index of the angle or the torsion (depending on   *
c     *       what we are biasing on) starting from iunit              *
c     * iufrom is the unit previous to iunit, this is used for the     *
c     *         torsion bias section and goes with the xyzcords        *
c     * xyzcord are the arrays with the relevent atom coordinates      *
c     *         *cord(1) is iunit                                      *
c     *         *cord(2) is iuother                                    *
c     *         *cord(3) is iufrom                                     *
c     *                                                                *
c     * originally written 01-23-2001 by M.G. Martin                   *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_DISTANCE
#define FUNCTION_EXPON
#define FUNCTION_FEBIAS
#define FUNCTION_GAUSSPROB
#define FUNCTION_VANGLE
#define FUNCTION_VBOND
#define FUNCTION_VTORSION
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure
      integer biastype,index,iunit,imolty,bxdflag,ibox,iufrom
      double precision mean,stddev,lolimit,hilimit
      double precision xcord,ycord,zcord
      dimension xcord(3),ycord(3),zcord(3)
c     --- parameters
      double precision minbias
      parameter (minbias = 1d-40)
c     --- local variables
c     --- logical scalars
      logical lfound,lovrlap
c     --- character string scalars
      character*50 twobbias_style,threebbias_style
c     --- integer scalars
      integer iuother,btone,bttwo,junit,kunit,iunext,iuthree
      integer goodind,ivib,ibend,iumiss,icord,iphi,torind,itor,itemp
      integer ijben2test,ijben3test,ijtor2test,ijtor3test,ijtor4test
      integer ijvibtest,inbend,intorsion,invibration,itvibtest
c     --- double precision scalars
      double precision lenone,lentwo,thetac,scale,theta
      double precision rxuij,ryuij,rzuij,length,dvalue
      double precision botlensq,midlen,prelen,nexlen,bigang,lowang
      double precision minang,maxang,minlen,maxlen,eqlength,eqangle
      double precision vbias,twopi,beta,a_len,b_len,c_len,r_len
      double precision phic,dg_mag_sq,bg_mag,duma,dumb,dumc,minboxdim
      double precision two_bond_bias_sdev_multiplier,mult
c     --- double precision arrays
      double precision bvect
      dimension bvect(3)
      double precision gvect
      dimension gvect(3)
      double precision dvect
      dimension dvect(3)
      double precision mvect
      dimension mvect(2,3)
      double precision xtilda
      dimension xtilda(3)
      double precision ytilda
      dimension ytilda(3)
      double precision phi
      dimension phi(2)
      double precision xatmp,yatmp,zatmp
      dimension xatmp(3),yatmp(3),zatmp(3)
      double precision xttmp,yttmp,zttmp
      dimension xttmp(4),yttmp(4),zttmp(4)
      double precision ddum3
      dimension ddum3(3)

c      write(6,*) 'starting FEBIAS'

c     --- initialize febias and lovrlap
      twh_febias = 1.0d0
      lovrlap = .false.
      theta = 0.0d0
c     --- retrieve constants
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_beta(GLB_GET,beta)
      call twh_minboxdim(GLB_GET,minboxdim)

      if ( biastype .eq. 1 ) then
c       --- two bond fixed endpoint biasing
         call twh_two_bond_fe_bias_style(GLB_GET,twobbias_style)
         if ( twobbias_style .eq. 'none' ) then
c           --- nothing to be done, return with the default bias
         elseif ( twobbias_style 
     &           .eq. 'analytic Boltzmann using angles' ) then
c           --- this is an angle bias based on the angle energy given 
c           --- the equilibrium bond lengths
            call twh_ibtben1(GLB_GET,imolty,iunit,index,btone)
            call twh_ibtben2(GLB_GET,imolty,iunit,index,bttwo)
            call twh_vibcoeff(GLB_GET,btone,0,lenone)
            call twh_vibcoeff(GLB_GET,bttwo,0,lentwo)
c           --- get the atom types for the rest of the angle
            call twh_ijben2(GLB_GET,imolty,iunit,index,junit)
            call twh_ijben3(GLB_GET,imolty,iunit,index,kunit)
c           --- compute the iunit-iuother length from the cords
            call twh_mimage(bxdflag,ibox,xcord(2),ycord(2),zcord(2)
     &           ,xcord(1),ycord(1),zcord(1),rxuij,ryuij,rzuij)
            length = twh_distance(rxuij,ryuij,rzuij)

c           --- intitialize vbias
            vbias = 0.0d0
            if ( length .gt. lenone+lentwo ) then
c              --- use an angle of Pi and compute stretching on the bond
c              --- add in a bit of noise to avoid exactly zero
               thetac = -1.0d0+1d-10
               scale = length/(lenone+lentwo)
               lenone = scale*lenone
               lentwo = scale*lentwo
               vbias = vbias + twh_vbond(lfailure,.false.,bxdflag,ibox
     &              ,btone,lenone,imolty,iunit,junit,lovrlap) 
     &              + twh_vbond(lfailure,.false.,bxdflag,ibox,bttwo
     &              ,lentwo,imolty,junit,kunit,lovrlap)
               if ( lfailure ) return
            else
c              --- compute the angle assuming the bond lengths are at 
c              --- their equilibrium values
               thetac = (lenone*lenone + lentwo*lentwo - length*length)
     &              /(2.0d0*lenone*lentwo)
            endif

c           --- get bias energies from vangle
            vbias = vbias + twh_vangle(lfailure,.false.,lovrlap
     &           ,bxdflag,ibox,imolty
     &           ,iunit,index,0,thetac,lenone,lentwo,ddum3,ddum3,ddum3)
            if ( lfailure ) return
c           --- convert the vbias energy into a probability
            if ( .not. lovrlap ) then
               twh_febias = twh_expon(-beta*vbias)
c              --- account for the phase space term
               theta = twh_arccos(thetac)
               twh_febias = twh_febias*dsin(theta)
            endif
         elseif ( twobbias_style .eq. 
     &           'analytic Boltzmann dihedral energy sum'
     &           ) then
c           --- we have one missing atom between the one we are growing
c           --- and our target.  
c           --- M.D. Rintoul helped out with this function
c           --- set the unit numbers, iunit and iufrom are already set
            call twh_ijben2(GLB_GET,imolty,iunit,index,iumiss)
            call twh_ijben3(GLB_GET,imolty,iunit,index,iuother)
c           --- compute the "a_len" ( |iunit-ifrom| )
            a_len = dsqrt((xcord(1)-xcord(3))**2 
     &           + (ycord(1)-ycord(3))**2 + (zcord(1)-zcord(3))**2)
c           --- set "b_len" to the equilibrium iunit-missing bond length
            call twh_ibtben1(GLB_GET,imolty,iunit,index,btone)
            call twh_vibcoeff(GLB_GET,btone,0,b_len)
c           --- set "c_len" to the equil. missing-iother bond length
            call twh_ibtben2(GLB_GET,imolty,iunit,index,bttwo)
            call twh_vibcoeff(GLB_GET,bttwo,0,c_len)

c           --- safety check that a reasonable solution is possible by
c           --- comparing iunit-iuother distance with b_len and c_len
            length = dsqrt( (xcord(2)-xcord(1))**2 
     &           + (ycord(2)-ycord(1))**2 + (zcord(2)-zcord(1))**2)
            if ( length .gt. b_len + c_len ) then
c              --- this is a pretty crummy conformation
               twh_febias = minbias
               return
            endif
c           --- set up the "beta" vector (iunit-iufrom)
            bvect(1) = xcord(1)-xcord(3)
            bvect(2) = ycord(1)-ycord(3)
            bvect(3) = zcord(1)-zcord(3)
c           --- set up the "gamma" vector (iother-iufrom)
            gvect(1) = xcord(2)-xcord(3)
            gvect(2) = ycord(2)-ycord(3)
            gvect(3) = zcord(2)-zcord(3)
c           --- set theta to the equilibrium ifrom--iunit--missing angle
            lfound = .false.
            ibend = 0
            call twh_inben(GLB_GET,imolty,iufrom,inbend)
            do while ( ibend .lt. inbend .and. (.not. lfound) )
               ibend = ibend + 1
               call twh_ijben2(GLB_GET,imolty,iufrom,ibend,ijben2test)
               call twh_ijben3(GLB_GET,imolty,iufrom,ibend,ijben3test)
               if ( ijben2test .eq. iunit .and. ijben3test .eq. iumiss 
     &              ) then
                  lfound = .true.
                  call twh_itben(GLB_GET,imolty,iufrom,ibend,itemp)
                  call twh_bencoeff(GLB_GET,itemp,0,theta)
               endif
            enddo
            if ( .not. lfound ) then
               write(6,*) 'FEBIAS: angle not found'
               lfailure = .true.
               return
            endif
c           --- compute "r_len", the radius of the cone
            r_len = b_len*dsin(theta)
c           --- compute the squared magnitude of "delta" - "gamma" 
            thetac = dcos(theta)
            dg_mag_sq = 0.0d0
            do icord = 1,3
               dvect(icord) = bvect(icord)*(1.0d0-(b_len/a_len)*thetac)
               dg_mag_sq = dg_mag_sq + (dvect(icord)- gvect(icord))**2
            enddo
c           --- compute the "ytilda" vector and magnitude as a cross 
c           --- of "beta" and "gamma"
            call twh_crossproduct(bvect(1),bvect(2),bvect(3),gvect(1)
     &           ,gvect(2),gvect(3),ytilda(1),ytilda(2),ytilda(3))
            bg_mag = 0.0d0
            do icord = 1,3
               bg_mag = bg_mag + ytilda(icord)*ytilda(icord)
            enddo
            bg_mag = dsqrt(bg_mag)
c           --- turn "ytilda" into a unit vector
            do icord = 1,3
               ytilda(icord) = ytilda(icord)/bg_mag
            enddo
c           --- compute "xtilda"
            call twh_crossproduct(bvect(1),bvect(2),bvect(3),ytilda(1)
     &           ,ytilda(2),ytilda(3),xtilda(1),xtilda(2),xtilda(3))
c           --- turn "xtilda" into a unit vector
            do icord = 1,3
               xtilda(icord) = xtilda(icord)/a_len
            enddo
c           --- compute phic
            phic = (c_len*c_len - r_len*r_len - dg_mag_sq)*a_len /
     &           (2.0d0*r_len*bg_mag)
c           --- compute the two potential values of phi
            phi(1) = twh_arccos(phic)
            phi(2) = twopi - phi(1)

            do iphi = 1,2
c              --- compute the two possible positions of iumiss
               do icord = 1,3
                  mvect(iphi,icord) = dvect(icord) 
     &                 + r_len*(xtilda(icord)*dcos(phi(iphi)) 
     &                 + ytilda(icord)*dsin(phi(iphi)))
               enddo
            enddo

c           --- store the atoms for the iunit-imiss-iother angle
            xatmp(1) = bvect(1)
            yatmp(1) = bvect(2)
            zatmp(1) = bvect(3)
            xatmp(3) = gvect(1)
            yatmp(3) = gvect(2)
            zatmp(3) = gvect(3)

c           --- store the iufrom-iunit,imiss-iother torsion atoms
            xttmp(1) = 0.0d0
            yttmp(1) = 0.0d0
            zttmp(1) = 0.0d0
            xttmp(2) = bvect(1)
            yttmp(2) = bvect(2)
            zttmp(2) = bvect(3)
            xttmp(4) = gvect(1)
            yttmp(4) = gvect(2)
            zttmp(4) = gvect(3)
c           --- compute the index for this torsion
            call twh_intor(GLB_GET,imolty,iufrom,intorsion)
            do itor = 1,intorsion
               call twh_ijtor2(GLB_GET,imolty,iufrom,itor,ijtor2test)
               call twh_ijtor3(GLB_GET,imolty,iufrom,itor,ijtor3test)
               call twh_ijtor4(GLB_GET,imolty,iufrom,itor,ijtor4test)
               if ( ijtor2test .eq. iunit .and.
     &              ijtor3test .eq. iumiss .and.
     &              ijtor4test .eq. iuother ) then
                  torind = itor
               endif
            enddo

            do iphi = 1,2
               vbias = 0.0d0
               xatmp(2) = mvect(iphi,1)
               yatmp(2) = mvect(iphi,2)
               zatmp(2) = mvect(iphi,3)
               vbias = vbias + twh_vangle(lfailure,.false.,lovrlap
     &              ,bxdflag,ibox,imolty
     &              ,iunit,index,1,thetac,dumb,dumc,xatmp,yatmp,zatmp)
               if ( lfailure ) return
c              --- compute the torsion energy
               xttmp(3) = mvect(iphi,1)
               yttmp(3) = mvect(iphi,2)
               zttmp(3) = mvect(iphi,3)
               vbias = vbias + twh_vtorsion(lfailure,.false.
     &              ,bxdflag,ibox,1
     &              ,imolty,iufrom,torind,xttmp,yttmp,zttmp,duma
     &              ,lovrlap)
               if ( lfailure ) return
            enddo

c           --- convert the vbias energy into a probability
            if ( .not. lovrlap ) then
               twh_febias = twh_expon(-beta*vbias)
            endif

         elseif ( twobbias_style .eq. 'autofit gaussian' .or.
     &           twobbias_style 
     &           .eq. 'self adapting gaussian using 1-3 distance' ) then
c           --- bias is based on a gaussian distribution on the 
c           --- iunit-iuother distance
c           --- compute the iunit-iuother length from the cords
            call twh_mimage(bxdflag,ibox,xcord(2),ycord(2),zcord(2)
     &           ,xcord(1),ycord(1),zcord(1),rxuij,ryuij,rzuij)
            length = twh_distance(rxuij,ryuij,rzuij)
c           --- set the biasing mean and standard deviation
            call twh_two_bond_bias_info(GLB_GET,imolty,iunit,index,1
     &           ,mean)
            call twh_two_bond_bias_sdev_multiplier(GLB_GET
     &           ,two_bond_bias_sdev_multiplier)
            call twh_two_bond_bias_info(GLB_GET,imolty,iunit,index,2
     &           ,dvalue)
            stddev = dvalue*two_bond_bias_sdev_multiplier
            lolimit = 0.0d0
            hilimit = minboxdim
c           --- get the gaussian probability
            twh_febias = twh_gaussprob(length,mean,stddev
     &           ,lolimit,hilimit)
         else
            write(6,*) 'FEBIAS: '
     &           ,'unknown two_bond_fe_bias_style'
            write(6,*) twobbias_style
            lfailure = .true.
            return
         endif
      elseif ( biastype .eq. 2 ) then
c        --- three bond fixed endpoint biasing
         call twh_three_bond_fe_bias_style(GLB_GET,threebbias_style)
         if ( threebbias_style .eq. 'none' ) then
c           --- nothing to do here
         elseif ( threebbias_style .eq. 
     &           'analytic using max and min 2-4 distance' 
     &           .or. threebbias_style .eq.
     &           'autofit gaussian using max and min 2-4 distance'
     &           ) then
c           --- find the rest of the atoms in the biasing torsion
            call twh_ijtor2(GLB_GET,imolty,iunit,index,iunext)
            call twh_ijtor3(GLB_GET,imolty,iunit,index,iuthree)
            call twh_ijtor4(GLB_GET,imolty,iunit,index,iuother)
c           --- get all the lengths we will need
c           --- botlen is the iuother - iufrom distance
            call twh_mimage(bxdflag,ibox,xcord(2),ycord(2),zcord(2)
     &           ,xcord(3),ycord(3),zcord(3),rxuij,ryuij,rzuij)
            botlensq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
c           --- midlen is the iuother - iunit distance
            call twh_mimage(bxdflag,ibox,xcord(2),ycord(2),zcord(2)
     &           ,xcord(1),ycord(1),zcord(1),rxuij,ryuij,rzuij)
            midlen = twh_distance(rxuij,ryuij,rzuij)
c           --- prelen is the iunit - iufrom distance
            call twh_mimage(bxdflag,ibox,xcord(3),ycord(3),zcord(3)
     &           ,xcord(1),ycord(1),zcord(1),rxuij,ryuij,rzuij)
            prelen = twh_distance(rxuij,ryuij,rzuij)
c           --- nexlen is the equilibrium bond length between 
c           --- iunit-iunext
            call twh_invib(GLB_GET,imolty,iunit,invibration)
            do ivib = 1,invibration
               call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
               if ( ijvibtest .eq. iunext ) then
                  call twh_itvib(GLB_GET,imolty,iunit,ivib,itvibtest)
                  call twh_vibcoeff(GLB_GET,itvibtest,0,nexlen)
                  goto 10
               endif
            enddo
            write(6,*) 'FEBIAS: bond length not found in febias'
            lfailure = .true.
            return
 10         continue
c           --- find the iufrom-iunit,iunext angle 
            call twh_inben(GLB_GET,imolty,iufrom,inbend)
            do ibend = 1,inbend
               call twh_ijben2(GLB_GET,imolty,iufrom,ibend,ijben2test)
               call twh_ijben3(GLB_GET,imolty,iufrom,ibend,ijben3test)
               if ( ijben2test .eq. iunit .and. ijben3test .eq. iunext 
     &              ) then
                  call twh_itben(GLB_GET,imolty,iufrom,ibend,itemp)
                  call twh_bencoeff(GLB_GET,itemp,0,bigang)
                  goto 20
               endif
            enddo
c           --- if we get here then something went wrong
            write(6,*) 'FEBIAS: angle not found'
            write(6,*) 'iufrom,iunit,iunext',iufrom,iunit,iunext
            write(6,*) 'iuthree,iuother',iuthree,iuother
            write(6,*) 'torsion index',index
            lfailure = .true.
            return
 20         continue
c           --- determine the iufrom-iunit-iuother angle using the law 
c           --- of cosines
            lowang = twh_arccos( (botlensq-midlen*midlen-prelen*prelen)
     &           /(-2.0d0*midlen*prelen) )
c           --- determine the minimum iunext-iunit-iuother angle
            minang = bigang-lowang
c           --- determine the minimum iuother-iunext distance using the
c           --- law of cosines
            minlen = dsqrt( midlen*midlen + nexlen*nexlen 
     &           - 2.0d0*midlen*nexlen*dcos(minang) )
c           --- determine the maximum iunext-iunit-iuother angle
            maxang = (twopi-bigang)-lowang
c           --- determine the maximum iuother-iunext distance using the
c           --- law of cosines
            maxlen = dsqrt( midlen*midlen + nexlen*nexlen 
     &           - 2.0d0*midlen*nexlen*dcos(maxang) )

c           --- want to bias using the iunext-iuthree-iuother angle
c           --- potential
            call twh_inben(GLB_GET,imolty,iunext,inbend)
            do ibend = 1,inbend
               call twh_ijben2(GLB_GET,imolty,iunext,ibend,ijben2test)
               call twh_ijben3(GLB_GET,imolty,iunext,ibend,ijben3test)
               if ( ijben2test .eq. iuthree .and.
     &              ijben3test .eq. iuother ) then
                  goodind = ibend
                  goto 30
               endif
            enddo
            write(6,*) 'FEBIAS: angle type not found'
            lfailure = .true.
            return
 30         continue

c           --- now we are finally to the part where we do the biasing
            if ( threebbias_style .eq. 
     &           'analytic using max and min 2-4 distance' ) then
c              --- compute equilibrium bond lengths for this angle
               call twh_ibtben1(GLB_GET,imolty,iunext,goodind,btone)
               call twh_ibtben2(GLB_GET,imolty,iunext,goodind,bttwo)
               call twh_vibcoeff(GLB_GET,btone,0,lenone)
               call twh_vibcoeff(GLB_GET,bttwo,0,lentwo)
c              --- get the equilibrium angle
               call twh_itben(GLB_GET,imolty,iunext,goodind,itemp)
               call twh_bencoeff(GLB_GET,itemp,0,eqangle)
c              --- compute the equilibrium length from law of cosine
               eqlength = dsqrt( lenone*lenone + lentwo*lentwo 
     &              - 2.0d0*lenone*lentwo*dcos(eqangle) )
c              --- intialize vbias
               vbias = 0.0d0
c              --- compute biasing potential
               if ( eqlength .lt. minlen ) then
                  if ( minlen .lt. lenone+lentwo ) then
c                    --- put the strain on the angle using the 
c                    --- equlibrium bond lengths
c                    --- compute strained angle using law of cosine
                     thetac = (minlen*minlen -lenone*lenone 
     &                    -lentwo*lentwo)/ (-2.0d0*lenone*lentwo)
                  else
c                    --- crank the angle to Pi, and stretch bond lengths
                     thetac = -1.0d0
c                    --- compute the scaled bond lengths
                     scale = minlen/(lenone+lentwo)
                     lenone = scale*lenone
                     lentwo = scale*lentwo
                     vbias = vbias
     &                    + twh_vbond(lfailure,.false.,bxdflag,ibox
     &                    ,btone,lenone,imolty,iunext,iuthree,lovrlap) 
     &                    + twh_vbond(lfailure,.false.,bxdflag,ibox
     &                    ,bttwo,lentwo,imolty,iuthree,iuother,lovrlap) 
                     if ( lfailure ) return
                  endif
c                 --- compute biasing energy from the angle
                  vbias = vbias + twh_vangle(lfailure,.false.,lovrlap
     &                 ,bxdflag,ibox,imolty,iunext,goodind,0,thetac
     &                 ,lenone,lentwo,ddum3,ddum3,ddum3)
                  if ( lfailure ) return
               elseif ( eqlength .le. maxlen ) then
c                 --- the equilbrium length for this angle falls in the 
c                 --- range accessible from this site.  No biasing
c                 --- energy
                  vbias = 0.0d0
               else
c                 --- the desired equilibrium length is larger than the 
c                 --- maximum posible length from this position. use a 
c                 --- bias energy from squeezing the angle down so it 
c                 --- gives maxlen
c                 --- compute strained angle using law of cosine
                  thetac = (maxlen*maxlen - lenone*lenone 
     &                 - lentwo*lentwo) / (-2.0d0*lenone*lentwo)
c                 --- compute biasing energy from the angle
                  vbias = vbias + twh_vangle(lfailure,.false.,lovrlap
     &                 ,bxdflag,ibox,imolty,iunext,goodind,0,thetac
     &                 ,lenone,lentwo,ddum3,ddum3,ddum3)
                  if ( lfailure ) return
               endif
c              --- compute febias from the biasing energy
               if ( .not. lovrlap ) then
                  twh_febias = twh_expon(-beta*vbias)
               endif
            elseif ( threebbias_style .eq.
     &              'autofit gaussian using max and min 2-4 distance'
     &              ) then
c              --- assign the gaussian parameters
               call twh_two_bond_bias_info(GLB_GET,imolty,iunext
     &              ,goodind,1,mean)
               call twh_two_bond_bias_info(GLB_GET,imolty,iunext
     &              ,goodind,2,dvalue)
               call twh_three_bond_bias_sdev_multiplier(GLB_GET,mult)
               stddev = dvalue*mult
               lolimit = 0.0d0
c              --- want the integrated probability distribution from 
c              --- minlen to maxlen on bounds [0,Inf)
               call twh_integratedgaussprob(lfailure,LIM_VALUE
     &              ,LIM_POSINF,minlen
     &              ,maxlen,mean,stddev,lolimit,hilimit,twh_febias)
               if ( lfailure ) return
            else
               write(6,*) 'FEBIAS: unknown three_bond_style'
               write(6,*) threebbias_style
               lfailure = .true.
               return
            endif
         elseif ( threebbias_style .eq.
     &           'self adapting gaussian using 1-4 distance' ) then
c           --- assign the gaussian parameters
            call twh_three_bond_bias_info(GLB_GET,imolty,iunit,index
     &           ,1,mean)
            call twh_three_bond_bias_info(GLB_GET,imolty,iunit,index
     &           ,2,dvalue)
            call twh_three_bond_bias_sdev_multiplier(GLB_GET,mult)
            stddev = dvalue*mult
            lolimit = 0.0d0
            hilimit = minboxdim
            call twh_mimage(bxdflag,ibox,xcord(2),ycord(2),zcord(2)
     &           ,xcord(1),ycord(1),zcord(1),rxuij,ryuij,rzuij)
            length = twh_distance(rxuij,ryuij,rzuij)
            twh_febias = twh_gaussprob(length,mean,stddev,lolimit
     &           ,hilimit)
         else
            write(6,*) 'FEBIAS: unknown '
     &           ,'three_bond_fixed_endpoint_bias_style'
            write(6,*) threebbias_style
            lfailure = .true.
            return
         endif
      else
         write(6,*) 'FEBIAS: biastype ',biastype,' not implemented'
         lfailure = .true.
      endif

      if ( lovrlap ) then
c        --- if we had an overlap set bias to the minbias
         twh_febias = minbias
      elseif ( twh_febias .lt. minbias ) then
c        --- make sure twh_febias is larger than minbias
         twh_febias = minbias
      endif

c      write(6,*) 'finished FEBIAS'
      return
      end
