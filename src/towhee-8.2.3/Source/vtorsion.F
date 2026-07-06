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
      function twh_vtorsion(lfailure,lreport,bxdflag,ibox,cordstyle
     &     ,imolty,unit_1,torind,xcord,ycord,zcord,phic,loverlap)
c     ******************************************************************
c     * Computes the torsional energy plus cross terms and 1-4 nonbond *
c     * for a particular torsion type.                                 *
c     * ibox: the box the molecule is in - needed for minimum image    *
c     * cordstyle: tells the function what kind of angle we want to    *
c     *            work with and whether it is passed to the function  *
c     *            directly, or should be computed from the coordinates*
c     *            -1: normal dihedral passed as an angle              *
c     *             0: special case, compute the normal dihedral angle *
c     *                from the coordinates and then just return that  *
c     *                angle instead of returning the torsional energy *
c     *                This is used for some of the non-uniform        *
c     *                sampling.                                       *
c     *             1: normal dihedral that must be computed from the  *
c     *                coordinates.  (1-2 x 2-3 * 2-3 x 3-4)           *
c     * imolty: molecule type                                          *
c     * unit_1: unit numbers of the starting atom of the torsion       *
c     * torind: the index of the torsion                               *
c     * x,y,z cord : coordinates of the four atoms in the torsion      *
c     * phic: cosine of the torsional angle                            *
c     * loverlap: true for hard overlap, note it is initialized only in*
c     *  the calling routine.  Do not initialize in any of the vxx     *
c     *  functions                                                     *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_DISTANCE
#define FUNCTION_EXPON
#define FUNCTION_DOTPRODUCT
#define FUNCTION_VTORSION
#define FUNCTION_VTWOBODY
#include "functions.h"
#define DEBUG_VTORSION 0
c     --- variables passed to and from function
      logical lfailure,lreport,loverlap
      integer bxdflag,ibox,imolty,torind,unit_1,cordstyle
      double precision phic
      double precision xcord,ycord,zcord
      dimension xcord(4),ycord(4),zcord(4)
c     --- local logicals
      logical lclasstwo
c     --- local character strings
      character*(FFNAMELEN) testname,torstring
      dimension testname(4),torstring(4)
c     --- local integer scalars
      integer itor,ntii,ntjj,bij,bjk,bkl,aijk,ajkl,torstyle,itest
      integer torffnum,itmp,testtype,invib_two,invib_thr,torffraw
      integer abcord,ivec,unit_2,unit_3,unit_4,index,torofcode
      integer ntorloop,icoeff
c     --- local integer arrays
      integer atype,btype
      dimension btype(3)
      dimension atype(2)
c     --- local double precision scalars
      double precision phi,rijsq,cosphi,cos2phi,cos3phi,rij
      double precision phi_eq,phi_diff
      double precision xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2
      double precision xcc,ycc,zcc,tcc
      double precision torscale,rxij,ryij,rzij
      double precision xone,xtwo,yone,ytwo,zone,ztwo,cosshiftphi
      double precision totbond,dv,onepi,vonefour
      double precision vbontor,vangtor,vaator,vbb
      double precision eqlen_bij,eqlen_bjk,eqlen_bkl
      double precision eqang_ijk,eqang_jkl
#if DEBUG_VTORSION
      double precision radtodeg
#endif
c     --- local double precision arrays
      double precision xvec,yvec,zvec
      dimension xvec(4),yvec(4),zvec(4)
      double precision length
      dimension length(4)
      double precision angle
      dimension angle(2)
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)
      double precision torcoeff
      dimension torcoeff(MINTORCOEFF:MAXTORCOEFF)
      
c     --- initialize vtorsion
      twh_vtorsion = 0.0d0
      lclasstwo = .false.
c     --- determine the torsion forfield number and style
      abcord = iabs(cordstyle)
      if ( abcord .eq. 0 ) then
c        --- just computing an angle and returning, set a dummy ff
         torffnum = 1
      elseif ( abcord .eq. 1 ) then
c        --- regular torsion - use ittor
         call twh_ittor(GLB_GET,imolty,unit_1,torind,torffraw)
         torffnum = abs(torffraw)
      else
         write(6,*) 'invalid cordstyle in VTORSION',cordstyle
         lfailure = .true.
         return
      endif
c     --- get the torsion style code
      call twh_torstyle(GLB_GET,torffnum,torstyle)
#if DEBUG_VTORSION
      write(6,*) 'imolty,unit1,torind:',imolty,unit_1,torind
      call twh_constant_radtodeg(GLB_GET,radtodeg)
#endif
      if ( cordstyle .lt. 0 ) then
c        --- use the angle that was passed to the fucntion
         if ( torstyle .eq. TOR_STYLE_COMP_X_AD ) then
            write(6,*) 'cannot pass angles to VTORSION for Class 2'
            lfailure = .true.
            return
         endif
      elseif ( cordstyle .eq. 0 ) then
c        --- want the standard angle, but need it over the full range
c        --- of -Pi to Pi so we can use this as an adjustment in the
c        --- non-uniform selection of torsional angles.
         do ivec = 1,3
            call twh_mimage(bxdflag,ibox,xcord(ivec+1),ycord(ivec+1)
     &           ,zcord(ivec+1),xcord(ivec),ycord(ivec),zcord(ivec)
     &           ,xvec(ivec),yvec(ivec),zvec(ivec))
         enddo
c        --- we also need a fourth vector for the triple product
         call twh_mimage(bxdflag,ibox,xcord(2),ycord(2),zcord(2)
     &        ,xcord(3),ycord(3),zcord(3),xvec(4),yvec(4),zvec(4))
c        --- calculate cross products 1-2 x 2-3 
         call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(2),yvec(2)
     &        ,zvec(2),xaa1,yaa1,zaa1)
c        --- calculate cross products 3-4 x 2-3
         call twh_crossproduct(xvec(3),yvec(3),zvec(3),xvec(2),yvec(2)
     &        ,zvec(2),xaa2,yaa2,zaa2)
c        --- calculate dot product of cross products
         phic = twh_dotproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2) 
         phic = twh_arccos(-phic)
c        --- calculate cross product of cross products
         call twh_crossproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2
     &        ,xcc,ycc,zcc)
c        --- calculate scalar triple product
         tcc = xcc*(xvec(4)) + ycc*(yvec(4)) + zcc*(zvec(4))
c        --- compute the angle over a range -pi to +pi
         if ( tcc .gt. 0.0d0 ) phic = -phic
         twh_vtorsion = phic
c        --- not actually computing any energies in this case
         return
      elseif ( cordstyle .eq. 1 ) then
c        --- compute the dihedral angle
c        --- defined so that a cis bond is at 0.0 phi
         do ivec = 1,3
            call twh_mimage(bxdflag,ibox,xcord(ivec+1),ycord(ivec+1)
     &           ,zcord(ivec+1),xcord(ivec),ycord(ivec),zcord(ivec)
     &           ,xvec(ivec),yvec(ivec),zvec(ivec))
         enddo
c        --- calculate cross products 1-2 x 2-3 
         call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(2),yvec(2)
     &        ,zvec(2),xaa1,yaa1,zaa1)
c        --- calculate cross products 3-4 x 2-3
         call twh_crossproduct(xvec(3),yvec(3),zvec(3),xvec(2),yvec(2)
     &        ,zvec(2),xaa2,yaa2,zaa2)
c        --- calculate dot product of cross products
         phic = -twh_dotproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2)
      endif
c     --- retrieve all torcoeff values
      do icoeff = MINTORCOEFF,MAXTORCOEFF
         call twh_torcoeff(GLB_GET,torffnum,icoeff,torcoeff(icoeff))
      enddo
      if ( torstyle .eq. TOR_STYLE_SIM_HARM ) then
c        --- simple harmonic potential 
c        --- used for TraPPE double bonds (cis vs trans) 
         phi = twh_arccos(phic)
         twh_vtorsion = torcoeff(0)
     &        *(phi-torcoeff(1))**2.0d0
      elseif ( torstyle .eq. TOR_STYLE_OPLS_COS  ) then
c        --- OPLS style cosine series 
c        --- used in the OPLS-AA force field
         phi = twh_arccos(phic)
         twh_vtorsion = torcoeff(1)*(1.d0 + phic) +
     &        torcoeff(2)*(1.d0 - dcos(2.0d0*phi)) +
     &        torcoeff(3)*(1.d0 + dcos(3.0d0*phi))
#if DEBUG_VTORSION
         write(6,*) 'phi,vtorsion:',phi*radtodeg
     &        ,twh_vtorsion/CON_KCALTOK
#endif
      elseif ( torstyle .eq. TOR_STYLE_GCA_COS ) then
c        --- Gromos/Charmm/Amber cosine potential
         phi = twh_arccos(phic)
         call twh_ntorloop(GLB_GET,torffnum,ntorloop)
         do itor = 1,ntorloop
            index = 3*(itor-1)
            twh_vtorsion = twh_vtorsion 
     &           + torcoeff(index+1)* 
     &           ( 1.0d0 + dcos(torcoeff(index+2)*phi
     &           - torcoeff(index+3))
     &           )
         enddo
      elseif ( torstyle .eq. TOR_STYLE_GCA_CUB ) then
c        --- Gromos/Charmm/Amber cosine potential
c        --- plus a harmonic term that Charmm would traditionally call
c        --- an improper torsion despite the fact that the bonding 
c        --- pattern is actually that of a regular torsion
         phi = twh_arccos(phic)
         call twh_ntorloop(GLB_GET,torffnum,ntorloop)
         do itor = 1,ntorloop
            index = 3*(itor-1)
            twh_vtorsion = twh_vtorsion 
     &           + torcoeff(index+1)* 
     &           ( 1.0d0 + dcos(torcoeff(index+2)*phi
     &           - torcoeff(index+3))
     &           )
         enddo
c        --- reset index
         call twh_ntorloop(GLB_GET,torffnum,ntorloop)
         index = 3*ntorloop
c        --- add in the harmonic term
         twh_vtorsion = twh_vtorsion 
     &        + torcoeff(index+1)
     &        *(phi-torcoeff(index+2))**2.0d0

      elseif ( torstyle .eq. TOR_STYLE_COMP_X_AD .or.
     &        torstyle .eq. TOR_STYLE_COMP_R_AD ) then
c        --- Compass cosine series using the clever atom-based terms 
c        --- to sort out all of the cross terms
c        --- torstyle X includes cross terms
c        --- torstyle R does not include cross terms
         phi = twh_arccos(phic)
         cosphi = phic
         cos2phi = dcos(2.0d0*phi)
         cos3phi = dcos(3.0d0*phi)

c        --- simple torsion part
         twh_vtorsion = torcoeff(0)*(1.0d0 - cosphi)
     &        + torcoeff(1)*(1.0d0 - cos2phi)
     &        + torcoeff(2)*(1.0d0 - cos3phi)

         if ( torstyle .eq. TOR_STYLE_COMP_X_AD ) then
            lclasstwo = .true.
c           --- need all of the unit numbers to determine the
c           --- order of the torsion
            call twh_ijtor2(GLB_GET,imolty,unit_1,torind,unit_2)
            call twh_ijtor3(GLB_GET,imolty,unit_1,torind,unit_3)
            call twh_ijtor4(GLB_GET,imolty,unit_1,torind,unit_4)

c           --- Compute the cross terms
c           --- Determine whether we are in normal or inverse order
            call twh_ntype(GLB_GET,imolty,unit_1,testtype)
            call twh_ctname(GLB_GET,testtype,testname(1))
            call twh_ntype(GLB_GET,imolty,unit_2,testtype)
            call twh_ctname(GLB_GET,testtype,testname(2))
            call twh_ntype(GLB_GET,imolty,unit_3,testtype)
            call twh_ctname(GLB_GET,testtype,testname(3))
            call twh_ntype(GLB_GET,imolty,unit_4,testtype)
            call twh_ctname(GLB_GET,testtype,testname(4))
            do itest = 1,4
               call twh_torstring(GLB_GET,torffnum,itest
     &              ,torstring(itest))
            enddo
            if ( testname(1) .eq. torstring(1) .and.
     &           testname(2) .eq. torstring(2) .and.
     &           testname(3) .eq. torstring(3) .and.
     &           testname(4) .eq. torstring(4) ) then
c              --- normal order
               bij = 1
               bjk = 2
               bkl = 3
               aijk = 1
               ajkl = 2
            elseif ( testname(4) .eq. torstring(1) .and.
     &              testname(3) .eq. torstring(2) .and.
     &              testname(2) .eq. torstring(3) .and.
     &              testname(1) .eq. torstring(4) ) then
c              --- inverse order
               bij = 3
               bjk = 2
               bkl = 1
               aijk = 2
               ajkl = 1
            else
               write(6,*) 'VTORSION: torffnum',torffnum
               write(6,*) 'ctname'
     &              ,testname(1),testname(2),testname(3),testname(4)
               do itest = 1,4
                  call twh_torstring(GLB_GET,torffnum,itest
     &                 ,torstring(itest))
               enddo
               write(6,*) 'torstring',(torstring(itest),itest=1,4)
               write(6,*) 'torsion does not match up in VTORSION'
               lfailure = .true.
               return
            endif

c           --- compute bond lengths using *vec (already mimaged)
            do itmp = 1,3
               rxij = xvec(itmp)
               ryij = yvec(itmp)
               rzij = zvec(itmp)
               length(itmp) = twh_distance(rxij,ryij,rzij)
            enddo
c           --- compute angles using *vec (already mimaged)
            do itmp = 1,2
               xone = -xvec(itmp)
               xtwo = xvec(itmp+1)
               yone = -yvec(itmp)
               ytwo = yvec(itmp+1)
               zone = -zvec(itmp)
               ztwo = zvec(itmp+1)
               angle(itmp) = twh_dotproduct(xone,yone,zone,xtwo,ytwo
     &              ,ztwo,daa1,daa2)
            enddo
c           --- store the bond and angle types
            call twh_ibttor1(GLB_GET,imolty,unit_1,torind,btype(1))
            call twh_ibttor2(GLB_GET,imolty,unit_1,torind,btype(2))
            call twh_ibttor3(GLB_GET,imolty,unit_1,torind,btype(3))
            call twh_iattor1(GLB_GET,imolty,unit_1,torind,atype(1))
            call twh_iattor2(GLB_GET,imolty,unit_1,torind,atype(2))
c           --- retrieve the equilibrium bond lengths
            call twh_vibcoeff(GLB_GET,btype(bij),0,eqlen_bij)
            call twh_vibcoeff(GLB_GET,btype(bjk),0,eqlen_bjk)
            call twh_vibcoeff(GLB_GET,btype(bkl),0,eqlen_bkl)
c           --- retrieve the equilibrium angles
            call twh_bencoeff(GLB_GET,atype(aijk),0,eqang_ijk)
            call twh_bencoeff(GLB_GET,atype(ajkl),0,eqang_jkl)
c           --- convert angles from cosine to normal
            angle(1) = twh_arccos(angle(1))
            angle(2) = twh_arccos(angle(2))
c           --- bond-torsion terms
            vbontor = 
c           --- bond(ij) - torsion cross
     &           (length(bij) - eqlen_bij)*
     &           (torcoeff(3)*cosphi 
     &           + torcoeff(4)*cos2phi
     &           + torcoeff(5)*cos3phi)
c           --- bond(jk) - torsion cross
     &           + (length(bjk) - eqlen_bjk)*
     &           (torcoeff(6)*cosphi 
     &           + torcoeff(7)*cos2phi
     &           + torcoeff(8)*cos3phi)
c           --- bond(kl) - torsion cross
     &           + (length(bkl) - eqlen_bkl)*
     &           (torcoeff(9)*cosphi 
     &           + torcoeff(10)*cos2phi
     &           + torcoeff(11)*cos3phi)
c           --- angle-torsion terms
            vangtor =
c           --- angle(ijk) - torsion cross
     &           (angle(aijk) - eqang_ijk)*
     &           (torcoeff(12)*cosphi 
     &           + torcoeff(13)*cos2phi
     &           + torcoeff(14)*cos3phi)
c           --- angle(jkl) - torsion cross
     &           + (angle(ajkl) - eqang_jkl)*
     &           (torcoeff(15)*cosphi 
     &           + torcoeff(16)*cos2phi
     &           + torcoeff(17)*cos3phi)
c           --- angle(ijk) - angle(jkl) - torsion cross
            vaator = torcoeff(18)*
     &           (angle(aijk) - eqang_ijk)*
     &           (angle(ajkl) - eqang_jkl)*cosphi
c           --- bond(ij) - bond(kl) cross
            vbb =  torcoeff(19)*
     &           (length(bij) - eqlen_bij)*(length(bkl) - eqlen_bkl)
         endif

      elseif ( torstyle .eq. TOR_STYLE_TRAPPE_COS ) then
c        --- TraPPE simple cosine function
c        --- trappe is defined so that trans is at 0.0
         call twh_constant_pi(GLB_GET,onepi)
         phi = twh_arccos(phic)-onepi
         twh_vtorsion = torcoeff(0)
     &        *(1.0d0 - dcos(2.0d0*phi + torcoeff(1)))

      elseif ( torstyle .eq. TOR_STYLE_NB ) then
c        --- nonbonded potential only -or absolutely no torsional energy
c        --- any scale factor is accounted for in ffnonbond
c        --- nothing computed here, all down below

      elseif ( torstyle .eq. TOR_STYLE_COMP_X_EX ) then
c        --- Compass cosine series plus cross terms using an explicit
c        --- parameter declaration and negative typing to determine 
c        --- the order
         phi = twh_arccos(phic)
         cosshiftphi = dcos(phi - torcoeff(3))
c        --- simple torsion part
         twh_vtorsion = 
     &        torcoeff(0)
     &        *(1.0d0 - cosshiftphi )
     &        + torcoeff(1)
     &        *(1.0d0 - dcos(2.0d0*phi - torcoeff(4)) )
     &        + torcoeff(2)
     &        *(1.0d0 - dcos(3.0d0*phi - torcoeff(5)) )

c        --- the unshifted values of the cosines are used in the cross
c        --- terms
         cosphi = dcos(phi)
         cos2phi = dcos(2.0d0*phi)
         cos3phi = dcos(3.0d0*phi)
c        --- compute bond lengths using *vec (already mimaged)
         do itmp = 1,3
            rxij = xvec(itmp)
            ryij = yvec(itmp)
            rzij = zvec(itmp)
            length(itmp) = twh_distance(rxij,ryij,rzij)
         enddo
c        --- compute angles using *vec (already mimaged)
         do itmp = 1,2
            xone = -xvec(itmp)
            xtwo = xvec(itmp+1)
            yone = -yvec(itmp)
            ytwo = yvec(itmp+1)
            zone = -zvec(itmp)
            ztwo = zvec(itmp+1)
            angle(itmp) = twh_arccos(twh_dotproduct(xone,yone,zone
     &           ,xtwo,ytwo,ztwo,daa1,daa2))
         enddo
c        --- Compute the cross terms
c        --- Determine whether we are in normal or inverse order
c        --- use the raw torff number that includes a sign
         if ( torffraw .gt. 0 ) then
c           --- normal order, don't change a thing
            bij = 1
            bjk = 2
            bkl = 3
            aijk = 1
            ajkl = 2
         else
c           --- inverse order, switch all of the bond lengths and the
c           --- angles
            bij = 3
            bjk = 2
            bkl = 1
            aijk = 2
            ajkl = 1
         endif
         lclasstwo = .true.
         vbontor =
c        --- bond(jk) - torsion cross
     &        (length(bjk) - torcoeff(9))*
     &        (torcoeff(6)*cosphi 
     &        + torcoeff(7)*cos2phi
     &        + torcoeff(8)*cos3phi)
c        --- bond(ij) - torsion cross
     &        + (length(bij) - torcoeff(13))*
     &        (torcoeff(10)*cosphi 
     &        + torcoeff(11)*cos2phi
     &        + torcoeff(12)*cos3phi)
c        --- bond(kl) - torsion cross
     &        + (length(bkl) - torcoeff(17))*
     &        (torcoeff(14)*cosphi 
     &        + torcoeff(15)*cos2phi
     &        + torcoeff(16)*cos3phi)
         vangtor =
c        --- angle(ijk) - torsion cross
     &        (angle(aijk) - torcoeff(24))*
     &        (torcoeff(18)*cosphi 
     &        + torcoeff(19)*cos2phi
     &        + torcoeff(20)*cos3phi)
c        --- angle(jkl) - torsion cross
     &        + (angle(ajkl) - torcoeff(25))*
     &        (torcoeff(21)*cosphi 
     &        + torcoeff(22)*cos2phi
     &        + torcoeff(23)*cos3phi)
c        --- angle(ijk) - angle(jkl) - torsion cross
         vaator = torcoeff(26)*
     &        (angle(aijk) - torcoeff(27))*
     &        (angle(ajkl) - torcoeff(28))*cosshiftphi
c        --- bond(ij) - bond(kl) cross
         vbb = torcoeff(29)*
     &        (length(bij) - torcoeff(30))*
     &        (length(bkl) - torcoeff(31))

      elseif ( torstyle .eq. TOR_STYLE_COS_POW ) then
c        --- cosine power series
         call twh_ntorloop(GLB_GET,torffnum,ntorloop)
         do itor = 0,ntorloop
            twh_vtorsion = twh_vtorsion 
     &           + torcoeff(itor)*(phic**itor)
         enddo
      elseif ( torstyle .eq. TOR_STYLE_O_OPLS_C ) then
c        --- old OPLS style cosine series 
c        --- used in the OPLS-UA force field
         phi = twh_arccos(phic)
         twh_vtorsion = torcoeff(0) + 
     &        torcoeff(1)*(1.d0 + phic) +
     &        torcoeff(2)*(1.d0 - dcos(2.0d0*phi)) +
     &        torcoeff(3)*(1.d0 + dcos(3.0d0*phi))
      elseif ( torstyle .eq. TOR_STYLE_SUM_COS ) then
c        --- Sum2003 style torsion potential
         phi = twh_arccos(phic)
         phi_eq = torcoeff(0)
         phi_diff = phi - phi_eq
         call twh_ntorloop(GLB_GET,torffnum,ntorloop)
         do itor = 0,ntorloop-1
            index = itor + 1
            twh_vtorsion = twh_vtorsion 
     &           + torcoeff(index)* 
     &           ( 1.0d0 - dcos(dble(itor)*phi_diff ))
         enddo
      elseif ( torstyle .eq. TOR_STYLE_O_OPLS_T ) then
c         two term torsional potential used in W.L. Jorgensen;
c         C.J. Swenson; J. Am. Chem. Soc. 107, 569-578 (1985).
         phi = twh_arccos(phic)
         twh_vtorsion = twh_vtorsion 
     &        + torcoeff(1)*(1.0d0 - phic)
     &        + torcoeff(2)*(1.0d0 - dcos(2.0d0*phi))
      elseif ( torstyle .eq. TOR_STYLE_UFF ) then
c        --- UFF style 1 - cosine potential
c        --- energy is divided by the total number of torsions that 
c        --- occur across these two atoms
         call twh_ijtor2(GLB_GET,imolty,unit_1,torind,unit_2)
         call twh_ijtor3(GLB_GET,imolty,unit_1,torind,unit_3)
         call twh_invib(GLB_GET,imolty,unit_2,invib_two)
         call twh_invib(GLB_GET,imolty,unit_3,invib_thr)
         totbond = dble((invib_two-1)*(invib_thr-1))
         phi = twh_arccos(phic)
         twh_vtorsion = torcoeff(1)
     &        *( 1.0d0 - torcoeff(2)
     &        *dcos( torcoeff(3)*phi) )/totbond
      elseif ( torstyle .eq. TOR_STYLE_DREIDING ) then
c         --- DREIDING k[ 1 - Cos ( n (phi - phi_0))]
         phi = twh_arccos(phic)
         twh_vtorsion = torcoeff(1)*(1.0d0 
     &        - dcos(torcoeff(2)
     &        *(phi - torcoeff(3))))
c        --- energy is divided by the total number of torsions that 
c        --- occur across these two atoms
         call twh_ijtor2(GLB_GET,imolty,unit_1,torind,unit_2)
         call twh_ijtor3(GLB_GET,imolty,unit_1,torind,unit_3)
         call twh_invib(GLB_GET,imolty,unit_2,invib_two)
         call twh_invib(GLB_GET,imolty,unit_3,invib_thr)
         totbond = dble((invib_two-1)*(invib_thr-1))
         twh_vtorsion = twh_vtorsion / totbond
      elseif ( torstyle .eq. TOR_STYLE_TF_COS ) then
c        --- 2-fold cosine
         phi = twh_arccos(phic)
         twh_vtorsion = torcoeff(1)
     &        *(1.0d0 - dcos(2.0d0*phi))
      elseif ( torstyle .eq. TOR_STYLE_TRAPPE_PLAN ) then
c        --- Potoff torsions from J. Phys. Chem. B 108 14130-14136 2004
c        --- note that there is a typo in equation 5 of that paper
c        --- the correct functional form is
c        U = c1[1 + cos(phi + f1)] + c2[1 - cos^2(phi)]
         phi = twh_arccos(phic)
         twh_vtorsion = torcoeff(1)
     &        *(1.0d0 + dcos(phi + torcoeff(3)) )
     &        + torcoeff(2)*(1.0d0 - phic*phic)
      elseif ( torstyle .eq. TOR_STYLE_SQUARE_WELL ) then
c        --- square well torsions originally implemented for Elli2002
c        --- can't just use nonbond below as has different parameters
c        --- compute the 1-4 distance squared
         call twh_mimage(bxdflag,ibox,xcord(4),ycord(4),zcord(4)
     &        ,xcord(1),ycord(1),zcord(1),rxij,ryij,rzij)
         rijsq = rxij*rxij + ryij*ryij + rzij*rzij
         rij = dsqrt(rijsq)
         if ( twh_cmp_lt(rij, torcoeff(1)) .or.
     &        twh_cmp_gt(rij, torcoeff(2)) ) then
c           --- hard overlap
            loverlap = .true.
         else
c           --- assign default energy
            twh_vtorsion = torcoeff(3)
         endif
      elseif ( torstyle .eq. TOR_STYLE_AMBER ) then
c        --- Gromos/Charmm/Amber cosine potential divided by the 
c        --- total number of torsions sharing these central atoms
         phi = twh_arccos(phic)
         call twh_ntorloop(GLB_GET,torffnum,ntorloop)
         do itor = 1,ntorloop
            index = 3*(itor-1)
            twh_vtorsion = twh_vtorsion 
     &           + torcoeff(index+1)* 
     &           ( 1.0d0 + dcos(torcoeff(index+2)*phi
     &           - torcoeff(index+3))
     &           )
         enddo
c        --- divide out by the total number of torsions on these atoms
         call twh_ijtor2(GLB_GET,imolty,unit_1,torind,unit_2)
         call twh_ijtor3(GLB_GET,imolty,unit_1,torind,unit_3)
         call twh_invib(GLB_GET,imolty,unit_2,invib_two)
         call twh_invib(GLB_GET,imolty,unit_3,invib_thr)
         totbond = dble((invib_two-1)*(invib_thr-1))
         twh_vtorsion = twh_vtorsion / totbond
      elseif ( torstyle .eq. TOR_STYLE_OPLS_FL  ) then
c        --- OPLS-AA Fluorocarbon 4 parameter cosine series 
c        --- used in the OPLS-AA force field.
c        --- J. Phys. Chem. A. v105, 4118-4125 (2001). 
         phi = twh_arccos(phic)
         twh_vtorsion = torcoeff(1)*(1.d0 + phic) +
     &        torcoeff(2)*(1.d0 - dcos(2.0d0*phi)) +
     &        torcoeff(3)*(1.d0 + dcos(3.0d0*phi)) +
     &        torcoeff(4)*(1.d0 - dcos(4.0d0*phi))
      elseif ( torstyle .eq. TOR_STYLE_MRD ) then
c        --- multiple rigid dihedral
         loverlap = .true.
         phi = twh_arccos(phic)
         call twh_ntorloop(GLB_GET,torffnum,ntorloop)
         do itor = 1,ntorloop
            if ( dabs( torcoeff(itor) - phi) .lt. 
     &           torcoeff(0) ) then
               loverlap = .false.
               twh_vtorsion = 0.0d0
            endif
         enddo
         call twh_constant_pi(GLB_GET,onepi)
      elseif ( torstyle .eq. TOR_STYLE_FL ) then
c        --- fluoroalkane from S.T. Cui; J.I. Siepmann; H.D. Cochran;
c        P.T. Cummings; "Intermolecular potentials and vapor-liquid
c        phase equilibria of perfluorinated alkanes",
c        Fluid Phase Equilib. 146 51-61 (1998)
         phi = twh_arccos(phic)
         twh_vtorsion = torcoeff(0)
     &        + torcoeff(1)*(1.0d0 - phic)
     &        + torcoeff(2)*(1.0d0 - dcos(3.0d0*phi))
     &        + torcoeff(3)*(1.0d0 - phic)**5
     &        + torcoeff(4)
     &        *twh_expon(-torcoeff(5)*phi**2)
      elseif ( torstyle .eq. TOR_STYLE_DRUKKER ) then
c        --- simple two term cosine originally from 
c        --- K. Drukker; G. Wu; G. C. Schatz;
c        --- "Model simulations of DNA denaturation dynamics"
c        --- J. Chem. Phys. 114 579-590 (2001)
         phi = twh_arccos(phi)
         twh_vtorsion = torcoeff(1)
     &        *(1.0d0 - dcos(phi - torcoeff(0)))
      else
         write(6,*) 'VTORSION: non-defined torsional type'
         lfailure = .true.
         return
      endif

      if ( lreport ) then
         call twh_report(GLB_INCR,REP_TOR_REG,twh_vtorsion)
      endif

      if ( lclasstwo ) then
         if ( lreport ) then
            call twh_report(GLB_INCR,REP_TOR_BT,vbontor)
            call twh_report(GLB_INCR,REP_TOR_AT,vangtor)
            call twh_report(GLB_INCR,REP_TOR_AAT,vaator)
            call twh_report(GLB_INCR,REP_TOR_BB,vbb)
         endif
         twh_vtorsion = twh_vtorsion + vbontor + vangtor + vaator + vbb
      endif

c     --- compute the 1-4 nonbonded terms, if they are used
      call twh_torofcode(GLB_GET,imolty,unit_1,torind,torofcode)
      if ( torofcode .ne. CODE_ZERO ) then
         if ( cordstyle .ne. 1 ) then
            write(6,*) 'cannot pass angle to VTORSION with 1-4'
            lfailure = .true.
            return
         endif
c        --- compute 1-4 vdw interactions
c        --- compute the 1-4 distance squared
         call twh_mimage(bxdflag,ibox,xcord(4),ycord(4),zcord(4)
     &        ,xcord(1),ycord(1),zcord(1),rxij,ryij,rzij)
         rijsq = rxij*rxij + ryij*ryij + rzij*rzij

c        --- need unit number of the opposite atom
         call twh_ijtor4(GLB_GET,imolty,unit_1,torind,unit_4)
         call twh_ntype(GLB_GET,imolty,unit_1,ntii)
         call twh_ntype(GLB_GET,imolty,unit_4,ntjj)

         vonefour = twh_vtwobody(lfailure,rijsq,ntii,ntjj,3,loverlap
     &        ,vf,dv )
         if ( torofcode .eq. CODE_RESCALE ) then
c           --- get the scaling factor and rescale
            call twh_itscale(GLB_GET,imolty,unit_1,torind,torscale)
            vonefour = torscale*vonefour
         endif
         if ( lreport ) then
            call twh_report(GLB_INCR,REP_TOR_VDW,vonefour)
         endif

c        --- add into the torsion energy
         twh_vtorsion = twh_vtorsion + vonefour
         if ( lfailure ) return
         if ( loverlap ) return

c        --- compute 1-4 coulombic interactions
         rij = dsqrt(rijsq)

c        --- note that we are only computing the real space portion
c        --- of the coulombic interaction here.  The correction for 
c        --- scaling factors due to the ewald sum is accounted for
c        --- in the eng** routines.  It is done this way for two reasons
c        --- 1) So we don't need any error function info in this step
c        --- 2) So I don't need any special logic to determine 1-4 
c        --- interactions in engatom.

         call twh_vcoulomb(0,bxdflag,ibox,rij,imolty,unit_1,imolty
     &        ,unit_4,vf,dv,loverlap,vonefour)
         vonefour = torcoeff(-1)*vonefour

         if ( torofcode .eq. CODE_RESCALE ) then
            vonefour = torscale*vonefour
         endif
         if ( lreport ) then
            call twh_report(GLB_INCR,REP_TOR_CHG,vonefour)
         endif
         twh_vtorsion = twh_vtorsion + vonefour
         if ( loverlap ) return
      endif

c      write(6,*) 'finished VTORSION'

      return
      end
