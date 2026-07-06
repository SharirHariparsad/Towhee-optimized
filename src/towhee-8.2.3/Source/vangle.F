#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
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
      function twh_vangle(lfailure,lreport,lovrlap,bxdflag,ibox
     &     ,imolty,iunit,ibend,iflag,angcos,distone,disttwo
     &     ,xcord,ycord,zcord)
c     ******************************************************************
c     * computes the bond bending energy given an angle type and the   *
c     * angle.                                                         *
c     * lreport: store all of the angle terms for a report             *
c     * ibox: the box the molecule is in (used for mimage)             *
c     * imolty: the molecule type                                      *
c     * iunit: the unit of the first atom in the bend                  *
c     * ibend: the index of the bend                                   *
c     * iflag: 0 if we should use the angcos,distone,disttwo info      *
c     *        1 if we should use the xcord,ycord,zcord info           *
c     *       -1 if we want to return the angle, instead of the        *
c     *             energy, from the xcord,ycord,zcord info            *
c     * xcord: xcoordinates of the three atoms                         *
c     * ycord: ycoordinates of the three atoms                         *
c     * zcord: zcoordinates of the three atoms                         *
c     * angcos: the cosine of the bond bending angle                   *
c     * distone: the distance between atoms 1-2                        *
c     * disttwo: the distance between atoms 2-3                        *
c     * lovrlap: returned as true if there is a hard overlap in the    *
c     *          nonbonded part of certain bending potentials          *
c     *          Note: we are a bit tricky as lovrlap is not intialized*
c     *          in this routine, should be intialized by the caller   *
c     *                                                                *
c     * rewritten from old sumup subroutines 01-21-2000 by M.G. Martin *
c     * last modified 03-30-2020 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_DOTPRODUCT
#define FUNCTION_LAWOFCOSINE
#define FUNCTION_VANGLE
#define FUNCTION_VTWOBODY
#include "functions.h"

#define DEBUG_VANGLE 0
c     --- variables passed to and from the function
      logical lfailure,lovrlap,lreport
      integer bxdflag,ibox,imolty,iunit,ibend,iflag
      double precision distone,disttwo,angcos
      double precision xcord,ycord,zcord
      dimension xcord(3),ycord(3),zcord(3)
c     --- local variables
c     --- logical scalars
      logical ldononbond,lbonang
c     --- integer scalars
      integer style,iutype,jutype,type,junit,iorder,bone,btwo,iallow
      integer bangtype,ibcoeff
c     --- double precision scalars
      double precision dist,distsq,vbondang,vvdw,vcharge
      double precision angle,differ,value,prefac,diffsquare,diffcube
      double precision xone,yone,zone,xtwo,ytwo,ztwo,dv
      double precision fac_one,fac_two,benprefact
      double precision len_bone,len_btwo
#if DEBUG_VANGLE
      double precision radtodeg
#endif
c     --- double precision parameters
      double precision tolerance
      parameter (tolerance = 1d-5)
c     --- double precision arrays
      double precision vf
      dimension vf(MAX_FOREIGN_LAMBDA)
      double precision bencoeff
      dimension bencoeff(MINBENDCOEFF:MAXBENDCOEFF)
      double precision bangcoeff
      dimension bangcoeff(2)
      
c     --- initialize
      twh_vangle = 0.0d0
      ldononbond = .false.
      if ( abs(iflag) .eq. 1 ) then
c        --- compute the bond lengths and angles from the coordinates
         call twh_mimage(bxdflag,ibox,xcord(1),ycord(1),zcord(1)
     &        ,xcord(2),ycord(2),zcord(2),xone,yone,zone )
         call twh_mimage(bxdflag,ibox,xcord(3),ycord(3),zcord(3)
     &        ,xcord(2),ycord(2),zcord(2),xtwo,ytwo,ztwo )
         angcos = twh_dotproduct(xone,yone,zone,xtwo,ytwo,ztwo,distone
     &        ,disttwo)
         if ( iflag .eq. -1 ) then
            twh_vangle = twh_arccos(angcos)
            return
         endif
      endif

#if DEBUG_VANGLE
      angle = twh_arccos(angcos)
      call twh_constant_radtodeg(GLB_GET,radtodeg)
      call twh_ijben2(GLB_GET,imolty,iunit,ibend,bone)
      call twh_ijben3(GLB_GET,imolty,iunit,ibend,btwo)
      write(6,*) 'imolty,units,angle,:',imolty,iunit,bone,btwo
     &     ,angle*radtodeg
#endif

c     --- determine type, style and order of the angle
      call twh_itben(GLB_GET,imolty,iunit,ibend,type)
      call twh_anglestyle(GLB_GET,type,style)
c     --- determine if there are bond-angle terms
      call twh_lbonang(GLB_GET,type,lbonang)
c     --- retrieve all bencoeff values for this type
      do ibcoeff = MINBENDCOEFF,MAXBENDCOEFF
         call twh_bencoeff(GLB_GET,type,ibcoeff,bencoeff(ibcoeff))
      enddo
      if ( style .eq. ANG_STYLE_RIGID ) then
c        --- rigid angle, make sure it is within tolerance, else overlap
c        --- local store bencoeffs
         angle = twh_arccos(angcos)
         value = (angle - bencoeff(0))**2
         if (twh_cmp_gt(value, bencoeff(1)) ) then 
c           --- out of tolerance, hard overlap
            lovrlap = .true.
#if DEBUG_VANGLE
            write(6,*) 'vangle angle:',angle*radtodeg
            call twh_constant_radtodeg(GLB_GET,radtodeg)
            write(6,*) 'coeff angle:',bencoeff(0)*radtodeg
#endif
         else
c           --- in tolerance, no energy
            twh_vangle = 0.0d0
         endif
      elseif ( style .eq. ANG_STYLE_HARMONIC ) then
c        --- standard harmonic bond bending potential
c        --- energy = k* (angle - equil)**2
c        --- local store bencoeffs
         angle = twh_arccos(angcos)
         twh_vangle = bencoeff(1)
     &        *(angle-bencoeff(0))**2
      elseif ( style .eq. ANG_STYLE_DREIDING ) then
c        --- Dreiding K[ 1 + Cos(angle)]
c        --- local store bencoeffs
         ibcoeff = 1
         twh_vangle = bencoeff(1)*(1.0d0 + angcos)
      elseif ( style .eq. ANG_STYLE_HARM_COS ) then
c        --- harmonic cosine potential
c        --- local store bencoeffs
c        --- energy = k* (cos(angle) - cos(equil))**2
         twh_vangle = bencoeff(1)
     &        *(angcos - bencoeff(-1))**2
      elseif ( style .eq. ANG_STYLE_COMPASS_AUTO ) then
c        --- Compass quartic angle with autodetection
         angle = twh_arccos(angcos)
         differ = angle - bencoeff(0)
         value = differ*differ
         twh_vangle = bencoeff(1)*(value)
     &        + bencoeff(2)*(value*differ)
     &        + bencoeff(3)*(value*value)
         if ( lbonang ) then
c           --- add in both of the bond-angle cross terms
            call twh_ibtben1(GLB_GET,imolty,iunit,ibend,bone)
            call twh_ibtben2(GLB_GET,imolty,iunit,ibend,btwo)
            call twh_orderben(GLB_GET,imolty,iunit,ibend,iorder)
            call twh_vibcoeff(GLB_GET,bone,0,len_bone)
            call twh_vibcoeff(GLB_GET,btwo,0,len_btwo)
            if ( iorder .eq. ORD_FORWARD ) then
c              --- force constants match up with bond order
               prefac = bencoeff(4)*(distone - len_bone)
     &              + bencoeff(5)*(disttwo - len_btwo)
            elseif ( iorder .eq. ORD_REVERSE ) then
c              --- force constants match with inverted bond order
               prefac = bencoeff(5)*(distone - len_bone)
     &              + bencoeff(4)*(disttwo - len_btwo)
            else
               write(6,*) 'ang type',type
               write(6,*) 'iorder',iorder
               write(6,*) 'VANGLE: bad type'
               lfailure = .true.
               return
            endif
            vbondang = prefac*differ
         endif
      elseif ( style .eq. ANG_STYLE_CHARMM_UB ) then
c         --- Charmm harmonic with urey-bradley
         angle = twh_arccos(angcos)
         twh_vangle = bencoeff(1)
     &        *((angle-bencoeff(0))**2)

c        --- compute the distance between the 1-3 atoms
         dist = twh_lawofcosine(distone,disttwo,angcos)
         twh_vangle = twh_vangle 
     &        + bencoeff(3)
     &        *((dist-bencoeff(2))**2)
      elseif ( style .eq. ANG_STYLE_NONBOND ) then
c        --- nonbonded interaction between the 1-3 atoms
         ldononbond = .true.
      elseif ( style .eq. ANG_STYLE_NB_HARM ) then
c        --- harmonic plus nonbond and coulombic between the 1-3 atoms
         angle = twh_arccos(angcos)
         twh_vangle = bencoeff(1)
     &        *(angle-bencoeff(0))**2
         ldononbond = .true.
      elseif ( style .eq. ANG_STYLE_COMPASS_EX ) then
c        --- Compass quartic angle with explicit ordering of terms
         angle = twh_arccos(angcos)
         differ = angle - bencoeff(0)
         value = differ*differ
         twh_vangle = bencoeff(1)*(value)
     &        + bencoeff(2)*(value*differ)
     &        + bencoeff(3)*(value*value)
         if ( lbonang ) then
c           --- add in both of the bond angle cross terms
            call twh_orderben(GLB_GET,imolty,iunit,ibend,iorder)
            if ( iorder .eq. ORD_FORWARD ) then
c              --- force constants match up with bond order
               prefac = bencoeff(4)
     &              *(distone - bencoeff(5))
     &              + bencoeff(6)
     &              *(disttwo - bencoeff(7))
            elseif ( iorder .eq. ORD_REVERSE ) then
c              --- force constants match with inverted bond order
               prefac = bencoeff(6)
     &              *(distone - bencoeff(7))
     &              + bencoeff(4)
     &              *(disttwo - bencoeff(5))
            else
               write(6,*) 'VANGLE: bad iorder value ',iorder
               lfailure = .true.
               return
            endif
            vbondang = prefac*differ
         endif
      elseif ( style .eq. ANG_STYLE_FOURIER ) then
c        --- Fourier expansion with constant minus single term
         angle = twh_arccos(angcos)
         call twh_benprefact(GLB_GET,imolty,iunit,ibend,benprefact)
         twh_vangle = benprefact*(1.0d0 + bencoeff(2)
     &        *dcos( bencoeff(1)*angle))
      elseif ( style .eq. ANG_STYLE_TTFOURIER ) then
c        --- 3 term Fourier expansion
         angle = twh_arccos(angcos)
         call twh_benprefact(GLB_GET,imolty,iunit,ibend,benprefact)
         twh_vangle = benprefact*(bencoeff(1)
     &        + bencoeff(2)*angcos
     &        + bencoeff(3)*dcos(2.0d0*angle) )
      elseif ( style .eq. ANG_STYLE_NONE ) then
c        --- no energy for this term
         twh_vangle = 0.0d0
      elseif ( style .eq. ANG_STYLE_SIXTH_AUTO ) then
c        --- MM2 sixth power angle with autodetection
         angle = twh_arccos(angcos)
         differ = angle - bencoeff(0)
         value = differ*differ
         twh_vangle = bencoeff(1)*value*
     &        (1.0d0 + bencoeff(2)*value*value)
         if ( lbonang ) then
c           --- add in both of the bond-angle cross terms
c           --- check based upon only the first letter for this term
            call twh_orderben(GLB_GET,imolty,iunit,ibend,iorder)
            if ( iorder .eq. ORD_FORWARD ) then
c              --- force constants match up with bond order
               call twh_ibtben1(GLB_GET,imolty,iunit,ibend,bone)
               call twh_ibtben2(GLB_GET,imolty,iunit,ibend,btwo)
            elseif ( iorder .eq. ORD_REVERSE ) then
c              --- force constants match with inverted bond order
               call twh_ibtben2(GLB_GET,imolty,iunit,ibend,bone)
               call twh_ibtben1(GLB_GET,imolty,iunit,ibend,btwo)
            else
               write(6,*) 'ang type',type
               write(6,*) 'orderben:',iorder
               write(6,*) 'VANGLE: bad type'
               lfailure = .true.
               return
            endif
            call twh_vibcoeff(GLB_GET,bone,0,len_bone)
            call twh_vibcoeff(GLB_GET,btwo,0,len_btwo)
            vbondang = bencoeff(4)*differ*(
     &           (distone - len_bone) + (disttwo - len_btwo) )
         endif
      elseif ( style .eq. ANG_STYLE_SQUARE_WELL ) then
c        --- infinite square well bending angle
c        --- compute the distance between the 1-3 atoms
         dist = twh_lawofcosine(distone,disttwo,angcos)
         if ( twh_cmp_lt(dist, bencoeff(1)) .or.
     &        twh_cmp_gt(dist, bencoeff(2)) ) then
            lovrlap = .true.
         else
            twh_vangle = bencoeff(3)
         endif
      elseif ( style .eq. ANG_STYLE_MULTI_RIGID ) then
c        --- rigid angle with two allowed angles
c        --- make sure it is within tolerance of one of the allowed
c        --- angles, else overlap
         lovrlap = .true.
         angle = twh_arccos(angcos)
         do iallow = 0,1
            value = (angle - bencoeff(iallow))**2
            if ( .not. twh_cmp_gt(value, bencoeff(2)) 
     &           ) then 
c              --- in tolerance, no energy
               lovrlap = .false.
               twh_vangle = 0.0d0
            endif
         enddo
      elseif ( style .eq. ANG_STYLE_MMFF ) then
c        --- MMFF cubic plus a bond-angle cross term
         angle = twh_arccos(angcos)
         differ = angle - bencoeff(0)
         diffsquare = differ*differ
         diffcube = differ*diffsquare
         twh_vangle = bencoeff(1)*diffsquare
     &        + bencoeff(2)*diffcube
c        --- cross terms
         if ( lbonang ) then
c           --- add in both of the bond-angle cross terms
            call twh_ibtben1(GLB_GET,imolty,iunit,ibend,bone)
            call twh_ibtben2(GLB_GET,imolty,iunit,ibend,btwo)
            call twh_orderben(GLB_GET,imolty,iunit,ibend,iorder)
            call twh_vibcoeff(GLB_GET,bone,0,len_bone)
            call twh_vibcoeff(GLB_GET,btwo,0,len_btwo)
            if ( iorder .eq. ORD_FORWARD ) then
c              --- force constants match up with bond order
               fac_one = bencoeff(3)*(distone - len_bone)
               fac_two = bencoeff(4)*(disttwo - len_btwo)
            elseif ( iorder .eq. ORD_REVERSE ) then
c              --- force constants match with inverted bond order
               fac_one = bencoeff(4)*(distone - len_bone)
               fac_two = bencoeff(3)*(disttwo - len_btwo)
            else
               write(6,*) 'ang type',type
               write(6,*) 'benorder:',iorder
               write(6,*) 'VANGLE: bad type'
               lfailure = .true.
               return
            endif
            vbondang = (fac_one + fac_two)*differ
#if DEBUG_VANGLE
            write(6,'(4(a,1x,f8.3,1x))')
     &           'angle:',angle*radtodeg
     &           ,'differ:',differ*radtodeg
     &           ,'delta one:',distone-len_bone
     &           ,'ba one:',fac_one*differ/CON_KCALTOK
            write(6,'(4(a,1x,f8.3,1x))')
     &           'angle:',angle*radtodeg
     &           ,'differ:',differ*radtodeg
     &           ,'delta two:',disttwo-len_btwo
     &           ,'ba two:',fac_two*differ/CON_KCALTOK
#endif            
         endif
      elseif ( style .eq. ANG_STYLE_HARM_COS_OT ) then
c        --- harmonic cosine plus nonbond between the 1-3 atoms
         twh_vangle = bencoeff(1)
     &        *(angcos - bencoeff(-1))**2
         ldononbond = .true.

      elseif ( style .eq. ANG_STYLE_MMFFC_ANG ) then
c        --- MMFF cubic plus a bond-angle cross term
         angle = twh_arccos(angcos)
         differ = angle - bencoeff(0)
         diffsquare = differ*differ
         diffcube = differ*diffsquare
         twh_vangle = bencoeff(1)*diffsquare
     &        + bencoeff(2)*diffcube
c        --- cross terms
         call twh_itbang(GLB_GET,imolty,iunit,ibend,bangtype)
         if ( bangtype .ne. 0 ) then
c           --- special case logic for bond-ang
            lbonang = .true.
c           --- add in both of the bond-angle cross terms
            call twh_ibtben1(GLB_GET,imolty,iunit,ibend,bone)
            call twh_ibtben2(GLB_GET,imolty,iunit,ibend,btwo)
            call twh_orderbang(GLB_GET,imolty,iunit,ibend,iorder)
            call twh_vibcoeff(GLB_GET,bone,0,len_bone)
            call twh_vibcoeff(GLB_GET,btwo,0,len_btwo)
            call twh_bencoeff(GLB_GET,bangtype,1,bangcoeff(1))
            call twh_bencoeff(GLB_GET,bangtype,2,bangcoeff(2))
            if ( iorder .eq. ORD_FORWARD ) then
c              --- force constants match up with bond order
               fac_one = bangcoeff(1)*(distone-len_bone)
               fac_two = bangcoeff(2)*(disttwo-len_btwo)
            elseif ( iorder .eq. ORD_REVERSE ) then
c              --- force constants match with inverted bond order
               fac_one = bangcoeff(2)*(distone-len_bone)
               fac_two = bangcoeff(1)*(disttwo-len_btwo)
            else
               write(6,*) 'iutype',iutype,'bang type',bangtype
               write(6,*) 'orderbang:',iorder
               write(6,*) 'VANGLE: bad type'
               lfailure = .true.
               return
            endif
            vbondang = (fac_one + fac_two)*differ
#if DEBUG_VANGLE
            write(6,'(4(a,1x,f8.3,1x))')
     &           'angle:',angle*radtodeg
     &           ,'differ:',differ*radtodeg
     &           ,'delta one:',distone-len_bone
     &           ,'ba one:',fac_one*differ/CON_KCALTOK
            write(6,'(4(a,1x,f8.3,1x))')
     &           'angle:',angle*radtodeg
     &           ,'differ:',differ*radtodeg
     &           ,'delta two:',disttwo-len_btwo
     &           ,'ba two:',fac_two*differ/CON_KCALTOK
#endif            
         else
            lbonang = .false.
         endif


      else
         write(6,*) 'VANGLE: non defined bond angle type'
         lfailure = .true.
         return
      endif

      if ( lreport ) then
         call twh_report(GLB_INCR,REP_ANG_REG,twh_vangle)
      endif

      if ( lbonang ) then
         if ( lreport ) then
            call twh_report(GLB_INCR,REP_ANG_BA,vbondang)
         endif
         twh_vangle = twh_vangle + vbondang
      endif

      if ( ldononbond ) then
c        --- compute the distance between the 1-3 atoms
         dist = twh_lawofcosine(distone,disttwo,angcos)
         distsq = dist*dist
c        --- get the nonbonded vdw energy 
         call twh_ijben3(GLB_GET,imolty,iunit,ibend,junit)
         call twh_ntype(GLB_GET,imolty,iunit,iutype)
         call twh_ntype(GLB_GET,imolty,junit,jutype)
         vvdw = twh_vtwobody(lfailure,distsq
     &        ,iutype,jutype,0,lovrlap,vf,dv)
         if ( lovrlap .or. lfailure ) return

         call twh_vcoulomb( 0,bxdflag,ibox,dist,imolty,iunit
     &     ,imolty,junit,vf, dv,lovrlap,vcharge)
         if ( lovrlap ) return

         if ( lreport ) then
            call twh_report(GLB_INCR,REP_ANG_VDW,vvdw)
            call twh_report(GLB_INCR,REP_ANG_CHG,vcharge)
         endif
         twh_vangle = twh_vangle + vvdw + vcharge
      endif

      return
      end
