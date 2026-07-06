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
      subroutine twh_writelammps(lfailure)
c     ******************************************************************
c     * creates an input file for Lammps (named lammps_data) which is  *
c     * suitable for reading into the LAMMPS massively parallel MD code*
c     * updated to work with the current C++ version of LAMMPS and not *
c     * the older versions                                             *
c     * also creates part of the lammps_input file that has to do with *
c     * the force field styles.  This is written to lammps_input.out   *
c     *                                                                *
c     * originally written 08-07-2001 by M.G. Martin                   *
c     * last modified 03-31-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_OLDGETATOMNUM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local parameters
      integer maxwarn
      parameter (maxwarn=30)
c     --- local variables
c     --- logical scalars
      logical lnew,lfound,lclass2angle,lclass2dihed,ldone,ladd,lshift
c     --- logical arrays
      logical lassign
      dimension lassign(TAAMAX)
      logical lwarnimp
      dimension lwarnimp(maxwarn)
      logical lwarnvibstyle
      dimension lwarnvibstyle(maxwarn)
      logical lwarnangstyle
      dimension lwarnangstyle(maxwarn)
      logical lwarndihedstyle
      dimension lwarndihedstyle(maxwarn)
c     --- character string scalars
      character*2 numbers
      character*30 filename
c     --- character string arrays
      character*(FFNAMELEN) ctname,torstring
      dimension ctname(4),torstring(4)
c     --- integer scalars
      integer iwarn,numboxes,cpnumber
      integer lam_atom,lam_bond,lam_angle,lam_dihedral,lam_improp
      integer lam_atom_types,lam_bond_types,lam_angle_types
      integer lam_angang,lam_angang_types,nunit,idim
      integer lam_dihedral_types,lam_improp_types,itemp
      integer index,nchain,ntorloop
      integer filenum,ibox,ivib,ibend,itor,imp,type,ichain,imolty,iunit
      integer itricky,count_atom,count_molecule,count_bond,count_types
      integer inpnum,offset,count_angle,iwinner,count_dihedral
      integer itotal,count_improp,icoeff,iaa,jaa,kaa,atoma,atomc,atomd
      integer test0,test1,test2,improperform,ntype
      integer unit_1,unit_2,unit_3,unit_4,ndiv_current,nmolty
      integer itest,iorder,jorder,jvib,nvib
      integer ijaa0,ijaa1,ijaa2,ijben3test,ijben2test
      integer ijimp2,ijimp3,ijimp4,ijtor2test,ijtor3test,ijtor4test
      integer ijvibtest,inaa,itaatemp,inbend,inimproper,intorsion
      integer invibration,invib_two,invib_thr,ittortest,testboxi
      integer torstyle,impstyle,vibstyle,anglestyle
c     --- integer arrays
      integer code_atom
      dimension code_atom(NNTYPE)
      integer code_bond
      dimension code_bond(TVIBMAX)
      integer code_angle
      dimension code_angle(TBENMAX)
      integer code_improp
      dimension code_improp(-TTORMAX:TIMPMAX)
      integer code_angang
#define WLAAMAX 50
      dimension code_angang(-WLAAMAX:WLAAMAX,-WLAAMAX:WLAAMAX
     &     ,-WLAAMAX:WLAAMAX)
      integer decode_atom
      dimension decode_atom(NNTYPE)
      integer decode_bond
      dimension decode_bond(TVIBMAX)
      integer decode_angle
      dimension decode_angle(TBENMAX)
      integer decode_dihedral
      dimension decode_dihedral(3*TTORMAX,2)
      integer decode_improp,itype
      dimension decode_improp(TTORMAX,0:2)
      integer triple
      dimension triple(3)
      integer div_dihedral
      dimension div_dihedral(3*TTORMAX)
c     --- double precision scalars
      double precision epsilon,sigma,epsilonof,sigmaof,ewald_prec
      double precision rcelect,hmsum,dvalue,rcut,aacoeff,mass
      double precision tolerance,scaletest,dtemp,temperature
      parameter( tolerance = 1d-5)
      double precision xtemp,ytemp,ztemp,radtodeg,qqatom
      double precision impcoeff_zero,impcoeff_one
c     --- double precision arrays
      double precision scale_dihedral
      dimension scale_dihedral(3*TTORMAX)
      double precision tmp
      dimension tmp(6)
      double precision vibcoeff
      dimension vibcoeff(MINVIBCOEFF:MAXVIBCOEFF)
      double precision bencoeff
      dimension bencoeff(MINBENDCOEFF:MAXBENDCOEFF)
      double precision torcoeff
      dimension torcoeff(MINTORCOEFF:MAXTORCOEFF)
      
c      write(6,*) 'starting writelammps'

c     --- retrieve constants
      call twh_constant_radtodeg(GLB_GET,radtodeg)
      call twh_numboxes(GLB_GET,numboxes)
      call twh_rcut(GLB_GET,rcut)

      do index = 1,3
         triple(index) = 0
      enddo

c     --- initialize lwarn
      do iwarn = 1,maxwarn
         lwarnimp(iwarn) = .false.
         lwarnvibstyle(iwarn) = .false.
         lwarnangstyle(iwarn) = .false.
         lwarndihedstyle(iwarn) = .false.
      enddo
c     --- get the classical potential
      call twh_cp_number(GLB_GET,cpnumber)

c     --- open lammps_input
      inpnum = 91
      filename = 'lammps_input.out'
      call twh_openfile(lfailure,inpnum,filename,'DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfailure ) return
      write(inpnum,*) 'All of the force field styles are outputted'
      write(inpnum,*) 'for every type.  You should check to make sure'
      write(inpnum,*) 'that they are all the same style, otherwise'
      write(inpnum,*) 'LAMMPS will have problems.'
      write(inpnum,*)

c     --- output the temperature
      call twh_temperature(GLB_GET,temperature)
      write(inpnum,'(a20,f10.4,a7)') 'create temp uniform ',temperature
     & ,' 123456'

      call twh_nchain(GLB_GET,nchain)
      do ibox = 1,numboxes
         filename = 'lammps_data'
         write(numbers,'(i2.2)') ibox
         filename(12:13) = numbers
         filenum = 90
         call twh_openfile(lfailure,filenum,filename,'DEFAULT','UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) return
c        --- record box number in input file
         write(inpnum,*) 'Partial input for simulation box ',ibox
         write(inpnum,*)

c        --- initialize atomcode
         do type = 1, NNTYPE
            code_atom(type) = 0
         enddo
         do type = 1,TVIBMAX
            code_bond(type) = 0
         enddo
         do type = 1,TBENMAX
            code_angle(type) = 0
         enddo
         do type = -TTORMAX,TIMPMAX
            code_improp(type) = 0
         enddo
c        --- initialize lam_*
         lam_atom = 0
         lam_bond = 0
         lam_angle = 0
         lam_angang = 0
         lam_dihedral = 0
         lam_improp = 0
c        --- initialize lam_*_types
         lam_atom_types = 0
         lam_bond_types = 0
         lam_angle_types = 0
         lam_angang_types = 0
         lam_dihedral_types = 0
         lam_improp_types = 0
c        --- initialize code_angang
         do iaa = -WLAAMAX,WLAAMAX
            do jaa = -WLAAMAX,WLAAMAX
               do kaa = -WLAAMAX,WLAAMAX
                  code_angang(iaa,jaa,kaa) = 0
               enddo
            enddo
         enddo
         
c        --- tally up atom,bonds,angles,dihedrals,impropers and types
         do ichain = 1,nchain
            call twh_nboxi(GLB_GET,ichain,testboxi)
            if ( testboxi .eq. ibox ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               do iunit = 1,nunit
c                 --- increment lam_atom
                  lam_atom = lam_atom + 1
                  call twh_ntype(GLB_GET,imolty,iunit,type)
                  if ( code_atom(type) .eq. 0 ) then
c                   --- we do not yet have an atom of this type, need 
c                   --- to add entry into code_atom
                     lam_atom_types = lam_atom_types + 1
                     code_atom(type) = lam_atom_types
                     decode_atom(lam_atom_types) = type
                  endif
                  call twh_invib(GLB_GET,imolty,iunit,invibration)
                  do ivib = 1,invibration
                     call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
                     if ( ijvibtest .gt. iunit ) then
                        lam_bond = lam_bond + 1
                        call twh_itvib(GLB_GET,imolty,iunit,ivib,type)
                        if ( code_bond(type) .eq. 0 ) then
c                          --- do not yet have a bond of this type
c                          --- add entry into code_bond
                           lam_bond_types = lam_bond_types + 1
                           code_bond(type) = lam_bond_types
                           decode_bond(lam_bond_types) = type
                        endif
                     endif
                  enddo
                  call twh_inben(GLB_GET,imolty,iunit,inbend)
                  do ibend = 1,inbend
                     call twh_itben(GLB_GET,imolty,iunit,ibend,type)
                     ladd = .false.
                     call twh_anglestyle(GLB_GET,type,anglestyle)
                     if ( anglestyle .eq. ANG_STYLE_COMPASS_EX ) then
c                       --- add this to the list in all cases as 
c                       --- the positive type number already selected
c                       --- this angle a single time
                        ladd = .true.
                     elseif ( anglestyle
     &                       .eq. ANG_STYLE_COMPASS_AUTO ) then
c                       --- only add in the case the constants are 
c                       --- in the correct order, or when the names
c                       --- are identical add based on iunit 
                        call twh_orderben(GLB_GET,imolty,iunit,ibend
     &                       ,iorder)
                        if ( iorder .eq. ORD_FORWARD ) then
c                          --- possible they are both forward, if that
c                          --- is the case then we resort to the unit
c                          --- comparison again
                           call twh_ijben3(GLB_GET,imolty,iunit,ibend
     &                          ,itest)
                           jvib = 0
                           call twh_inben(GLB_GET,imolty,itest,nvib)
                           lfound = .false.
                           do while ( jvib .lt. nvib
     &                          .and. .not. lfound )
                              jvib = jvib + 1
                              call twh_ijben3(GLB_GET,imolty,itest,jvib
     &                             ,ijben3test)
                              if ( ijben3test .eq. iunit ) then
c                                 --- grab the inverse order
                                 call twh_orderben(GLB_GET,imolty,itest
     &                                ,jvib,jorder)
                                 lfound = .true.
                              endif
                           enddo
                           if ( .not. lfound ) then
                              write(6,*) 'WRITELAMMPS: problem reverse'
     &                             ,' in the bending terms'
                              return
                           endif
                           if ( iorder .eq. jorder ) then
c                             --- only add if ijben3 > iunit
                              call twh_ijben3(GLB_GET,imolty,iunit
     &                             ,ibend,ijben3test)
                              if ( ijben3test .gt. iunit ) then
                                 ladd = .true.
c                                --- stash the vibcoeffs in bencoeffs
                                 call twh_ibtben1(GLB_GET,imolty,iunit
     &                                ,ibend,itemp)
                                 call twh_vibcoeff(GLB_GET,itemp,0
     &                                ,dtemp)
                                 call twh_bencoeff(GLB_SET,type,9,dtemp)
                                 call twh_ibtben2(GLB_GET,imolty,iunit
     &                                ,ibend,itemp)
                                 call twh_vibcoeff(GLB_GET,itemp,0
     &                                ,dtemp)
                                 call twh_bencoeff(GLB_SET,type,10
     &                                ,dtemp)
                              endif
                           else
c                             --- other string is reverse so always add
                              ladd = .true.
                           endif
                        endif
                     else
c                       --- add based solely on ijben3 > iunit
                        call twh_ijben3(GLB_GET,imolty,iunit,ibend
     &                       ,ijben3test)
                        if ( ijben3test .gt. iunit ) ladd = .true.
                     endif
                     if ( ladd ) then
                        lam_angle = lam_angle + 1
                        if ( code_angle(type) .eq. 0 ) then
c                          --- do not yet have an angle of this type
c                          --- add entry into code_angle
                           lam_angle_types = lam_angle_types + 1
                           code_angle(type) = lam_angle_types
                           decode_angle(lam_angle_types) = type
                        endif
                     endif
                  enddo
                  call twh_inaa(GLB_GET,imolty,iunit,inaa)
                  if ( inaa .ne. 0 ) then
c                    --- unfortunately the angle-angle terms are stored 
c                    --- quite a bit differently in LAMMPS than Towhee
c                    --- LAMMPS combines 3 Towhee angle-angle terms into
c                    --- a single term.  Need to search through and find
c                    --- appropriate triples
c                    --- reinitialize tripple array
                     do iaa = 1,inaa
                        lassign(iaa) = .false.
                     enddo
                     ldone = .false.
                     do while (.not. ldone ) 
c                       --- find the first unassigned angang
                        triple(1) = 0
                        iaa = 0
                        do while ( triple(1) .eq. 0 ) 
                           iaa = iaa + 1
                           if ( .not. lassign(iaa) ) then
c                             --- take this as the first in the triple
                              call twh_itaa(GLB_GET,imolty,iunit,iaa
     &                             ,triple(1))
                              lassign(iaa) = .true.
                           endif
                        enddo
c                       --- extract the first triple
                        call twh_ijaa0(GLB_GET,imolty,iunit,iaa,atomc)
                        call twh_ijaa1(GLB_GET,imolty,iunit,iaa,atoma)
                        call twh_ijaa2(GLB_GET,imolty,iunit,iaa,atomd)
                        do jaa = 1,inaa
                           if ( .not. lassign(jaa) ) then
c                             --- see if some permutation of these atoms
c                             --- matches the first atoms
                              call twh_ijaa0(GLB_GET,imolty,iunit,jaa
     &                             ,test0)
                              call twh_ijaa1(GLB_GET,imolty,iunit,jaa
     &                             ,test1)
                              call twh_ijaa2(GLB_GET,imolty,iunit,jaa
     &                             ,test2)
                              if ( atoma .eq. test0 
     &                             .and. atomc .eq. test1
     &                             .and. atomd .eq. test2 ) then
c                                 --- this fills the second slot in the 
c                                 --- lammps angle-angle
                                 call twh_itaa(GLB_GET,imolty,iunit,jaa
     &                                ,triple(2))
                                 lassign(jaa) = .true.
                              elseif (atoma .eq. test0 
     &                                .and. atomc .eq. test2
     &                                .and. atomd .eq. test1 ) then
c                                 --- this fills the second slot in the 
c                                 --- lammps angle-angle
                                 call twh_itaa(GLB_GET,imolty,iunit,jaa
     &                                ,itaatemp)
                                 triple(2) = -itaatemp
                                 lassign(jaa) = .true.
                              elseif ( atomd .eq. test0 
     &                             .and. atoma .eq. test1
     &                             .and. atomc .eq. test2 ) then
c                                 --- this fills the third slot in the 
c                                 --- lammps angle-angle
                                 call twh_itaa(GLB_GET,imolty,iunit,jaa
     &                                ,triple(3))
                                 lassign(jaa) = .true.
                              elseif ( atomd .eq. test0 
     &                             .and. atoma .eq. test2
     &                             .and. atomc .eq. test1 ) then
                                 call twh_itaa(GLB_GET,imolty,iunit,jaa
     &                                ,itaatemp)
                                 triple(3) = -itaatemp
                                 lassign(jaa) = .true.
                              endif
                           endif
                        enddo
c                       --- safety check
                        if ( triple(1)*triple(2)*triple(3) .eq. 0 ) then
c                          --- something went wrong
                           write(6,*) 'WRITELAMMPS: could not map '
     &                          ,'angle-angle terms'
                           goto 666
                        endif
c                       --- increase lam_angang
                        lam_angang = lam_angang + 1
c                       --- see if this type has been added to the list
                        if ( code_angang(triple(1),triple(2),triple(3)) 
     &                       .eq. 0 ) then
c                          --- this is a new type of lammps angang
                           lam_angang_types = lam_angang_types + 1
                           code_angang(triple(1),triple(2),triple(3)) 
     &                          = lam_angang_types
                        endif
c                       --- determine if we are done
                        ldone = .true.
                        do iaa = 1,inaa
                           if ( .not. lassign(iaa) ) ldone = .false.
                        enddo
                     enddo
                  endif
                  call twh_intor(GLB_GET,imolty,iunit,intorsion)
                  do itor = 1,intorsion
                     call twh_ittor(GLB_GET,imolty,iunit,itor,type)
                     unit_1 = iunit
                     call twh_ijtor2(GLB_GET,imolty,iunit,itor,unit_2)
                     call twh_ijtor3(GLB_GET,imolty,iunit,itor,unit_3)
                     call twh_ijtor4(GLB_GET,imolty,iunit,itor,unit_4)
                     if ( type .gt. 0 ) then
                        ladd = .false.
                        call twh_torstyle(GLB_GET,type,torstyle)
                        if ( torstyle .eq. TOR_STYLE_COMP_X_EX ) then
                           ladd = .true.
                        elseif ( torstyle
     &                          .eq. TOR_STYLE_COMP_X_AD ) then
                           call twh_ntype(GLB_GET,imolty,unit_1,ntype)
                           call twh_ctname(GLB_GET,ntype,ctname(1))
                           call twh_ntype(GLB_GET,imolty,unit_2,ntype)
                           call twh_ctname(GLB_GET,ntype,ctname(2))
                           call twh_ntype(GLB_GET,imolty,unit_3,ntype)
                           call twh_ctname(GLB_GET,ntype,ctname(3))
                           call twh_ntype(GLB_GET,imolty,unit_4,ntype)
                           call twh_ctname(GLB_GET,ntype,ctname(4))
                           if ( ( ctname(1) .eq. ctname(4) )
     &                          .and. ( ctname(2) .eq. ctname(3) )
     &                          ) then
c                             --- reversable, only store if ij4 > iunit
                              call twh_ijtor4(GLB_GET,imolty,iunit,itor
     &                             ,ijtor4test)
                              if ( ijtor4test .gt. iunit ) ladd = .true.
                           else
c                             --- add if the atoms match in order
                              do itest = 1,4
                                 call twh_torstring(GLB_GET,type
     &                                ,itest,torstring(itest))
                              enddo
                              if ( ctname(1) .eq. torstring(1) .and.
     &                             ctname(2) .eq. torstring(2) .and.
     &                             ctname(3) .eq. torstring(3) .and.
     &                             ctname(4) .eq. torstring(4) ) then
                                 ladd = .true.
                              endif
                           endif
                           if ( ladd ) then
c                             --- stash all of the parameters
                              call twh_ibttor1(GLB_GET,imolty,unit_1
     &                             ,itor,itemp)
                              call twh_vibcoeff(GLB_GET,itemp,0,dtemp)
                              call twh_torcoeff(GLB_SET,type,20,dtemp)
                              call twh_ibttor2(GLB_GET,imolty,unit_1
     &                             ,itor,itemp)
                              call twh_vibcoeff(GLB_GET,itemp,0,dtemp)
                              call twh_torcoeff(GLB_SET,type,21,dtemp)
                              call twh_ibttor3(GLB_GET,imolty,unit_1
     &                             ,itor,itemp)
                              call twh_vibcoeff(GLB_GET,itemp,0,dtemp)
                              call twh_torcoeff(GLB_SET,type,22,dtemp)
                              call twh_iattor1(GLB_GET,imolty,unit_1
     &                             ,itor,itemp)
                              call twh_bencoeff(GLB_GET,itemp,0,dtemp)
                              call twh_torcoeff(GLB_SET,type,23,dtemp)
                              call twh_iattor2(GLB_GET,imolty,unit_1
     &                             ,itor,itemp)
                              call twh_bencoeff(GLB_GET,itemp,0,dtemp)
                              call twh_torcoeff(GLB_SET,type,24,dtemp)
                           endif
                        else
c                          --- any other case just compar ijt4 > iunit
                           call twh_ijtor4(GLB_GET,imolty,iunit,itor
     &                          ,ijtor4test)
                           if ( ijtor4test .gt. iunit ) ladd = .true.
                        endif
                        if ( torstyle .eq. TOR_STYLE_AMBER ) then
c                          --- need the number of torsions across the
c                          --- two central atoms as this changes the 
c                          --- scaling of the energies
                           call twh_invib(GLB_GET,imolty,unit_2
     &                          ,invib_two)
                           call twh_invib(GLB_GET,imolty,unit_3
     &                          ,invib_thr)
                           ndiv_current = (invib_two -1)*(invib_thr -1)
                        else
c                          --- not used, just set to 0
                           ndiv_current = 0
                        endif
                        if ( ladd ) then
c                          --- converting from Towhee torsions to 
c                          --- LAMMPS torsions is tricky if ntorloop
c                          --- is something other than 0 or 1.
                           call twh_ntorloop(GLB_GET,type,ntorloop)
                           if ( ntorloop .eq. 0 ) then
                              lam_dihedral = lam_dihedral + 1
                           else
                              lam_dihedral = lam_dihedral + ntorloop
                           endif
                           lnew = .true.
                           do count_types = 1,lam_dihedral_types
c                             --- see if we have a complete match
                              if ( type .eq. 
     &                             decode_dihedral(count_types,1)) then
c                                --- check scaling factor
                                 if ( ntorloop .eq. 0 ) then
                                    call twh_itscale(GLB_GET,imolty
     &                                   ,iunit,itor,scaletest)
                                 else
                                    call twh_itscale(GLB_GET,imolty
     &                                   ,iunit,itor,scaletest)
                                    scaletest = scaletest/dble(ntorloop)
                                 endif
                                 if ( (scaletest 
     &                                - scale_dihedral(count_types))**2
     &                                .lt. tolerance ) then
c                                   --- check div
                                    if ( ndiv_current .eq. 
     &                                   div_dihedral(count_types) )then
                                       lnew = .false.
                                    endif
                                 endif
                              endif
                           enddo

                           if ( lnew ) then
c                             ---do not yet have a dihedral of this type
c                             --- add entry into decode_dihedral
                              if ( ntorloop .eq. 0 ) then
c                                --- easy case
                                 lam_dihedral_types 
     &                                = lam_dihedral_types+1
                                 decode_dihedral(lam_dihedral_types,1)
     &                                = type
                                 decode_dihedral(lam_dihedral_types,2)
     &                                = 0
                                 call twh_itscale(GLB_GET,imolty,iunit
     &                                ,itor,scale_dihedral(
     &                                lam_dihedral_types))
                                 div_dihedral(lam_dihedral_types) =
     &                                ndiv_current
                              else
c                                --- need to record this as multiple 
c                                --- torsions and torsion types
                                 do itricky = 1,ntorloop
                                    lam_dihedral_types 
     &                                   = lam_dihedral_types + 1
                                    decode_dihedral(lam_dihedral_types
     &                                   ,1) = type
                                    decode_dihedral(lam_dihedral_types
     &                                   ,2) = itricky
                                    call twh_itscale(GLB_GET,imolty
     &                                   ,iunit,itor,dvalue)
                                    scale_dihedral(lam_dihedral_types) =
     &                                   dvalue/dble(ntorloop)
                                    div_dihedral(lam_dihedral_types) =
     &                                   ndiv_current
                                 enddo
                              endif
                           endif
c                          --- see if we also need to split this into an
c                          --- improper
                           call twh_torstyle(GLB_GET,type,torstyle)
                           if ( torstyle .eq. TOR_STYLE_GCA_CUB ) then
c                             --- the last two terms of this are treated
c                             --- as a LAMMPS improper
                              lam_improp = lam_improp + 1
                              if ( code_improp(-type) .eq. 0 ) then
c                                --- do not yet have this entry
c                                --- add entry into code_improp
                                 lam_improp_types = lam_improp_types + 1
                                 code_improp(-type) = lam_improp_types
                                 decode_improp(lam_improp_types,1) 
     &                                = type
                                 decode_improp(lam_improp_types,2) = 0
                              endif
                           endif
                        endif
                     endif
                  enddo
c                 --- check the true impropers
                  call twh_inimprop(GLB_GET,imolty,iunit,inimproper)
                  do imp = 1,inimproper
                     call twh_itimprop(GLB_GET,imolty,iunit,imp,type)
                     call twh_impform(GLB_GET,type,improperform)
                     if ( improperform .ne. 2 
     &                    .and. improperform .ne. 3 ) then
                        if ( .not. lwarnimp(improperform)
     &                       ) then
                           write(6,*) 'WRITELAMMPS: cannot output'
     &                          ,' impform ',improperform,' to LAMMPS'
                           lwarnimp(improperform) = .true.
                        endif
                     else
                        lam_improp = lam_improp + 1
                        if ( code_improp(type) .eq. 0 ) then
c                          --- do not yet have an improper of this type
c                          --- add entry into code_improp
                           lam_improp_types = lam_improp_types + 1
                           code_improp(type) = lam_improp_types
                           decode_improp(lam_improp_types,1) = type
                           decode_improp(lam_improp_types,2) = 1
                        endif
                     endif
                  enddo
               enddo
            endif
         enddo

c        --- output information to data_lammps
         write(filenum,*) 'LAMMPS data file generated by the Towhee'
     &        ,' Monte Carlo code'
         write(filenum,*)
         write(filenum,*) lam_atom,' atoms'
         write(filenum,*) lam_bond,' bonds'
         write(filenum,*) lam_angle,' angles'
         write(filenum,*) lam_dihedral,' dihedrals'
c        --- combine the impropers and the angle-angle
         write(filenum,*) lam_improp+lam_angang,' impropers'
         write(filenum,*)
         write(filenum,*) lam_atom_types,' atom types'
         if ( lam_bond_types .ne. 0 ) 
     &        write(filenum,*) lam_bond_types,' bond types'
         if ( lam_angle_types .ne. 0 ) 
     &        write(filenum,*) lam_angle_types,' angle types'
         if ( lam_dihedral_types .ne. 0 )
     &        write(filenum,*) lam_dihedral_types,' dihedral types'
c        --- combine the improp with the angang types
         if ( lam_improp_types + lam_angang_types .ne. 0 ) 
     &        write(filenum,*) lam_improp_types+lam_angang_types
     &        ,' improper types'
         write(filenum,*)
c        --- box dimensions
c        --- lammps has to use rectangular parallelipipid boxes so 
c        --- just add all of the components together and mush this into
c        --- that kind of box
         hmsum = 0.0d0
         do idim = 1,3
            call twh_hmatrix(GLB_GET,BXD_REAL,ibox,idim,1,dvalue)
            hmsum = hmsum + dabs(dvalue)
         enddo
         write(filenum,*) '0.0 ',hmsum,' xlo xhi'
         hmsum = 0.0d0
         do idim = 1,3
            call twh_hmatrix(GLB_GET,BXD_REAL,ibox,idim,2,dvalue)
            hmsum = hmsum + dabs(dvalue)
         enddo
         write(filenum,*) '0.0 ',hmsum,' ylo yhi'
         hmsum = 0.0d0
         do idim = 1,3
            call twh_hmatrix(GLB_GET,BXD_REAL,ibox,idim,3,dvalue)
            hmsum = hmsum + dabs(dvalue)
         enddo
         write(filenum,*) '0.0 ',hmsum,' zlo zhi'
         write(filenum,*)

         if ( lam_atom .ne. 0 ) then
c           --- output masses
            write(filenum,3)
 3          format('Masses')
            write(filenum,*)
            do count_types = 1,lam_atom_types
               type = decode_atom(count_types)
               call twh_mass(GLB_GET,type,mass)
               write(filenum,*) count_types,mass
            enddo
            write(filenum,*)
         endif

         if ( lam_atom_types .ne. 0 ) then
c           --- output pair coeffs
            if ( cpnumber .ne. CP_EMBEDDED_ATOM_METHOD ) then
               write(filenum,'(a11)') 'Pair Coeffs'
               write(filenum,*)
            endif
            do count_types = 1,lam_atom_types
               type = decode_atom(count_types)
               if ( cpnumber .eq. CP_LENNARD_JONES .or. 
     &              cpnumber .eq. CP_12_6_PLUS_SOLVATION .or.
     &              cpnumber .eq. CP_12_6_PLUS_12_10_H_BOND .or.
     &              cpnumber .eq. CP_UFF_12_6
     &              ) then
c                 --- need to output potential style to input
                  write(inpnum,'(a30,1x,2(f7.2,1x))') 
     &                 'pair_style lj/charmm/coul/long',rcut,rcut+0.1d0
                  write(inpnum,'(a24)') 'kspace_style pppm 1.0e-4'
c                 --- Lennard-Jones 12-6
                  type = twh_getnbtype(type,type)
                  call twh_nbcoeff(GLB_GET,2,type,epsilon)
                  call twh_nbcoeff(GLB_GET,1,type,sigma)
                  call twh_nbcoeff(GLB_GET,4,type,epsilonof)
                  call twh_nbcoeff(GLB_GET,3,type,sigmaof)
                  write(filenum,*) count_types,epsilon/CON_KCALTOK
     &                 ,sigma,epsilonof/CON_KCALTOK,sigmaof
                  call twh_nb_to_hbtype(GLB_GET,type,itest)
                  if ( itest .ne. 0 ) then
                     write(6,*) 'WRITELAMMPS: cannot output the 12-10'
     &                    ,' portion for lammps'
                  endif
               elseif ( cpnumber .eq. CP_9_6 ) then
c                 --- still old style as 9-6 not (yet) in new LAMMPS
c                 --- LJ 9-6
                  call twh_lshift(GLB_GET,lshift)
                  if ( lshift ) then
                     offset = 1
                  else
                     offset = 0
                  endif
c                 --- need to output potential style to input
                  write(inpnum,'(a36,1x,f7.2,1x,i1)') 
     &                 'nonbond style          class2/cutoff'
     &                 ,rcut,offset
                  call twh_ewald_prec(GLB_GET,ewald_prec)
                  call twh_rcelect(GLB_GET,BXD_REAL,ibox,rcelect)
                  write(inpnum,*) 'coulomb style               ewald'
     &                 ,rcelect,ewald_prec
c                 --- Lennard-Jones 9-6
                  type = twh_getnbtype(type,type)
                  call twh_nbcoeff(GLB_GET,2,type,epsilon)
                  call twh_nbcoeff(GLB_GET,1,type,sigma)
                  sigma = (sigma)**(1.0/3.0)
                  write(filenum,*) count_types,epsilon/CON_KCALTOK
     &                 ,sigma
               elseif ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c                 --- need to output potential style to input
                  write(inpnum,'(a14)') 'pair_style eam'
                  call twh_nmolty(GLB_GET,nmolty)
                  if ( nmolty .eq. 1 ) then
c                    --- we can handle this case with the writedynamo
                     write(inpnum,'(a21)') 'pair_coeff dynamo_out'
                  endif
               else
                  write(6,*) 'cpnumber: ',cpnumber
                  write(6,*) 'invalid classical potential number for'
     &                 ,' outputting LAMMPS'
                  goto 666
               endif
            enddo
            if ( cpnumber .ne. CP_EMBEDDED_ATOM_METHOD ) then
               write(filenum,*)
            endif
         endif

         if ( lam_bond_types .ne. 0 ) then
c           --- output bond coeffs
            write(filenum,5)
 5          format('Bond Coeffs')
            write(filenum,*)
            do count_types = 1,lam_bond_types
               type = decode_bond(count_types)
               call twh_vibstyle(GLB_GET,type,vibstyle)
               if ( vibstyle .eq. BOND_STYLE_HARMONIC ) then
c                 --- output bond style
                  write(inpnum,11)
 11               format('bond style             harmonic')
c                 --- output coeffs
                  call twh_vibcoeff(GLB_GET,type,0,vibcoeff(0))
                  call twh_vibcoeff(GLB_GET,type,1,vibcoeff(1))
                  write(filenum,*) count_types,vibcoeff(1)/CON_KCALTOK
     &                 ,vibcoeff(0)
               elseif ( vibstyle .eq. BOND_STYLE_COMPASS ) then
c                 --- compass quartic class2 potential
                  write(inpnum,'(a27)') 'bond style           class2'
c                 --- output coeffs
                  do icoeff = 0,3
                     call twh_vibcoeff(GLB_GET,type,icoeff
     &                    ,vibcoeff(icoeff))
                  enddo
                  write(filenum,*) count_types,vibcoeff(0)
     &                 ,(vibcoeff(icoeff)/CON_KCALTOK,icoeff=1,3)
               else
                  if ( .not. lwarnvibstyle(vibstyle) ) then
                     write(6,*) 'WRITELAMMPS: invalid vibration style '
     &                    ,vibstyle,' has nothing similar in lammps'
                     lwarnvibstyle(vibstyle) = .true.
                  endif
c                 --- fake entry
                  write(filenum,*) count_types,(0.0d0,icoeff=1,4)
               endif
            enddo
            write(filenum,*)
         endif
         
         lclass2angle = .false.
         if ( lam_angle_types .ne. 0 ) then
c           --- output angle coeffs
            write(filenum,6)
 6          format('Angle Coeffs')
            write(filenum,*)
            do count_types = 1,lam_angle_types
               type = decode_angle(count_types)
c              --- retrieve all bencoeffs
               do icoeff = MINBENDCOEFF,MAXBENDCOEFF
                  call twh_bencoeff(GLB_GET,type,icoeff
     &                 ,bencoeff(icoeff))
               enddo
               call twh_anglestyle(GLB_GET,type,anglestyle)
               if ( anglestyle .eq. ANG_STYLE_HARMONIC ) then
c                 --- standard harmonic, output style as charmm for
c                 --- compatibility
                  write(inpnum,13)
c                 --- output coeffs
                  write(filenum,*) count_types
     &                 ,bencoeff(1)/CON_KCALTOK
     &                 ,bencoeff(0)*radtodeg,0.0d0,0.0d0
               elseif ( anglestyle .eq. ANG_STYLE_COMPASS_AUTO ) then
c                 --- compass with autodection of terms
c                 --- only output the regular bond terms here, but set
c                 --- logical to output the rest later
                  lclass2angle = .true.
                  write(inpnum,'(a22)') 'angle style     class2'
                  write(filenum,*) count_types
     &                 ,bencoeff(0)*radtodeg
     &                 ,(bencoeff(icoeff)/CON_KCALTOK,icoeff=1,3)
               elseif ( anglestyle .eq. ANG_STYLE_CHARMM_UB ) then
c                 --- harmonic with urey-bradley
                  write(inpnum,13)
 13               format('angle style            charmm')
c                 --- output coeffs
                  write(filenum,*) count_types,bencoeff(1)/CON_KCALTOK
     &                 ,bencoeff(0)*radtodeg,bencoeff(3)/CON_KCALTOK
     &                 ,bencoeff(2)

               elseif ( anglestyle .eq. ANG_STYLE_COMPASS_EX ) then
c                 --- compass with explicit declaration of terms
c                 --- only output the regular bond terms here, but set
c                 --- logical to output the rest later
                  lclass2angle = .true.
                  write(inpnum,'(a22)') 'angle style     class2'
                  write(filenum,*) count_types
     &                 ,bencoeff(0)*radtodeg
     &                 ,(bencoeff(icoeff)/CON_KCALTOK,icoeff=1,3)
               else
c                 --- don't know how to output this for lammps
                  if ( .not. lwarnangstyle(anglestyle) ) then
                     write(6,*) 'WRITELAMMPS: invalid angle style '
     &                    ,anglestyle,' has nothing similar in lammps'
                     lwarnangstyle(anglestyle) = .true.
                  endif
c                 --- output a dummy entry
                  write(filenum,*) count_types
     &                 ,bencoeff(0)*radtodeg,(0.0d0,icoeff=1,3)
               endif
            enddo
            write(filenum,*)
         endif

         if ( lclass2angle ) then
c           --- need to output the bondbond and the bondangle coeffs
c           --- output bondbond coeffs
            write(filenum,'(a15)') 'BondBond Coeffs'
            write(filenum,*)
            do count_types = 1,lam_angle_types
               type = decode_angle(count_types)
c              --- retrieve all bencoeffs
               do icoeff = MINBENDCOEFF,MAXBENDCOEFF
                  call twh_bencoeff(GLB_GET,type,icoeff
     &                 ,bencoeff(icoeff))
               enddo
               call twh_anglestyle(GLB_GET,type,anglestyle)
               if ( anglestyle .eq. ANG_STYLE_COMPASS_AUTO ) then
c                 --- output the bondbond coefficients in the approp
c                 --- order ( slightly different from style 8 )
                  write(filenum,*) count_types
     &                 ,bencoeff(6)/CON_KCALTOK
     &                 ,(bencoeff(icoeff),icoeff=9,10)
               elseif ( anglestyle .eq. ANG_STYLE_COMPASS_EX ) then
c                 --- output the bondbond coefficients
                  write(filenum,*) count_types
     &                 ,bencoeff(8)/CON_KCALTOK
     &                 ,(bencoeff(icoeff),icoeff=9,10)
               else
c                 --- these terms don't exist for this angle style
c                 --- can't really mix class2 with others in LAMMPS
c                 --- don't know how to output this for lammps
                  write(6,*) 'invalid bondbond style ',anglestyle
                  write(6,*) 'LAMMPS cannot mix class2 with others'
                  goto 666
               endif
            enddo
            write(filenum,*)

c           --- output bondangle coeffs
            write(filenum,'(a16)') 'BondAngle Coeffs'
            write(filenum,*)
            do count_types = 1,lam_angle_types
               type = decode_angle(count_types)
               do icoeff = MINBENDCOEFF,MAXBENDCOEFF
                  call twh_bencoeff(GLB_GET,type,icoeff
     &                 ,bencoeff(icoeff))
               enddo
               call twh_anglestyle(GLB_GET,type,anglestyle)
               if ( anglestyle .eq. ANG_STYLE_COMPASS_AUTO ) then
c                 --- output the bondangle coefficients (slightly 
c                 --- different from style 8 )
                  write(filenum,*) count_types
     &                 ,bencoeff(4)/CON_KCALTOK
     &                 ,bencoeff(5)/CON_KCALTOK
     &                 ,bencoeff(9),bencoeff(10)
               elseif ( anglestyle .eq. ANG_STYLE_COMPASS_EX ) then
c                 --- output the bondangle coefficients
                  write(filenum,*) count_types
     &                 ,bencoeff(4)/CON_KCALTOK
     &                 ,bencoeff(6)/CON_KCALTOK
     &                 ,bencoeff(5),bencoeff(7)
               else
c                 --- these terms don't exist for this angle style
c                 --- can't really mix class2 with others in LAMMPS
c                 --- don't know how to output this for lammps
                  write(6,*) 'invalid bondangle style ',anglestyle
                  write(6,*) 'LAMMPS cannot mix class2 with others'
                  goto 666
               endif
            enddo
            write(filenum,*)
         endif

         lclass2dihed = .false.
         if ( lam_dihedral_types .ne. 0 ) then
c           --- output dihedral coeffs
            write(filenum,7)
 7          format('Dihedral Coeffs')
            write(filenum,*)
            do count_types = 1,lam_dihedral_types
               type = decode_dihedral(count_types,1)
               call twh_torstyle(GLB_GET,type,torstyle)
               do icoeff = MINTORCOEFF,MAXTORCOEFF
                  call twh_torcoeff(GLB_GET,type,icoeff
     &                 ,torcoeff(icoeff))
               enddo
               if ( torstyle .eq. TOR_STYLE_GCA_COS
     &              .or. torstyle .eq. TOR_STYLE_GCA_CUB ) then
c                 --- Amber/Charmm/Gromos cosine series
                  write(inpnum,14)
 14               format('dihedral style         charmm')
c                 --- output dihedral coeffs
                  itricky = decode_dihedral(count_types,2)
                  if ( itricky .eq. 0 ) then
c                    --- this is a null torsion, needed for 1-4
                     write(filenum,*) count_types,' 0.0 1 0.0 1.0'
                  else
c                    --- normal charmm torsion style
c                    --- convert itricky into needed notation
                     itricky = 3*(itricky-1)
                     write(filenum,*) count_types
     &                    ,torcoeff(itricky+1)/CON_KCALTOK
     &                    ,nint(torcoeff(itricky+2))
     &                    ,torcoeff(itricky+3)*radtodeg
     &                    ,' 1.0'
                  endif
               elseif ( torstyle .eq. TOR_STYLE_COMP_X_AD ) then
c                 --- compass with autodetected cross terms
                  lclass2dihed = .true.
                  write(inpnum,'(a22)') 'dihedral style  class2'
                  write(filenum,*) count_types
     &                 ,torcoeff(0)/CON_KCALTOK,0.0d0
     &                 ,torcoeff(1)/CON_KCALTOK,0.0d0
     &                 ,torcoeff(2)/CON_KCALTOK,0.0d0
               elseif ( torstyle .eq. TOR_STYLE_COMP_R_AD ) then
c                 --- compass without cross terms
                  write(inpnum,'(a22)') 'dihedral style  class2'
                  write(filenum,*) count_types
     &                 ,torcoeff(0)/CON_KCALTOK,0.0d0
     &                 ,torcoeff(1)/CON_KCALTOK,0.0d0
     &                 ,torcoeff(2)/CON_KCALTOK,0.0d0
               elseif ( torstyle .eq. TOR_STYLE_COMP_X_EX ) then
c                 --- compass with explicitly declared cross terms
                  lclass2dihed = .true.
                  write(inpnum,'(a22)') 'dihedral style  class2'
                  write(filenum,*) count_types
     &                 ,torcoeff(0)/CON_KCALTOK
     &                 ,torcoeff(3)*radtodeg
     &                 ,torcoeff(1)/CON_KCALTOK
     &                 ,torcoeff(4)*radtodeg
     &                 ,torcoeff(2)/CON_KCALTOK
     &                 ,torcoeff(5)*radtodeg
               elseif ( torstyle .eq. TOR_STYLE_COS_POW ) then
c                 --- cosine power series - LAMMPS always does a 5 term
c                 --- series so if we have more then that we cannot
c                 --- output
                  call twh_ntorloop(GLB_GET,type,ntorloop)
                  if ( ntorloop .le. 5 ) then
                     write(filenum,*) count_types
     &                    ,(torcoeff(itype)/CON_KCALTOK,itype=1,5)
                  else
                     write(6,*) 'cannot output dihedral style 10'
     &                    ,' with ntorloop > 5'
                     goto 666
                  endif
               elseif ( torstyle .eq. TOR_STYLE_AMBER ) then
c                 --- Amber/Charmm/Gromos cosine series / ndiv
                  write(inpnum,14)
c                 --- output dihedral coeffs
                  itricky = decode_dihedral(count_types,2)
                  if ( itricky .eq. 0 ) then
c                    --- this is a null torsion, needed for 1-4
                     write(filenum,*) count_types,' 0.0 1 0.0 1.0'
                  else
c                    --- normal charmm torsion style
c                    --- convert itricky into needed notation
                     itricky = 3*(itricky-1)
                     write(filenum,*) count_types
     &                    ,torcoeff(itricky+1)
     &                    /(dble(div_dihedral(count_types))
     &                    *CON_KCALTOK)
     &                    ,nint(torcoeff(itricky+2))
     &                    ,torcoeff(itricky+3)*radtodeg,' 1.0'
                  endif
               else
c                 --- invalid style
                  if ( .not. lwarndihedstyle(torstyle) ) then
                     write(6,*) 'WRITELAMMPS: invalid dihedral style '
     &                    ,torstyle,' has nothing similar in lammps'
                     lwarndihedstyle(torstyle) = .true.
                  endif
c                 --- output dummy entry
                  write(filenum,*) count_types,' 0.0 1 0.0 '
     &                 ,torcoeff(-1)*scale_dihedral(count_types)
               endif
            enddo
            write(filenum,*)
         endif

         if ( lclass2dihed ) then
c           --- output all of the torsion cross terms
c           --- output middlebondtorsion coeff
            write(filenum,'(a24)') 'MiddleBondTorsion Coeffs'
            write(filenum,*)
            do count_types = 1,lam_dihedral_types
               type = decode_dihedral(count_types,1)
               do icoeff = MINTORCOEFF,MAXTORCOEFF
                  call twh_torcoeff(GLB_GET,type,icoeff
     &                 ,torcoeff(icoeff))
               enddo
               call twh_torstyle(GLB_GET,type,torstyle)
               if ( torstyle .eq. TOR_STYLE_COMP_X_AD ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(3)/CON_KCALTOK
     &                 ,torcoeff(4)/CON_KCALTOK
     &                 ,torcoeff(5)/CON_KCALTOK
     &                 ,torcoeff(20)
               elseif ( torstyle .eq. TOR_STYLE_COMP_R_AD ) then
c                  --- nothing to do here, but comes up as this mixes
c                  --- in the same forcefield as those that have cross
c                  --- terms
               elseif ( torstyle .eq. TOR_STYLE_COMP_X_EX ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(6)/CON_KCALTOK
     &                 ,torcoeff(7)/CON_KCALTOK
     &                 ,torcoeff(8)/CON_KCALTOK
     &                 ,torcoeff(9)
               else
c                 --- invalid style
                  write(6,*) 'middlebondtorsion style ',torstyle
                  write(6,*) 'cannot output this to LAMMPS'
                  goto 666
               endif
            enddo
            write(filenum,*)

c           --- output endbondtorsion coeff
            write(filenum,'(a21)') 'EndBondTorsion Coeffs'
            write(filenum,*)
            do count_types = 1,lam_dihedral_types
               type = decode_dihedral(count_types,1)
               call twh_torstyle(GLB_GET,type,torstyle)
               if ( torstyle .eq. TOR_STYLE_COMP_X_AD ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(6)/CON_KCALTOK
     &                 ,torcoeff(7)/CON_KCALTOK
     &                 ,torcoeff(8)/CON_KCALTOK
     &                 ,torcoeff(9)/CON_KCALTOK
     &                 ,torcoeff(10)/CON_KCALTOK
     &                 ,torcoeff(11)/CON_KCALTOK
     &                 ,torcoeff(21),torcoeff(22)
               elseif ( torstyle .eq. TOR_STYLE_COMP_R_AD ) then
c                  --- nothing to do here, but comes up as this mixes
c                  --- in the same forcefield as those that have cross
c                  --- terms
               elseif ( torstyle .eq. TOR_STYLE_COMP_X_EX ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(10)/CON_KCALTOK
     &                 ,torcoeff(11)/CON_KCALTOK
     &                 ,torcoeff(12)/CON_KCALTOK
     &                 ,torcoeff(14)/CON_KCALTOK
     &                 ,torcoeff(15)/CON_KCALTOK
     &                 ,torcoeff(16)/CON_KCALTOK
     &                 ,torcoeff(13),torcoeff(17)
               else
c                 --- invalid style
                  write(6,*) 'endbondtorsion style ',torstyle
                  write(6,*) 'cannot output this to LAMMPS'
                  goto 666
               endif
            enddo
            write(filenum,*)

c           --- output angletorsion coeff
            write(filenum,'(a19)') 'AngleTorsion Coeffs'
            write(filenum,*)
            do count_types = 1,lam_dihedral_types
               type = decode_dihedral(count_types,1)
               do icoeff = MINTORCOEFF,MAXTORCOEFF
                  call twh_torcoeff(GLB_GET,type,icoeff
     &                 ,torcoeff(icoeff))
               enddo
               call twh_torstyle(GLB_GET,type,torstyle)
               if ( torstyle .eq. TOR_STYLE_COMP_X_AD ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(12)/CON_KCALTOK
     &                 ,torcoeff(13)/CON_KCALTOK
     &                 ,torcoeff(14)/CON_KCALTOK
     &                 ,torcoeff(15)/CON_KCALTOK
     &                 ,torcoeff(16)/CON_KCALTOK
     &                 ,torcoeff(17)/CON_KCALTOK
     &                 ,torcoeff(23)*radtodeg
     &                 ,torcoeff(24)*radtodeg
               elseif ( torstyle .eq. TOR_STYLE_COMP_R_AD ) then
c                  --- nothing to do here, but comes up as this mixes
c                  --- in the same forcefield as those that have cross
c                  --- terms
               elseif ( torstyle .eq. TOR_STYLE_COMP_X_EX ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(18)/CON_KCALTOK
     &                 ,torcoeff(19)/CON_KCALTOK
     &                 ,torcoeff(20)/CON_KCALTOK
     &                 ,torcoeff(21)/CON_KCALTOK
     &                 ,torcoeff(22)/CON_KCALTOK
     &                 ,torcoeff(23)/CON_KCALTOK
     &                 ,torcoeff(24)*radtodeg
     &                 ,torcoeff(25)*radtodeg
               else
c                 --- invalid style
                  write(6,*) 'angletorsion style ',torstyle
                  write(6,*) 'cannot output this to LAMMPS'
                  goto 666
               endif
            enddo
            write(filenum,*)

c           --- output angleangletorsion coeff
            write(filenum,'(a24)') 'AngleAngleTorsion Coeffs'
            write(filenum,*)
            do count_types = 1,lam_dihedral_types
               type = decode_dihedral(count_types,1)
               do icoeff = MINTORCOEFF,MAXTORCOEFF
                  call twh_torcoeff(GLB_GET,type,icoeff
     &                 ,torcoeff(icoeff))
               enddo
               call twh_torstyle(GLB_GET,type,torstyle)
               if ( torstyle .eq. TOR_STYLE_COMP_X_AD ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(18)/CON_KCALTOK
     &                 ,torcoeff(23)*radtodeg
     &                 ,torcoeff(24)*radtodeg
               elseif ( torstyle .eq. TOR_STYLE_COMP_R_AD ) then
c                  --- nothing to do here, but comes up as this mixes
c                  --- in the same forcefield as those that have cross
c                  --- terms
               elseif ( torstyle .eq. TOR_STYLE_COMP_X_EX ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(26)/CON_KCALTOK
     &                 ,torcoeff(27)*radtodeg
     &                 ,torcoeff(28)*radtodeg
               else
c                 --- invalid style
                  write(6,*) 'angleangletorsion style ',torstyle
                  write(6,*) 'cannot output this to LAMMPS'
                  goto 666
               endif
            enddo
            write(filenum,*)

c           --- bondbond13 coeff
            write(filenum,'(a17)') 'BondBond13 Coeffs'
            write(filenum,*)
            do count_types = 1,lam_dihedral_types
               type = decode_dihedral(count_types,1)
               do icoeff = MINTORCOEFF,MAXTORCOEFF
                  call twh_torcoeff(GLB_GET,type,icoeff
     &                 ,torcoeff(icoeff))
               enddo
               call twh_torstyle(GLB_GET,type,torstyle)
               if ( torstyle .eq. TOR_STYLE_COMP_X_AD ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(19)/CON_KCALTOK
     &                 ,torcoeff(20),torcoeff(22)
               elseif ( torstyle .eq. TOR_STYLE_COMP_R_AD ) then
c                  --- nothing to do here, but comes up as this mixes
c                  --- in the same forcefield as those that have cross
c                  --- terms
               elseif ( torstyle .eq. TOR_STYLE_COMP_X_EX ) then
c                 --- compass with explicitly declared cross terms
                  write(filenum,*) count_types
     &                 ,torcoeff(29)/CON_KCALTOK
     &                 ,torcoeff(30),torcoeff(31)
               else
c                 --- invalid style
                  write(6,*) 'bondbond13 style ',torstyle
                  write(6,*) 'cannot output this to LAMMPS'
                  goto 666
               endif
            enddo
            write(filenum,*)
         endif

         if ( lam_improp_types .ne. 0 .or. lam_angang_types .ne. 0) then
c           --- output improper coeffs
c           --- note that for class 2 LAMMPS has this wacky rule that 
c           --- improper types and angle-angle types are combined.  Fill
c           --- with zeros for any angle-angle types
            write(filenum,8)
 8          format('Improper Coeffs')
            write(filenum,*)
            do count_types = 1,lam_improp_types
               type = decode_improp(count_types,1)
               do icoeff = MINTORCOEFF,MAXTORCOEFF
                  call twh_torcoeff(GLB_GET,type,icoeff
     &                 ,torcoeff(icoeff))
               enddo
               if ( decode_improp(count_types,2) .eq. 0 ) then
                  call twh_torstyle(GLB_GET,type,torstyle)
c                 --- out of plane version of improper stored in towhee
c                 --- as a regular torsion
                  if ( torstyle .eq. TOR_STYLE_GCA_CUB ) then
c                    --- take the harmonic part of this torsion
                     write(inpnum,15)
                     call twh_ntorloop(GLB_GET,type,ntorloop)
                     index = 3*ntorloop +1
                     write(filenum,*) count_types
     &                    ,torcoeff(index+0)/CON_KCALTOK
     &                    ,torcoeff(index+1)*radtodeg
                  else
                     write(6,*) 'improper from regular tor style '
     &                    ,torstyle
                     write(6,*) 'cannot output this for LAMMPS'
                     goto 666
                  endif
               elseif ( decode_improp(count_types,2) .eq. 1 ) then
                  call twh_impstyle(GLB_GET,type,impstyle)
c                 --- improper stored as an impropers
                  call twh_impcoeff(GLB_GET,type,0,impcoeff_zero)
                  call twh_impcoeff(GLB_GET,type,1,impcoeff_one)
                  if ( impstyle .eq. 1 ) then
c                    --- harmonic
                     write(inpnum,15)
 15                  format('improper style         harmonic')
                     write(filenum,*) count_types
     &                    ,impcoeff_zero/CON_KCALTOK
     &                    ,impcoeff_one*radtodeg
                  elseif ( impstyle .eq. 5 ) then
c                    --- wilson oop class 2 harmonic
                     write(inpnum,'(a29)') 
     &                    'improper style         class2'
                     write(filenum,*) count_types
     &                    ,impcoeff_zero/CON_KCALTOK
     &                    ,impcoeff_one*radtodeg
                  else
                     write(6,*) 'oop improper style ',impstyle
                     write(6,*) 'cannot output this for LAMMPS'
                     goto 666
                  endif
               else
c                 --- stereo improper type
                  write(6,*) 'cannot output stereo improper to LAMMPS'
                  goto 666
               endif
            enddo
c           --- here we fill in for any angle-angle types
            do count_types = lam_improp_types+1
     &           ,lam_improp_types+lam_angang_types
               write(filenum,*) count_types,0.0d0,0.0d0
            enddo
            write(filenum,*)
         endif

         if ( lam_angang .ne. 0 ) then
c           --- output angle-angle coeficients
            write(filenum,'(a17)') 'AngleAngle Coeffs'
            write(filenum,*)
c           --- first we have to fill in for any impropers
            do icoeff = 1,6
               tmp(icoeff) = 0.0d0
            enddo
            do count_types = 1,lam_improp_types
               write(filenum,*) count_types,(tmp(icoeff),icoeff=1,6)
            enddo
            do count_types = 1,lam_angang_types
c              --- types are stored a bit differently  
c              --- search through code_angang to find this type
               do iaa = -WLAAMAX,WLAAMAX
                  do jaa = -WLAAMAX,WLAAMAX
                     do kaa = -WLAAMAX,WLAAMAX
                        if ( code_angang(iaa,jaa,kaa) .eq. count_types 
     &                       ) then
c                          --- determine the coeficients
                           call twh_aacoeff(GLB_GET,iaa,0,aacoeff)
                           tmp(1) = aacoeff/CON_KCALTOK
                           call twh_aacoeff(GLB_GET,abs(jaa),0
     &                          ,aacoeff)
                           tmp(2) = aacoeff/CON_KCALTOK
                           call twh_aacoeff(GLB_GET,abs(kaa),0,aacoeff)
                           tmp(3) = aacoeff/CON_KCALTOK
                           call twh_aacoeff(GLB_GET,iaa,1,aacoeff)
                           tmp(4) = aacoeff*radtodeg
                           call twh_aacoeff(GLB_GET,iaa,2,aacoeff)
                           tmp(5) = aacoeff*radtodeg
                           if ( jaa .gt. 0 ) then
                              call twh_aacoeff(GLB_GET,jaa,2,aacoeff)
                              tmp(6) = aacoeff*radtodeg
                           else
                              call twh_aacoeff(GLB_GET,abs(jaa),1
     &                             ,aacoeff)
                              tmp(6) = aacoeff*radtodeg
                           endif
                           write(filenum,*) count_types+lam_improp_types
     &                          ,(tmp(icoeff),icoeff=1,6)
                        endif
                     enddo
                  enddo
               enddo
            enddo
            write(filenum,*)
         endif

         if ( lam_atom .ne. 0 ) then
c           --- output atom positions and types for each atom
            write(filenum,'(a5)') 'Atoms'
            write(filenum,*)
            count_molecule = 0
            count_atom = 0
            do ichain = 1,nchain
               call twh_nboxi(GLB_GET,ichain,testboxi)
               if ( testboxi .eq. ibox) then
                  count_molecule = count_molecule + 1
                  call twh_moltyp(GLB_GET,ichain,imolty)
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1,nunit
                     count_atom = count_atom + 1
                     call twh_ntype(GLB_GET,imolty,iunit,type)
                     call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &                    ,xtemp,ytemp,ztemp)
                     if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c                       --- special format without charges or molecule
c                       --- types
                        write(filenum,'(1x,i6,1x,i3,3(1x,f9.5))') 
     &                       count_atom,code_atom(type)
     &                       ,xtemp,ytemp,ztemp
                     else
c                       --- default format
                        call twh_qqatom(GLB_GET,imolty,iunit,qqatom)
                        write(filenum,'(1x,i6,1x,i6,1x,i3,2x,f6.3
     &                       ,3(1x,f9.5),3(1x,i4))') 
     &                       count_atom,count_molecule
     &                       ,code_atom(type),qqatom
     &                       ,xtemp,ytemp,ztemp,0,0,0
                     endif
                  enddo
               endif
            enddo
            if ( lam_atom .ne. count_atom ) then
               write(6,*) 'Some sort of error in writelammps'
               write(6,*) 'lam_atom ',lam_atom
               write(6,*) 'count_atom ',count_atom
            endif
            write(filenum,*)
         endif

         if ( lam_bond .ne. 0 ) then
            write(filenum,2)
 2          format('Bonds')
            write(filenum,*)
            count_bond = 0
            do ichain = 1,nchain
               call twh_nboxi(GLB_GET,ichain,testboxi)
               if ( testboxi .eq. ibox ) then
                  call twh_moltyp(GLB_GET,ichain,imolty)
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1,nunit
                     call twh_invib(GLB_GET,imolty,iunit,invibration)
                     do ivib = 1,invibration
                        call twh_ijvib(GLB_GET,imolty,iunit,ivib
     &                       ,ijvibtest)
                        if ( ijvibtest .gt. iunit ) then
                           count_bond = count_bond + 1
                           call twh_itvib(GLB_GET,imolty,iunit,ivib
     &                          ,type)
                           write(filenum,*) count_bond
     &                          ,code_bond(type)
     &                          ,twh_oldgetatomnum(ichain,iunit)
     &                          ,twh_oldgetatomnum(ichain,ijvibtest)
                        endif
                     enddo
                  enddo
               endif
            enddo
            if ( lam_bond .ne. count_bond ) then
               write(6,*) 'Error in writelammps'
               write(6,*) 'lam_bond ',lam_bond
               write(6,*) 'count_bond ',count_bond
            endif
            write(filenum,*)
         endif

         if ( lam_angle .ne. 0 ) then
            write(filenum,9)
 9          format('Angles')
            write(filenum,*)
            count_angle = 0
            do ichain = 1,nchain
               call twh_nboxi(GLB_GET,ichain,testboxi)
               if ( testboxi .eq. ibox ) then
                  call twh_moltyp(GLB_GET,ichain,imolty)
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1,nunit
                     call twh_inben(GLB_GET,imolty,iunit,inbend)
                     do ibend = 1,inbend
                        call twh_itben(GLB_GET,imolty,iunit,ibend,type)
                        call twh_anglestyle(GLB_GET,type,anglestyle)
                        call twh_ijben3(GLB_GET,imolty,iunit,ibend
     &                       ,ijben3test)
                        if ( anglestyle .eq. ANG_STYLE_COMPASS_EX .or. 
     &                       ijben3test .gt.iunit ) then
                           call twh_ijben2(GLB_GET,imolty,iunit
     &                          ,ibend,ijben2test)
                           count_angle = count_angle + 1
                           write(filenum,*) count_angle
     &                          ,code_angle(type)
     &                          ,twh_oldgetatomnum(ichain,iunit)
     &                          ,twh_oldgetatomnum(ichain,ijben2test)
     &                          ,twh_oldgetatomnum(ichain,ijben3test)
                        endif
                     enddo
                  enddo
               endif
            enddo
            if ( lam_angle .ne. count_angle ) then
               write(6,*) 'problem with angles in writelammps'
               write(6,*) 'lam_angle ',lam_angle
               write(6,*) 'count_angle ',count_angle
            endif
            write(filenum,*)
         endif

         if ( lam_dihedral .ne. 0 ) then
            write(filenum,16)
 16         format('Dihedrals')
            write(filenum,*)
            count_dihedral = 0
            do ichain = 1,nchain
               call twh_nboxi(GLB_GET,ichain,testboxi)
               if ( testboxi .eq. ibox ) then
                  call twh_moltyp(GLB_GET,ichain,imolty)
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1,nunit
                     call twh_intor(GLB_GET,imolty,iunit,intorsion)
                     do itor = 1,intorsion
                        call twh_ittor(GLB_GET,imolty,iunit,itor,type)
                        if ( type .gt. 0 ) then
                           call twh_ijtor4(GLB_GET,imolty,iunit,itor
     &                          ,ijtor4test)
                           call twh_torstyle(GLB_GET,type,torstyle)
                           if ( torstyle .eq. TOR_STYLE_COMP_X_EX .or. 
     &                          ijtor4test .gt. iunit ) then
                              call twh_ntorloop(GLB_GET,type,ntorloop)
                              if ( ntorloop .eq. 0 ) then
                                 itotal = 1
                              else
                                 itotal = ntorloop
                              endif
                              do itricky = 1,itotal
                                 count_dihedral = count_dihedral + 1
c                                --- have to search to find the lammps 
c                                --- dihedral types
                                 lfound = .false.
                                 do count_types = 1,lam_dihedral_types
c                                 --- see if we have a complete match
                                    if ( type .eq. 
     &                                   decode_dihedral(count_types,1)
     &                                   ) then
                                       if ( ntorloop .eq. 0 .or. 
     &                                      itricky .eq. 
     &                                      decode_dihedral(count_types
     &                                      ,2) ) then
c                                         --- check scaling factor
                                          if ( ntorloop .eq. 0 ) then
                                             call twh_itscale(GLB_GET
     &                                            ,imolty,iunit,itor
     &                                            ,scaletest)
                                          else
                                             call twh_itscale(GLB_GET
     &                                            ,imolty,iunit,itor
     &                                            ,dvalue)
                                             scaletest = dvalue /
     &                                            dble(ntorloop)
                                          endif
                                          if ( (scaletest - 
     &                                         scale_dihedral(
     &                                         count_types))**2
     &                                         .lt.  tolerance ) then
c                                          --- match found
                                             lfound = .true.
                                             iwinner = count_types
                                          endif
                                       endif
                                    endif
                                 enddo
                                 if ( .not. lfound ) then
                                    write(6,*) 'WRITELAMMPS'
                                    write(6,*) 'dihedral matchup failed'
                                    goto 666
                                 endif
                                 call twh_ijtor2(GLB_GET,imolty,iunit
     &                                ,itor,ijtor2test)
                                 call twh_ijtor3(GLB_GET,imolty,iunit
     &                                ,itor,ijtor3test)
                                 call twh_ijtor4(GLB_GET,imolty,iunit
     &                                ,itor,ijtor4test)
                                 write(filenum,*) count_dihedral,iwinner
     &                                ,twh_oldgetatomnum(ichain,iunit)
     &                                ,twh_oldgetatomnum(ichain
     &                                ,ijtor2test)
     &                                ,twh_oldgetatomnum(ichain
     &                                ,ijtor3test)
     &                                ,twh_oldgetatomnum(ichain
     &                                ,ijtor4test)
                              enddo
                           endif
                        endif
                     enddo
                  enddo
               endif
            enddo
            write(filenum,*)
         endif
         
         if ( lam_improp + lam_angang .ne. 0 ) then
            write(filenum,17)
 17         format('Impropers')
            write(filenum,*)
            count_improp = 0
            do ichain = 1,nchain
               call twh_nboxi(GLB_GET,ichain,testboxi)
               if ( testboxi .eq. ibox ) then
                  call twh_moltyp(GLB_GET,ichain,imolty)
                  call twh_nunit(GLB_GET,imolty,nunit)
                  do iunit = 1,nunit
c                    --- retrieve inaa
                     call twh_inaa(GLB_GET,imolty,iunit,inaa)
c                    --- lammps impropers that towhee stores as improper
                     call twh_inimprop(GLB_GET,imolty,iunit,inimproper)
                     do imp = 1,inimproper
                        call twh_itimprop(GLB_GET,imolty,iunit,imp,type)
                        call twh_impform(GLB_GET,type,improperform)
                        if ( improperform .eq. 2
     &                       .or. improperform .eq. 3 ) then
                           count_improp = count_improp + 1
                           call twh_ijimprop2(GLB_GET,imolty,iunit,imp
     &                          ,ijimp2)
                           call twh_ijimprop3(GLB_GET,imolty,iunit,imp
     &                          ,ijimp3)
                           call twh_ijimprop4(GLB_GET,imolty,iunit,imp
     &                          ,ijimp4)
                           if ( improperform .eq. 2 ) then
c                             --- traditional oop improper ordering
                              write(filenum,*) count_improp
     &                             ,code_improp(type)
     &                             ,twh_oldgetatomnum(ichain,iunit)
     &                             ,twh_oldgetatomnum(ichain,ijimp2)
     &                             ,twh_oldgetatomnum(ichain,ijimp3)
     &                             ,twh_oldgetatomnum(ichain,ijimp4)
                           elseif ( improperform .eq. 3 ) then
c                          --- wilson oop improper ordering
                              write(filenum,*) count_improp
     &                             ,code_improp(type)
     &                             ,twh_oldgetatomnum(ichain,ijimp2)
     &                             ,twh_oldgetatomnum(ichain,iunit)
     &                             ,twh_oldgetatomnum(ichain,ijimp3)
     &                             ,twh_oldgetatomnum(ichain,ijimp4)
                           endif
                        endif
                     enddo
c                    --- lammps impropers that towhee stores as torsions
                     call twh_intor(GLB_GET,imolty,iunit,intorsion)
                     do itor = 1,intorsion
                        call twh_ijtor4(GLB_GET,imolty,iunit,itor
     &                       ,ijtor4test)
                        call twh_ittor(GLB_GET,imolty,iunit,itor
     &                       ,ittortest)
                        call twh_torstyle(GLB_GET,abs(ittortest)
     &                       ,torstyle)
                        if ( ijtor4test .gt. iunit .and.
     &                       torstyle .eq. TOR_STYLE_GCA_CUB ) then
c                           --- towhee regular torsion that is a lmp imp
                           count_improp = count_improp + 1
                           type = abs(ittortest)
                           call twh_ijtor2(GLB_GET,imolty,iunit,itor
     &                          ,ijtor2test)
                           call twh_ijtor3(GLB_GET,imolty,iunit,itor
     &                          ,ijtor3test)
                           call twh_ijtor4(GLB_GET,imolty,iunit,itor
     &                          ,ijtor4test)
                           write(filenum,*) count_improp
     &                          ,code_improp(-type)
     &                          ,twh_oldgetatomnum(ichain,iunit)
     &                          ,twh_oldgetatomnum(ichain,ijtor2test)
     &                          ,twh_oldgetatomnum(ichain,ijtor3test)
     &                          ,twh_oldgetatomnum(ichain,ijtor4test)
                        endif
                     enddo
c                    --- lammps impropers that towhee stores as angang
                     do iaa = 1,inaa
                        lassign(iaa) = .false.
                     enddo
                     ldone = .false.
                     do while (.not. ldone .and. inaa .ne. 0 ) 
c                       --- find the first unassigned angang
                        triple(1) = 0
                        iaa = 0
                        do while ( triple(1) .eq. 0 ) 
                           iaa = iaa + 1
                           if ( .not. lassign(iaa) ) then
c                             --- take this as the first in the triple
                              call twh_itaa(GLB_GET,imolty,iunit,iaa
     &                             ,triple(1))
                              lassign(iaa) = .true.
                           endif
                        enddo
c                       --- extract the first triple
                        call twh_ijaa0(GLB_GET,imolty,iunit,iaa,atomc)
                        call twh_ijaa1(GLB_GET,imolty,iunit,iaa,atoma)
                        call twh_ijaa2(GLB_GET,imolty,iunit,iaa,atomd)
                        do jaa = 1,inaa
                           if ( .not. lassign(jaa) ) then
c                             --- see if some permutation of these atoms
c                             --- matches the first atoms
                              call twh_ijaa0(GLB_GET,imolty,iunit,jaa
     &                             ,test0)
                              call twh_ijaa1(GLB_GET,imolty,iunit,jaa
     &                             ,test1)
                              call twh_ijaa2(GLB_GET,imolty,iunit,jaa
     &                             ,test2)
                              if ( atoma .eq. test0 
     &                             .and. atomc .eq. test1
     &                             .and. atomd .eq. test2 ) then
c                                 --- this fills the second slot in the 
c                                 --- lammps angle-angle
                                 call twh_itaa(GLB_GET,imolty,iunit,jaa
     &                                ,triple(2))
                                 lassign(jaa) = .true.
                              elseif (atoma .eq. test0 
     &                                .and. atomc .eq. test2
     &                                .and. atomd .eq. test1 ) then
c                                 --- this fills the second slot in the 
c                                 --- lammps angle-angle
                                 call twh_itaa(GLB_GET,imolty,iunit,jaa
     &                                ,itaatemp)
                                 triple(2) = -itaatemp
                                 lassign(jaa) = .true.
                              elseif ( atomd .eq. test0 
     &                             .and. atoma .eq. test1
     &                             .and. atomc .eq. test2 ) then
c                                 --- this fills the third slot in the 
c                                 --- lammps angle-angle
                                 call twh_itaa(GLB_GET,imolty,iunit,jaa
     &                                ,triple(3))
                                 lassign(jaa) = .true.
                              elseif ( atomd .eq. test0 
     &                             .and. atoma .eq. test2
     &                             .and. atomc .eq. test1 ) then
                                 call twh_itaa(GLB_GET,imolty,iunit,jaa
     &                                ,itaatemp)
                                 triple(3) = -itaatemp
                                 lassign(jaa) = .true.
                              endif
                           endif
                        enddo
c                       --- safety check
                        if ( triple(1)*triple(2)*triple(3) .eq. 0 ) then
c                          --- something went wrong
                           write(6,*) 'WRITELAMMPS: could not map '
     &                          ,'angle-angle terms'
                           goto 666
                        endif
c                       --- get the type
                        type = code_angang(triple(1),triple(2)
     &                       ,triple(3)) + lam_improp_types
c                       --- output the atom numbers and type
                        count_improp = count_improp + 1
                        call twh_ijaa0(GLB_GET,imolty,iunit,iaa,ijaa0)
                        call twh_ijaa1(GLB_GET,imolty,iunit,iaa,ijaa1)
                        call twh_ijaa2(GLB_GET,imolty,iunit,iaa,ijaa2)
                        write(filenum,*) count_improp,type
     &                       ,twh_oldgetatomnum(ichain,ijaa1)
     &                       ,twh_oldgetatomnum(ichain,iunit)
     &                       ,twh_oldgetatomnum(ichain,ijaa0)
     &                       ,twh_oldgetatomnum(ichain,ijaa2)
c                       --- determine if we are done
                        ldone = .true.
                        do iaa = 1,inaa
                           if ( .not. lassign(iaa) ) ldone = .false.
                        enddo
                     enddo
                  enddo
               endif
            enddo
            write(filenum,*)
            if ( count_improp .ne. lam_improp + lam_angang ) then
               write(6,*) 'WRITELAMMPS: improper types did not add up'
            endif
         endif

 666     close(filenum)
      enddo
      close(inpnum)
c      write(6,*) 'done with writelammps'

      return
      end
