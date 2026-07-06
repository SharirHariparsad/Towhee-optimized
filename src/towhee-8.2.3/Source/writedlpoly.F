#include "preproc.h"
c     *****************************************************************
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
      subroutine twh_writedlpoly(lfailure)
c     ******************************************************************
c     * creates input files CONFIG and FIELD for DL_POLY, a molecular  *
c     * dynamics simulation code available from Daresbury labs, UK.    *
c     * The Dl_POLY code can be obtained at                            *
c     * http://www.cse.clrc.ac.uk/msi/software/DL_POLY                 *
c     * Note that user still needs to prepare the CONTROL file         *
c     *                                                                *
c     * originally written 08-05-2006 by O. Yazaydin                   *
c     * modified 12-17-2007 by C.M. Tenney                             *
c     *     angle types 2 (DREIDING) and 3 (harmonic cosine)           *
c     *     dihedral types 15 (DREIDING) and (not) 8 (nonbond)         *
c     *     improper form 4 type 7 (DREIDING 1-cos)                    *
c     *     fixed unit cell vector output and changed to triclinic pbc *
c     * last modified 03-31-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_LHERE
#define FUNCTION_OLDGETATOMNUM
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- local logicals
      logical ladd,lconstraints
c     --- local character strings
      character*30 filename
      character*(FFNAMELEN) nbnamei,nbnamej
c     --- local integer scalars
      integer dlp_bond,dlp_angle,dlp_dihedral,cpnumber
      integer nuni,testmolty,reset,itype,jtype,ij
      integer ibox,ivib,ibend,itor,type,ichain,imolty,iunit
      integer count_atom,count_bond,nchain,nunit
      integer style,inpnum,count_angle,count_dihedral,nmolty
      integer dlp_improper, count_improper, ncmt, icmt,numboxes
      integer ijben2test,ijben3test,ijimp2,ijimp3,ijimp4
      integer ijtor2test,ijtor3test,ijtor4test,ijvibtest,inbend
      integer inimproper,intorsion,invibration,invib_two,invib_thr
      integer testboxi,icoeff
c     --- local double precision
      double precision xtemp,ytemp,ztemp,radtodeg,charge,impcoeff
      double precision totbond,hm1,hm2,hm3,mass,bc_one,bc_zero
      double precision nbcoeff_13,nbcoeff_14
c     --- local double precision arrays
      double precision vibcoeff
      dimension vibcoeff(0:1)
      double precision torcoeff
      dimension torcoeff(0:3)
      
c     --- retrieve constants
      call twh_constant_radtodeg(GLB_GET,radtodeg)
      call twh_nchain(GLB_GET,nchain)
      call twh_numboxes(GLB_GET,numboxes)

      do ibox=1, numboxes
c      --- start of FIELD file
c      --- open FIELD
       inpnum = 23
       write (filename, '(a6,i1)') 'FIELD_',ibox
       call twh_openfile(lfailure,inpnum,filename,'DEFAULT','UNKNOWN'
     &      ,'DEFAULT')
       if ( lfailure ) return
       write(inpnum,'(a36)') 'DL_POLY FIELD file written by Towhee'
       write(inpnum,'(a8)') 'UNITS kj'
       call twh_nmolty(GLB_GET,nmolty)
       write(inpnum,'(a10,i2)') 'MOLECULES ', nmolty
       do imolty=1,nmolty
          call twh_nunit(GLB_GET,imolty,nunit)
          write(inpnum,'(a17,i2)') 'Name of Molecule ', imolty
          call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
          write(inpnum,'(a9,i6)') 'NUMMOLS  ',ncmt
          write(inpnum,'(a7,i6)') 'ATOMS  ', nunit
          do iunit = 1,nunit
             call twh_ntype(GLB_GET,imolty,iunit,type)
             call twh_qqatom(GLB_GET,imolty,iunit,charge)
             call twh_mass(GLB_GET,type,mass)
             call twh_nbname(GLB_GET,type,nbnamei)
             write(inpnum,'(a8,2f12.5,2a5)')
     &            nbnamei,mass,charge,'    1','    0'
          enddo
c         ---- number of bonds,angles,torsions
          dlp_bond = 0
          dlp_angle = 0
          dlp_dihedral = 0
          dlp_improper = 0
          do ichain = 1,nchain
             call twh_nboxi(GLB_GET,ichain,testboxi)
             if ( testboxi .eq. ibox ) then
                call twh_moltyp(GLB_GET,ichain,testmolty)
                if ( imolty .eq. testmolty ) then
                   do iunit = 1,nunit
c                     --- number of bonds
                      call twh_invib(GLB_GET,imolty,iunit,invibration)
                      do ivib = 1,invibration
                         call twh_ijvib(GLB_GET,imolty,iunit,ivib
     &                        ,ijvibtest)
                         if ( ijvibtest .gt. iunit ) then
                            dlp_bond = dlp_bond + 1
                         endif
                      enddo
c                     --- number of angles
                      call twh_inben(GLB_GET,imolty,iunit,inbend)
                      do ibend = 1,inbend
                         call twh_itben(GLB_GET,imolty,iunit,ibend,type)
                         ladd = .false.
c                        --- add based solely on ijben3 > iunit
                         call twh_ijben3(GLB_GET,imolty,iunit,ibend
     &                        ,ijben3test)
                         if ( ijben3test .gt. iunit ) ladd = .true.
                         if ( ladd ) then
                            dlp_angle = dlp_angle + 1
                         endif
                      enddo
c                     --- number of torsions
                      call twh_intor(GLB_GET,imolty,iunit,intorsion)
                      do itor = 1,intorsion
                         call twh_ittor(GLB_GET,imolty,iunit,itor
     &                        ,type)
                         if ( type .gt. 0 ) then
                            ladd = .false.
c                           --- any other case compare ijt4 > iunit
                            call twh_ijtor4(GLB_GET,imolty,iunit,itor
     &                           ,ijtor4test)
                            if ( ijtor4test .gt. iunit ) then
                               call twh_torstyle(GLB_GET,type,style)
                               if ( style .ne. 8) then
                                  ladd = .true.
                               endif
                            endif
                            if ( ladd ) then
                               dlp_dihedral = dlp_dihedral + 1
                            endif
                         endif
                      enddo
c                     --- end number of torsions

c                     --- number of impropers
                      call twh_inimprop(GLB_GET,imolty,iunit,inimproper)
                      do itor = 1,inimproper
                         call twh_itimprop(GLB_GET,imolty,iunit,itor
     &                        ,type)
                         if ( type .gt. 0 ) then
c c                         --- compare ijimprop4 > iunit
                            dlp_improper = dlp_improper + 1
                         endif
                      enddo
c                     --- end number of impropers
                   enddo
                   goto 111
                endif
             endif
          enddo
 111      continue
c         ---  output bond information
          if ( dlp_bond .ne. 0 ) then
             count_bond = 0
c            --- search for CONSTRAINTS
c            --- DL_POLY treats fixed bonds with SHAKE
c            --- and lists under a seperate section called CONSTRAINTS
             lconstraints = .false.
             do ichain = 1,nchain
                call twh_nboxi(GLB_GET,ichain,testboxi)
                if ( testboxi .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      call twh_oldatomnum(GLB_GET,ichain,1,reset)
                      reset = reset-1
                      do iunit = 1,nunit
                         call twh_invib(GLB_GET,imolty,iunit
     &                        ,invibration)
                         do ivib = 1,invibration
                            call twh_itvib(GLB_GET,imolty,iunit,ivib
     &                           ,type)
                            if (type .gt.0) then
                               call twh_ijvib(GLB_GET,imolty,iunit,ivib
     &                              ,ijvibtest)
                               if ( ijvibtest .gt. iunit ) then
                                  call twh_vibstyle(GLB_GET,type
     &                                 ,style)
                                  if ( style .eq. BOND_STYLE_FIXED) then
                                     count_bond = count_bond + 1
                                     if (count_bond .eq. 1) then
                                        lconstraints = .true.
                                        write(inpnum,'(a13,i5)')
     &                                       'CONSTRAINTS  ', dlp_bond
                                     endif
                                     call twh_vibcoeff(GLB_GET,type,0
     &                                    ,vibcoeff(0))
                                     write(inpnum,'(2i5,f12.5)')
     &                                    twh_oldgetatomnum(ichain
     &                                    ,iunit) - reset
     &                                    ,twh_oldgetatomnum(ichain
     &                                    ,ijvibtest) - reset
     &                                    ,vibcoeff(0)
                                  endif
                               endif
                            endif
                         enddo
                      enddo
                      goto 222
                   endif
                endif
             enddo
 222         continue
c            --- look for flexible bonds
c            --- subtract # of fixed bonds from flexible bonds
             if (count_bond .gt. 0) dlp_bond = dlp_bond - count_bond
c            --- if there are no flexible bonds continue with angles
             if (dlp_bond .eq. 0) goto 444
             write(inpnum,'(a8,i5)') 'BONDS   ', dlp_bond
             do ichain = 1,nchain
                call twh_nboxi(GLB_GET,ichain,testboxi)
                if ( testboxi .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      reset = twh_oldgetatomnum(ichain,1)-1
                      do iunit = 1,nunit
                         call twh_invib(GLB_GET,imolty,iunit
     &                        ,invibration)
                         do ivib = 1,invibration
                            call twh_itvib(GLB_GET,imolty,iunit,ivib
     &                           ,type)
                            if (type .gt.0) then
                               call twh_ijvib(GLB_GET,imolty,iunit,ivib
     &                              ,ijvibtest)
                               if ( ijvibtest .gt. iunit ) then
                                  count_bond = count_bond + 1
                                  call twh_vibstyle(GLB_GET,type
     &                                 ,style)
                                  if ( style .eq. BOND_STYLE_HARMONIC
     &                                 ) then
                                     call twh_vibcoeff(GLB_GET,type,0
     &                                    ,vibcoeff(0))
                                     call twh_vibcoeff(GLB_GET,type,1
     &                                    ,vibcoeff(1))
                                     write(inpnum,'(a5,2i5,2f12.5)')
     &                                    'harm '
     &                                    ,twh_oldgetatomnum(ichain
     &                                    ,iunit) - reset
     &                                    ,twh_oldgetatomnum(ichain
     &                                    ,ijvibtest) - reset
     &                                    ,vibcoeff(1)
     &                                    *2.0d0/CON_KJTOK
     &                                    ,vibcoeff(0)
                                  else
                                     if ( .not. lconstraints ) then
                                        write(6,*) 'WRITEDLPOLY: ',
     &                                       'invalid bond style '
     &                                       ,style,' has nothing '
     &                                       ,'similar in dlpoly'
                                     endif
                                  endif
                               endif
                            endif
                         enddo
                      enddo
                      goto 333
                   endif
                endif
             enddo
 333         continue
             if ( dlp_bond .ne. count_bond ) then
                write(6,*) 'Error in writedlpoly'
                write(6,*) 'dlp_bond ',dlp_bond
                write(6,*) 'count_bond ',count_bond
             endif
          endif
 444      continue
c         ----  output angle information
          write(inpnum,'(a8,i5)') 'ANGLES   ', dlp_angle
          if ( dlp_angle .ne. 0 ) then
             count_angle = 0
             do ichain = 1,nchain
                call twh_nboxi(GLB_GET,ichain,testboxi)
                if ( testboxi .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      reset = twh_oldgetatomnum(ichain,1)-1
                      do iunit = 1,nunit
                         call twh_inben(GLB_GET,imolty,iunit,inbend)
                         do ibend = 1,inbend
                            call twh_itben(GLB_GET,imolty,iunit,ibend
     &                           ,type)
                            call twh_ijben3(GLB_GET,imolty,iunit,ibend
     &                           ,ijben3test)
                            if ( ijben3test .gt.iunit ) then
                               count_angle = count_angle + 1
                               call twh_ijben2(GLB_GET,imolty,iunit
     &                              ,ibend,ijben2test)
                               call twh_ijben3(GLB_GET,imolty,iunit
     &                              ,ibend,ijben3test)
                               call twh_anglestyle(GLB_GET,type,style)
                               if ( style .eq. ANG_STYLE_HARMONIC ) then
                                  call twh_bencoeff(GLB_GET,type,0
     &                                 ,bc_zero)
                                  call twh_bencoeff(GLB_GET,type,1
     &                                 ,bc_one)
                                  write(inpnum,'(a5,3i5,2f12.5)')
     &                                 'harm '
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,iunit) - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,ijben2test) - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,ijben3test) - reset
     &                                 ,bc_one*2.0d0/CON_KJTOK
     &                                 ,bc_zero*radtodeg
                               elseif ( style .eq. ANG_STYLE_DREIDING
     &                                 ) then
                                  call twh_bencoeff(GLB_GET,type,1
     &                                 ,bc_one)
                                  write(inpnum,'(a5,3i5,3f12.5)')
     &                                 'cos  '
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,iunit) - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,ijben2test) - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,ijben3test)- reset
     &                                 ,bc_one*2.0d0/CON_KJTOK
     &                                 , 0.0, 1.0
                               elseif ( style .eq. ANG_STYLE_HARM_COS
     &                                 ) then
                                  call twh_bencoeff(GLB_GET,type,0
     &                                 ,bc_zero)
                                  call twh_bencoeff(GLB_GET,type,1
     &                                 ,bc_one)
                                  write(inpnum,'(a5,3i5,2f12.5)')
     &                                 'hcos '
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,iunit) - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,ijben2test) - reset
     &                                 ,twh_oldgetatomnum(ichain
     &                                 ,ijben3test) - reset
     &                                 ,bc_one*2.0d0/CON_KJTOK
     &                                 ,bc_zero*radtodeg
                               else
                                  write(6,*) 'WRITEDLPOLY: ',
     &                                 'invalid angle style ',style, 
     &                                 ' has nothing similar '
     &                                 ,'in dlpoly'
                               endif
                            endif
                         enddo
                      enddo
                      goto 555
                   endif
                endif
             end do
 555         continue
             if ( dlp_angle .ne. count_angle ) then
                write(6,*) 'problem with angles in writelammps'
                write(6,*) 'dlp_angle ',dlp_angle
                write(6,*) 'count_angle ',count_angle
             endif
          endif
c         ---- output torsion information
          write(inpnum,'(a11,i5)') 'DIHEDRALS  ', dlp_dihedral
          if ( dlp_dihedral .ne. 0 ) then
             count_dihedral = 0
             do ichain = 1,nchain
                call twh_nboxi(GLB_GET,ichain,testboxi)
                if ( testboxi .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      reset = twh_oldgetatomnum(ichain,1)-1
                      do iunit = 1,nunit
                         call twh_intor(GLB_GET,imolty,iunit,intorsion)
                         do itor = 1,intorsion
                            call twh_ittor(GLB_GET,imolty,iunit,itor
     &                           ,type)
                            if ( type .gt. 0 ) then
                               call twh_ijtor4(GLB_GET,imolty,iunit
     &                              ,itor,ijtor4test)
                               if ( ijtor4test .gt. iunit ) then
c                                 --- ignore style 8 as it has
c                                 --- (one-four interactions only)
                                  call twh_torstyle(GLB_GET,type,style)
                                  if ( style .ne. 8) then
                                     count_dihedral = count_dihedral
     &                                    + 1
                                     call twh_ijtor2(GLB_GET,imolty
     &                                    ,iunit,itor,ijtor2test)
                                     call twh_ijtor3(GLB_GET,imolty
     &                                    ,iunit,itor,ijtor3test)
                                     do icoeff = 0,3
                                        call twh_torcoeff(GLB_GET,type
     &                                       ,icoeff,torcoeff(icoeff))
                                     enddo
                                     if ( style .eq. 2) then
                                        write(inpnum
     &                                       ,'(a5,4i5,3f12.5,2a8)')
     &                                       'cos3 '
     &                                       ,twh_oldgetatomnum(
     &                                       ichain,iunit) - reset
     &                                       ,twh_oldgetatomnum(
     &                                       ichain,ijtor2test) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,ijtor3test) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,ijtor4test) - reset
     &                                       ,torcoeff(1)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,torcoeff(2)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,torcoeff(3)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,'     0.5', '     0.5'
                                     else if ( style .eq. 11) then
                                        write(inpnum
     &                                       ,'(a5,4i5,3f12.5,2a8
     &,f12.5)') 'opls '
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,iunit) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,ijtor2test) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,ijtor3test) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,ijtor4test) - reset
     &                                       ,torcoeff(0)
     &                                       /CON_KJTOK
     &                                       ,torcoeff(1)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,torcoeff(2)
     &                                       *2.0d0/CON_KJTOK
     &                                       ,'     0.5', '     0.5'
     &                                       ,torcoeff(3)
     &                                       *2.0d0/CON_KJTOK
                                     else if ( style .eq. 15) then
                                        call twh_invib(GLB_GET,imolty
     &                                       ,ijtor2test,invib_two)
                                        call twh_invib(GLB_GET,imolty
     &                                       ,ijtor3test,invib_thr)
                                        totbond = dble((invib_two-1)
     &                                       *(invib_thr-1))
                                        write(inpnum,'(a5,4i5,3f12.5,2a8
     &,f12.5)')
     &                                       'cos '
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,iunit) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,ijtor2test) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,ijtor3test) - reset
     &                                       ,twh_oldgetatomnum(ichain
     &                                       ,ijtor4test) - reset
     &                                       ,torcoeff(1)
     &                                       /CON_KJTOK/totbond
     &                                       ,torcoeff(2)
     &                                       *torcoeff(3)
     &                                       *radtodeg-180.0d0
     &                                       ,torcoeff(2)
     &                                       ,'     1.0', '     1.0'
                                     else
                                        write(6,*) 'WRITEDLPOLY: ',
     &                                       'invalid torsion style '
     &                                       ,style,
     &                                       ' has nothing similar'
     &                                       ,' in dlpoly'
                                     endif
                                  endif
                               endif
                            endif
                         enddo
                      enddo
                      goto 666
                   endif
                endif
             enddo
 666         continue
             if ( dlp_dihedral .ne. count_dihedral ) then
                write(6,*) 'problem with angles in writedlpoly'
                write(6,*) 'dlp_dihedral ',dlp_dihedral
                write(6,*) 'count_dihedral ',count_dihedral
             endif
          endif

c         ---- output improper torsion (inversion) information
          write(inpnum,'(a11,i5)') 'INVERSIONS  ', dlp_improper
          if ( dlp_improper .ne. 0 ) then
             count_improper = 0
             do ichain = 1,nchain
                call twh_nboxi(GLB_GET,ichain,testboxi)
                if ( testboxi .eq. ibox ) then
                   call twh_moltyp(GLB_GET,ichain,imolty)
                   if (imolty .eq. testmolty ) then
                      reset = twh_oldgetatomnum(ichain,1)-1
                      do iunit = 1,nunit
                         call twh_inimprop(GLB_GET,imolty,iunit
     &                        ,inimproper)
                         do itor = 1,inimproper
                            call twh_itimprop(GLB_GET,imolty,iunit,itor
     &                           ,type)
                            if ( type .gt. 0 ) then
                               count_improper = count_improper + 1
                               call twh_impstyle(GLB_GET,type,style)
                               if (style .eq. 7) then
                                  call twh_ijimprop2(GLB_GET,imolty
     &                                 ,iunit,itor,ijimp2)
                                  call twh_ijimprop3(GLB_GET,imolty
     &                                 ,iunit,itor,ijimp3)
                                  call twh_ijimprop4(GLB_GET,imolty
     &                                 ,iunit,itor,ijimp4)
                                  call twh_impcoeff(GLB_GET,type,1
     &                                 ,impcoeff)
                                  write(inpnum,'(a5,4i5,f12.5)')
     &                                 'plan '
     &                                 ,twh_oldgetatomnum(ichain,iunit)
     &                                 - reset
     &                                 ,twh_oldgetatomnum(ichain,ijimp2)
     &                                 - reset
     &                                 ,twh_oldgetatomnum(ichain,ijimp3)
     &                                 - reset
     &                                 ,twh_oldgetatomnum(ichain,ijimp4)
     &                                 - reset
     &                                 ,impcoeff/CON_KJTOK
                               else
                                  write(6,*) 'WRITEDLPOLY: ',
     &                                 'invalid improper style ',style,
     &                                 ' has nothing similar in dlpoly'
                               endif
                            endif
                         enddo
                      enddo
                      goto 777
                   endif
                endif
             enddo
 777         continue
             if ( dlp_improper .ne. count_improper ) then
                write(6,*) 'problem with angles in writedlpoly'
                write(6,*) 'dlp_improper ',dlp_improper
                write(6,*) 'count_improper ',count_improper
             endif
          endif

          write(inpnum,'(a6)') 'FINISH'
       enddo
c      --- ouptput nonbonded parameters
c      --- get the classical potential
       call twh_cp_number(GLB_GET,cpnumber)
c      --- determine the number of unique nonbonded interactions (nuni)
       nuni=0
       do itype = 1,NNTYPE
          do jtype = itype,NNTYPE
             if ( twh_get_lhere(itype)
     &            .and. twh_get_lhere(jtype) ) then
                nuni=nuni+1
             endif
          enddo
       enddo
       write(inpnum,'(a6,i5)') 'VDW   ', nuni
       if ( cpnumber .eq. CP_LENNARD_JONES ) then
          do itype = 1,NNTYPE
             do jtype = itype,NNTYPE
                if ( twh_get_lhere(itype)
     &               .and. twh_get_lhere(jtype) ) then
                   ij = twh_getnbtype(itype,jtype)
                   call twh_nbname(GLB_GET,itype,nbnamei)
                   call twh_nbname(GLB_GET,jtype,nbnamej)
                   call twh_nbcoeff(GLB_GET,13,ij,nbcoeff_13)
                   call twh_nbcoeff(GLB_GET,14,ij,nbcoeff_14)
                   write(inpnum,'(2a8,1x,a4,2f12.5)')
     &                  nbnamei,nbnamej,'lj  '
     &                  ,nbcoeff_14/(4.0d0*CON_KJTOK)
     &                  ,(nbcoeff_13)**(1.0/6.0)
                endif
             enddo
          enddo
       end if
       write(inpnum,'(a5)') 'CLOSE'
       close (inpnum)
      enddo
c     --- end of FIELD file
c     --- start of CONFIG file
      do ibox=1,numboxes
         inpnum = 33
         write (filename, '(a7,i1)') 'CONFIG_',ibox
c        --- filename = 'CONFIG'
         call twh_openfile(lfailure,inpnum,filename,'DEFAULT','UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) return
         write(inpnum,'(a37)') 'DL_POLY CONFIG file written by Towhee'
         write(inpnum,'(2a10)') '         0','         3'
c        --- output box dimensions
         call twh_hmatrix(GLB_GET,BXD_REAL,ibox,1,1,hm1)
         call twh_hmatrix(GLB_GET,BXD_REAL,ibox,1,2,hm2)
         call twh_hmatrix(GLB_GET,BXD_REAL,ibox,1,3,hm3)
         write(inpnum,'(3f20.12)') hm1,hm2,hm3
         call twh_hmatrix(GLB_GET,BXD_REAL,ibox,2,1,hm1)
         call twh_hmatrix(GLB_GET,BXD_REAL,ibox,2,2,hm2)
         call twh_hmatrix(GLB_GET,BXD_REAL,ibox,2,3,hm3)
         write(inpnum,'(3f20.12)') hm1,hm2,hm3
         call twh_hmatrix(GLB_GET,BXD_REAL,ibox,3,1,hm1)
         call twh_hmatrix(GLB_GET,BXD_REAL,ibox,3,2,hm2)
         call twh_hmatrix(GLB_GET,BXD_REAL,ibox,3,3,hm3)
         write(inpnum,'(3f20.12)') hm1,hm2,hm3
c        --- output coordinates
         count_atom = 0
         do imolty = 1, nmolty
            call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
            call twh_nunit(GLB_GET,imolty,nunit)
            do icmt = 1, ncmt
               call twh_chainlist(GLB_GET,icmt,ibox,imolty,ichain)
               do iunit = 1,nunit
                  count_atom = count_atom + 1
                  call twh_ntype(GLB_GET,imolty,iunit,type)
                  call twh_coordinates(GLB_GET,CRD_REAL,ichain
     &                 ,iunit,xtemp,ytemp,ztemp)
c                 --- default format
                  call twh_nbname(GLB_GET,type,nbnamei)
                  write(inpnum,'(a8,i10)') nbnamei,count_atom
                  write(inpnum,'(3f20.12)') xtemp,ytemp,ztemp
               enddo
            enddo
         enddo
      enddo
c     --- end of CONFIG file

      return
      end
