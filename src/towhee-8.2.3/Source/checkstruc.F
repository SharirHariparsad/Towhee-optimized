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
      subroutine twh_checkstruc(lfailure)
c     ******************************************************************
c     * this subroutine runs through the input file structure and makes*
c     * sure that all of the bonds, bends, and torsions are set up in a*
c     * consistent way.                                                *
c     * it then sets up some arrays of information that is useful when *
c     * computing bending angles, and torsions                         *
c     *                                                                *
c     * Originally written 12-15-2000 by M.G. Martin                   *
c     * Last modified 06-03-2021 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
#define DEBUG_CHECKSTRUC 0
c     --- variables passed to/from the subroutine
      logical lfailure
c     --- local variables
c     --- local logicals
      logical lfound,lbnbn,lintraterm,loutdft,lonefourtor,lbonbon
c     --- local characters
      character*(BONPATLEN) bondpatt
c     --- integer scalars
      integer imolty,ivib,jvib,iben,jben,itor,jtor,iunit,iaa,jaa,itype
      integer iimp,jimp,itest,itestb,itestc,junit,count,iof,jof,itemp
      integer aunit,bunit,cunit,kunit,atype,ii,kk,ktype,matom,type
      integer iorder,nmolty,nunit,nstep,numtor
      integer itbentest,jtbentest,anglestyle,cpnumber
      integer ijaa0i,ijaa0j,ijaa1i,ijaa1j,ijaa2i,ijaa2j
      integer ijben2test,ijben3test,ijben2j,ijben3j
      integer ijimp2test,ijimp3test,ijimp4test,ijimp2j,ijimp3j,ijimp4j
      integer ijof5test,ijvibtest,ijvibj
      integer ijtor2i,ijtor2j,ijtor3i,ijtor3j,ijtor4i,ijtor4j
      integer inaaunit,invibunit,invibtest
      integer inbend,inbenda,inbendb,inbendt,inimproper,itimproper
      integer inonefivei,inonefivej,intorsion,intorsionc,itaatest
      integer ittortest,ittori,ittorj,itvibi,itvibj
c     --- local double precision arrays
      integer torcross
      dimension torcross(NUMAX)
c     --- local dp scalars
      double precision rij,rjk,rik,zi,zk,theta,dvalue,thetatwo

      write(6,*) 'Verifying input structures are consistent'

c     --- initialize flag to tell if we have any bond-bond interactions
      lbnbn = .false.
c     --- determine whether we need to stop if we observe a non-valid
c     --- intrmolecular setup.  we might not want to stop if all we are
c     --- doing is converting from lammps into tramonto and not doing
c     --- any towhee simulation steps
      call twh_nstep(GLB_GET,nstep)
      call twh_loutdft(GLB_GET,loutdft)
      if ( loutdft .and. nstep .eq. 0 ) then
         lintraterm = .false.
      else
         lintraterm = .true.
      endif

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1, nmolty
c        --- check that each atom has self-consistent and valid values
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1, nunit
            call twh_invib(GLB_GET,imolty,iunit,invibunit)
c           --- check vibrations
            do ivib = 1, invibunit
c              --- check that vibration atoms are between 1 and nunit
               call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
               if ( ijvibtest .lt. 1 .or. ijvibtest .gt. nunit )then
                  write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                 ,iunit,ivib
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,ijvibtest
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that vibration type is between 1 and tvibmax
               call twh_itvib(GLB_GET,imolty,iunit,ivib,itvibi)
               if ( itvibi .lt. 1 .or. itvibi .gt. TVIBMAX )then
                  write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                 ,iunit,ivib
                  write(6,*) 'Vibration type is nonvalid value ',itvibi
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
c           --- check for multiple vibrations across the same atoms
            do ivib = 1, invibunit-1
               call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
               do jvib = ivib+1, invibunit
                  call twh_ijvib(GLB_GET,imolty,iunit,jvib,ijvibj)
                  if ( ijvibtest .eq. ijvibj ) then
                     write(6,*) 'Problem: imolty,iunit,ivib,jvib',
     &                    imolty,iunit,ivib,jvib
                     write(6,*) 'Are both connected to atom '
     &                    ,ijvibtest
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo
c           --- check that each pair of atoms that is bonded to this
c           --- atom also have a bending interaction
            do ivib = 1, invibunit-1
               call twh_ijvib(GLB_GET,imolty,iunit,ivib,aunit)
               do jvib = ivib+1,invibunit
                  call twh_ijvib(GLB_GET,imolty,iunit,jvib,bunit)
                  lfound = .false.
                  call twh_inben(GLB_GET,imolty,aunit,inbenda)
                  do iben = 1,inbenda
                     call twh_ijben2(GLB_GET,imolty,aunit,iben
     &                    ,ijben2test)
                     call twh_ijben3(GLB_GET,imolty,aunit,iben
     &                    ,ijben3test)
                     if ( ijben2test .eq. iunit 
     &                    .and. ijben3test .eq. bunit ) then
                        lfound = .true.
                     endif
                  enddo
                  if ( .not. lfound ) then
                     write(6,*) 'imolty,iunit,ivib,jvib'
     &                    ,imolty,iunit,ivib,jvib
                     write(6,*) 'atoms ',aunit,bunit,' are connected to'
     &                    ,' iunit, but have no bend'
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo

c           --- check bending
            call twh_inben(GLB_GET,imolty,iunit,inbend)
            do iben = 1, inbend
c              --- check that bending atoms are between 1 and nunit
               call twh_ijben2(GLB_GET,imolty,iunit,iben,ijben2test)
               call twh_ijben3(GLB_GET,imolty,iunit,iben,ijben3test)
               if ( ijben2test .lt. 1 .or. ijben2test .gt. nunit .or.
     &              ijben3test .lt. 1 .or. ijben3test .gt. nunit )then
                  write(6,*) 'Problem: imolty,iunit,iben',imolty
     &                 ,iunit,iben
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,ijben2test,ijben3test
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that bending type is between 1 and tbenmax
               call twh_itben(GLB_GET,imolty,iunit,iben,itbentest)
               if ( itbentest .lt. 1 .or. itbentest .gt. TBENMAX )then
                  write(6,*) 'Problem: imolty,iunit,iben',imolty
     &                 ,iunit,iben
                  write(6,*) 'Bending type is nonvalid value '
     &                 ,itbentest
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
c           --- check for multiple bends across the same atoms
            do iben = 1, inbend-1
               call twh_ijben2(GLB_GET,imolty,iunit,iben,ijben2test)
               call twh_ijben3(GLB_GET,imolty,iunit,iben,ijben3test)
               do jben = iben+1, inbend
                  call twh_ijben2(GLB_GET,imolty,iunit,jben,ijben2j)
                  call twh_ijben3(GLB_GET,imolty,iunit,jben,ijben3j)
                  if ( ijben2test .eq. ijben2j .and. 
     &                 ijben3test .eq. ijben3j ) then
                     write(6,*) 'Problem: imolty,iunit,iben,jben',
     &                    imolty,iunit,iben,jben
                     write(6,*) 'Are both connected to atoms '
     &                    ,ijben2test,ijben3test
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo
c           --- check that all bends which are also bonded to another
c           --- atom have a torsion interaction.
            do iben = 1, inbend
               call twh_ijben3(GLB_GET,imolty,iunit,iben,aunit)
               call twh_ijben2(GLB_GET,imolty,iunit,iben,bunit)
               do ivib = 1, invibunit
                  call twh_ijvib(GLB_GET,imolty,iunit,ivib,cunit)
                  if ( cunit .ne. bunit .and. cunit .ne. aunit ) then
                     lfound = .false.
                     call twh_intor(GLB_GET,imolty,cunit,numtor)
                     do itor = 1,numtor
                        call twh_ijtor2(GLB_GET,imolty,cunit,itor
     &                       ,ijtor2i)
                        call twh_ijtor3(GLB_GET,imolty,cunit,itor
     &                       ,ijtor3i)
                        call twh_ijtor4(GLB_GET,imolty,cunit,itor
     &                       ,ijtor4i)
                        if ( ijtor2i .eq. iunit
     &                       .and. ijtor3i .eq. bunit
     &                       .and. ijtor4i .eq. aunit ) then
                           lfound = .true.
                        endif
                     enddo
                     if ( .not. lfound ) then
                        write(6,*) 'imolty,iunit,iben,ivib',imolty,iunit
     &                       ,iben,ivib
                        write(6,*) 'bonded to an atom and an angle that'
     &                       ,' do not form a torsion'
                        if ( lintraterm ) then
                           lfailure = .true.
                           return
                        endif
                     endif
                  endif
               enddo
            enddo

c           --- check torsions
            call twh_intor(GLB_GET,imolty,iunit,numtor)
            do itor = 1, numtor
               call twh_ijtor2(GLB_GET,imolty,iunit,itor,ijtor2i)
               call twh_ijtor3(GLB_GET,imolty,iunit,itor,ijtor3i)
               call twh_ijtor4(GLB_GET,imolty,iunit,itor,ijtor4i)
c              --- check that torsion atoms are between 1 and nunit
               if ( ijtor2i .lt. 1 .or. ijtor2i .gt. nunit
     &              .or. ijtor3i .lt. 1 .or. ijtor3i .gt. nunit
     &              .or. ijtor4i .lt. 1 .or. ijtor4i .gt. nunit )then
                  write(6,*) 'Problem: imolty,iunit,itor',imolty
     &                 ,iunit,itor
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,ijtor2i,ijtor3i,ijtor4i
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif

c              --- check that torsion type is between 1 and ttormax
               call twh_ittor(GLB_GET,imolty,iunit,itor,ittortest)
               if ( abs(ittortest) .eq. 0 .or.
     &              abs(ittortest) .gt. TTORMAX )then
                  write(6,*) 'Problem: imolty,iunit,itor',imolty
     &                 ,iunit,itor
                  write(6,*) 'Torsion type is nonvalid value '
     &                 ,ittortest
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo

c           --- check for multiple torsions across the same atoms
            call twh_intor(GLB_GET,imolty,iunit,numtor)
            do itor = 1, numtor-1
               call twh_ijtor2(GLB_GET,imolty,iunit,itor,ijtor2i)
               call twh_ijtor3(GLB_GET,imolty,iunit,itor,ijtor3i)
               call twh_ijtor4(GLB_GET,imolty,iunit,itor,ijtor4i)
               do jtor = itor+1,numtor
                  call twh_ijtor2(GLB_GET,imolty,iunit,jtor,ijtor2j)
                  call twh_ijtor3(GLB_GET,imolty,iunit,jtor,ijtor3j)
                  call twh_ijtor4(GLB_GET,imolty,iunit,jtor,ijtor4j)
                  if ( ijtor2i .eq. ijtor2j .and. 
     &                 ijtor3i .eq. ijtor3j .and.
     &                 ijtor4i .eq. ijtor4j ) then
                     write(6,*) 'Problem: imolty,iunit,itor,jtor',
     &                    imolty,iunit,itor,jtor
                     write(6,*) 'Are both connected to atoms '
     &                    ,ijtor2i,ijtor3i,ijtor4i
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo

c           --- check angle-angle terms
            call twh_inaa(GLB_GET,imolty,iunit,inaaunit)
            do iaa = 1, inaaunit
c              --- check that ang-ang atoms are between 1 and nunit
               call twh_ijaa0(GLB_GET,imolty,iunit,iaa,ijaa0i)
               call twh_ijaa1(GLB_GET,imolty,iunit,iaa,ijaa1i)
               call twh_ijaa2(GLB_GET,imolty,iunit,iaa,ijaa2i)
               if ( ijaa0i .lt. 1 .or. ijaa0i .gt. nunit
     &              .or. ijaa1i .lt. 1 .or. ijaa1i .gt. nunit
     &              .or. ijaa2i .lt. 1 .or. ijaa2i .gt. nunit
     &              )then
                  write(6,*) 'Problem: imolty,iunit,iaa',imolty
     &                 ,iunit,iaa
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,ijaa0i,ijaa1i,ijaa2i
                  if ( lintraterm )  then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that ang-ang type is between 1 and taamax
               call twh_itaa(GLB_GET,imolty,iunit,iaa,itaatest)
               if ( itaatest .lt. 1 .or. itaatest .gt. TAAMAX )then
                  write(6,*) 'Problem: imolty,iunit,iaa',imolty
     &                 ,iunit,iaa
                  write(6,*) 'Ang-ang type is nonvalid value '
     &                 ,itaatest
                  if ( lintraterm )  then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
c           --- check for multiple ang-ang across the same atoms
            do iaa = 1, inaaunit-1
               do jaa = iaa+1, inaaunit
                  call twh_ijaa0(GLB_GET,imolty,iunit,iaa,ijaa0i)
                  call twh_ijaa0(GLB_GET,imolty,iunit,jaa,ijaa0j)
                  call twh_ijaa1(GLB_GET,imolty,iunit,iaa,ijaa1i)
                  call twh_ijaa1(GLB_GET,imolty,iunit,jaa,ijaa1j)
                  call twh_ijaa2(GLB_GET,imolty,iunit,iaa,ijaa2i)
                  call twh_ijaa2(GLB_GET,imolty,iunit,jaa,ijaa2j)
                  if (  ijaa0i .eq. ijaa0j .and. 
     &                 ijaa1i .eq. ijaa1j .and.
     &                 ijaa2i .eq. ijaa2j ) then
                     write(6,*) 'Problem: imolty,iunit,iaa,jaa',
     &                    imolty,iunit,iaa,jaa
                     write(6,*) 'Are both connected to atoms '
     &                    ,ijaa0i,ijaa1i,ijaa2i
                     if ( lintraterm )  then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo

c           --- check improper torsions
            call twh_inimprop(GLB_GET,imolty,iunit,inimproper)
            do iimp = 1,inimproper 
c              --- check that atoms are between 1 and nunit
               call twh_ijimprop2(GLB_GET,imolty,iunit,iimp,ijimp2test)
               call twh_ijimprop3(GLB_GET,imolty,iunit,iimp,ijimp3test)
               call twh_ijimprop4(GLB_GET,imolty,iunit,iimp,ijimp4test)
               if ( ijimp2test .lt. 1 .or. ijimp2test .gt. nunit .or.
     &              ijimp3test .lt. 1 .or. ijimp3test .gt. nunit .or.
     &              ijimp4test .lt. 1 .or. ijimp4test .gt. nunit) then
                  write(6,*) 'Problem: imolty,iunit,iimp',imolty
     &                 ,iunit,iimp
                  write(6,*) 'Atom is connected to nonvalid atom '
     &                 ,ijimp2test,ijimp3test,ijimp4test
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that imp torsion type is between 1 and timpmax
               call twh_itimprop(GLB_GET,imolty,iunit,iimp,itimproper)
               if ( itimproper .lt. 1 .or. itimproper .gt. TIMPMAX )then
                  write(6,*) 'Problem: imolty,iunit,iimp',imolty
     &                 ,iunit,iimp
                  write(6,*) 'Improper torsion type is nonvalid value '
     &                 ,itimproper
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
c           --- check for multiple imp torsions across the same atoms
            call twh_inimprop(GLB_GET,imolty,iunit,inimproper)
            do iimp = 1,inimproper-1
               call twh_ijimprop2(GLB_GET,imolty,iunit,iimp,ijimp2test)
               call twh_ijimprop3(GLB_GET,imolty,iunit,iimp,ijimp3test)
               call twh_ijimprop4(GLB_GET,imolty,iunit,iimp,ijimp4test)
               do jimp = iimp+1,inimproper
                  call twh_ijimprop2(GLB_GET,imolty,iunit,jimp,ijimp2j)
                  call twh_ijimprop3(GLB_GET,imolty,iunit,jimp,ijimp3j)
                  call twh_ijimprop4(GLB_GET,imolty,iunit,jimp,ijimp4j)
                  if ( ijimp2test .eq. ijimp2j .and. 
     &                 ijimp3test .eq. ijimp3j .and.
     &                 ijimp4test .eq. ijimp4j ) then
                     write(6,*) 'Problem: imolty,iunit,iimp,jimp',
     &                    imolty,iunit,iimp,jimp
                     write(6,*) 'Are both connected to atoms '
     &                    ,ijimp2test,ijimp3test,ijimp4test
                     if ( lintraterm ) then
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo
c           --- check one-five interactions
            call twh_inof(GLB_GET,imolty,iunit,inonefivei)
            do iof = 1,inonefivei
               call twh_ijof5(GLB_GET,imolty,iunit,iof,matom)
               if ( matom .lt. 1 .or. matom .gt. nunit ) then
c                 --- ijof5 out of bounds
                  write(6,*) 'CHECKSTRUC: atom out of allowed range'
                  write(6,*) 'imolty,iunit,ijof5',imolty,iunit,matom
                  lfailure = .true.
                  return
               endif
               call twh_itof(GLB_GET,imolty,iunit,iof,type)
               if ( type .lt. 1 .or. type .gt. TOFMAX ) then
c                 --- itof out of bound
                  write(6,*) 'CHECKSTRUC: type out of allowed range'
                  write(6,*) 'imolty,iunit,itof',imolty,iunit,type
                  lfailure = .true.
                  return
               endif
c              --- check multiples
               call twh_inof(GLB_GET,imolty,iunit,inonefivej)
               do jof = 1,inonefivej
                  if ( iof .ne. jof ) then
                     call twh_ijof5(GLB_GET,imolty,iunit,jof,ijof5test)
                     if ( ijof5test .eq. matom ) then
                        write(6,*) 'CHECKSTRUC: multiple identical 1-5'
                        write(6,*) 'imolty,iunit,iof,ijof5',imolty,iunit
     &                       ,iof,matom
                        write(6,*) 'imolty,iunit,jof,ijof5',imolty,iunit
     &                       ,jof,ijof5test
                        lfailure = .true.
                        return
                     endif
                  endif
               enddo
            enddo
         enddo

c        --- check the terms are consistent from one atom to the next
c        --- for the terms that should be symmetric
         do iunit = 1,nunit
c           --- check vibrations
            call twh_invib(GLB_GET,imolty,iunit,invibunit)
            do ivib = 1, invibunit
c              --- check that atom is not connected to itself
               call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
               if ( ijvibtest .eq. iunit )then
                  write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                 ,iunit,ivib
                  write(6,*) 'Atom is connected to itself ',ijvibtest
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that each other atom also contains this 
c              --- interaction
               call twh_ijvib(GLB_GET,imolty,iunit,ivib,itest)
               lfound = .false.
               call twh_invib(GLB_GET,imolty,itest,invibtest)
               do jvib = 1, invibtest
                  call twh_ijvib(GLB_GET,imolty,itest,jvib,ijvibj)
                  if ( ijvibj .eq. iunit ) then
                     call twh_itvib(GLB_GET,imolty,itest,jvib,itvibi)
                     call twh_itvib(GLB_GET,imolty,iunit,ivib,itvibj)
                     if ( itvibi .ne. itvibj ) then
                        write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                       ,iunit,ivib
                        write(6,*) 'Has nonsymmetric type: itest,jvib'
     &                       ,itest,jvib
                        write(6,*) 'iunit ivib type ',itvibi
                        write(6,*) 'itest jvib type ',itvibj
                        if ( lintraterm ) then
                           lfailure = .true.
                           return
                        endif
                     endif
                     lfound = .true.
                  endif
               enddo
               if ( .not. lfound ) then 
                  write(6,*) 'Problem: imolty,iunit,ivib',imolty
     &                 ,iunit,ivib
                  write(6,*) 'There is no symmetric vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
         enddo

         do iunit = 1,nunit
c           --- check bends
            call twh_inben(GLB_GET,imolty,iunit,inbend)
            do iben = 1, inbend
               call twh_ijben2(GLB_GET,imolty,iunit,iben,itest)
               call twh_ijben3(GLB_GET,imolty,iunit,iben,itestb)
c              --- check that atom #2 is bonded to this atom
               lfound = .false.
               call twh_invib(GLB_GET,imolty,iunit,invibunit)
               do ivib = 1,invibunit
                  call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
                  if ( ijvibtest .eq. itest ) then
                     lfound = .true.
c                    --- store bond type in this angle
                     call twh_itvib(GLB_GET,imolty,iunit,ivib,itemp)
                     call twh_ibtben1(GLB_SET,imolty,iunit,iben,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijben2,ijben3,iben'
     &                 ,imolty,iunit,itest,itestb,iben
                  write(6,*) 'Atoms 1-2 connected in a bend are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #3 is bonded to atom #2
               lfound = .false.
               call twh_invib(GLB_GET,imolty,itest,invibtest)
               do ivib = 1,invibtest
                  call twh_ijvib(GLB_GET,imolty,itest,ivib,ijvibtest)
                  if ( ijvibtest .eq. itestb ) then
                     lfound = .true.
c                    --- store bond type in this angle
                     call twh_itvib(GLB_GET,imolty,itest,ivib,itemp)
                     call twh_ibtben2(GLB_SET,imolty,iunit,iben,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijben2,ijben3,iben'
     &                 ,imolty,iunit,itest,itestb,iben
                  write(6,*) 'Atoms 2-3 connected in a bend are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif

c              --- check that each other atom also contains this 
c              --- interaction
               lfound = .false.
               call twh_inben(GLB_GET,imolty,itestb,inbendb)
               do jben = 1, inbendb
                  call twh_ijben2(GLB_GET,imolty,itestb,jben,ijben2test)
                  call twh_ijben3(GLB_GET,imolty,itestb,jben,ijben3test)
                  if ( ijben3test .eq. iunit .and. 
     &                 ijben2test .eq. itest ) then
                     call twh_itben(GLB_GET,imolty,iunit,iben,itbentest)
                     call twh_itben(GLB_GET,imolty,itestb,jben
     &                    ,jtbentest)
                     if ( jtbentest .ne. itbentest ) then
                        write(6,*) 'Problem: imolty,iunit,iben',imolty
     &                       ,iunit,iben
                        write(6,*) 'Has nonsymmetric type: itestb,itest'
     &                       ,',jben',itestb,itest,jben
                        write(6,*) 'itest iben type ',itbentest
                        write(6,*) 'itestb jben type ',jtbentest
                        if ( lintraterm ) then
                           lfailure = .true.
                           return
                        endif
                     endif
                     lfound = .true.
                  endif
               enddo
               if ( .not. lfound ) then 
                  write(6,*) 'Problem: imolty,iunit,iben',imolty
     &                 ,iunit,iben
                  write(6,*) 'There is no symmetric bend'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
               call twh_itben(GLB_GET,imolty,iunit,iben,itbentest)
               call twh_lbonbon(GLB_GET,itbentest,lbonbon)
               if ( lbonbon ) then
c                 --- we have at least on bond-bond interaction
                  lbnbn = .true.
c                 --- do not double count as we are on the center atom
                  if ( iunit .lt. itestb ) then
                     call twh_inbnbn(GLB_INCR,imolty,itest,1)
                     call twh_inbnbn(GLB_GET,imolty,itest,count)
                     call twh_itben(GLB_GET,imolty,iunit,iben,itemp)
                     call twh_orderben(GLB_GET,imolty,iunit,iben,iorder)
                     if ( iorder .eq. ORD_REVERSE ) then
                        itemp = -itemp
                     endif
                     call twh_itbnbn(GLB_SET,imolty,itest,count,itemp)
                     call twh_ijbnbn1(GLB_SET,imolty,itest,count,iunit)
                     call twh_ijbnbn2(GLB_SET,imolty,itest,count,itestb)
                     call twh_ibtben1(GLB_GET,imolty,iunit,iben,itemp)
                     call twh_ibtbnbn1(GLB_SET,imolty,itest,count,itemp)
                     call twh_ibtben2(GLB_GET,imolty,iunit,iben,itemp)
                     call twh_ibtbnbn2(GLB_SET,imolty,itest,count,itemp)
                  endif
               endif
               call twh_itben(GLB_GET,imolty,iunit,iben,atype)
               call twh_anglestyle(GLB_GET,atype,anglestyle)
               if ( anglestyle .eq. ANG_STYLE_FOURIER
     &              .or. anglestyle .eq. ANG_STYLE_TTFOURIER ) then
c                 --- this only works with the UFF 12-6 classical
c                 --- potential
                  call twh_cp_number(GLB_GET,cpnumber)
                  if ( cpnumber .ne. CP_UFF_12_6 ) then
                     write(6,*) 'CHECKSTRUC: cannot combine anglestyle'
     &                    ,anglestyle
                     write(6,*) 'with anything other than a'
     &                    ,' classical_potential of UFF 12-6'
                     lfailure = .true.
                     return
                  endif
c                 --- using the UFF angle terms means we need to 
c                 --- setup the prefactor here to save time later on
c                 --- set zi
                  call twh_ntype(GLB_GET,imolty,iunit,itype)
                  ii = twh_getnbtype(itype,itype)
                  call twh_nbcoeff(GLB_GET,10,ii,zi)
c                 --- set zk
                  call twh_ijben3(GLB_GET,imolty,iunit,iben,kunit)
                  call twh_ntype(GLB_GET,imolty,kunit,ktype)
                  kk = twh_getnbtype(ktype,ktype)
                  call twh_nbcoeff(GLB_GET,10,kk,zk)
c                 --- set theta naught
                  call twh_itben(GLB_GET,imolty,iunit,iben,itemp)
                  call twh_bencoeff(GLB_GET,itemp,0,theta)
c                 --- get the equilibrium bond lengths
                  call twh_ibtben1(GLB_GET,imolty,iunit,iben,itemp)
                  call twh_vibcoeff(GLB_GET,itemp,0,rij)
                  call twh_ibtben2(GLB_GET,imolty,iunit,iben,itemp)
                  call twh_vibcoeff(GLB_GET,itemp,0,rjk)
c                 --- determine the rik squared value
                  rik = dsqrt(rij*rij + rjk*rjk 
     &                 - 2.0d0*rij*rjk*dcos(theta))
c                 --- set benprefact
                  dvalue = CON_KCALTOK*(664.12d0)
     &                 *(zi*zk/(rik**5))*( 
     &                 3.0d0*rij*rjk*(1.0d0 - dcos(theta)**2)
     &                 - rik*rik*dcos(theta))
                  if ( anglestyle .eq. ANG_STYLE_FOURIER ) then
c                    --- divide out by n^2
                     call twh_itben(GLB_GET,imolty,iunit,iben,itemp)
                     call twh_bencoeff(GLB_GET,itemp,1,thetatwo)
                     dvalue = dvalue/(thetatwo**2)
                  endif
                  call twh_benprefact(GLB_SET,imolty,iunit,iben,dvalue)
               endif
            enddo
         enddo

         do iunit = 1,nunit
c           --- check torsions
            call twh_intor(GLB_GET,imolty,iunit,intorsion)
            do itor = 1,intorsion
               call twh_ijtor2(GLB_GET,imolty,iunit,itor,itest)
               call twh_ijtor3(GLB_GET,imolty,iunit,itor,itestb)
               call twh_ijtor4(GLB_GET,imolty,iunit,itor,itestc)
c              --- check that atom #2 is bonded to this atom
               lfound = .false.
               call twh_invib(GLB_GET,imolty,iunit,invibunit)
               do ivib = 1,invibunit
                  call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
                  if ( ijvibtest .eq. itest ) then
                     lfound = .true.
                     call twh_itvib(GLB_GET,imolty,iunit,ivib,itemp)
                     call twh_ibttor1(GLB_SET,imolty,iunit,itor,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 1-2 connected in torsion are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #3 is bonded to atom #2
               lfound = .false.
               call twh_invib(GLB_GET,imolty,itest,invibtest)
               do ivib = 1,invibtest
                  call twh_ijvib(GLB_GET,imolty,itest,ivib,ijvibtest)
                  if ( ijvibtest .eq. itestb ) then
                     lfound = .true.
                     call twh_itvib(GLB_GET,imolty,itest,ivib,itemp)
                     call twh_ibttor2(GLB_SET,imolty,iunit,itor,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 2-3 connected in torsion are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #4 is bonded to atom #3
               lfound = .false.
               call twh_invib(GLB_GET,imolty,itestb,invibtest)
               do ivib = 1,invibtest
                  call twh_ijvib(GLB_GET,imolty,itestb,ivib,ijvibtest)
                  if ( ijvibtest .eq. itestc ) then
                     lfound = .true.
                     call twh_itvib(GLB_GET,imolty,itestb,ivib,itemp)
                     call twh_ibttor3(GLB_SET,imolty,iunit,itor,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 3-4 connected in torsion are not '
     &                 ,'connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that the 1-2-3 bend exists
               lfound = .false.
               call twh_inben(GLB_GET,imolty,iunit,inbend)
               do iben = 1, inbend
                  call twh_ijben2(GLB_GET,imolty,iunit,iben,ijben2test)
                  call twh_ijben3(GLB_GET,imolty,iunit,iben,ijben3test)
                  if ( ijben2test .eq. itest .and.
     &                 ijben3test .eq. itestb ) then
                     lfound = .true.
                     call twh_itben(GLB_GET,imolty,iunit,iben,itemp)
                     call twh_iattor1(GLB_SET,imolty,iunit,itor,itemp) 
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 1-2-3 connected in torsion are not '
     &                 ,'connected by a bend'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that the 2-3-4 bend exists
               lfound = .false.
               call twh_inben(GLB_GET,imolty,itest,inbend)
               do iben = 1, inbend
                  call twh_ijben2(GLB_GET,imolty,itest,iben,ijben2test)
                  call twh_ijben3(GLB_GET,imolty,itest,iben,ijben3test)
                  if ( ijben2test .eq. itestb .and.
     &                 ijben3test .eq. itestc ) then
                     lfound = .true.
                     call twh_itben(GLB_GET,imolty,itest,iben,itemp)
                     call twh_iattor2(GLB_SET,imolty,iunit,itor,itemp)
                  endif
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijtor2,ijtor3'
     &                 ,',ijtor4,itor',imolty,iunit,itest,itestb,itestc
     &                 ,itor
                  write(6,*) 'Atoms 2-3-4 connected in torsion are not '
     &                 ,'connected by a bend'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that each other atom also contains this 
c              --- interaction
               lfound = .false.
               call twh_intor(GLB_GET,imolty,itestc,intorsionc)
               do jtor = 1,intorsionc
                  call twh_ijtor4(GLB_GET,imolty,itestc,jtor,ijtor4j)
                  call twh_ijtor3(GLB_GET,imolty,itestc,jtor,ijtor3j)
                  call twh_ijtor2(GLB_GET,imolty,itestc,jtor,ijtor2j)
                  if ( ijtor4j .eq. iunit .and. ijtor3j .eq. itest
     &                 .and. ijtor2j .eq. itestb ) then
                     call twh_ittor(GLB_GET,imolty,itestc,jtor,ittorj)
                     call twh_ittor(GLB_GET,imolty,iunit,itor,ittori)
                     if ( abs(ittorj) .ne. abs(ittori) ) then
                        write(6,*) 'Problem: imolty,iunit,itor',imolty
     &                       ,iunit,itor
                        write(6,*) 'forward atoms ',iunit,itest,itestb
     &                       ,itestc
                        write(6,*) 'reverse atoms ',itestc
     &                       ,ijtor2j,ijtor3j,ijtor4j
                        write(6,*) 'iunit itor forward type ',ittori
                        write(6,*) 'itestb jtor reverse type ',ittorj
                        if ( lintraterm ) then
                           lfailure = .true.
                           return
                        endif
                     endif
                     lfound = .true.
                  endif
               enddo
               if ( .not. lfound ) then 
                  write(6,*) 'Problem: imolty,iunit,itor',imolty
     &                 ,iunit,itor
                  write(6,*) 'There is no symmetric torsion'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
         enddo

c        --- setup the angletypes for the angle-angle terms
         do iunit = 1,nunit
            call twh_inaa(GLB_GET,imolty,iunit,inaaunit)
            do iaa = 1,inaaunit
               call twh_ijaa0(GLB_GET,imolty,iunit,iaa,itest)
               call twh_ijaa1(GLB_GET,imolty,iunit,iaa,itestb)
               call twh_ijaa2(GLB_GET,imolty,iunit,iaa,itestc)
               call twh_inben(GLB_GET,imolty,itest,inbendt)
               do iben = 1,inbendt
                  call twh_ijben2(GLB_GET,imolty,itest,iben,ijben2test)
                  if ( iunit .eq. ijben2test ) then
                     call twh_ijben3(GLB_GET,imolty,itest,iben
     &                    ,ijben3test)
                     if (itestb .eq. ijben3test ) then
                        call twh_itben(GLB_GET,imolty,itest,iben,itemp)
                        call twh_iataa1(GLB_SET,imolty,iunit,iaa,itemp)
                     elseif( itestc .eq. ijben3test ) then
                        call twh_itben(GLB_GET,imolty,itest,iben,itemp)
                        call twh_iataa2(GLB_SET,imolty,iunit,iaa,itemp)
                     endif
                  endif
               enddo
            enddo
         enddo

c        --- check improper torsions
         do iunit = 1,nunit
            call twh_inimprop(GLB_GET,imolty,iunit,inimproper)
            do iimp = 1, inimproper
c              --- check that atom #2 is bonded to first atom
               call twh_ijimprop2(GLB_GET,imolty,iunit,iimp,itest)
               lfound = .false.
               call twh_invib(GLB_GET,imolty,iunit,invibunit)
               do ivib = 1,invibunit
                  call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
                  if ( ijvibtest .eq. itest ) lfound = .true.
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijimprop2,iimp'
     &                 ,imolty,iunit,itest,iimp
                  write(6,*) 'Atoms 1-2 connected in improper are '
     &                 ,'not connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #3 is bonded to first atom
               call twh_ijimprop3(GLB_GET,imolty,iunit,iimp,itestb)
               lfound = .false.
               do ivib = 1,invibunit
                  call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
                  if ( ijvibtest .eq. itestb ) lfound = .true.
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijimprop3,iimp'
     &                 ,imolty,iunit,itestb,iimp
                  write(6,*) 'Atoms 1-3 connected in improper are '
     &                 ,'not connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
c              --- check that atom #4 is bonded to first atom
               call twh_ijimprop4(GLB_GET,imolty,iunit,iimp,itestc)
               lfound = .false.
               do ivib = 1,invibunit
                  call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
                  if ( ijvibtest .eq. itestc ) lfound = .true.
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'Problem: imolty,iunit,ijimprop4,iimp'
     &                 ,imolty,iunit,itestc,iimp
                  write(6,*) 'Atoms 1-4 connected in improper are '
     &                 ,'not connected by a vibration'
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            enddo
         enddo
c        --- check special one-five
         do iunit = 1,nunit
            call twh_inof(GLB_GET,imolty,iunit,inonefivei)
            do iof = 1,inonefivei
               call twh_ijof5(GLB_GET,imolty,iunit,iof,matom)
               lfound = .false.
               call twh_inof(GLB_GET,imolty,matom,inonefivej)
               do jof = 1,inonefivej
                  call twh_ijof5(GLB_GET,imolty,matom,jof,ijof5test)
                  if ( ijof5test .eq. iunit ) lfound = .true.
               enddo
               if ( .not. lfound ) then
                  write(6,*) 'CHECKSTRUC: non-symetric one-five'
                  write(6,*) 'imolty,iunit,iof',imolty,iunit,iof
                  write(6,*) 'has no return interaction from ijof5'
     &                 ,matom
                  lfailure = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- need some special logic for torsions that have 1-4 nonbonded
c     --- interactions.  In certain cases (6-membered rings) an atom
c     --- can have multiple torsions which end at the same 4th bead
c     --- In this case we need to scale the nonbonded interactions to 
c     --- avoid counting the nonbond terms multiple times
c     --- for other ring structures it is possible that the 1-4 is also
c     --- a 1-3 or 1-2 (smaller rings) and we want to disable the
c     --- nonbonded calculation
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
c           --- initialize the array of cross terms
            do junit = 1,nunit
               torcross(junit) = 0
            enddo
            call twh_intor(GLB_GET,imolty,iunit,numtor)
            do itor = 1,numtor
c              --- check to see if there are nonbonded terms for this
c              --- torsion
               call twh_ittor(GLB_GET,imolty,iunit,itor,itype)
               call twh_loftor(GLB_GET,iabs(itype),lonefourtor)
               if ( lonefourtor ) then
                  call twh_ijtor4(GLB_GET,imolty,iunit,itor,ijtor4i)
                  torcross(ijtor4i) = torcross(ijtor4i) + 1
               endif
            enddo
c           --- now make sure we don't have bonded or angle interaction
c           --- with that 4th bead (possible in rings)
            do itor = 1,numtor
               call twh_ijtor4(GLB_GET,imolty,iunit,itor,junit)
               lfound = .false.
               call twh_invib(GLB_GET,imolty,iunit,invibunit)
               do ivib = 1,invibunit
                  call twh_ijvib(GLB_GET,imolty,iunit,ivib,ijvibtest)
                  if ( junit .eq. ijvibtest ) lfound = .true.
               enddo
               call twh_inben(GLB_GET,imolty,iunit,inbend)
               do iben = 1,inbend
                  call twh_ijben3(GLB_GET,imolty,iunit,iben,ijben3test)
                  if ( junit .eq.ijben3test ) lfound = .true.
               enddo
c              --- store the scaling factor for each torsion
               if ( lfound .or. torcross(junit) .eq. 0 ) then
c                 --- we do not want to count these nonbonded interact
                  call twh_torofcode(GLB_SET,imolty,iunit,itor
     &                 ,CODE_ZERO)
               else
c                 --- use the torcross to compute scale factor
                  if ( torcross(junit) .gt. 1 ) then
                     dvalue = 1.0d0 / dble(torcross(junit))
                     call twh_torofcode(GLB_SET,imolty,iunit,itor
     &                    ,CODE_RESCALE)
                     call twh_itscale(GLB_SET,imolty,iunit,itor,dvalue)
                  endif
               endif
            enddo
         enddo
      enddo

c     --- check that the hybridization of the atoms is consistent with
c     --- the number of bonds.  Note that since united-atom force fields
c     --- are possible there may not be an exact match, this is more of
c     --- a less than or equal to senario
c     --- we mostly just print out warnings, but if we have too many 
c     --- bonds to an atom for its hybridization then we are in trouble
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_invib(GLB_GET,imolty,iunit,ivib)
            call twh_ntype(GLB_GET,imolty,iunit,itype)
            call twh_bondpatt(GLB_GET,itype,bondpatt)
            if ( bondpatt.eq. 'sp3' ) then
               if ( ivib .gt. 4 ) then
                  write(6,*) 'too many bonds to an sp3 atom'
                  write(6,*) 'imolty,iunit,ivib',imolty,iunit,ivib
                  write(6,*) 'itype,bondpatt',itype,bondpatt
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            elseif ( bondpatt .eq. 'sp2' ) then
               if ( ivib .gt. 3 ) then
                  write(6,*) 'too many bonds to an sp2 atom'
                  write(6,*) 'imolty,iunit,ivib',imolty,iunit,ivib
                  write(6,*) 'itype,bondpatt',itype,bondpatt
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            elseif ( bondpatt .eq. 'sp' ) then
               if ( ivib .gt. 2 ) then
                  write(6,*) 'too many bonds to an sp atom'
                  write(6,*) 'imolty,iunit,ivib',imolty,iunit,ivib
                  write(6,*) 'itype,bondpatt',itype,bondpatt
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            elseif ( bondpatt .eq. 's' ) then
               if ( ivib .gt. 1 ) then
                  write(6,*) 'too many bonds to an s atom'
                  write(6,*) 'imolty,iunit,ivib',imolty,iunit,ivib
                  write(6,*) 'itype,bondpatt',itype,bondpatt
                  if ( lintraterm ) then
                     lfailure = .true.
                     return
                  endif
               endif
            endif
         enddo
      enddo

c     --- check that the molecule is fully connected and there is no
c     --- subset of the molecule that is not connected to the rest.
      do imolty = 1, nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         call twh_checkconnect(imolty,nunit)
         if ( lfailure ) return
      enddo

      return
      end

      subroutine twh_checkconnect(imolty,nunit)
c     ******************************************************************
c     * checks the bonding pattern of a molecule to make sure it is    *
c     * fully connected.  Stores this information in lbondgraph        *
c     *                                                                *
c     * originally written 05-30-2013 by M.G. Martin                   *
c     * last modified 07-15-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer imolty,nunit
c     --- local variables
      logical lexist,lbondgraph
      integer workcount,ivib,numvib,iother,icount,icurrent,iunit
      integer worklist
      dimension worklist(nunit*NNBOND)

c     --- initialise lexist to all false
      call twh_lexist(GLB_INIT,nunit,.false.)
c     --- start with atom 1 and progress through the bonded stucture
      workcount = 1
      worklist(workcount) = 1

      do while ( workcount .gt. 0 ) 
         icurrent = worklist(1)
         call twh_lexist(GLB_GET,icurrent,lexist)
         if ( lexist ) then
c           --- this unit has already been found to exist and tested
c           --- nothing else to do here
         else
c           --- add all of the atoms bonded to icurrent to the worklist
            call twh_invib(GLB_GET,imolty,icurrent,numvib)
            do ivib = 1,numvib
               call twh_ijvib(GLB_GET,imolty,icurrent,ivib,iother)
               workcount = workcount + 1
               worklist(workcount) = iother
            enddo
c           --- set this atom to exist
            call twh_lexist(GLB_SET,icurrent,.true.)
         endif
c        --- decrement and shift the worklist
         do icount = 1,workcount - 1
            worklist(icount) = worklist(icount+1)
         enddo
         workcount = workcount - 1
      enddo

      lbondgraph = .true.
      do iunit = 1,nunit
         call twh_lexist(GLB_GET,iunit,lexist)
         if ( .not. lexist ) then
c           --- this atom was never found in the pattern
c           --- warn on first occurance
            if ( lbondgraph ) then
               write(6,'(a,1x,i4,1x,a,1x,i6,1x,a)')
     &              'Molecule type',imolty,'unit'
     &              ,iunit,' not connected via bonding graph to atom 1'
               write(6,'(2x,a,1x,a)')
     &              'therefore Monte Carlo moves will not be allowed'
     &              ,'to alter atom positions for this molecule type'
               lbondgraph = .false.
            endif
         endif
      enddo
c     --- set lbondgraph
      call twh_lbondgraph(GLB_SET,imolty,lbondgraph)

      return
      end

