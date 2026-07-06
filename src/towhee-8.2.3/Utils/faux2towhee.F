c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004 Marcus G. Martin                            *
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
      program faux2towhee
c     *****************************************************************
c     * converts the files from David Faux into a format for towhee   *
c     * originally written 02-09-2001 by M.G. Martin                  *
c     * last modified 06-28-2004 by M.G. Martin                       *
c     *****************************************************************
      implicit none
      integer natoms,nangles,nnbond,nntype
      parameter (nnbond = 4, natoms = 672, nangles = 1152,nntype=10)

      logical lfound

      integer invib,ijvib,inpstyle
      dimension invib(natoms),ijvib(natoms,nnbond)
      integer inben,ijben1,ijben2
      dimension inben(natoms)
      dimension ijben1(natoms,nnbond**2),ijben2(natoms,nnbond**2)
      integer type
      dimension type(natoms)

      integer ivib,ibend,iuone,iutwo,iuthree
      double precision qqu,xpos,ypos,zpos,xadjust,yadjust,zadjust
      dimension qqu(natoms),xpos(natoms),ypos(natoms),zpos(natoms)

      integer count
      dimension count(4)
      integer iatom,iangle,iend,idum,ispot,triplet
      dimension triplet(nangles,3)

      character*4 cdum
      character*6 name
      dimension name(nntype)

c     --- here we set the inpstyle that we wish to generate
      inpstyle = 2

      open(20,file='faux_4a.atoms',form='formatted')
      do iatom = 1,natoms
         read(20,*) idum,idum,type(iatom),qqu(iatom),xpos(iatom)
     &        ,ypos(iatom),zpos(iatom)
      enddo
      close(20)

      open(21,file='faux_4a.ang',form='formatted')
c     --- skip the first 10 lines
      do iangle = 1,10
         read(21,*)
      enddo
      do iangle = 1,nangles
         read(21,*) cdum,(triplet(iangle,ispot),ispot=1,3)
      enddo
      close(21)
      
c     --- initialize arrays
      do iatom = 1,natoms
         invib(iatom) = 0
         inben(iatom) = 0
      enddo

c     --- find vibrations and bends
      do iangle = 1,nangles
         do iend = 1,2
            if ( iend .eq. 1 ) then
               iuone = triplet(iangle,1)
               iuthree = triplet(iangle,3)
            else
               iuone = triplet(iangle,3)
               iuthree = triplet(iangle,1)
            endif
            iutwo = triplet(iangle,2)

c           --- check the iuone list for bonds
            lfound = .false.
            do ivib = 1,invib(iuone)
               if ( ijvib(iuone,ivib) .eq. iutwo ) lfound = .true.
            enddo
            if ( .not. lfound ) then
c              --- add this to the vibration list 
               invib(iuone) = invib(iuone) + 1
               ijvib(iuone,invib(iuone)) = iutwo
            endif

c           --- check the iutwo list for bonds
            lfound = .false.
            do ivib = 1,invib(iutwo)
               if ( ijvib(iutwo,ivib) .eq. iuone ) lfound = .true.
            enddo
            if ( .not. lfound ) then
c              --- add this to the vibration list 
               invib(iutwo) = invib(iutwo) + 1
               ijvib(iutwo,invib(iutwo)) = iuone
            endif

c           --- check the iuone list for angles
            lfound = .false.
            do ibend = 1,inben(iuone)
               if ( ijben1(iuone,ibend) .eq. iutwo .and.
     &              ijben2(iuone,ibend) .eq. iuthree ) lfound=.true.
            enddo
            if ( .not. lfound ) then
c              --- add this to the angle list
               inben(iuone) = inben(iuone) + 1
               ijben1(iuone,inben(iuone)) = iutwo
               ijben2(iuone,inben(iuone)) = iuthree
            endif

c           --- check the iuthree list for angles
            lfound = .false.
            do ibend = 1,inben(iuthree)
               if ( ijben1(iuthree,ibend) .eq. iutwo .and.
     &              ijben2(iuthree,ibend) .eq. iuone ) lfound=.true.
            enddo
            if ( .not. lfound ) then
c              --- add this to the angle list
               inben(iuthree) = inben(iuthree) + 1
               ijben1(iuthree,inben(iuthree)) = iutwo
               ijben2(iuthree,inben(iuthree)) = iuone
            endif
         enddo
      enddo

c     --- initialize counter
      do iend = 1,4
         count(iend) = 0
      enddo

c     --- setup type converter
      name(1) = 'OZ'
      name(2) = 'Si'
      name(3) = 'Al'
      name(4) = 'Na'

c     --- convert atom types
      do iatom = 1,natoms
         if ( type(iatom) .eq. 10 ) then
            type(iatom) = 1 
            count(1) = count(1) + 1
        elseif ( type(iatom) .eq. 8 ) then
            type(iatom) = 2
            count(2) = count(2) + 1
         elseif ( type(iatom) .eq. 9 ) then
            type(iatom) = 3
            count(3) = count(3) + 1
         elseif ( type(iatom) .eq. 11 ) then
            type(iatom) = 4
            count(4) = count(4) + 1
         else
            write(6,*) 'unknown atom type ',iatom,type(iatom)
         endif
      enddo

      do iend = 1,4
         write(6,*) 'atom count',iend,count(iend)
      enddo

c     --- output all of the information for towhee
      open(3,file='towhee_altinp',form='formatted')
      write(3,*) 'inpstyle'
      write(3,*) inpstyle
      write(3,*) 'nunit'
      write(3,*) natoms
      write(3,*) 'nmaxcbmc'
      write(3,*) natoms
      write(3,*) 'lelect'
      write(3,*) 'T'
      if ( inpstyle .eq. 0 ) then
         write(3,*) 'lpdb'
         write(3,*) 'F'
         do iatom = 1,natoms
            write(3,*) 'unit ntype qqatom'
            write(3,*) iatom,type(iatom),qqu(iatom)
            write(3,*) 'vibration'
            write(3,*) invib(iatom)
            do ivib = 1,invib(iatom)
               write(3,*) ijvib(iatom,ivib),49
            enddo
            write(3,*) 'bending'
            write(3,*) inben(iatom)
            do ibend = 1,inben(iatom)
               write(3,*) ijben1(iatom,ibend),ijben2(iatom,ibend),48
            enddo
            write(3,*) 'torsion'
            write(3,*) 0
            write(3,*) 'angle-angle'
            write(3,*) 0
            write(3,*) 'improper torsion'
            write(3,*) 0
         enddo
      elseif ( inpstyle .eq. 2 ) then
         write(3,*) 'forcefield'
         write(3,*) 'Faux'
         do iatom = 1,natoms
            write(3,*) 'unit ntype qqatom'
            write(3,*) iatom,' ',char(39),name(type(iatom)),char(39)
     &           ,' ',qqu(iatom)
            write(3,*) 'vibration'
            write(3,*) invib(iatom)
            if ( invib(iatom) .ne. 0 ) then
               write(3,*) (ijvib(iatom,ivib),ivib=1,invib(iatom))
            endif
            write(3,*) 'improper torsion'
            write(3,*) 0
         enddo
      else
         write(6,*) 'invalid inpstyle ',inpstyle
         stop
      endif
c     --- close the output file
      close(3)

      xadjust = 24.555d0/2.0d0
      yadjust = 24.555d0/2.0d0
      zadjust = 24.555d0/2.0d0
c     --- output all of the information for towhee_coords
      open(9,file='towhee_coords',form='formatted')
      do iatom = 1,natoms
         write(9,*) xpos(iatom)+xadjust
     &        ,ypos(iatom)+yadjust,zpos(iatom)+zadjust,qqu(iatom)
      enddo
      close(9)

      end
      


