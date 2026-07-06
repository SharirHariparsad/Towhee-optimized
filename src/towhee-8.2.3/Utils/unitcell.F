c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
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
      program unitcell
c     ******************************************************************
c     * takes a unit cell and creates a towhee_coords file that has    *
c     * as many images of the unit cell in each dimension that the user*
c     * specifies                                                      *
c     *                                                                *
c     * originally 05-15-2002 by M.G. Martin                           *
c     * last modified 06-29-2004 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- parameters
      integer maxatom
      parameter(maxatom=300)

c     --- local variables
      character*10 stringid
      dimension stringid(maxatom)
      integer iatom,natom,xdup,ydup,zdup
      integer dup
      dimension dup(3)
      double precision boxlength
      dimension boxlength(3)
      double precision xdisp,ydisp,zdisp
      double precision apos,bpos,cpos
      dimension apos(maxatom),bpos(maxatom),cpos(maxatom)

      write(6,*) 'This program only works for retangular boxes'
      write(6,*) 'Enter the original unit cell x-dimension'
      read(5,*) boxlength(1)
      write(6,*) 'Enter the original unit cell y-dimension'
      read(5,*) boxlength(2)
      write(6,*) 'Enter the original unit cell z-dimension'
      read(5,*) boxlength(3)

      write(6,*) 'Enter the number of atoms in the original unit cell'
      read(5,*) natom
      if ( natom .gt. maxatom ) then
         write(6,*) 'natom ',natom,' exceeds maxatom',maxatom
         stop
      endif
      write(6,*) 'Enter the number of new duplicates in the x-dimension'
      read(5,*) dup(1)
      write(6,*) 'Enter the number of new duplicates in the y-dimension'
      read(5,*) dup(2)
      write(6,*) 'Enter the number of new duplicates in the z-dimension'
      read(5,*) dup(3)

c     --- read the old positions
      open(10,file='towhee_coords',form='formatted')
      do iatom = 1,natom
         read(10,*) apos(iatom),bpos(iatom),cpos(iatom),stringid(iatom)
      enddo
      close(10)

c     --- write the new positions
      open(12,file='towhee_newcoords',form='formatted')
      do iatom = 1,natom
         do xdup = 1,dup(1)
            xdisp = dble(xdup-1)*boxlength(1)
            do ydup = 1,dup(2)
               ydisp = dble(ydup-1)*boxlength(2)
               do zdup = 1,dup(3)
                  zdisp = dble(zdup-1)*boxlength(3)
                  write(12,*) apos(iatom)+xdisp
     &                 ,bpos(iatom)+ydisp
     &                 ,cpos(iatom)+zdisp,' ',stringid(iatom)
               enddo
            enddo
         enddo
      enddo
      close(12)
      write(6,*) 'New box dimensions'
     &     ,dble(dup(1))*boxlength(1)
     &     ,dble(dup(2))*boxlength(2)
     &     ,dble(dup(3))*boxlength(3)

      end

