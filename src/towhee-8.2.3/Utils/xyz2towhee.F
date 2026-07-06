c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2016 Marcus G. Martin                       *
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
      program xyz2towhee
c     ******************************************************************
c     * converts from xyz to towhee format                             *
c     * this file is rewritten to convert various files into Towhee    *
c     *                                                                *
c     * originally written 03-27-2003 by M.G. Martin                   *
c     * last modified 08-25-2016 by M.G. Martin                        *
c     ******************************************************************
      implicit none

      character*2 atomname
      character*30 fileinput,fileoutput
      integer natoms,iatom
      double precision xxx,yyy,zzz

c     --- filename passed in when running the code
      write(6,*) 'Please enter an input filename'
      read(5,*) fileinput
      open(unit=10,file=fileinput,form="formatted")

      write(6,*) 'Please enter an output filename'
      read(5,*) fileoutput
      open(unit=20,file=fileoutput,form="formatted")

c     --- get the number of atoms
      read(10,*) natoms

c     --- get the atom positions and types
      read(10,*)
      do iatom = 1,natoms
         read(10,*) atomname,xxx,yyy,zzz
         write(20,*) xxx,yyy,zzz,atomname
      enddo
      close(10)
      close(20)

      end
