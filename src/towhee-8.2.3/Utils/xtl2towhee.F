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
      program xtl2towhee
c     ******************************************************************
c     * converts files from the xtl format into the Towhee format      *
c     * currently this is not especially general, but should prove     *
c     * useful for most cases with a bit of tinkering                  *
c     *                                                                *
c     * originally written 09-02-2004 by M.G. Martin                   *
c     * last modified 09-16-2004 by M.G. Martin                        *
c     ******************************************************************
      implicit none
      integer MAXFILE
      parameter (MAXFILE=20)
c     ---- local variables
      character*1 testatom,letter
      dimension testatom(MAXFILE)
      character*10 name
      character*20 filename
      integer natoms,iatom
      integer iskip,nskip,filetot,myfile,index
      integer file_O,file_H,file_OH,file_OHH

      logical lfound,lmerge

      double precision xdim,ydim,zdim,xmin,ymin,zmin,xmax,ymax,zmax
      double precision xcord,ycord,zcord,charge

      write(6,*) 'Input the filename you wish to convert'
      read(5,'(a20)') filename

      open(20,file=filename,form='formatted')

c     --- skip header lines
      nskip = 3
      do iskip = 1,nskip
         read(20,*)
      enddo
c     --- read in the box dimensions
      read(20,*) xdim,ydim,zdim
c     --- skip some more lines
      nskip = 5
      do iskip = 1,nskip
         read(20,*)
      enddo

      write(6,*) 'Please enter the number of atoms in this file'
      read(5,*) natoms

c     --- run through the file once to determine the maximum and
c     --- minimum box dimensions
      do iatom = 1,natoms
         read(20,*) name,xcord,ycord,zcord,charge
c        --- test for min and max
         if ( iatom .eq. 1 ) then
c           --- just set min and max to these values
            xmin = xcord
            ymin = ycord
            zmin = zcord
            xmax = xcord
            ymax = ycord
            zmax = zcord
         else
c           --- test to see if these values exceed the current
c           --- min and max values
            if ( xcord .lt. xmin ) xmin = xcord
            if ( ycord .lt. ymin ) ymin = ycord
            if ( zcord .lt. zmin ) zmin = zcord
            
            if ( xcord .gt. xmax ) xmax = xcord
            if ( ycord .gt. ymax ) ymax = ycord
            if ( zcord .gt. zmax ) zmax = zcord
         endif
      enddo

c     --- output the original minimums and maximums
      write(6,*) 'original xmin,xmax',xmin*xdim,xmax*xdim
      write(6,*) 'original ymin,ymax',ymin*ydim,ymax*ydim
      write(6,*) 'original zmin,zmax',zmin*zdim,zmax*zdim
c     --- output the final box lengths
      write(6,*) 'x boxlength',xdim*(xmax-xmin)
      write(6,*) 'y boxlength',ydim*(ymax-ymin)
      write(6,*) 'z boxlength',zdim*(zmax-zmin)

c     --- go through the file again and output the adjust coords

c     --- reinitialize filename
      filename = 'output.'
      filetot = 0
      do myfile = 1,MAXFILE
         testatom(myfile) = ''
      enddo

c     --- rewind the file
      rewind(20)

c     --- skip header lines
      nskip = 3
      do iskip = 1,nskip
         read(20,*)
      enddo
c     --- read in the box dimensions
      read(20,*) xdim,ydim,zdim
c     --- skip some more lines
      nskip = 5
      do iskip = 1,nskip
         read(20,*)
      enddo

      do iatom = 1,natoms
         read(20,*) name,xcord,ycord,zcord,charge
c        --- see if we already have found this atom type
         lfound = .false.
         letter = name(1:1)
         do myfile = 1,filetot
            if ( letter .eq. testatom(myfile) ) then
c               --- found a match
               lfound = .true.
               index = 20+myfile
            endif
         enddo
         if ( .not. lfound ) then
c            --- add this to the list
            filetot = filetot + 1
c           --- bounds check filetot
            if ( filetot .gt. MAXFILE ) then
               write(6,*) 'filetot',filetot
               write(6,*) 'Error: you need to increase MAXFILE'
               stop
            endif
c           --- set testatom
            testatom(filetot) = letter
            index = 20+filetot
            filename(8:8) = letter
            write(6,*) 'name =',name,'='
            write(6,*) 'letter =',letter,'='
            write(6,*) 'filename =',filename,'='
c           --- open the file
            open(index,file=filename,form='formatted')
         endif
c        --- output the information to the appropriate file
         write(index,*) xcord*xdim,ycord*ydim
     &        ,zcord*zdim,charge,' ',name
      enddo

c     --- close files
      close(20)
      do myfile = 1,filetot
         index = 20+myfile
         close(index)
         rewind(index)
      enddo

      write(6,*) 'would you like to merge the O and H files?'
      read(5,*) lmerge

      if ( lmerge ) then
c       --- complete and total hack for a gibbsite file I was 
c       --- working with

c        --- open the O file
         file_O = 10
         filename = 'output.O'
         open(file_O,file=filename,form='formatted')
c        --- open the H file
         file_H = 11
         filename = 'output.H'
         open(file_H,file=filename,form='formatted')
c        --- open the OH file
         file_OH = 12
         filename = 'combo.OH'
         open(file_OH,file=filename,form='formatted')
c        --- open the OHH file
         file_OHH = 13
         filename = 'combo.OHH'
         open(file_OHH,file=filename,form='formatted')

         do iatom = 1,576
            read(file_O,*) xcord,ycord,zcord,charge,name
            write(file_OH,*) xcord,ycord,zcord,charge,' ',name
            read(file_H,*) xcord,ycord,zcord,charge,name
            write(file_OH,*) xcord,ycord,zcord,charge,' ',name
         enddo

         do iatom = 1,517
            read(file_O,*) xcord,ycord,zcord,charge,name
            write(file_OHH,*) xcord,ycord,zcord,charge,' ',name
            read(file_H,*) xcord,ycord,zcord,charge,name
            write(file_OHH,*) xcord,ycord,zcord,charge,' ',name
            read(file_H,*) xcord,ycord,zcord,charge,name
            write(file_OHH,*) xcord,ycord,zcord,charge,' ',name
         enddo
        
         close(file_O)
         close(file_H)
         close(file_OH)
         close(file_OHH)
        
      endif

      end
      


