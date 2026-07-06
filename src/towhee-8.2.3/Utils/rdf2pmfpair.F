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
      program rdf2pmfpair
c     ******************************************************************
c     * this program converts a radial distribution function into a    *
c     * pair-aproximiation of the potential of mean force              *
c     * designed for use with the tabulated pair potentype             *
c     *                                                                *
c     * originally written 09-03-2004 by M.G. Martin                   *
c     * last modified 09-23-2004 by M.G. Martin                        *
c     ******************************************************************
      implicit none
      integer npoint,ipoint
      double precision temperature,rdf,pmf,position,large,rcut
      character*50 filename

      write(6,*) 'Please input the rdf filename'
      read(5,*) filename

      open(20,file=filename,form='formatted')

      write(6,*) 'Please input the number of data lines'
      read(5,*) npoint

      write(6,*) 'Please input the temperature [K]'
      read(5,*) temperature

      write(6,*) 'Please input the cutoff'
     &     ,' (will be set to zero at cutoff)'
      read(5,*) rcut

c     --- set large number
      large = -temperature*dlog(0.00000001d0)

c     --- open the output filename
      filename = 'towhee_pmf'
      open(30,file=filename,form='formatted')

c     --- output some header information
      write(30,'(a11)') 'pmf version'
      write(30,*) 1
      write(30,'(a10)') 'datapoints'
      write(30,*) npoint+1
      write(30,'(a4)') 'data'

      do ipoint = 1,npoint
         read(20,*) position,rdf
c        --- need to be careful not to take the log of zero
         if ( rdf .le. 0.0d0 ) then
c           --- just set to a really large number
            pmf = large
         else
c           --- compute the potential of mean force pair term
            pmf = -temperature*dlog(rdf)
         endif
         write(30,*) position,pmf
      enddo

c     --- make sure the cutoff is larger than the final point
      if ( rcut .lt. position ) then
         write(6,*) 'ERROR: final point is larger than cutoff'
         write(6,*) 'final point:',position
         write(6,*) 'cutoff',rcut
         stop
      endif

c     --- output the final point
      write(30,*) rcut,0.0d0

      end
