#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2020 Marcus G. Martin                       *
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
      function twh_vthreebody(lfailure,rij,rjk,thetac,ntii,ntjj,ntkk )
c     ******************************************************************
c     * computes the three-body term between non-bonded atoms in       *
c     * certain force fields                                           *
c     *                                                                *
c     * originally written 06-26-2002 by  M.G. Martin                  *
c     * last update 03-31-2020 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_GETNBTYPE
#define FUNCTION_VTHREEBODY
#include "functions.h"
c     --- local defines for logic
#define LOCAL_IJ 1
#define LOCAL_JJ 2
#define LOCAL_JK 3
c     --- variables passed to/from the function
      logical lfailure
      integer ntii,ntjj,ntkk
      double precision rij,rjk,thetac
c     --- local variables
c     --- integer scalars
      integer ij,jk,jj,cpnumber,icoeff
c     --- double precision arrays
      double precision nbcoeff
      dimension nbcoeff(MINNBCOEFF:MAXNBCOEFF,3)
      
c     --- initialize vthreebody
      twh_vthreebody = 0.0d0
      call twh_cp_number(GLB_GET,cpnumber)
      if ( cpnumber .eq. CP_STILLINGER_WEBER ) then
c        --- set cross types
         ij = twh_getnbtype(ntii,ntjj)
         do icoeff = 1,9
            call twh_nbcoeff(GLB_GET,icoeff,ij,nbcoeff(icoeff,LOCAL_IJ))
         enddo
         jk = twh_getnbtype(ntjj,ntkk)
         do icoeff = 1,9
            call twh_nbcoeff(GLB_GET,icoeff,jk,nbcoeff(icoeff,LOCAL_JK))
         enddo
c        --- set self type of the center atom for use with some of the
c        --- ambiguous constants
         jj = twh_getnbtype(ntjj,ntjj)
         do icoeff = 1,9
            call twh_nbcoeff(GLB_GET,icoeff,jj,nbcoeff(icoeff,LOCAL_JJ))
         enddo
c        --- Stillinger-Weber three body potential
c        --- nb 1, 7, and 9 use the self terms of the central atom
c        --- all other constants depend on the cross type
         twh_vthreebody = nbcoeff(1,LOCAL_JJ)*nbcoeff(7,LOCAL_JJ)
     &        *(twh_expon(nbcoeff(8,LOCAL_IJ)
     &        /(rij/nbcoeff(4,LOCAL_IJ) - nbcoeff(6,LOCAL_IJ))
     &        + nbcoeff(8,LOCAL_JK)/(rjk/nbcoeff(4,LOCAL_JK)
     &        - nbcoeff(6,LOCAL_JK)))
     &        * (thetac - nbcoeff(9,LOCAL_JJ))**2)
      else
         write(6,*) 'VTHREEBODY: invalid classical_potential'
     &        ,cpnumber
         lfailure = .true.
      endif

c      write(6,*) 'finished VTHREEBODY'

      return
      end
