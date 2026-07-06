#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2016 Marcus G. Martin,                      *
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
      subroutine twh_putcartinbox (bxdflag,ibox,rx,ry,rz)
c     ******************************************************************
c     * puts the atom inside the central simulation box                *
c     *                                                                *
c     * originally written 12-11-2001 by M.G. Martin                   *
c     * last modified 07-25-2016 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      integer bxdflag,ibox
      double precision rx,ry,rz
c     --- local variables
      integer shiftmax
      parameter (shiftmax = 10)
      double precision arba,arbb,arbc

c     --- transform into arbitrary coordinate system
      call twh_carttoarb(bxdflag,ibox,rx,ry,rz,arba,arbb,arbc)

c     --- put back into central box
      call twh_putarbinbox(arba,arbb,arbc)

c     --- transform back to cartesian coordinates
      call twh_arbtocart(bxdflag,ibox,arba,arbb,arbc,rx,ry,rz)

      return
      end
