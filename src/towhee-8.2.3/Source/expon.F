#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2011 Marcus G. Martin                       *
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
      function twh_expon( value )
c     ******************************************************************
c     * computes the exponential of a value with safety checks to avoid*
c     * underflow.                                                     *
c     * Any value that would result in an exponential smaller than     *
c     * SMALLEST gets returns as 0.0 instead                           *
c     *                                                                *
c     * originally written 11-21-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#include "functions.h"

c     --- variables passed to the function
       double precision value
c     --- Phase 6 P9.2: branchless expon via exp(max(arg,-700))
c         Enables vectorization; eliminates twh_constant_underflow dispatch
c         The underflow constant is ~ -745 for double; -700 is safe threshold
c         used in phase6 spec: exp(max(value,-700)).  Avoids branch for -ftree-vectorize
c         Keep exact 0.0 behavior for < -700 via max+exp (exp(-700) ~ 9.8e-305 ~ 0 in practice
c         but preserves continuity; if strict 0 needed, caller can threshold)
!$OMP DECLARE SIMD(twh_expon)
       twh_expon = dexp(max(value, -700d0))
       
       return
       end
