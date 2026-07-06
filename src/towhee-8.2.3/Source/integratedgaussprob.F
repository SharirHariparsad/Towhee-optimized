#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2016 Marcus G. Martin                       *
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
      subroutine twh_integratedgaussprob(lfailure,lolimstyle,hilimstyle
     &     ,lovalue,hivalue,mean,stddev,lolimit,hilimit,probability)
c     ******************************************************************
c     * returns the integrated probability that a value falls in the   *
c     * range of lovalue to hivalue inside of a Gaussian distribution  *
c     * with the mean, standard deviation, and limits passed into this *
c     * function                                                       *
c     *                                                                *
c     * originally written 10-06-2005 by M.G. Martin                   *
c     * last modified 07-11-2016 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DERF
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure
      integer lolimstyle,hilimstyle
      double precision lovalue,hivalue,mean,stddev,lolimit,hilimit
c     --- local variables
      double precision hilimitscale,lolimitscale,scale
      double precision lovaluescale,hivaluescale,probability
      double precision hilimint,lolimint

c     --- set the scale
      scale = 1.0d0/(dsqrt(2.0d0)*stddev)
c     --- determine integrated lo limit
      if ( lolimstyle .eq. LIM_VALUE ) then
c        --- check lo range
         if ( lolimit .gt. lovalue ) then
            write(6,*) 'INTEGRAGEDGAUSSPROB: invalid lo range'
            write(6,*) 'lolimit,lovalue:',lolimit,lovalue
            lfailure = .true.
            return
         endif
c        --- determine integrated probabilty of the low limit
         lolimitscale = (lolimit-mean)*scale
         lolimint = twh_derf(lolimitscale)
      elseif ( lolimstyle .eq. LIM_NEGINF ) then
c        --- error function of negative infinity is zero
         lolimint = 0.0d0
      else
         write(6,*) 'INTEGRATEDGAUSSPROB: unknown lolimstyle:'
     &        ,lolimstyle
         lfailure = .true.
         return
      endif

      if ( hilimstyle .eq. LIM_VALUE ) then
c        --- check hi range
         if ( hilimit .lt. hivalue ) then
            write(6,*) 'INTEGRAGEDGAUSSPROB: invalid range'
            write(6,*) 'hivalue,hilimit:',hivalue,hilimit
            lfailure = .true.
            return
         endif
c        --- determine integrated probability of the hi limit
         hilimitscale = (hilimit-mean)*scale
         hilimint = twh_derf(hilimitscale)
      elseif ( hilimstyle .eq. LIM_POSINF ) then
c        --- error function of positive infinity is one
         hilimint = 1.0d0
      else
         write(6,*) 'INTEGRATEDGAUSSPROB: unknown hilimstyle'
     &        ,hilimstyle
         lfailure = .true.
         return
      endif

      hivaluescale = (hivalue-mean)*scale
      lovaluescale = (lovalue-mean)*scale
      probability = 
     &     (twh_derf(hivaluescale) -twh_derf(lovaluescale))
     &     /(hilimint  - lolimint)

      return
      end
