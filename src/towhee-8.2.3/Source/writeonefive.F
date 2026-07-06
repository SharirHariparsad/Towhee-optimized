#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2020 Marcus G. Martin                       *
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
      subroutine twh_writeonefive(lfailure,index)
c     ******************************************************************
c     * writes all of the information about the one-five interactions  *
c     * to the file index                                              *
c     *                                                                *
c     * originally written 12-15-2004 by M.G. Martin                   *
c     * last modified 03-30-2020 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LOFHERE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer index
c     --- local variables
      logical lfound
      character*(FFNAMELEN) ofnamesone,ofnamesfiv
      integer type,ofstyle,icoeff
      double precision ofcoeff
      dimension ofcoeff(MINOFCOEFF:MAXOFCOEFF)

c     --- see if there are any special one-five types
      lfound = .false.
      do type = 1,TOFMAX
         if ( twh_get_lofhere(type) )  then
            if ( .not. lfound ) then
c              --- this is the first found interaction - write the intro
               write(index,'(a22)') 'Special One-Five Types'
               lfound = .true.
            endif
            call twh_ofnames(GLB_GET,type,1,ofnamesone)
            call twh_ofnames(GLB_GET,type,5,ofnamesfiv)
            call twh_ofstyle(GLB_GET,type,ofstyle)
            do icoeff = MINOFCOEFF,MAXOFCOEFF
               call twh_ofcoeff(GLB_GET,type,icoeff,ofcoeff(icoeff))
            enddo
            if ( ofstyle .eq. 1 ) then
               write(index,10) type,ofnamesone,ofnamesfiv
     &              ,ofcoeff(1),ofcoeff(2),ofcoeff(0)
            elseif ( ofstyle .eq. 2 ) then
               write(index,20) type,ofnamesone,ofnamesfiv,ofcoeff(1)
            else
               write(index,*) 'WRITEONEFIVE: unknown ofstyle '
     &              ,ofstyle,' for type ',type
               lfailure = .true.
               return
            endif
         endif
      enddo

 10   format('Type:',1x,i3,1x,'Atom(i):',1x,a10,1x,'Atom(j):',1x,a10,1x
     &     ,'Sigma:',1x,f7.3,1x,'Epsilon:',1x,f7.3,1x,'Shift:',1x,f7.4)
 20   format('Type:',1x,i3,1x,'Atom(i):',1x,a10,1x,'Atom(j):',1x,a10,1x
     &     ,'A/r^12:',1x,e10.4)

      return
      end
    
