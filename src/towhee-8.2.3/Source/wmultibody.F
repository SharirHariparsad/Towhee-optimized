#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2013 Marcus G. Martin                       *
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
      subroutine twh_wmultibody(lfailure,rijsq,ntii,ntjj,ichain,jchain
     &     ,wmulti)
c     ******************************************************************
c     * Computes the intermolecular pair virial function divided by r^2*
c     * as defined by equation 2.60 in Allen and Tildesley             *
c     * w(r) = r dv(r)/dr                                              *
c     * w(r)/r^2 = (dv(r)/dr) /r                                       *
c     * does this for the embedded atom multibody portion using the    *
c     * chain rule d/dr F(rho(r)) = dF/drho * drho/dr                  *
c     * have to do this twice, once for i and once for -j              *
c     *                                                                *
c     * originally written 03-25-2005 by M.G. Martin                   *
c     * last update 08-17-2013 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_LT
#define FUNCTION_EXPON
#define FUNCTION_EXTRACTDENS
#define FUNCTION_EXTRACTEMBED
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure
      integer ntii,ntjj,ichain,jchain
      double precision rijsq,wmulti
c     --- local variables
c     --- logical scalars
      logical lcontinue
c     --- character string scalars
      character*20 eamembedstyle,eamdensstyle
c     --- integer scalars
      integer itrial,icompute,icoeff,trialchain,iother,cpnumber
      integer eamindex
c     --- double precision scalars
      double precision x,rij,density,eamdensone,eamdenstwo
      double precision rhoone,rhotwo,maxdist,densdiff,eamalpha,eambeta
      double precision rhozero,rhomin,bela,belb,belc,belm,beln
      double precision densratio
c     --- double precision arrays
      double precision embedone,embedtwo,embedthr
      double precision embedprime,rhoprime
      dimension embedprime(2),rhoprime(2)
      double precision belrho
      dimension belrho(0:9)

c      write(6,*) 'starting WMULTIBODY'

c     --- initialize wmultibody
      wmulti = 0.0d0

      call twh_cp_number(GLB_GET,cpnumber)
      if ( cpnumber .eq. CP_EMBEDDED_ATOM_METHOD ) then
c        --- compute rij
         rij = dsqrt(rijsq)
c        --- the embedding and density parts
         do itrial = 1,2
            rhoprime(itrial) = 0.0d0
            embedprime(itrial) = 0.0d0
            if ( itrial .eq. 1 ) then
               icompute = ntii
               iother = ntjj
               trialchain = ichain
            else
               icompute = ntjj
               iother = ntii
               trialchain = jchain
            endif
c           --- density derivative
            call twh_eam_dens_style(GLB_GET,iother,icompute
     &           ,eamdensstyle)
            if ( eamdensstyle .eq. 'Ackland cubic sum' )then
c              --- cubic sum with heaviside
               call twh_eam_ndens(GLB_GET,iother,icompute,icoeff)
               if ( icoeff .gt. 0 ) then
                  lcontinue = .true.
               else
                  lcontinue = .false.
               endif
               do while ( lcontinue )
                  call twh_eam_dens(GLB_GET,1,iother,icompute,icoeff
     &                 ,eamdensone)
                  x = eamdensone - rij
                  call twh_eam_dens(GLB_GET,2,iother,icompute,icoeff
     &                 ,eamdenstwo)
                  rhoprime(itrial) = rhoprime(itrial) +
     &                 eamdenstwo*(-3.0d0)*x*x
                  icoeff = icoeff - 1
                  if ( icoeff .eq. 0 ) then
                     lcontinue = .false.
                  else
                     call twh_eam_dens(GLB_GET,1,iother,icompute,icoeff
     &                    ,eamdensone)
                     if ( eamdensone .le. rij ) lcontinue = .false.
                  endif
               enddo
            elseif ( eamdensstyle .eq. 'exponential' ) then
c              --- first check max density cutoff
               call twh_eam_dens(GLB_GET,1,iother,icompute,2,maxdist)
               if ( rij .le. maxdist ) then
c                 --- pone*ptwo*exp(ptwo*r)
                  call twh_eam_dens(GLB_GET,2,iother,icompute,1,rhoone)
                  call twh_eam_dens(GLB_GET,2,iother,icompute,2,rhotwo)
                  rhoprime(itrial) = rhoone*rhotwo*twh_expon(rhotwo*rij)
               endif
            elseif ( eamdensstyle .eq. 'table' ) then
               rhoprime(itrial) = twh_extractdens(lfailure,INTERP_DERIV
     &              ,iother,icompute,rij)
               if ( lfailure ) return
            elseif ( eamdensstyle .eq. 'none' ) then
c              --- no density contribution from this pair interaction
               rhoprime(itrial) = 0.0d0
            else
               write(6,*) 'WMULTIBODY: unknown eam_dens_style'
               write(6,*) eamdensstyle
               lfailure = .true.
               return
            endif
c           --- embedding derivative
c           --- set the density
            call twh_eam_rho(GLB_GET,CRD_REAL,trialchain,density)
            call twh_eam_embed_style(GLB_GET,icompute,eamembedstyle)
            if ( eamembedstyle .eq. 'power 0.5 and 2' ) then
c              --- power series with 0.5, and 2
               if ( density .ne. 0.0d0 ) then
                  call twh_eam_embed(GLB_GET,2,icompute,1,embedone)
                  call twh_eam_embed(GLB_GET,2,icompute,2,embedtwo)
                  embedprime(itrial) = (0.5d0)*embedone/dsqrt(density) 
     &                 +(2.0d0)*embedtwo*density
               endif
            elseif ( eamembedstyle .eq. 'power 0.5, 2, and 4' ) then
c              --- power series with 0.5, 2, and 4 powers
               if ( density .ne. 0.0d0 ) then
                  call twh_eam_embed(GLB_GET,2,icompute,1,embedone)
                  call twh_eam_embed(GLB_GET,2,icompute,2,embedtwo)
                  call twh_eam_embed(GLB_GET,2,icompute,3,embedthr)
                  embedprime(itrial) =
     &                 (0.5d0)*embedone/dsqrt(density)
     &                 +(2.0d0)*embedtwo*density
     &                 +(4.0d0)*embedthr*(density**3)
               endif
            elseif ( eamembedstyle .eq. 'Belashchenko' ) then
c              --- depends on the value of rho
               call twh_eam_embed(GLB_GET,1,icompute,5,rhomin)
               if ( density .lt. rhomin ) then
c                 --- rho < rhomin alpha/(2*sqrt(rho)) + beta
                  call twh_eam_embed(GLB_GET,1,ntii,4,eamalpha)
                  call twh_eam_embed(GLB_GET,2,ntii,4,eambeta)
                  embedprime(itrial) = eamalpha/(2.0d0*dsqrt(density))
     &                 + eambeta
               else
c                 --- 2atwo(rho-rho0) + 3athree(rho-rho0)^2
                  call twh_eam_embed(GLB_GET,2,ntii,5,rhozero)
                  call twh_eam_embed(GLB_GET,2,ntii,2,embedtwo)
                  call twh_eam_embed(GLB_GET,2,ntii,3,embedthr)
                  densdiff = density - rhozero
                  embedprime(itrial) = (2.0d0)*embedtwo*(densdiff)
     &                 + (3.0d0)*embedthr*densdiff*densdiff
               endif
            elseif ( eamembedstyle .eq. 'Belashchenko 10 rho' ) then
c              --- Belaschenko piecewise potential based on the density
               call twh_eam_embed(GLB_GET,1,ntii,11,belrho(0))
               call twh_eam_embed(GLB_GET,2,ntii,11,belrho(1))
               call twh_eam_embed(GLB_GET,3,ntii,11,belrho(2))
               call twh_eam_embed(GLB_GET,1,ntii,12,belrho(3))
               call twh_eam_embed(GLB_GET,2,ntii,12,belrho(4))
               call twh_eam_embed(GLB_GET,3,ntii,12,belrho(5))
               call twh_eam_embed(GLB_GET,1,ntii,13,belrho(6))
               call twh_eam_embed(GLB_GET,2,ntii,13,belrho(7))
               call twh_eam_embed(GLB_GET,3,ntii,13,belrho(8))
               call twh_eam_embed(GLB_GET,1,ntii,14,belrho(9))
c              --- go small to large density
               if ( twh_cmp_lt(density,belrho(7)) ) then
c                 --- 0 to belrho(7)
                  call twh_eam_embed(GLB_GET,1,ntii,8,bela)
                  call twh_eam_embed(GLB_GET,2,ntii,8,belb)
                  call twh_eam_embed(GLB_GET,3,ntii,8,belc)
                  densdiff = density-belrho(7)
                  densratio = density/belrho(7)
c                 --- chain rule, 2 parts
                  embedprime(itrial) = ((belb + 2.0d0*belc*densdiff)
     &                 *(2.0d0*densratio - densratio*densratio))
     &                 + ((bela + belb*densdiff
     &                 + belc*densdiff*densdiff)
     &                 *(2.0d0/belrho(7) - 2.0d0*density/belrho(7)**2)
     &                 )
               elseif ( twh_cmp_lt(density,belrho(1)) ) then
c                 --- somewhere in belrho(7) to belrho(1)
c                 --- this group uses the same function, find the
c                 --- subinterval
                  if ( twh_cmp_lt(density,belrho(6)) ) then
                     eamindex = 7
                  elseif ( twh_cmp_lt(density,belrho(5)) ) then
                     eamindex = 6
                  elseif ( twh_cmp_lt(density,belrho(4)) ) then
                     eamindex = 5
                  elseif ( twh_cmp_lt(density,belrho(3)) ) then
                     eamindex = 4
                  elseif ( twh_cmp_lt(density,belrho(2)) ) then
                     eamindex = 3
                  else
                     eamindex = 2
                  endif
                  densdiff = density - belrho(eamindex-1)
                  call twh_eam_embed(GLB_GET,2,ntii,eamindex,belb)
                  call twh_eam_embed(GLB_GET,3,ntii,eamindex,belc)
                  embedprime(itrial) = belb + 2.0d0*belc*densdiff
               elseif ( twh_cmp_lt(density,belrho(8)) ) then
c                 --- this interval doesn't need a b and refers to
c                 --- another central rho that is between belrho(1)
c                 ---- and belrho(8)
                  call twh_eam_embed(GLB_GET,3,ntii,1,belc)
                  densdiff = density - belrho(0)
                  embedprime(itrial) = 2.0d0*belc*densdiff
               elseif ( twh_cmp_lt(density,belrho(9)) ) then
c                 --- between belrho(8) and belrho(9), getting to high
c                 --- density
                  call twh_eam_embed(GLB_GET,2,ntii,9,belb)
                  call twh_eam_embed(GLB_GET,3,ntii,9,belc)
                  call twh_eam_embed(GLB_GET,2,ntii,14,belm)
                  densdiff = density - belrho(8)
                  embedprime(itrial) = belb
     &                 + belm*belc*(densdiff**(belm-1.0d0))
               else
c                 --- very high density, beyond belrho(9)
                  call twh_eam_embed(GLB_GET,2,ntii,10,belb)
                  call twh_eam_embed(GLB_GET,3,ntii,10,belc)
                  call twh_eam_embed(GLB_GET,3,ntii,14,beln)
                  densdiff = density - belrho(9)
                  embedprime(itrial) = belb 
     &                 + beln*belc*(densdiff**(beln-1.0d0))
               endif
            elseif ( eamembedstyle .eq. 'table' ) then
               embedprime(itrial) = twh_extractembed(lfailure
     &              ,INTERP_DERIV,icompute,density)
               if ( lfailure ) return
            else
               write(6,*) 'WMULTIBODY: unknown eam_embed_style'
               write(6,*) eamembedstyle
               lfailure = .true.
               return
            endif
         enddo
c        --- compute wmultibody (i terms + jterms)/r
         wmulti = ( rhoprime(1)*embedprime(1) 
     &        + rhoprime(2)*embedprime(2) )/rij
      endif

      return
      end
