#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2020 Marcus G. Martin                       *
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
      subroutine twh_radialpressure(lfailure,pressure,ibox)
c     ******************************************************************
c     * computes the pressure from a snapshot by using the radial      *
c     * distribution funtion definition of pressure for hard sphere or *
c     * square well systems.  Does this in a non-traditional manner by *
c     * looking at the current radial distribution function and using  *
c     * a bin width to get an estimate of the various rdf values       *
c     *                                                                *
c     * The algorithm implemented here is correct in the limit that    *
c     * delta goes to zero and the molecules are monatomic.  Anyone's  *
c     * guess as to how well it works away from those limits, but I    *
c     * implemented it as a curiosity for the discontinuous well       *
c     * potentials (MGM 06-27-2005)                                    *
c     *                                                                *
c     * originally written 06-27-2005 by M.G. Martin                   *
c     * last modified 03-31-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ibox
      double precision pressure
c     --- parameters
      integer nwellmax
      parameter (nwellmax=5)
c     --- local variables
c     --- logical scalars
      logical lelect
c     --- integer scalars
      integer ichain,jchain,ilist,jlist,jstart,ivalue,iwell,icoeff
      integer ii,jj,ntii,ntjj,ntij,imolty,jmolty,nmolty
      integer numiunit,numjunit,ncmti,ncmtj
      integer nwells,ensemble,cpnumber
c     --- double precision scalars
      double precision rxui,ryui,rzui,rxuj,ryuj,rzuj,rx_ab,ry_ab,rz_ab
      double precision rijsq,rdfconst,convertfact,boxvolume
      double precision densfact,pcfactor,onepi,temperature
      double precision radial_pressure_delta
c     --- double precision arrays
      double precision countwell
      dimension countwell(2,nwellmax,CROSSTYPEMAX)
      double precision wellbounds
      dimension wellbounds(3,nwellmax,CROSSTYPEMAX)
      double precision atcount
      dimension atcount(NNTYPE)
      double precision gofr
      dimension gofr(2)
      double precision nbcoeff
      dimension nbcoeff(MINNBCOEFF:MAXNBCOEFF)

c     --- initialize pressure
      pressure = 0.0d0
c     --- set the converstion factor to go from K molecule / A^3 to kPa
      convertfact = ((1.0d30*8.314d0)/(6.02d23*1000.0d0))
c     --- retrieve constants
      call twh_radial_pressure_delta(GLB_GET,radial_pressure_delta)

c     --- safety check to make sure we don't have coulomb
      call twh_lelect(GLB_GET,lelect)
      if ( lelect ) then
         write(6,*) 'RADIALPRESSURE: not functional for combination of'
         write(6,*) 'hard core potentials with coulombic interactions'
         write(6,*) 'please disable the pressure calculation for this'
     &        ,' system'
         return
      endif

c     --- no pressure in the ideal box of grand canonical
      call twh_ensemble(GLB_GET,ensemble)
      if ( ensemble .eq. ENS_muVT .and. ibox .eq. 0 ) return

c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_temperature(GLB_GET,temperature)
      call twh_boxvolume(GLB_GET,BXD_REAL,ibox,boxvolume)

c     --- initialize the countwell array
      do ivalue = 1,2
         do iwell = 1,nwellmax
            do ntij = 1,CROSSTYPEMAX
               countwell(ivalue,iwell,ntij) = 0.0d0
            enddo
         enddo
      enddo
c     --- initialize the atomcount array
      do ntii = 1,NNTYPE
         atcount(ntii) = 0.0d0
      enddo

c     --- set nwell and well boundaries depending upon the classical
c     --- potential
      call twh_cp_number(GLB_GET,cpnumber)
      if ( cpnumber .eq. CP_HARD_SPHERE ) then
         iwell = 0
c        --- set the hard sphere boundary
         iwell = iwell + 1
         do ntij = 1,CROSSTYPEMAX
            do icoeff = 1,2
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            wellbounds(1,iwell,ntij) = ( nbcoeff(1) 
     &           - radial_pressure_delta )**2
            wellbounds(2,iwell,ntij) = nbcoeff(2)
            wellbounds(3,iwell,ntij) = ( nbcoeff(1) 
     &           + radial_pressure_delta )**2
         enddo
c        --- set nwells
         nwells = iwell
      elseif ( cpnumber .eq. CP_SQUARE_WELL ) then
         iwell = 0
c        --- set the hard sphere boundary
         iwell = iwell + 1
         do ntij = 1,CROSSTYPEMAX
            do icoeff = 1,2
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            wellbounds(1,iwell,ntij) = ( dsqrt(nbcoeff(1))
     &           - radial_pressure_delta )**2
            wellbounds(2,iwell,ntij) = nbcoeff(1)
            wellbounds(3,iwell,ntij) = ( dsqrt(nbcoeff(1)) 
     &           + radial_pressure_delta )**2
         enddo
c        --- set the attractive well boundary
         iwell = iwell + 1
         do ntij = 1,CROSSTYPEMAX
            call twh_nbcoeff(GLB_GET,2,ntij,nbcoeff(2))
            wellbounds(1,iwell,ntij) = ( dsqrt(nbcoeff(2)) 
     &           - radial_pressure_delta )**2
            wellbounds(2,iwell,ntij) = nbcoeff(2)
            wellbounds(3,iwell,ntij) = ( dsqrt(nbcoeff(2)) 
     &           + radial_pressure_delta )**2
         enddo
c        --- set nwells
         nwells = iwell
      elseif ( cpnumber .eq. CP_HARD_2580_MULTISTEP ) then
c        --- hard wall plus 4 wells
         do ntij = 1,CROSSTYPEMAX
            do icoeff = 1,11
               call twh_nbcoeff(GLB_GET,icoeff,ntij,nbcoeff(icoeff))
            enddo
            iwell = 0
c           --- hard boundary
            iwell = iwell + 1
            wellbounds(1,iwell,ntij) = ( nbcoeff(1) 
     &           - radial_pressure_delta )**2
            wellbounds(2,iwell,ntij) = nbcoeff(7)
            wellbounds(3,iwell,ntij) = ( nbcoeff(1) 
     &           + radial_pressure_delta )**2
c           --- 1.2 sigma boundary
            iwell = iwell + 1
            wellbounds(1,iwell,ntij) = ( (1.2d0)*nbcoeff(1)
     &           - radial_pressure_delta )**2
            wellbounds(2,iwell,ntij) = nbcoeff(8)
            wellbounds(3,iwell,ntij) = ( (1.2d0)*nbcoeff(1)
     &           + radial_pressure_delta )**2
c           --- 1.5 sigma boundary
            iwell = iwell + 1
            wellbounds(1,iwell,ntij) = ( (1.5d0)*nbcoeff(1)
     &           - radial_pressure_delta )**2
            wellbounds(2,iwell,ntij) = nbcoeff(9)
            wellbounds(3,iwell,ntij) = ( (1.5d0)*nbcoeff(1)
     &           + radial_pressure_delta )**2
c           --- 1.8 sigma boundary
            iwell = iwell + 1
            wellbounds(1,iwell,ntij) = ( (1.8d0)*nbcoeff(1)
     &           - radial_pressure_delta )**2
            wellbounds(2,iwell,ntij) = nbcoeff(10)
            wellbounds(3,iwell,ntij) = ( (1.8d0)*nbcoeff(1)
     &           + radial_pressure_delta )**2
c           --- 2.0 sigma boundary
            iwell = iwell + 1
            wellbounds(1,iwell,ntij) = ( (2.0d0)*nbcoeff(1)
     &           - radial_pressure_delta )**2
            wellbounds(2,iwell,ntij) = nbcoeff(11)
            wellbounds(3,iwell,ntij) = ( (2.0d0)*nbcoeff(1)
     &           + radial_pressure_delta )**2
            nwells = iwell
         enddo
      else
         write(6,*) 'RADIALPRESSURE: unknown classical_potential'
         lfailure = .true.
         return
      endif

c     --- loop over all molecules that are in this box
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,numiunit)
         call twh_ncmt(GLB_GET,ibox,imolty,ncmti)
         do ilist = 1,ncmti
            call twh_chainlist(GLB_GET,ilist,ibox,imolty,ichain)
c           --- loop over all beads of ichain  
            do ii = 1, numiunit
               call twh_ntype(GLB_GET,imolty,ii,ntii)
c              --- increase atom type count
               atcount(ntii) = atcount(ntii) + 1.0d0
               call twh_coordinates(GLB_GET,CRD_REAL,ichain,ii,rxui,ryui
     &              ,rzui)
c              --- loop over all other possible molecule types that 
c              --- have an index equal to or higher than the current 
c              --- molecule type
               do jmolty = imolty,nmolty
                  call twh_nunit(GLB_GET,jmolty,numjunit)
                  if ( imolty .eq. jmolty ) then
c                    --- avoid counting twice
                     jstart = ilist + 1
                  else
c                    --- different molecule types so cross is avoided at
c                    --- that higher level
                     jstart = 1
                  endif
                  call twh_ncmt(GLB_GET,ibox,jmolty,ncmtj)
                  do jlist = jstart,ncmtj
                     call twh_chainlist(GLB_GET,jlist,ibox,jmolty
     &                    ,jchain)
c                    --- loop over all beads of jchain  
                     do jj = 1, numjunit
                        call twh_ntype(GLB_GET,jmolty,jj,ntjj)
c                       --- set cross term type
                        ntij = twh_getnbtype(ntii,ntjj)
                        call twh_coordinates(GLB_GET,CRD_REAL,jchain,jj
     &                       ,rxuj,ryuj,rzuj)
c                       --- minimum image the pair separations
                        call twh_mimage(BXD_REAL,ibox,rxui,ryui,rzui
     &                       ,rxuj,ryuj,rzuj,rx_ab,ry_ab,rz_ab )
                        rijsq = rx_ab*rx_ab+ry_ab*ry_ab+rz_ab*rz_ab
                        do iwell = 1,nwells
c                          --- check lower boundary
                           if( wellbounds(1,iwell,ntij) .lt. rijsq .and.
     &                          rijsq .lt. wellbounds(2,iwell,ntij))then
                              countwell(1,iwell,ntij) 
     &                             = countwell(1,iwell,ntij) + 1.0d0
                           endif
c                          --- check upper boundary
                           if( wellbounds(2,iwell,ntij) .lt. rijsq .and.
     &                          rijsq .lt. wellbounds(3,iwell,ntij))then
                              countwell(2,iwell,ntij) 
     &                             = countwell(2,iwell,ntij) + 1.0d0
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

c     --- compute the ideal pressure from the molecule number density
c     --- add this ideal pressure term (kT N/V)
      do imolty = 1,nmolty
         call twh_ncmt(GLB_GET,ibox,imolty,ncmti)
         pressure = pressure + (temperature*dble(ncmti)/boxvolume)
     &        *convertfact
      enddo

      do ntii = 1,NNTYPE
         if ( atcount(ntii) .ne. 0 ) then
            do ntjj = ntii,NNTYPE
               if ( atcount(ntjj) .ne. 0 ) then
                  ntij = twh_getnbtype(ntii,ntjj)
                  if ( ntii .eq. ntjj ) then
c                    --- identical types, total const (N*M-1)/2
                     rdfconst = (atcount(ntii)*atcount(ntjj)-1.0d0)
     &                    /2.0d0
c                    --- just an i-i so pressure count factor is 1.0
                     pcfactor = 1.0d0
                  else
c                    --- different types so total is N*M
                     rdfconst = atcount(ntii)*atcount(ntjj)
c                    --- need to account for i-j and j-i so double
                     pcfactor = 2.0d0
                  endif
c                 --- include the 4/3 Pi portion
                  rdfconst = rdfconst*4.0d0*onepi/3.0d0
                  do iwell=1,nwells
                     do ivalue = 1,2
                        gofr(ivalue) = boxvolume*
     &                       countwell(ivalue,iwell,ntij)/(rdfconst
     &                       *(wellbounds(ivalue+1,iwell,ntij)**(1.5d0)
     &                       -wellbounds(ivalue,iwell,ntij)**(1.5d0)))
                     enddo
c                    --- densities factor
                     densfact = atcount(ntii)*atcount(ntjj)
     &                    /(boxvolume**2)
c                    --- add the non-ideal pressure term
                     pressure = pressure - pcfactor*(2.0d0*onepi/3.0d0)
     &                    *wellbounds(2,iwell,ntij)**(1.5d0)
     &                    *(gofr(1)-gofr(2))*densfact*convertfact
                  enddo
               endif
            enddo
         endif
      enddo

      return
      end
