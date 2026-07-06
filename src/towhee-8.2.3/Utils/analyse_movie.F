c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * see the file license.gpl for the full license information      *
c     * Copyright (C) 2001-2005 Marcus G. Martin                       *
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
      program analyse_movie
c     ******************************************************************
c     * This program analyses the towhee_movie file and outputs        *
c     * 2-dimensional atomic density profiles.                         *
c     * Written by J.A. Greathouse on 03-23-2010.                      *
c     * Based on analyse_movie.F                                       *
c     ******************************************************************
      implicit none
c     --- values passed to the main routine
      integer nframe,nchain,nmolty,nbox,numax,tor_max,nnbond
c     --- local variables
      logical lanalyse
      integer imolty,nunit,intor,movieversion,iunit,itor,nvib,nbin
      integer tor_sum,tor_smax

      write(6,*) 'Opening towhee_movie'
c     --- open the movie file
      open(10,file='towhee_movie',form='formatted')

c     --- read in the movie version and check it
      read(10,*) movieversion
      write(6,*) 'towhee_movie input version ',movieversion
      if ( movieversion .ne. 2 ) then
         write(6,*) 'invalid movieversion'
         write(6,*) 'you are probably trying to read in a version '
         write(6,*) 'generated prior to towhee version 4.7.4'
         write(6,*) 'That old version is no longer supported'
         stop
      endif

c     --- scan the file to determine the array dimensions that are 
c     --- required
      read(10,*) nframe,nchain,nmolty,nbox

c     --- scan to find the maximum value of nunit
      numax = 0
      nnbond = 0
      tor_max = 0
      tor_sum = 0
      tor_smax = 0
      do imolty = 1,nmolty
         read(10,*) nunit
         if ( nunit .gt. numax ) numax = nunit
c        --- set nnbond
         do iunit = 1,nunit
            read(10,*) nvib
            if ( nvib .gt. nnbond ) nnbond = nvib
         enddo
c        --- set tor_max
         do iunit = 1,nunit
            read(10,*) intor
            if ( intor .gt. tor_max ) tor_max = intor
            tor_sum = tor_sum + intor
            do itor = 1,intor
               read(10,*) 
            enddo
         enddo
         if ( tor_sum .gt. tor_smax ) tor_smax = tor_sum
      enddo

      call reset_movie()

c     --- radial distribution functions
      write(6,*) 'Please answer the questions T (true) or F (false)'
      write(6,*) 'Compute intermolecular'
     &     ,' radial distribution functions? (T/F)'
      read(5,*) lanalyse
      if ( lanalyse ) then
         write(6,*) 'Input number of bins for rdf'
         read(5,*) nbin
         if ( nbin .le. 0 ) then
            write(6,*) 'number of bins must be positive'
            stop
         endif
         call analyse_rdf(nframe,nchain,nbox,nmolty,numax,nbin
     &        ,nnbond,tor_max)
         call reset_movie()
      endif

      write(6,*) 'Compute atom distribution profiles? (T/F)'
      read(5,*) lanalyse
      if ( lanalyse ) then
         write(6,*) 'Input number of bins for atom dist profile'
         read(5,*) nbin
         if ( nbin .le. 0 ) then
            write(6,*) 'number of bins must be positive'
            stop
         endif
         call analyse_profile(nframe,nchain,nbox,nmolty,numax,nbin
     &        ,nnbond,tor_max)
         call reset_movie()
      endif


      write(6,*) 'Compute the vibration distribution? (T/F)'
      read(5,*) lanalyse
      if ( lanalyse ) then
         write(6,*) 'Input number of bins for vibration profile'
         read(5,*) nbin
         if ( nbin .le. 0 ) then
            write(6,*) 'number of bins must be positive'
            stop
         endif
         call analyse_vib(nframe,nchain,nbox,nmolty,numax,nbin
     &        ,nnbond,tor_max)
         call reset_movie()
      endif

      write(6,*) 'Compute the bond angle distribution? (T/F)'
      read(5,*) lanalyse
      if ( lanalyse ) then
         write(6,*) 'Input number of bins for bending profile'
         read(5,*) nbin
         if ( nbin .le. 0 ) then
            write(6,*) 'number of bins must be positive'
            stop
         endif
         call analyse_bend(nframe,nchain,nbox,nmolty,numax,nbin
     &        ,nnbond,tor_max)
         call reset_movie()
      endif

      write(6,*) 'Compute the dihedral angle distribution? (T/F)'
      read(5,*) lanalyse
      if ( lanalyse ) then
         write(6,*) 'Input number of bins for dihedral profile'
         read(5,*) nbin
         if ( nbin .le. 0 ) then
            write(6,*) 'number of bins must be positive'
            stop
         endif
         call analyse_dihed(nframe,nchain,nbox,nmolty,numax,nbin
     &        ,nnbond,tor_max,tor_smax)
         call reset_movie()
      endif

      end

c     --- some subroutines stolen from the main code.  these
c     --- are modified to work with this analysis program
      subroutine carttoarb (nbox,ibox,hinverse
     &     ,rxuij,ryuij,rzuij,arba,arbb,arbc)
c     ******************************************************************
c     * transforms from the cartesian coordinate system to the         *
c     * arbitrary coordinate system (on 0,1)                           *
c     *                                                                *
c     * originally written 12-14-2001 by M.G. Martin                   *
c     * last modified 06-07-2002 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      integer nbox,ibox
      double precision hinverse
      dimension hinverse(nbox,3,3)
      double precision arba,arbb,arbc,rxuij,ryuij,rzuij
c     --- local variables

      arba = rxuij*hinverse(ibox,1,1) 
     &     + ryuij*hinverse(ibox,2,1)
     &     + rzuij*hinverse(ibox,3,1)

      arbb = rxuij*hinverse(ibox,1,2) 
     &     + ryuij*hinverse(ibox,2,2)
     &     + rzuij*hinverse(ibox,3,2)

      arbc = rxuij*hinverse(ibox,1,3) 
     &     + ryuij*hinverse(ibox,2,3)
     &     + rzuij*hinverse(ibox,3,3)

      return
      end


      subroutine arbtocart (nbox,ibox,hmatrix
     &     ,arba,arbb,arbc,rxuij,ryuij,rzuij)
c     ******************************************************************
c     * transforms from the arbitrary coordinate system (on 0,1) to the*
c     * cartesian coordinate system                                    *
c     *                                                                *
c     * originally written 12-13-2001 by M.G. Martin                   *
c     * last modified 06-07-2002 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      integer nbox,ibox
      double precision hmatrix
      dimension hmatrix(nbox,3,3)
      double precision arba,arbb,arbc,rxuij,ryuij,rzuij

      rxuij = hmatrix(ibox,1,1)*arba 
     &     + hmatrix(ibox,2,1)*arbb
     &     + hmatrix(ibox,3,1)*arbc
      ryuij = hmatrix(ibox,1,2)*arba
     &     + hmatrix(ibox,2,2)*arbb
     &     + hmatrix(ibox,3,2)*arbc
      rzuij = hmatrix(ibox,1,3)*arba
     &     + hmatrix(ibox,2,3)*arbb
     &     + hmatrix(ibox,3,3)*arbc

      return
      end


      subroutine mimage (nbox,hmatrix,hinverse,ibox,rxuij,ryuij,rzuij)
c     ******************************************************************
c     * this subroutine determines the minimum distance between any two*
c     * representations of the coordinates through all of the periodic *
c     * images                                                         *
c     *                                                                *
c     * completely rewritten 02-14-2001 by M.G. Martin                 *
c     * last modified 06-07-2002 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      integer ibox,nbox
      double precision hmatrix,hinverse
      dimension hmatrix(nbox,3,3),hinverse(nbox,3,3)
      double precision rxuij,ryuij,rzuij
c     --- local variables
      integer idim
      double precision arbcord
      dimension arbcord(3)

c     --- shift the coordinate system
      call carttoarb(nbox,ibox,hinverse
     &     ,rxuij,ryuij,rzuij,arbcord(1),arbcord(2),arbcord(3))

      do idim = 1,3
c        --- find the minimum image
         if ( arbcord(idim) .gt. 0.5d0 ) then
            arbcord(idim) = arbcord(idim) - 1.0d0
         elseif ( arbcord(idim) .lt. -0.5d0 ) then
            arbcord(idim) = arbcord(idim) + 1.0d0
         endif
      enddo
c     --- transform back to cartesian coordinates
      call arbtocart(nbox,ibox,hmatrix
     &     ,arbcord(1),arbcord(2),arbcord(3),rxuij,ryuij,rzuij)

      return
      end

c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2000 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
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
c     *                                                                *
c     * See the file towhee.F for more information about the code      *
c     ******************************************************************
      subroutine inverthmatrix (nbox,ibox,hmatrix,hinverse,boxvolume)
c     ******************************************************************
c     * inverts the hmatrix for a particular box                       *
c     * also computes the boxvolume and stores it                      *
c     *                                                                *
c     * originally written 12-13-2001 by M.G. Martin                   *
c     * last modified 07-01-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      integer nbox,ibox
      double precision hmatrix,hinverse,boxvolume
      dimension hmatrix(nbox,3,3),hinverse(nbox,3,3)
      dimension boxvolume(nbox)

c     --- local variables
      double precision tolerance,determ
      parameter(tolerance=1.0d-7)

c     --- invert the hmatrix
      hinverse(ibox,1,1) = hmatrix(ibox,2,2)*hmatrix(ibox,3,3) 
     &     - hmatrix(ibox,2,3)*hmatrix(ibox,3,2)
      hinverse(ibox,2,1) = hmatrix(ibox,3,1)*hmatrix(ibox,2,3) 
     &     - hmatrix(ibox,3,3)*hmatrix(ibox,2,1)
      hinverse(ibox,3,1) = hmatrix(ibox,2,1)*hmatrix(ibox,3,2) 
     &     - hmatrix(ibox,2,2)*hmatrix(ibox,3,1)
      hinverse(ibox,1,2) = hmatrix(ibox,3,2)*hmatrix(ibox,1,3) 
     &     - hmatrix(ibox,3,3)*hmatrix(ibox,1,2)
      hinverse(ibox,2,2) = hmatrix(ibox,1,1)*hmatrix(ibox,3,3) 
     &     - hmatrix(ibox,1,3)*hmatrix(ibox,3,1)
      hinverse(ibox,3,2) = hmatrix(ibox,3,1)*hmatrix(ibox,1,2) 
     &     - hmatrix(ibox,3,2)*hmatrix(ibox,1,1)
      hinverse(ibox,1,3) = hmatrix(ibox,1,2)*hmatrix(ibox,2,3) 
     &     - hmatrix(ibox,1,3)*hmatrix(ibox,2,2)
      hinverse(ibox,2,3) = hmatrix(ibox,2,1)*hmatrix(ibox,1,3) 
     &     - hmatrix(ibox,2,3)*hmatrix(ibox,1,1)
      hinverse(ibox,3,3) = hmatrix(ibox,1,1)*hmatrix(ibox,2,2) 
     &     - hmatrix(ibox,1,2)*hmatrix(ibox,2,1)
      
c     --- compute determinant (boxvolume) and make sure it is non-zero
      determ = hmatrix(ibox,1,1)*hinverse(ibox,1,1) 
     &     + hmatrix(ibox,1,2)*hinverse(ibox,2,1)
     &     + hmatrix(ibox,1,3)*hinverse(ibox,3,1)
      boxvolume(ibox) = determ
      if ( dabs(determ) .lt. tolerance ) then
         write(6,*) 'INVERTHMATRIX: zero boxvolume ',determ
         stop
      endif
      
      hinverse(ibox,1,1) = hinverse(ibox,1,1) / determ
      hinverse(ibox,1,2) = hinverse(ibox,1,2) / determ
      hinverse(ibox,1,3) = hinverse(ibox,1,3) / determ
      hinverse(ibox,2,1) = hinverse(ibox,2,1) / determ
      hinverse(ibox,2,2) = hinverse(ibox,2,2) / determ
      hinverse(ibox,2,3) = hinverse(ibox,2,3) / determ
      hinverse(ibox,3,1) = hinverse(ibox,3,1) / determ
      hinverse(ibox,3,2) = hinverse(ibox,3,2) / determ
      hinverse(ibox,3,3) = hinverse(ibox,3,3) / determ

      return
      end


      subroutine analyse_rdf(nframe,nchain,nbox,nmolty,numax,nbin
     &     ,nnbond,tor_max)
c     ******************************************************************
c     * computes the radial distribution functions                     *
c     *                                                                *
c     * originally split out of the main routine 05-02-2005 M.G. Martin*
c     * last modified 08-25-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer nframe,nchain,nbox,nmolty,numax,nbin,nnbond,tor_max
c     --- functions
      character*3 numbers
c     --- local variables used in other subroutines
      integer invib,ijvib,intor,ijtor2,ijtor3,ijtor4
      dimension invib(nmolty,numax),ijvib(nmolty,numax,nnbond)
      dimension intor(nmolty,numax)
     &     ,ijtor2(nmolty,numax,tor_max)
     &     ,ijtor3(nmolty,numax,tor_max)
     &     ,ijtor4(nmolty,numax,tor_max)
      double precision rcut,addist
      dimension rcut(nbox)
      integer ncmt
      dimension ncmt(nbox,numax)
      double precision hmatrix,hinverse
      dimension hmatrix(nbox,3,3),hinverse(nbox,3,3)
      double precision volume
      dimension volume(nbox)
      double precision rxu,ryu,rzu,qqu
      dimension rxu(nchain,numax),ryu(nchain,numax),rzu(nchain,numax)
     &     ,qqu(nchain,numax)
      integer nunit,moltyp,nboxi
      dimension nunit(nmolty)
      dimension moltyp(nchain),nboxi(nchain)
c     --- truly local variables
      character*40 filename
      logical lchange
      integer ibox,imolty,jmolty,xx,yy,bin,ichain,iframe
      integer ii,jchain,jj,numxx,numyy
      double precision hist,bighist
      dimension hist(nmolty,numax,nmolty,numax,nbin)
      dimension bighist(nbox,nmolty,numax,nmolty,numax,nbin)
      double precision const,rlower,rupper,ideal,rcutsq,binstep
      double precision rdfcount
      dimension rdfcount(nbox,nmolty,nmolty)
      double precision rxui,ryui,rzui,rxuij,ryuij,rzuij,ruijsq,ruij
      double precision aframe,value

c     --- open the movie file and get the generic header info
      call getheader(nmolty,numax,nbox,nnbond,tor_max,nunit
     &     ,invib,ijvib,intor,ijtor2,ijtor3,ijtor4,rcut)

      write(6,*) 'Amount of additional displacement (0.0 for none).'
      read(5,*) addist

      write(6,*) 'The current cut-offs used to compute rdf in each'
      write(6,*) 'box are ',(rcut(ibox),ibox=1,nbox)
      write(6,*) 'Would you like to enter different cutoff values?'
     &     ,' (T/F)'
      read(5,*) lchange
      if ( lchange ) then
         do ibox=1,nbox
            write(6,*) 'Enter a new cut-off value for box ',ibox
            read(5,*) rcut(ibox)
         enddo
      endif

c     --- initialize some global radial distribution function vars
      do ibox = 1,nbox
         do imolty = 1,nmolty
            do jmolty = 1,nmolty
c              --- initialize rdfcount
               rdfcount(ibox,imolty,jmolty) = 0.0d0
               do xx = 1,numax
                  do yy = 1,numax
                     do bin = 1,nbin
c                       --- initialize bighist
                        bighist(ibox,imolty,xx,jmolty,yy,bin) = 0.0d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do iframe = 1, nframe
c        --- read in the frame step information
         call readframestep(nbox,nchain,nmolty,numax,ncmt
     &        ,hmatrix,hinverse
     &        ,volume,nunit,nboxi,moltyp,rxu,ryu,rzu,qqu)
c        --- compute the radial distribution histogram for this frame
         do ibox = 1,nbox
            rcutsq = rcut(ibox)*rcut(ibox)
            binstep = rcut(ibox)/dble(nbin)
c           --- initiallize hist
            do imolty = 1,nmolty
               do xx = 1,numax
                  do jmolty = 1,nmolty
                     do yy = 1,numax
                        do bin = 1,nbin
                           hist(imolty,xx,jmolty,yy,bin) = 0.0d0
                        enddo
                     enddo
                  enddo
               enddo
            enddo
c           --- loop through all of the molecules
            do ichain = 1, nchain
c              --- check if i is in relevant box
               if ( nboxi(ichain) .eq. ibox ) then
                  imolty = moltyp(ichain)
c                 --- loop over all beads ii of chain i 
                  do ii = 1, nunit(imolty)
                     rxui = rxu(ichain,ii)
                     ryui = ryu(ichain,ii)
                     rzui = rzu(ichain,ii)
c                    --- loop over all chains j with j>i as we only
c                    --- want intermolecular interactions
                     do jchain = ichain+1, nchain
c                       --- check for simulation box
                        if ( nboxi(jchain) .eq. ibox ) then
                           jmolty = moltyp(jchain)
c                          --- loop over all beads jj of chain jchain 
                           do jj = 1, nunit(jmolty)
                              rxuij = rxui - rxu(jchain,jj)
                              ryuij = ryui - ryu(jchain,jj)
                              rzuij = rzui - rzu(jchain,jj)
c                             --- minimum image the pair separations 
                              call mimage(nbox,hmatrix,hinverse
     &                             ,ibox,rxuij,ryuij,rzuij)
                              ruijsq = rxuij*rxuij + ryuij*ryuij 
     &                             + rzuij*rzuij
                              if ( ruijsq .lt. rcutsq ) then
                                 ruij = dsqrt(ruijsq)
                                 bin = dint(ruij/binstep) + 1 
                                 if ( imolty .eq. jmolty .and.
     &                                ii .eq. jj ) then
c                                   --- don't double count
                                    hist(imolty,ii,jmolty,jj,bin) = 
     &                                   hist(imolty,ii
     &                                   ,jmolty,jj,bin) + 1.0d0
                                 else
c                                   --- add a tally for forward
                                    hist(imolty,ii,jmolty,jj,bin) = 
     &                                   hist(imolty,ii
     &                                   ,jmolty,jj,bin) + 1.0d0
c                                   --- tally for reverse
                                    hist(jmolty,jj,imolty,ii,bin) = 
     &                                   hist(jmolty,jj
     &                                   ,imolty,ii,bin) + 1.0d0
                                 endif
                              endif
                           enddo
                        endif
                     enddo
                  enddo
               endif
            enddo

c           --- normalize the histogram and add it to the bighist
            do imolty = 1,nmolty
               numxx = ncmt(ibox,imolty)
               do jmolty = 1,nmolty
                  numyy = ncmt(ibox,jmolty)
                  if ( numxx * numyy .ne. 0 ) then
c                    --- increment rdfcount
                     rdfcount(ibox,imolty,jmolty) = 
     &                    rdfcount(ibox,imolty,jmolty) + 1.0d0
                     do xx = 1,nunit(imolty)
                        do yy = 1,nunit(jmolty)
c                          --- set const
                           const = (4.0d0*3.1415d0)
     &                          / (3.0d0 * volume(ibox) )
c                          --- determine whether any entries exist
                           if ( imolty .eq. jmolty 
     &                          .and. xx .eq. yy ) then
c                             --- identical types
c                             --- total interactions are (N*M-1)/2
                              const = const*dble(numxx*numyy-1)/2.0d0
                           else
c                             --- different types
c                             --- total interactions are N*M
                              const = const*numxx*numyy
                           endif
                           do bin = 1,nbin
c                             --- set rdf variables
                              rlower = dble(bin-1)*binstep
                              rupper = rlower + binstep
                              ideal = rupper**3 -rlower**3
c                             --- normalize hist
                              hist(imolty,xx,jmolty,yy,bin) =
     &                             hist(imolty,xx,jmolty,yy,bin)
     &                             /(const*ideal)
c                             --- add into bighist
                              bighist(ibox,imolty,xx,jmolty,yy,bin) =
     &                             bighist(ibox,imolty,xx,jmolty,yy,bin)
     &                             + hist(imolty,xx,jmolty,yy,bin)
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
      enddo

c     --- now that we are finished looping through all of the frames
c     --- we are ready for final processing and output
      write(6,*)
      write(6,*) 'type-type radial distribution functions in '
     &     ,'rdf_box_molecule(i)_unit(ii)_molecule(j)_unit(jj)'
      do ibox = 1,nbox
         binstep = rcut(ibox)/dble(nbin)
         do imolty = 1,nmolty
            do jmolty = 1,nmolty
               aframe = rdfcount(ibox,imolty,jmolty)
               do xx = 1,nunit(imolty)
                  do yy = 1,nunit(jmolty)
                     if ( aframe .gt. 0.5d0 ) then
c                       --- open up the output files
                        filename = 'rdf_b'
                        filename(6:8) = numbers(ibox)
c                       --- molecule type (i)
                        filename(9:10) = '_m'
                        filename(11:13) = numbers(imolty)
c                       --- atom 
                        filename(14:15) = '_u'
                        filename(16:18) = numbers(xx)
c                       --- molecule type (j)
                        filename(19:20) = '_m'
                        filename(21:23) = numbers(jmolty)
c                       --- atom (j)
                        filename(24:25) = '_u'
                        filename(26:28) = numbers(yy)
                        open(30,file=filename,form='formatted')
c                       --- output the rdf information
                        write(30,50) 0.0d0,addist,imolty,xx
     &                       ,jmolty,yy 
 50                     format(2f7.2,4i5)
                        do bin = 1,nbin
                           value = bighist(ibox,imolty,xx
     &                          ,jmolty,yy,bin) / aframe
                           rxuij =  binstep*(dble(bin)-0.5d0)
                           write(30,*) rxuij,value+addist
                        enddo
                        close(30)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

      subroutine getheader(nmolty,numax,nbox,nnbond,tor_max,nunit
     &     ,invib,ijvib,intor,ijtor2,ijtor3,ijtor4,rcut)
c     ******************************************************************
c     * opens towhee_movie and reads in the generic header information *
c     *                                                                *
c     * split out of mainprog 05-02-2005 by M.G. Martin                *
c     * last modified 05-03-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nmolty,numax,nbox,nnbond,tor_max
      integer invib,ijvib,intor,ijtor2,ijtor3,ijtor4
      dimension invib(nmolty,numax),ijvib(nmolty,numax,nnbond)
      dimension intor(nmolty,numax)
     &     ,ijtor2(nmolty,numax,tor_max)
     &     ,ijtor3(nmolty,numax,tor_max)
     &     ,ijtor4(nmolty,numax,tor_max)
      double precision rcut
      dimension rcut(nbox)
      integer nunit
      dimension nunit(nmolty)
c     --- local variables
      integer imolty,idum,jdum,ibox,ii,j,ivib

c     --- open the movie file
      open(10,file='towhee_movie',form='formatted')

c     --- skip the movie version
      read(10,*)
c     --- skip this line except for rcut
      read(10,*) (idum,jdum=1,4),(rcut(ibox),ibox=1,nbox)

c     --- read in all of the connectivity information
      do imolty = 1,nmolty
         read(10,*) nunit(imolty)
         if ( nunit(imolty) .gt. numax ) then
            write(6,*) 'moltype ',imolty,' units ',nunit(imolty)
     &           ,' is > numax ',numax
            stop
         endif
c        --- read bond connectivity information
         do ii = 1,nunit(imolty)
            read(10,*) invib(imolty,ii)
     &           ,(ijvib(imolty,ii,ivib),ivib=1,invib(imolty,ii))
         enddo

c        --- read torsional connectivity information
         do j = 1,nunit(imolty)
            read(10,*) intor(imolty,j)
            do ii = 1,intor(imolty,j)
               read(10,*) ijtor2(imolty,j,ii)
     &              ,ijtor3(imolty,j,ii),ijtor4(imolty,j,ii)
            enddo
         enddo
      enddo

      return
      end

      function numbers(index)
c     ******************************************************************
c     * turns an integer into a 3 character string                     *
c     *                                                                *
c     * originally written 05-03-2005 by M.G. Martin                   *
c     * last modified 05-09-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the function
      integer index
      character*3 numbers
c     --- local variables
      integer irun,ihundred,iten
      character*1 table
      dimension table(0:9)
      data (table(irun),irun=0,9) /
     &     '0','1','2','3','4','5','6','7','8','9'/
      if ( index .gt. 999 ) then
         write(6,*) 'NUMBERS: index out of bounds',index
         stop
      endif
c     --- assign hundreds digit
      irun = index
      ihundred = irun/100
      numbers(1:1) = table(ihundred)
c     --- assign tens digit
      irun = irun -(ihundred*100)
      iten = irun/10
      numbers(2:2) = table(iten)
c     --- assign ones digit
      irun = irun -(iten*10)
      numbers(3:3) = table(irun)

      return
      end

      subroutine readframestep(nbox,nchain,nmolty,numax,ncmt
     &     ,hmatrix,hinverse
     &     ,volume,nunit,nboxi,moltyp,rxu,ryu,rzu,qqu)
c     ******************************************************************
c     * reads in a single frame step of information for the analyse    *
c     * routines                                                       *
c     *                                                                *
c     * originially written 05-03-2005 by M.G. Martin                  *
c     * last modified 05-04-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer nbox,nchain,nmolty,numax
      integer ncmt
      dimension ncmt(nbox,nmolty)
      integer nunit
      dimension nunit(nmolty)
      integer nboxi,moltyp
      dimension nboxi(nchain),moltyp(nchain)
      double precision hmatrix,hinverse,volume
      dimension hmatrix(nbox,3,3),hinverse(nbox,3,3)
      dimension volume(nbox)
      double precision rxu,ryu,rzu,qqu
      dimension rxu(nchain,numax),ryu(nchain,numax),rzu(nchain,numax)
      dimension qqu(nchain,numax)
c     --- local variables
      integer ibox,idim,jdim,ichain,imolty,chnum,iunit

c     --- skip frame step number
      read(10,*)
c     --- read in ncmt for the boxes
      do ibox = 1,nbox
         read(10,*) (ncmt(ibox,imolty),imolty=1,nmolty)
c        --- read in the hmatrix and the boxvolume
         do idim = 1,3
            read(10,*) (hmatrix(ibox,idim,jdim),jdim=1,3)
         enddo
c        --- invert the matrix and get the box volume
         call inverthmatrix(nbox,ibox,hmatrix,hinverse,volume)
      enddo

      do ichain = 1,nchain
         read(10,*) chnum,imolty,nunit(imolty),nboxi(chnum)
         moltyp(chnum) = imolty
         if ( ichain .ne. chnum ) then
            stop 'ichain ne chnum'
         endif
         do iunit = 1, nunit(imolty)
            read(10,*) rxu(chnum,iunit),ryu(chnum,iunit)
     &           ,rzu(chnum,iunit),qqu(chnum,iunit)
         enddo
      enddo

      return
      end

      subroutine reset_movie()
c     ******************************************************************
c     * closes and rewinds the movie file                              *
c     *                                                                *
c     * originally written 05-03-2005 by M.G. Martin                   *
c     * last modified 05-03-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- close towhee_movie
      close(10)
c     --- rewind towhee_movie
      rewind(10)
      return
      end

      subroutine analyse_profile(nframe,nchain,nbox,nmolty,numax,nbin
     &     ,nnbond,tor_max)
c     ******************************************************************
c     * computes the atom distribution functions                       *
c     *                                                                *
c     * originally split out of the main routine 05-03-2005 M.G. Martin*
c     * last modified 07-26-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer nframe,nchain,nbox,nmolty,numax,nbin,nnbond,tor_max
c     --- functions
      character*3 numbers
c     --- local variables used in other subroutines
      integer invib,ijvib,intor,ijtor2,ijtor3,ijtor4
      dimension invib(nmolty,numax),ijvib(nmolty,numax,nnbond)
      dimension intor(nmolty,numax)
     &     ,ijtor2(nmolty,numax,tor_max)
     &     ,ijtor3(nmolty,numax,tor_max)
     &     ,ijtor4(nmolty,numax,tor_max)
      double precision rcut
      dimension rcut(nbox)
      integer nunit
      dimension nunit(nmolty)
      integer ncmt
      dimension ncmt(nbox,nmolty)
      double precision hmatrix,hinverse,volume
      dimension hmatrix(nbox,3,3),hinverse(nbox,3,3)
      dimension volume(nbox)
      integer nboxi,moltyp
      dimension nboxi(nchain),moltyp(nchain)
      double precision rxu,ryu,rzu,qqu
      dimension rxu(nchain,numax),ryu(nchain,numax),rzu(nchain,numax)
      dimension qqu(nchain,numax)
c     --- truly local variables
      character*40 filename
      integer ibox,imolty,iunit,ibin,jbin,kbin,iframe,ichain
      double precision xyprof,xzprof,yzprof,rx,ry,rz,dmidi,dmidj
      dimension xyprof(nbox,nmolty,numax,nbin,nbin)
      dimension xzprof(nbox,nmolty,numax,nbin,nbin)
      dimension yzprof(nbox,nmolty,numax,nbin,nbin)
      double precision pcount
      dimension pcount(nbox,nmolty)
      double precision arba,arbb,arbc,avgden

c     --- open the movie file and get the generic header info
      call getheader(nmolty,numax,nbox,nnbond,tor_max,nunit
     &     ,invib,ijvib,intor,ijtor2,ijtor3,ijtor4,rcut)

      do ibox = 1,nbox
         do imolty = 1,nmolty
            pcount(ibox,imolty) = 0.0d0
            do iunit = 1,numax
               do ibin = 1,nbin
                 do jbin = 1,nbin
                   xyprof(ibox,imolty,iunit,ibin,jbin) = 0.0d0
                   xzprof(ibox,imolty,iunit,ibin,jbin) = 0.0d0
                   yzprof(ibox,imolty,iunit,ibin,jbin) = 0.0d0
                 enddo
               enddo
            enddo
         enddo
      enddo

      do iframe = 1,nframe
         call readframestep(nbox,nchain,nmolty,numax,ncmt
     &        ,hmatrix,hinverse
     &        ,volume,nunit,nboxi,moltyp,rxu,ryu,rzu,qqu)

c        --- compute the atom distribution profiles for this frame
         do ichain = 1,nchain
            ibox = nboxi(ichain)
            if ( ibox .ne. 0 ) then
               imolty = moltyp(ichain)
c              --- increment pcount
               pcount(ibox,imolty) = pcount(ibox,imolty) + 1.0d0
               do iunit = 1,nunit(imolty)
                  rx = rxu(ichain,iunit)
                  ry = ryu(ichain,iunit)
                  rz = rzu(ichain,iunit)
c                 --- transform into arbitrary coordinates
                  call carttoarb(nbox,ibox,hinverse,rx,ry,rz,arba
     &                 ,arbb,arbc)
c                 --- put this in the box
                  if ( arba .lt. 0.0d0 ) arba = arba + 1.0d0
                  if ( arba .gt. 1.0d0 ) arba = arba - 1.0d0
                  if ( arbb .lt. 0.0d0 ) arbb = arbb + 1.0d0
                  if ( arbb .gt. 1.0d0 ) arbb = arbb - 1.0d0
                  if ( arbc .lt. 0.0d0 ) arbc = arbc + 1.0d0
                  if ( arbc .gt. 1.0d0 ) arbc = arbc - 1.0d0
                  ibin = dint(dble(nbin)*arba) + 1
                  jbin = dint(dble(nbin)*arbb) + 1
                  kbin = dint(dble(nbin)*arbc) + 1
                  xyprof(ibox,imolty,iunit,ibin,jbin) = 
     &                 xyprof(ibox,imolty,iunit,ibin,jbin) + 1.0d0
                  xzprof(ibox,imolty,iunit,ibin,kbin) = 
     &                 xzprof(ibox,imolty,iunit,ibin,kbin) + 1.0d0
                  yzprof(ibox,imolty,iunit,jbin,kbin) = 
     &                 yzprof(ibox,imolty,iunit,jbin,kbin) + 1.0d0
               enddo
            endif
         enddo
      enddo

c     --- output profile information
      write(6,*) 'xy profile in xypr_b#_t#_u#'
      write(6,*) 'xz profile in xzpr_b#_t#_u#'
      write(6,*) 'yz profile in yzpr_b#_t#_u#'
 
      do ibox = 1,nbox
         do imolty = 1,nmolty
            if ( pcount(ibox,imolty) .gt. 0.5d0 ) then
c              --- initialize filename
               filename = '--pr_b'
c              --- set the identical portions of the filename
               filename(7:9) = numbers(ibox)
               filename(10:11) = '_t'
               filename(12:14) = numbers(imolty)
               filename(15:16) = '_u'
c              --- loop through each unit of the molecule
               do iunit = 1,nunit(imolty)
                  filename(17:19) = numbers(iunit)
c                 --- open the three files
                  filename(1:2) = 'xy'
                  open(45,file=filename,form='formatted')
                  filename(1:2) = 'xz'
                  open(46,file=filename,form='formatted')
                  filename(1:2) = 'yz'
                  open(47,file=filename,form='formatted')
c                 --- output the molecule type and unit
                  write(45,*) imolty,iunit
                  write(46,*) imolty,iunit 
                  write(47,*) imolty,iunit 
c                 --- output the profiles
                  do ibin = 1,nbin
                    dmidi = (dble(ibin) - 0.5d0)/dble(nbin)
                    do jbin = 1,nbin
                       dmidj = (dble(jbin) - 0.5d0)/dble(nbin)
                       write(45,*) dmidi,dmidj,dble(nbin)
     &                    *xyprof(ibox,imolty,iunit,ibin,jbin)
     &                    /pcount(ibox,imolty)
                       write(46,*) dmidi,dmidj,dble(nbin)
     &                    *xzprof(ibox,imolty,iunit,ibin,jbin) 
     &                    /pcount(ibox,imolty)
                       write(47,*) dmidi,dmidj,dble(nbin)
     &                    *yzprof(ibox,imolty,iunit,ibin,jbin) 
     &                    /pcount(ibox,imolty)
                    enddo
                  enddo
c                 --- close the files
                  close(45)
                  close(46)
                  close(47)
               enddo
            endif
         enddo
      enddo

      return
      end

      subroutine analyse_vib(nframe,nchain,nbox,nmolty,numax,nbin
     &     ,nnbond,tor_max)
c     ******************************************************************
c     * computes the vibration distribution functions                  *
c     *                                                                *
c     * originally split out of the main routine 05-03-2005 M.G. Martin*
c     * last modified 07-27-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer nframe,nchain,nbox,nmolty,numax,nbin,nnbond,tor_max
c     --- functions
      character*3 numbers
c     --- local variables used in other subroutines
      integer invib,ijvib,intor,ijtor2,ijtor3,ijtor4
      dimension invib(nmolty,numax),ijvib(nmolty,numax,nnbond)
      dimension intor(nmolty,numax)
     &     ,ijtor2(nmolty,numax,tor_max)
     &     ,ijtor3(nmolty,numax,tor_max)
     &     ,ijtor4(nmolty,numax,tor_max)
      double precision rcut
      dimension rcut(nbox)
      integer nunit
      dimension nunit(nmolty)
      integer ncmt
      dimension ncmt(nbox,nmolty)
      double precision hmatrix,hinverse,volume
      dimension hmatrix(nbox,3,3),hinverse(nbox,3,3)
      dimension volume(nbox)
      integer nboxi,moltyp
      dimension nboxi(nchain),moltyp(nchain)
      double precision rxu,ryu,rzu,qqu
      dimension rxu(nchain,numax),ryu(nchain,numax),rzu(nchain,numax)
      dimension qqu(nchain,numax)
c     --- truly local variables
      character*40 filename
      integer bond,vib_1,vib_2,vib_num
      double precision vib_bin_size,maxvib,length
      parameter( maxvib = 2.0d0)
      dimension vib_1(nmolty,numax*nnbond),vib_2(nmolty,numax*nnbond)
      dimension vib_num(nmolty)
      double precision vib_bin,vib_tot,vib_avg
      dimension vib_bin(nbox,nmolty,numax*nnbond,nbin)
      dimension vib_avg(nbox,nmolty,numax*nnbond)
      dimension vib_tot(nbox,nmolty,numax*nnbond)
      integer ibin,imolty,ii,iv,iuvib,ibox,ichain,j,ip1,iframe
      double precision xtmp,ytmp,ztmp
      dimension xtmp(4),ytmp(4),ztmp(4)
      double precision total,value,sdev

c     --- open the movie file and get the generic header info
      call getheader(nmolty,numax,nbox,nnbond,tor_max,nunit
     &     ,invib,ijvib,intor,ijtor2,ijtor3,ijtor4,rcut)

c     --- compute vib_bin_size
      vib_bin_size = maxvib/dble(nbin)
c     --- figure out how many bonds exist and assign them numbers
      do imolty = 1,nmolty
         bond = 0
         do ii = 1,nunit(imolty)
            do iv = 1,invib(imolty,ii)
               iuvib = ijvib(imolty,ii,iv)
c              --- determine whether the beads connected to this unit
c              --- are of higher index than ii
               if ( iuvib .gt. ii ) then
                  bond = bond + 1
                  vib_1(imolty,bond) = ii
                  vib_2(imolty,bond) = iuvib
               endif
            enddo
         enddo
         vib_num(imolty) = bond
      enddo

      do ibox = 1,nbox
         do imolty = 1,nmolty
            do bond = 1,vib_num(imolty)
               do ibin = 1,nbin
                  vib_bin(ibox,imolty,bond,ibin) = 0.0d0
                  vib_tot(ibox,imolty,bond) = 0.0d0
                  vib_avg(ibox,imolty,bond) = 0.0d0
               enddo
            enddo
         enddo
      enddo

      do iframe = 1,nframe
         call readframestep(nbox,nchain,nmolty,numax,ncmt
     &        ,hmatrix,hinverse
     &        ,volume,nunit,nboxi,moltyp,rxu,ryu,rzu,qqu)
         do ibox=1,nbox
            do ichain = 1, nchain
               imolty = moltyp(ichain)
c              --- check if i is in relevant box 
               if ( nboxi(ichain) .eq. ibox ) then
c                 --- vibration
                  do bond = 1,vib_num(imolty)
                     j = vib_1(imolty,bond)
                     ip1 = vib_2(imolty,bond)
                     xtmp(1) = rxu(ichain,j) - rxu(ichain,ip1)
                     ytmp(1) = ryu(ichain,j) - ryu(ichain,ip1)
                     ztmp(1) = rzu(ichain,j) - rzu(ichain,ip1)
                     call mimage(nbox,hmatrix,hinverse,ibox
     &                    ,xtmp(1),ytmp(1),ztmp(1))
                     length = dsqrt( xtmp(1)**2 + ytmp(1)**2
     &                    + ztmp(1)**2 )
                     vib_avg(ibox,imolty,bond) = 
     &                    vib_avg(ibox,imolty,bond) + length
                     ibin = dint( length /vib_bin_size)+1
                     vib_bin(ibox,imolty,bond,ibin) = 
     &                    vib_bin(ibox,imolty,bond,ibin) + 1.0d0
                     vib_tot(ibox,imolty,bond) = 
     &                    vib_tot(ibox,imolty,bond)+ 1.0d0
                  enddo
               endif
            enddo
         enddo
      enddo

c     --- output the vibration distributions for each bond in the 
c     --- molecule
      write(6,*) 'vibration distributions output to vib_b#_t#_u#_u#'
      do ibox = 1,nbox
         do imolty = 1,nmolty
            do bond = 1,vib_num(imolty)
               filename = 'vib_b'
               filename(6:8) = numbers(ibox)
               filename(9:10) = '_t'
               filename(11:13) = numbers(imolty)
               filename(14:15) = '_u'
               filename(16:18) = numbers(vib_1(imolty,bond))
               filename(19:20) = '_u'
               filename(21:23) = numbers(vib_2(imolty,bond))
               open(67,file=filename,form='formatted')
               total = vib_tot(ibox,imolty,bond)*vib_bin_size
               if ( total .gt. 0.0d0 ) then
c                 --- compute the average bond length
                  vib_avg(ibox,imolty,bond) 
     &                 = vib_avg(ibox,imolty,bond)
     &                 /vib_tot(ibox,imolty,bond)
c                 --- output the distribution 
                  value = -vib_bin_size/2.0d0
                  sdev = 0.0d0
                  do ibin = 1,nbin
                     value = value + vib_bin_size
                     write(67,*) value
     &                    ,vib_bin(ibox,imolty,bond,ibin)/total
c                    --- compute standard deviation
                     sdev = sdev + vib_bin(ibox,imolty,bond,ibin)
     &                    *(value-vib_avg(ibox,imolty,bond))**2
                  enddo
                  write(67,*)
c                 --- output the average and sdev to the screen
                  write(6,*) 
     &                 'Moltyp ',imolty,' bond '
     &                 ,vib_1(imolty,bond),vib_2(imolty,bond)
     &                 ,' average ',vib_avg(ibox,imolty,bond)
     &                 ,' sdev '
     &                 ,dsqrt(sdev/vib_tot(ibox,imolty,bond))
               endif
               close(67)
            enddo
         enddo
      enddo                  

      return
      end

      subroutine analyse_bend(nframe,nchain,nbox,nmolty,numax,nbin
     &     ,nnbond,tor_max)
c     ******************************************************************
c     * computes the bending distribution functions                    *
c     *                                                                *
c     * originally split out of the main routine 05-03-2005 M.G. Martin*
c     * last modified 07-27-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer nframe,nchain,nbox,nmolty,numax,nbin,nnbond,tor_max
c     --- functions
      character*3 numbers
      double precision arccos
c     --- local variables used in other subroutines
      integer invib,ijvib,intor,ijtor2,ijtor3,ijtor4
      dimension invib(nmolty,numax),ijvib(nmolty,numax,nnbond)
      dimension intor(nmolty,numax)
     &     ,ijtor2(nmolty,numax,tor_max)
     &     ,ijtor3(nmolty,numax,tor_max)
     &     ,ijtor4(nmolty,numax,tor_max)
      double precision rcut
      dimension rcut(nbox)
      integer nunit
      dimension nunit(nmolty)
      integer ncmt
      dimension ncmt(nbox,nmolty)
      double precision hmatrix,hinverse,volume
      dimension hmatrix(nbox,3,3),hinverse(nbox,3,3)
      dimension volume(nbox)
      integer nboxi,moltyp
      dimension nboxi(nchain),moltyp(nchain)
      double precision rxu,ryu,rzu,qqu
      dimension rxu(nchain,numax),ryu(nchain,numax),rzu(nchain,numax)
      dimension qqu(nchain,numax)
c     --- truly local variables
      character*40 filename
      integer bend,iv,iuvib,iuv,iutest
      double precision angle_bin,angle_tot,onepi
     &     ,ang_bin_size,value,total,degree,angle_avg
      integer angle_1,angle_2,angle_3,angle_num
      dimension angle_num(nmolty)
      dimension angle_1(nmolty,numax*nnbond*nnbond)
      dimension angle_2(nmolty,numax*nnbond*nnbond)
      dimension angle_3(nmolty,numax*nnbond*nnbond)
      dimension angle_tot(nbox,nmolty,numax*nnbond*nnbond)
      dimension angle_avg(nbox,nmolty,numax*nnbond*nnbond)
      dimension angle_bin(nbox,nmolty,numax*nnbond*nnbond,nbin)
      double precision distone,disttwo
      integer imolty,ii,ibox,bin,iframe,ichain,j,ip1,ip2
      double precision xtmp,ytmp,ztmp
      dimension xtmp(4),ytmp(4),ztmp(4)
      double precision thetac,theta,sdev

c     --- open the movie file and get the generic header info
      call getheader(nmolty,numax,nbox,nnbond,tor_max,nunit
     &     ,invib,ijvib,intor,ijtor2,ijtor3,ijtor4,rcut)
c     --- set onepi
      onepi = arccos(-1.0d0)

c     --- compute ang_bin_size
      ang_bin_size = onepi/dble(nbin)
c     --- figure out how many angles exist and assign them numbers
      do imolty = 1,nmolty
         bend = 0
         do ii = 1,nunit(imolty)
            do iv = 1,invib(imolty,ii)
               iuvib = ijvib(imolty,ii,iv)
c              --- determine whether the beads connected to this unit
c              --- are of higher index than ii
               do iuv = 1,invib(imolty,iuvib)
                  iutest = ijvib(imolty,iuvib,iuv)
                  if ( iutest .gt. ii ) then
                     bend = bend + 1
                     angle_1(imolty,bend) = ii
                     angle_2(imolty,bend) = iuvib
                     angle_3(imolty,bend) = iutest
                  endif
               enddo
            enddo
         enddo
         angle_num(imolty) = bend
      enddo

      do ibox = 1,nbox
         do imolty = 1,nmolty
            do bend = 1,angle_num(imolty)
               do bin = 1,nbin
                  angle_bin(ibox,imolty,bend,bin) = 0.0d0
                  angle_tot(ibox,imolty,bend) = 0.0d0
                  angle_avg(ibox,imolty,bend) = 0.0d0
               enddo
            enddo
         enddo
      enddo

      do iframe = 1,nframe
         call readframestep(nbox,nchain,nmolty,numax,ncmt
     &        ,hmatrix,hinverse
     &        ,volume,nunit,nboxi,moltyp,rxu,ryu,rzu,qqu)
         do ibox=1,nbox
            do ichain = 1, nchain
               imolty = moltyp(ichain)
c              --- check if i is in relevant box 
               if ( nboxi(ichain) .eq. ibox ) then
                  do bend = 1,angle_num(imolty)
                     j = angle_1(imolty,bend)
                     ip1 = angle_2(imolty,bend)
                     ip2 = angle_3(imolty,bend)
                     xtmp(1) = rxu(ichain,j) - rxu(ichain,ip1)
                     ytmp(1) = ryu(ichain,j) - ryu(ichain,ip1)
                     ztmp(1) = rzu(ichain,j) - rzu(ichain,ip1)
                     call mimage(nbox,hmatrix,hinverse,ibox
     &                    ,xtmp(1),ytmp(1),ztmp(1))
                     distone = dsqrt( xtmp(1)**2 + ytmp(1)**2
     &                    + ztmp(1)**2 )
                     xtmp(2) = rxu(ichain,ip2) - rxu(ichain,ip1)
                     ytmp(2) = ryu(ichain,ip2) - ryu(ichain,ip1)
                     ztmp(2) = rzu(ichain,ip2) - rzu(ichain,ip1)
                     call mimage(nbox,hmatrix,hinverse,ibox
     &                    ,xtmp(2),ytmp(2),ztmp(2))
                     disttwo = dsqrt( xtmp(2)**2 + ytmp(2)**2
     &                    + ztmp(2)**2 )
                     thetac = ( xtmp(1)*xtmp(2) + ytmp(1)*ytmp(2)
     +                    + ztmp(1)*ztmp(2) )/(distone*disttwo)
                     theta = arccos(thetac)
                     angle_avg(ibox,imolty,bend) = 
     &                    angle_avg(ibox,imolty,bend) + theta
                     bin = dint(theta/ang_bin_size)+1
                     angle_bin(ibox,imolty,bend,bin) = 
     &                    angle_bin(ibox,imolty,bend,bin) + 1.0d0
                     angle_tot(ibox,imolty,bend) = 
     &                    angle_tot(ibox,imolty,bend)+ 1.0d0
                  enddo
               endif
            enddo
         enddo
      enddo

c     --- output the bending angle distributions for each angle 
c     --- in the molecule
      write(6,*) 'Bond bending distributions in ben_b#_t#_u#_u#_u#'
      degree = 180.0d0/onepi
      do ibox = 1,nbox
         do imolty = 1,nmolty
            do bend = 1,angle_num(imolty)
               total = angle_tot(ibox,imolty,bend)*ang_bin_size
     &              *degree
               if ( total .gt. 0.0d0 ) then
                  filename = 'ben_b'
                  filename(6:8) = numbers(ibox)
                  filename(9:10) = '_t'
                  filename(11:13) = numbers(imolty)
                  filename(14:15) = '_u'
                  filename(16:18) = numbers(angle_1(imolty,bend))
                  filename(19:20) = '_u'
                  filename(21:23) = numbers(angle_2(imolty,bend))
                  filename(24:25) = '_u'
                  filename(26:28) = numbers(angle_3(imolty,bend))
                  open(50,file=filename,form='formatted')
c                 --- compute the average angle
                  angle_avg(ibox,imolty,bend) = 
     &                 degree*angle_avg(ibox,imolty,bend)
     &                 /angle_tot(ibox,imolty,bend)
c                 --- output the entire distribution
                  value = -degree*ang_bin_size/2.0d0
                  sdev = 0.0d0
                  do bin = 1,nbin
                     value = value + degree*ang_bin_size
                     write(50,*) value
     &                    ,angle_bin(ibox,imolty,bend,bin)/total
c                    --- compute standard deviation
                     sdev = sdev + angle_bin(ibox,imolty,bend,bin)
     &                    *(value - angle_avg(ibox,imolty,bend))**2
                  enddo
c                 --- output average and standard deviation
                  write(6,*) 
     &                 'Moltyp ',imolty,' angle '
     &                 ,angle_1(imolty,bend),angle_2(imolty,bend)
     &                 ,angle_3(imolty,bend),' average '
     &                 ,angle_avg(ibox,imolty,bend),' sdev '
     &                 ,dsqrt(sdev/angle_tot(ibox,imolty,bend))
                  close(50)
               endif
            enddo
         enddo
      enddo                  

      return
      end

      subroutine analyse_dihed(nframe,nchain,nbox,nmolty,numax,nbin
     &     ,nnbond,tor_max,tor_smax)
c     ******************************************************************
c     * computes the bending distribution functions                    *
c     *                                                                *
c     * originally split out of the main routine 05-03-2005 M.G. Martin*
c     * last modified 07-01-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer nframe,nchain,nbox,nmolty,numax,nbin,nnbond,tor_max
      integer tor_smax
c     --- functions
      character*3 numbers
      double precision arccos
c     --- local variables used in other subroutines
      integer invib,ijvib,intor,ijtor2,ijtor3,ijtor4
      dimension invib(nmolty,numax),ijvib(nmolty,numax,nnbond)
      dimension intor(nmolty,numax)
     &     ,ijtor2(nmolty,numax,tor_max)
     &     ,ijtor3(nmolty,numax,tor_max)
     &     ,ijtor4(nmolty,numax,tor_max)
      double precision rcut
      dimension rcut(nbox)
      integer nunit
      dimension nunit(nmolty)
      integer ncmt
      dimension ncmt(nbox,nmolty)
      double precision hmatrix,hinverse,volume
      dimension hmatrix(nbox,3,3),hinverse(nbox,3,3)
      dimension volume(nbox)
      integer nboxi,moltyp
      dimension nboxi(nchain),moltyp(nchain)
      double precision rxu,ryu,rzu,qqu
      dimension rxu(nchain,numax),ryu(nchain,numax),rzu(nchain,numax)
      dimension qqu(nchain,numax)
c     --- truly local variables
      character*40 filename
      integer torsion,tor_1,tor_2,tor_3,tor_4,tor_num,itor
      dimension tor_1(nmolty,tor_smax),tor_2(nmolty,tor_smax)
     &     ,tor_3(nmolty,tor_smax),tor_4(nmolty,tor_smax)
      dimension tor_num(nmolty)
      integer iutor
      double precision tor_bin,tor_tot,xcc,ycc,zcc,tcc,tor_bin_size
      double precision xaa1,yaa1,zaa1,xa1a2,ya1a2,za1a2,daa1,da1a2
      double precision dot,thetac,theta
      dimension tor_bin(nbox,nmolty,tor_smax,nbin)
      dimension tor_tot(nbox,nmolty,tor_smax)
      integer imolty,ii,ibox,bin,iframe,ichain,gaudef,dum
      integer atom1,atom2,atom3,atom4
      double precision onepi,xtmp,ytmp,ztmp,total,value
      dimension xtmp(4),ytmp(4),ztmp(4)

c     --- open the movie file and get the generic header info
      call getheader(nmolty,numax,nbox,nnbond,tor_max,nunit
     &     ,invib,ijvib,intor,ijtor2,ijtor3,ijtor4,rcut)
c     --- set onepi
      onepi = arccos(-1.0d0)

c     --- compute tor_bin_size
      tor_bin_size = 360.0d0/dble(nbin)
c     --- figure out how many torsions exist and assign them numbers
      do imolty = 1,nmolty
         torsion = 0
         do ii = 1,nunit(imolty)
            do itor = 1,intor(imolty,ii)
               iutor = ijtor4(imolty,ii,itor)
c              --- determine whether final bead connected to this unit
c              --- is of higher index than ii
               if ( iutor .gt. ii ) then
                  torsion = torsion + 1
                  tor_1(imolty,torsion) = ii
                  tor_2(imolty,torsion) = ijtor2(imolty,ii,itor)
                  tor_3(imolty,torsion) = ijtor3(imolty,ii,itor)
                  tor_4(imolty,torsion) = iutor
               endif
            enddo
         enddo
         tor_num(imolty) = torsion
      enddo

c     --- initialize the arrays
      do ibox = 1,nbox
         do imolty = 1,nmolty
            do torsion = 1,tor_num(imolty)
               do bin = 1,nbin
                  tor_bin(ibox,imolty,torsion,bin) = 0.0d0
                  tor_tot(ibox,imolty,torsion) = 0.0d0
               enddo
            enddo
         enddo
      enddo

c     --- analyse the bonded terms distribution
      do iframe = 1,nframe
         call readframestep(nbox,nchain,nmolty,numax,ncmt
     &        ,hmatrix,hinverse
     &        ,volume,nunit,nboxi,moltyp,rxu,ryu,rzu,qqu)
         do ibox=1,nbox
            do ichain = 1, nchain
               imolty = moltyp(ichain)
c              --- check if i is in relevant box 
               if ( nboxi(ichain) .eq. ibox ) then
c                 --- molecule with dihedral potenials
                  gaudef = 1
                  dum = 0
                  do torsion = 1, tor_num(imolty)
                     atom1 = tor_1(imolty,torsion)
                     atom2 = tor_2(imolty,torsion)
                     atom3 = tor_3(imolty,torsion)
                     atom4 = tor_4(imolty,torsion)
                     xtmp(1) = rxu(ichain,atom2) - rxu(ichain,atom1)
                     ytmp(1) = ryu(ichain,atom2) - ryu(ichain,atom1)
                     ztmp(1) = rzu(ichain,atom2) - rzu(ichain,atom1)
                     call mimage(nbox,hmatrix,hinverse,ibox
     &                    ,xtmp(1),ytmp(1),ztmp(1))
                     xtmp(2) = rxu(ichain,atom3) - rxu(ichain,atom2)
                     ytmp(2) = ryu(ichain,atom3) - ryu(ichain,atom2)
                     ztmp(2) = rzu(ichain,atom3) - rzu(ichain,atom2)
                     call mimage(nbox,hmatrix,hinverse,ibox
     &                    ,xtmp(2),ytmp(2),ztmp(2))
                     xtmp(3) = rxu(ichain,atom4) - rxu(ichain,atom3)
                     ytmp(3) = ryu(ichain,atom4) - ryu(ichain,atom3)
                     ztmp(3) = rzu(ichain,atom4) - rzu(ichain,atom3)
                     call mimage(nbox,hmatrix,hinverse,ibox
     &                    ,xtmp(3),ytmp(3),ztmp(3))
c                    --- calculate cross product d_a x d_a-1
                     xaa1 = ytmp(1)*ztmp(2) - ztmp(1)*ytmp(2) 
                     yaa1 = ztmp(1)*xtmp(2) - xtmp(1)*ztmp(2)
                     zaa1 = xtmp(1)*ytmp(2) - ytmp(1)*xtmp(2)
c                    --- calculate cross product d_a-1 x d_a-2
                     xa1a2 = ytmp(2)*ztmp(3) - ztmp(2)*ytmp(3)
                     ya1a2 = ztmp(2)*xtmp(3) - xtmp(2)*ztmp(3)
                     za1a2 = xtmp(2)*ytmp(3) - ytmp(2)*xtmp(3)
c                    --- calculate lengths of cross products
                     daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                     da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
c                    --- calculate dot product of cross products
                     dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                     thetac = - dot / ( daa1 * da1a2 )
                     theta = arccos(thetac)
c                    ---calculate cross product of cross products 
                     xcc = yaa1*za1a2 - zaa1*ya1a2
                     ycc = zaa1*xa1a2 - xaa1*za1a2
                     zcc = xaa1*ya1a2 - yaa1*xa1a2
c                    --- calculate scalar triple product 
                     tcc = xcc*xtmp(2) + ycc*ytmp(2) + zcc*ztmp(2)
                     if ( tcc .lt. 0.0d0 ) theta = - theta
c                    --- convert theta to degrees 
                     theta = theta*(180.0d0/onepi)
c                    --- shift by 180 to align with the cis=0 convention
                     theta = theta + 180.0d0
c                    --- bin the torsion
                     bin = dint((theta)/tor_bin_size)+1
c                    --- possible that we have an angle right on the 
c                    --- boundary
                     if ( bin .gt. nbin ) then
c                       --- this is OK if right on the boundary
                        if ( theta .eq. 360.0d0 ) then
c                          --- just put this in the maximum bin
                           bin = nbin
                        else
c                          --- this should not happen
                           write(6,*) 'Bin out of range in torsions'
                           write(6,*) 'bin,nbin',bin,nbin
                           write(6,*) 'theta',theta
                           stop
                        endif
                     elseif ( bin .lt. 1 ) then
c                       --- should not happen
                        write(6,*) 'theta,bin,nbin',theta,bin,nbin
                        stop
                     endif
                     tor_bin(ibox,imolty,torsion,bin) = 
     &                    tor_bin(ibox,imolty,torsion,bin)+1.0d0
                     tor_tot(ibox,imolty,torsion) = 
     &                    tor_tot(ibox,imolty,torsion) + 1.0d0
                  enddo
               endif
            enddo
         enddo
      enddo

      write(6,*)
      write(6,*) 'torsion angle distribution in tor_b#_t#_u#_u#_u#_u#'
      do ibox  = 1,nbox
         do imolty = 1,nmolty
            do torsion = 1,tor_num(imolty)
c              --- compute total for the pdf
               total = tor_tot(ibox,imolty,torsion)*tor_bin_size
               if ( total .gt. 0.1d-10 ) then
c                 --- open the torsion distribution file
                  filename = 'tor_b'
                  filename(6:8) = numbers(ibox)
                  filename(9:10) = '_t'
                  filename(11:13) = numbers(imolty)
                  filename(14:15) = '_u'
                  filename(16:18) = numbers(tor_1(imolty,torsion))
                  filename(19:20) = '_u'
                  filename(21:23) = numbers(tor_2(imolty,torsion))
                  filename(24:25) = '_u'
                  filename(26:28) = numbers(tor_3(imolty,torsion))
                  filename(29:30) = '_u'
                  filename(31:33) = numbers(tor_4(imolty,torsion))
                  open(95,file=filename,form='formatted')
C                 --- output the torsion probablility distributions
                  value = 0.0d0 - (0.5d0*tor_bin_size)
                  do bin = 1,nbin
                     value = value + tor_bin_size
                     write(95,*) value
     &                    ,tor_bin(ibox,imolty,torsion,bin)/total
                  enddo
                  close(95)
               endif
            enddo
         enddo
      enddo

      return
      end

      function arccos( value )
c     ******************************************************************
c     * computes the arc cosine of a value, with safety checks to make *
c     * sure that value is between -1.0 and 1.0.  If value is outside  *
c     * of that range, it is set to the nearest extreme of the range   *
c     *                                                                *
c     * originally written 02-12-2001 by M.G. Martin                   *
c     * last modified 05-09-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the function
      double precision value,arccos

      if ( value .gt. 1.0d0 ) then
         arccos = 0.0d0
      elseif ( value .lt. -1.0d0 ) then
c        --- set to pi
         arccos = dacos(-1.0d0)
      else
         arccos = dacos(value)
      endif

      return
      end
