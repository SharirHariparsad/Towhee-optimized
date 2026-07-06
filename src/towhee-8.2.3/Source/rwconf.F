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
      subroutine twh_rwconf(lfailure, nchain, rwflag, ncycle, flag)
c     ******************************************************************
c     * wrapper around twh_rwconf_main that gets a few constants used  *
c     * to dimension arrays and then call the main routine             *
c     *                                                                *
c     * split into a wrapper 01-16-2018 by M.G. Martin                 *
c     * last modified 01-16-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nchain,rwflag,ncycle,flag
c     --- local variables
      integer nmolty,numboxes,nboxpair,nvmmax

      call twh_nmolty(GLB_GET,nmolty)
      call twh_numboxes(GLB_GET,numboxes)
      call twh_nboxpair(GLB_GET,nboxpair)
      call twh_nvmmax(GLB_GET,nvmmax)

      call twh_rwconf_main(lfailure,nchain,rwflag,ncycle,flag
     &     ,nmolty,numboxes,nboxpair,nvmmax)

      return
      end

      subroutine twh_rwconf_main(lfailure, nchain, rwflag, ncycle, flag
     &     ,nmolty,numboxes,nboxpair,nvmmax)
c     ******************************************************************
c     * this program reads from or writes to a files that contains all *
c     * of the information needed to use as an initial conformation    *
c     * nchain is the number of molecules to read or write             *
c     * rwflag: is either READ_FLAG or WRITE_FLAG                      *
c     * ncycle: is the total number of cycles that have passed         *
c     * flag: if rwflag is write then this determines the file         *
c     *           RWCONF_FINAL then write to towhee_final              *
c     *           RWCONF_BACKUP then write to towhee_backup            *
c     *           RWCONF_RESTART then write to towhee_restart_xxxxx    *
c     *              which is a file just like towhee_backup but which *
c     *              does not get overwritten every time               *
c     *           if rwflag is RW_READ then 'flag'is the atomcount     *
c     *                                                                *
c     * rewritten 04-13-2004 by M.G. Martin from writeconf & readtowhee*
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer nchain,rwflag,ncycle,flag,nmolty,numboxes,nboxpair,nvmmax
c     --- local variables
c     --- character string scalars
      character*14 stepchar
      character*30 filename
c     --- logical scalars
      logical lvalue,lsaverng,lnegative,lduplicate
c     --- integer scalars
      integer ibox,imolty,iunit,ichain,iversion,idim,jdim
      integer nnstep,itest,ncount,istart,iend,ntest,ipart,nboxi,jtest
      integer activechain,iadjust,currentversion,nmtype
      integer ivalue,jvalue,nunit,nchbox,ncmt
      integer globalpos,itemp,ensemble,rng_code
      integer input_code,ivolmove,ipair
      integer ivolstart,ivolend,testboxi
c     --- integer arrays
      integer nu_temp
      dimension nu_temp(nmolty)
      integer tempstore
      dimension tempstore(nchain)
      integer localncmt
      dimension localncmt(MAXBOX)
c     --- 32-bit integer arrays
      integer*4 nrng, irng
      dimension irng(1600)
c     --- double precision scalars
      double precision xtemp,ytemp,ztemp,rmtraa,rmtrac,rmrot,boxvolume
      double precision rmvol
c     --- double precision arrays
      double precision dtemp
      dimension dtemp(nvmmax)
      double precision fatmp
      dimension fatmp(3)

c     --- retrieve constants
      call twh_ensemble(GLB_GET,ensemble)

      if ( rwflag .eq. WRITE_FLAG ) then
         if ( flag .eq. RWCONF_FINAL ) then
            filename = 'towhee_final'
         elseif ( flag .eq. RWCONF_BACKUP ) then
            filename = 'towhee_backup'
         elseif ( flag .eq. RWCONF_RESTART ) then
c           --- turn ncycle into a character string
            write(stepchar,'(i14.14)') ncycle
            filename = 'towhee_restart_'//stepchar
         else
            write(6,*) 'Unknown flag passed to rwconf: ', flag
            stop
         endif
c        --- open the file
         call twh_openfile(lfailure,8,filename,'DEFAULT','UNKNOWN'
     &        ,'DEFAULT')
         if ( lfailure ) return
      elseif ( rwflag .eq. READ_FLAG ) then
c        --- open the file
         filename = 'towhee_initial'
         call twh_openfile(lfailure,8,filename,'DEFAULT','OLD'
     &        ,'DEFAULT')
         if ( lfailure ) return
      else
         write(6,*) 'RWCONF: unknown rwflag:',rwflag
         lfailure = .true.
         return
      endif

c     --- set the currentversion
      currentversion = 8

      lnegative = .false.
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- current version number
         iversion = currentversion
         write(8,*) iversion
      else
         write(6,*) 'Reading in initial conformation from towhee_inital'
c        --- read in the version number of this file
         read(8,*) iversion
         write(6,'(a16,i5)') 'Initial version:',iversion
c        --- make sure we have a legal iversion
         if ( iversion .lt. 0 ) then
            lnegative = .true.
            iversion = abs(iversion)
         endif
         if ( iversion .gt. currentversion ) then
            write(6,*) 'RWCONF: unknown version in towhee_initial'
            lfailure = .true.
            return
         endif
      endif

c     --- random number generator information
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- output the integer code for the random number generator
         call twh_random_code(GLB_GET,rng_code)
         write(8,*) rng_code
c        --- get the sequence of nrng 32-bit integers required for a 
c        --- clean restart of the random number generator
         if ( rng_code .eq. RNG_RANLUX ) then
            call twh_random_luxlevel(GLB_GET,irng(1))
         endif
         call rngsave(rng_code,nrng,irng)
         write(8,*) nrng
         write(8,*) (irng(itest),itest=1,nrng)
      else
         if ( iversion .ge. 5 ) then
            if ( iversion .eq. 5 .or. iversion .eq. 6 ) then
c              --- only had the lux random number generator
c              --- read in the sequence of 32-bit integers
               rng_code = RNG_RANLUX
            else
c              --- read in the random number generator code
               read(8,*) rng_code
            endif
c           --- check this against the rng_code in towhee_input
            call twh_random_code(GLB_GET,input_code)
            if ( rng_code .ne. input_code ) then
c              --- mismatch of input and initial rng_code
               write(6,'(a,1x,i5,1x,a,1x,i5)') 
     &              'towhee_input random number generator'
     &              ,input_code
     &              ,'overrides towhee_initial random number generator'
     &              ,rng_code
               lsaverng = .false.
c              --- initialize using the towhee_input rng
               call twh_init_random_single(lfailure)
               if ( lfailure ) return
            else
               lsaverng = .true.
            endif
c           --- read in the sequence of 32-bit integers
            if ( iversion .eq. 6 .and. lnegative ) then
c              --- special version that was output just for 6.2.16
c              --- through 6.2.19
               read(8,*) nrng
               istart = 1
            elseif ( iversion .eq. 5 .or. iversion .eq. 6 ) then
c              --- always used 26 integers for ranlux
               nrng = 26
c              --- store the current luxlevel in first position
               call twh_random_luxlevel(GLB_GET,irng(1))
               istart = 2
            else
c              --- starting with version 7 read from file
               read(8,*) nrng
               istart = 1
            endif
c           --- read the integer restart array
            read(8,*) (irng(itest),itest=istart,nrng)
            if ( iversion .eq. 6 .or. lnegative ) then
c              --- luxlevel was stored in the first position
               call twh_random_luxlevel(GLB_SET,irng(1))
c              --- shift down
               do istart = 1,nrng-1
                  irng(istart) = irng(istart+1)
               enddo
c              --- shift down nrng
               nrng = nrng - 1
            endif

c           --- if restarts are allowed and we have a matching RNG code
c           --- then initialize the RNG
            call twh_random_allow_restart(GLB_GET,lvalue)
            if ( lvalue .and. lsaverng ) then
               call twh_restart_random(lfailure,rng_code,nrng,irng)
               if ( lfailure ) return
            endif
         else
c           --- there are not any old values to restart with, use the
c           --- random_seed approach if restarts were allowed
            call twh_random_allow_restart(GLB_GET,lvalue)
            if ( lvalue ) then
               call twh_init_random_single(lfailure)
               if ( lfailure ) return
            endif
         endif
      endif

      if ( rwflag .eq. WRITE_FLAG ) then
         write(8,*) ncycle,numboxes,nmolty
      else
c        --- nnstep is in the initial file just for use in seeing how 
c        --- far a backup file in case of a crash
         if ( iversion .le. 5 ) then
c           --- older versions did not have nmolty here
            read(8,*) nnstep,itest
            jtest = nmolty
         else
            read(8,*) nnstep,itest,jtest
         endif
         if ( itest .ne. numboxes ) then
            write(6,*) 'RWCONF: boxes in towhee_inital does not'
     &           ,' match numbox'
            lfailure = .true.
            return
         endif
         if ( jtest .ne. nmolty ) then
            write(6,*)
     &           'RWCONF: conflicting information in towhee_initial'
     &           ,' and towhee_input'
            write(6,*) 'nmolty towhee_input',nmolty
            write(6,*) 'nmolty towhee_initial',jtest
            lfailure = .true.
            return
         endif

      endif

c     --- translational and rotational maximum displacements
      if ( rwflag .eq. WRITE_FLAG ) then
         do ibox = 1, numboxes
            do imolty = 1, nmolty
               call twh_rmtraa(GLB_GET,imolty,ibox,rmtraa)
               write(8,*) rmtraa
               call twh_rmtrac(GLB_GET,imolty,ibox,rmtrac)
               write(8,*) rmtrac
               call twh_rmrot(GLB_GET,imolty,ibox,rmrot)
               write(8,*) rmrot
            enddo
         enddo
      else
c        --- read in the maximum displacements
         do ibox = 1,numboxes
            do imolty = 1,nmolty
c              --- older version has multiple entries, just take the
c              --- first one
               read(8,*) rmtraa
               call twh_rmtraa(GLB_SET,imolty,ibox,rmtraa)
               read(8,*) rmtrac
               call twh_rmtrac(GLB_SET,imolty,ibox,rmtrac)
               read(8,*) rmrot
               call twh_rmrot(GLB_SET,imolty,ibox,rmrot)
            enddo
         enddo
         write(6,*) 'new maximum displacements read from towhee_initial'
         do ibox = 1, numboxes
            write(6,'(a4,1x,i5)')'box:',ibox
            do imolty = 1, nmolty
               write(6,'(a14,1x,i5)') 'molecule type:',imolty
               call twh_rmtraa(GLB_GET,imolty,ibox,rmtraa)
               write(6,'(1x,a36,1x,f10.6)')
     &              'Max displacement for Atom translate:',rmtraa
               call twh_rmtrac(GLB_GET,imolty,ibox,rmtrac)
               write(6,'(1x,a35,1x,f10.6)')
     &              'Max displacement for COM translate:',rmtrac
               call twh_rmrot(GLB_GET,imolty,ibox,rmrot)
               write(6,'(1x,a30,1x,f10.6)')
     &              'Max displacement for rotation:',rmrot
            enddo
         enddo
      endif

c     --- volume moves maximum displacements
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- output isotropic volume move maximum displacements
         do ivolmove = 1,nvmmax
            call twh_rmvol(GLB_GET,ivolmove,dtemp(ivolmove))
         enddo
         write(8,*) (dtemp(ivolmove),ivolmove=1,nvmmax)
      else
c        --- read in the standard volume move max displacements
         if ( iversion .lt. 8 ) then
c           --- older versions used either the number of boxes or
c           --- the number of pairs of boxes depending on the ensemble
            if ( ensemble .eq. ENS_NPT ) then
               ivolstart = 1
               ivolend = numboxes
            else
               ivolstart = numboxes + 1
               ivolend = numboxes + nboxpair
            endif
            lduplicate = .true.
         else
c           --- version 8 and later uses the sum of boxes and pairs
            ivolstart = 1
            ivolend = nvmmax
            lduplicate = .false.
         endif

         read(8,*) ( dtemp(ivolmove),ivolmove=ivolstart,ivolend)
         do ivolmove = ivolstart,ivolend
            call twh_rmvol(GLB_SET,ivolmove,dtemp(ivolmove))
            if ( lduplicate ) then
               if ( ivolmove .le. numboxes ) then
c                 --- also store in the pair portion if room
                  if ( ivolmove .le. nboxpair ) then
                     iadjust = ivolmove + numboxes
                     call twh_rmvol(GLB_SET,iadjust,dtemp(ivolmove))
                  endif
               else
c                 --- also store in the numboxes portion, if room
                  if ( ivolmove - numboxes .le. numboxes ) then
                     iadjust = ivolmove - numboxes
                     call twh_rmvol(GLB_SET,iadjust,dtemp(ivolmove))
                  endif
               endif
            endif
         enddo

c        --- output new maximum displacements
         write(6,'(a54)') 
     &        'Max displacement for single box isotropic volume moves'
         do ivolmove = 1,numboxes
            call twh_rmvol(GLB_GET,ivolmove,rmvol)
            write(6,'(1x,a4,1x,i2,1x,a11,1x,e12.4)')
     &           'Box:',ivolmove,'rmvol(box):',rmvol
         enddo
         write(6,'(a53)') 
     &        'Max displacement for box pairs isotropic volume moves'
         do ivolmove = numboxes+1,numboxes+nboxpair
            call twh_rmvol(GLB_GET,ivolmove,rmvol)
            write(6,'(1x,a5,1x,i2,1x,a11,1x,e12.4)')
     &           'Pair:',ivolmove-numboxes,'rmvol(box):',rmvol
         enddo
      endif

      if ( rwflag .eq. WRITE_FLAG ) then
c        --- output unit cell move maximum displacements
         do ivolmove = 1,nvmmax
            do idim = 1,3
               call twh_rmcell(GLB_GET,ivolmove,idim,1,xtemp)
               call twh_rmcell(GLB_GET,ivolmove,idim,2,ytemp)
               call twh_rmcell(GLB_GET,ivolmove,idim,3,ztemp)
               write(8,*) xtemp,ytemp,ztemp
            enddo
         enddo
      else
         lduplicate = .false.
c        --- account for previous differences in the rmcell
         if ( iversion .le. 1 ) then
c           --- these displacments did not exist until version 2
c           --- skip reading them, and leave them set to the defaults
c           --- which were set up in the pmcell section
            ivolstart = 1
            ivolend = 0
         elseif ( iversion .eq. 2 ) then
c           --- version 2 outputted these for either the number of
c           --- pairs of boxes (NVT) or the number of boxes (other).
c           --- this led to problems when switching between nvt and npt
            if ( ensemble .eq. ENS_NVT ) then
c              --- only pairs of boxes were stored
               ivolstart = numboxes + 1
               ivolend = numboxes + nboxpair
            else
c              --- only single box information was stored
               ivolstart = 1
               ivolend = numboxes
            endif
            lduplicate = .true.
         elseif ( iversion .le. 7 ) then
c           --- version 3-7 uses the larger of numboxes or nboxpair so
c           --- it is possible to switch between nvt and npt
            if ( numboxes .gt. nboxpair ) then
               ivolstart = 1
               ivolend = numboxes
            else
               ivolstart = numboxes + 1
               ivolend = numboxes + nboxpair
            endif
            lduplicate = .true.
         else
c           --- version 8 introduced a complete set of volume
c           --- displacements for each box and each box pair
c           --- and changed the definition of nvmmax to be the sum
c           --- of the number of boxes plus the number of box pairs
            ivolstart = 1
            ivolend = nvmmax
            lduplicate = .false.
         endif

c        --- read in the unit cell maximum displacements
         do ivolmove = ivolstart,ivolend
            do idim = 1,3
               read(8,*) xtemp,ytemp,ztemp
               call twh_rmcell(GLB_SET,ivolmove,idim,1,xtemp)
               call twh_rmcell(GLB_SET,ivolmove,idim,2,ytemp)
               call twh_rmcell(GLB_SET,ivolmove,idim,3,ztemp)
               if ( lduplicate ) then
                  if ( ivolmove .le. numboxes ) then
                     iadjust = ivolmove + numboxes
                  else
                     iadjust = ivolmove - numboxes
                  endif
                  if ( iadjust .gt. 0 .and. iadjust .lt. nvmmax ) then
                     call twh_rmcell(GLB_SET,ivolmove,idim,1,xtemp)
                     call twh_rmcell(GLB_SET,ivolmove,idim,2,ytemp)
                     call twh_rmcell(GLB_SET,ivolmove,idim,3,ztemp)
                  endif
               endif
            enddo
         enddo

c        --- standard output the new (or current) max unit cell disp
         write(6,'(a36)') 'Max disp. for unit cell perturbation'
         do ivolmove = 1,nvmmax
            do idim = 1,3
               call twh_rmcell(GLB_GET,ivolmove,idim,1,xtemp)
               call twh_rmcell(GLB_GET,ivolmove,idim,2,ytemp)
               call twh_rmcell(GLB_GET,ivolmove,idim,3,ztemp)
               if ( ivolmove .le. numboxes ) then
c                 --- individual box values
                  ibox = ivolmove
                  write(6,1001) ibox,idim,xtemp,ytemp,ztemp
               else
c                 --- pair box values
                  ipair = ivolmove - numboxes
                  call twh_pairbox(GLB_GET,ipair,1,ivalue)
                  call twh_pairbox(GLB_GET,ipair,2,jvalue)
                  write(6,1002) ivalue,jvalue,idim,xtemp,ytemp,ztemp
               endif
            enddo
         enddo
 1001    format(1x,' Box ',i3,' idim ',i1,' rmcell: ',3e12.4)
 1002    format(1x,' Boxes ',i3,' and ',i3,' idim ',i1
     &        ,' rmcell: ',3e12.4)
      endif

c     --- hmatrices
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- output box hmatrices
         do ibox = 1,numboxes
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_GET,BXD_REAL,ibox,idim,jdim
     &                 ,fatmp(jdim))
               enddo
               write(8,*) (fatmp(jdim),jdim=1,3)
            enddo
         enddo
      else
         write(6,*)
c        --- initialize hmatrix
         call twh_hmatrix(GLB_INIT,MAXBXD,numboxes,3,3,0.0d0)
         do ibox = 1,numboxes
            if ( iversion .le. 1 ) then
c              --- use the old boxlength notation
c              --- just read in the diagonal terms of hmatrix
               read(8,*) (fatmp(idim),idim=1,3)
               do idim = 1,3
                  jdim = idim
                  call twh_hmatrix(GLB_SET,BXD_REAL,ibox,idim,jdim
     &                 ,fatmp(idim))
               enddo
            else
c              --- read in the hmatrix
               do idim = 1,3
                  read(8,*) (fatmp(jdim),jdim=1,3)
                  do jdim = 1,3
                     call twh_hmatrix(GLB_SET,BXD_REAL,ibox,idim,jdim
     &                    ,fatmp(jdim))
                  enddo
               enddo
            endif
         enddo
         write(6,*) 'new box dimensions read from towhee_initial'
         do ibox = 1,numboxes
            do idim = 1,3
               do jdim = 1,3
                  call twh_hmatrix(GLB_GET,BXD_REAL,ibox,idim,jdim
     &                 ,fatmp(jdim))
               enddo
               write (6,1104) ibox,idim,(fatmp(jdim),jdim=1,3)
            enddo
c           --- error check box dimensions
            call twh_checkhmatrix(BXD_REAL,ibox,.true.,lfailure)
            if ( lfailure ) return
c           --- invert this matrix
            call twh_inverthmatrix(BXD_REAL,ibox)
c           --- make sure the volumes are positive
            call twh_boxvolume(GLB_GET,BXD_REAL,ibox,boxvolume)
            if ( boxvolume .lt. 0.0d0 ) then
               write(6,*) 'RWCONF: Box ',ibox
     &              ,' has a negative volume of ',boxvolume
               write(6,*) 'this means you have set up a left-handed '
     &              ,'coordinate system'
               lfailure = .true.
               return
            endif
         enddo
 1104    format(' Box ',i2,' hmatrix(',i1,',x): ',3f14.5)
         write(6,*)
      endif

c     --- number of chains
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- we only want to output the number of chains that are 
c        --- not in the ideal chain box
         activechain = 0
         do ibox = 1,numboxes
            call twh_nchbox(GLB_GET,ibox,nchbox)
            activechain = activechain + nchbox
         enddo
         write(8,*) activechain
      else
c        --- check that number of particles in towhee_input and 
c        --- towhee_initial agree 
         read(8,*) activechain
         if ( ( ensemble .ne. ENS_muVT .and. activechain .ne. nchain ) 
     &        .or.
     &        ( ensemble .eq. ENS_muVT .and. activechain .gt. nchain)
     &        ) then
            write(6,*)
     &           'RWCONF: conflicting information in towhee_initial'
     &           ,' and towhee_input'
            write(6,*) 'nchain towhee_input',nchain
            write(6,*) 'nchain towhee_initial',itest
            lfailure = .true.
            return
         endif
      endif

c     --- number of molecule types
c     --- compatability with older versions
      if ( iversion .le. 5 .and. rwflag .eq. READ_FLAG ) then
c        --- check that the number of molecule types in towhee_input
c        --- and towhee_initial agree
         read(8,*) jtest
         if ( jtest .ne. nmolty ) then
            write(6,*)
     &           'RWCONF: conflicting information in towhee_initial'
     &           ,' and towhee_input'
            write(6,*) 'nmolty towhee_input',nmolty
            write(6,*) 'nmolty towhee_initial',jtest
            lfailure = .true.
            return
         endif
      endif

c     --- number of units per chain for each molecule type
      if ( rwflag .eq. WRITE_FLAG ) then
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nu_temp(imolty))
         enddo
         write(8,*) (nu_temp(imolty),imolty=1,nmolty)
      else
         read(8,*) (nu_temp(imolty),imolty=1,nmolty)
         do imolty = 1,nmolty
            call twh_nunit(GLB_GET,imolty,nunit)
            if ( nu_temp(imolty) .ne. nunit ) then
               write(6,*) 'RWCONF: moltype,towhee_initial nunit '
     &              ,' towhee_input nunit'
     &              ,imolty,nu_temp(imolty),nunit
               write(6,*) 'Error: these values do not match'
               lfailure = .true.
               return
            endif
         enddo
      endif

c     --- molecule type for each chain in the data structure
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- only want to output the non-ideal box moltyp values
         itest = 0
         do ichain = 1,nchain
            call twh_nboxi(GLB_GET,ichain,testboxi)
            if ( testboxi .ne. 0 ) then
               itest = itest + 1
               call twh_moltyp(GLB_GET,ichain,tempstore(itest))
            endif
         enddo
         write(8,*) (tempstore(ichain),ichain=1,activechain)
      else
         read(8,*) (tempstore(ichain),ichain=1,activechain)
         do ichain = 1,activechain
            call twh_moltyp(GLB_SET,ichain,tempstore(ichain))
         enddo
      endif

c     --- box each molecule is in
      if ( rwflag .eq. WRITE_FLAG ) then
c        --- only want to output the non-ideal box locations
         itest = 0
         do ichain = 1,nchain
            call twh_nboxi(GLB_GET,ichain,nboxi)
            if ( nboxi .ne. 0 ) then
               itest = itest + 1
               tempstore(itest) = nboxi
            endif
         enddo
         write(8,*) (tempstore(ichain),ichain=1,activechain)
      else
         read(8,*) (tempstore(ichain),ichain=1,activechain)
         do ichain = 1,activechain
            call twh_nboxi(GLB_SET,ichain,tempstore(ichain))
         enddo
      endif

c     --- setup some of the arrays in the code
      if ( rwflag .eq. READ_FLAG ) then
c        --- obtain ncmt values
         do ibox = 0,numboxes
            call twh_nchbox(GLB_SET,ibox,0)
            do imolty = 1, nmolty
               call twh_ncmt(GLB_SET,ibox,imolty,0)
            enddo
         enddo

         do ichain = 1, activechain
c           --- make sure that a sensible box was generated
            call twh_nboxi(GLB_GET,ichain,nboxi)
            if ( nboxi .le. 0 .or. nboxi .gt. numboxes ) then
               write(6,*) 'RWCONF:Particle found in ill-defined box'
               lfailure = .true.
               return
            endif
c           --- increment the book-keeping arrays for each molecule
            call twh_nboxi(GLB_GET,ichain,ibox)
            call twh_nchbox(GLB_INCR,ibox,1)
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_ncmt(GLB_INCR,ibox,imolty,1)
         enddo

c        --- check that number of particles of each type is consistent
         ncount = activechain
         do imolty = 1, nmolty
            itest = 0
            do ibox = 1,numboxes
               call twh_ncmt(GLB_GET,ibox,imolty,ncmt)
               itest = itest + ncmt
            enddo

            call twh_nmolectyp(GLB_GET,imolty,nmtype)
            if ( ensemble .eq. ENS_muVT ) then
c              --- first check that we have a positive surplus
               if ( itest .gt. nmtype ) then
                  write(6,*) 'RWCONF: more molecules in towhee_initial'
                  write(6,*) 'than in towhee_input for type',imolty
                  lfailure = .true.
                  return
               endif
c              --- place any surplus molecules into ideal box 0
               ibox = 0
               do ntest = 1,nmtype - itest
                  iadjust = ncount + 1
                  call twh_nboxi(GLB_SET,iadjust,ibox)
                  call twh_moltyp(GLB_SET,iadjust,imolty)
                  call twh_ncmt(GLB_INCR,ibox,imolty,1)
                  call twh_nchbox(GLB_INCR,ibox,1)
                  ncount = iadjust
               enddo
            elseif ( itest .ne. nmtype ) then
c              --- testing for non uvt ensembles then it must match
               write(6,*) 'RWCONF: Particle type # inconsistency'
               write(6,*) 'type',imolty
               do ibox=1,numboxes
                  call twh_ncmt(GLB_GET,ibox,imolty,localncmt(ibox))
               enddo
               write(6,*) 'ncmt',(localncmt(ibox),ibox=1,numboxes)
               write(6,*) 'nmolectyp',nmtype
               lfailure = .true.
               return
            endif
         enddo

c        --- make sure ncount and nchain now match
         if ( ncount .ne. nchain ) then
            write(6,*) 'RWCONF: ncount and nchain do not match'
            write(6,*) 'ncount,nchain',ncount,nchain
            lfailure = .true.
            return
         endif

c        --- setup the globalpos array (two parts)
         flag = 0
         call twh_globalpos(GLB_SET,1,1)
         do ipart = 1,2
            if ( ipart .eq. 1 ) then
c              --- first part
               istart = 1
               iend = activechain
            else
               istart = activechain + 1
               iend = nchain
            endif
            do ichain = istart,iend
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               flag = flag + nunit
c              --- set the next globalpos
               itemp = ichain+1
               globalpos = 3*flag + 1
               call twh_globalpos(GLB_SET,itemp,globalpos)
            enddo
         enddo
      endif

c     --- coordinates
      if ( rwflag .eq. WRITE_FLAG ) then 
         do ichain = 1, nchain
c           --- only output the non-ideal coordinates
            call twh_nboxi(GLB_GET,ichain,testboxi)
            if ( testboxi .ne. 0 ) then
               call twh_moltyp(GLB_GET,ichain,imolty)
               call twh_nunit(GLB_GET,imolty,nunit)
               do iunit = 1, nunit
                  call twh_coordinates(GLB_GET,CRD_REAL,ichain,iunit
     &                 ,xtemp,ytemp,ztemp)
                  write (8,*) xtemp,ytemp,ztemp
               enddo
            endif
         enddo
      else
         do ichain = 1, activechain
            call twh_moltyp(GLB_GET,ichain,imolty)
            call twh_nunit(GLB_GET,imolty,nunit)
            do iunit = 1, nunit
               if ( iversion .eq. 0 ) then
c                 --- version 0 included charges on each atom
c                 --- this is no longer needed
                  read(8,*) xtemp,ytemp,ztemp,itest
               else
c                 --- later versions do not include charges
                  read(8,*) xtemp,ytemp,ztemp
               endif
c              --- assign the positions
               call twh_coordinates(GLB_SET,CRD_REAL,ichain,iunit
     &              ,xtemp,ytemp,ztemp)
            enddo
         enddo
      endif

c     --- rewind the file (only needed for backup but doesn't hurt)
      rewind(8)
      close(8)

      return
      end
