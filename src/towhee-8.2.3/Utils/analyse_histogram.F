      program reweight
c     ******************************************************************
c     * this program determines the weights used for patching multiple *
c     * histograms together. The Ferrenberg/Swendsen method outlined   *
c     * in Phys, Rev. Lett. 63,1195 (1989) is used to determine the    *
c     * weights.  This code has been modified from a program           *
c     * originally written by Jeff Potoff on 04/29/1998.               *
c     *                                                                *
c     * originally written 04-19-2004 by C.D. Lorenz                   *
c     * last modified 06-28-2004 by C.D. Lorenz                        *
c     ******************************************************************

      implicit none

c     --- parameters
      integer MAXCOMP,MAXVECTOR,MAXN,MAXFILES,MAXENTRY,MAXTEMP,MAXISING
      integer MAXCRIT,NUMBIN,MAXITER
      double precision MAXTOL,MAXPRESSTOL,CONV2BAR

c     MAXCOMP == the maximum number of components in a system
c     MAXVECTOR == the maximum number of optimizing dimensions used
c                  to determine critical properties (2*MAXCOMP)+1
c     MAXN == the maximum number of molecules of any component
c     MAXFILES == the maximum number of independent histogram files to be
c                 patched together
      parameter(MAXCOMP=10,MAXVECTOR=21,MAXN=3000,MAXFILES=30)
c     MAXENTRY == the maximum number of total entries in all files combined
c     MAXTOL == the tolerance for convergence of the weights
c     MAXPRESSTOL == the tolerance for the convergence of the pressures
      parameter(MAXENTRY=600000,MAXTOL=5.0e-5,MAXPRESSTOL=1.0e-3)
c     MAXTEMP == the maximum number of temperature entries in the phase 
c                coexistence data file 
c     MAXISING == the maximum number of data points in the ising.dat file
c     MAXCRIT == the maximum number of bins in the discretization of the 
c                field mixing parameter
c     NUMBIN == the number of energy bins used in critical point 
c               distribution calculations 
      parameter(MAXTEMP=100,MAXISING=1000,MAXCRIT=1000,NUMBIN=1001)
c     MAXITER == the maximum number of iterations allowed during simulation
c     CONV2BAR == the conversion to pressure in bar 
      parameter(MAXITER=5000,CONV2BAR=138.0d0)
  
c      common /constant1/ MAXCOMP,MAXVECTOR
c      common /constant1a/ MAXN,MAXFILES
c      common /constant2/ MAXENTRY,MAXTEMP
c      common /constant2a/ MAXISING,MAXCRIT
c      common /constant3/ MAXITER,NUMBIN
c      common /constant4/ MAXTOL,MAXPRESSTOL,CONV2BAR

   
c     --- global variables for reweight.f
      integer nmol,nentry,ntotentry,ntotfiles,ncomp,nising,num_part
      dimension nmol(MAXCOMP,MAXENTRY),nentry(MAXFILES)
      dimension num_part(2,MAXCOMP)

      double precision beta,energy,chempot,slope,oldweight,weight
      double precision chempotorig
      double precision forder,xising,yising,dens,avgnum,nmol_liq
      double precision nmol_gas,midpt
      dimension oldweight(MAXFILES),beta(MAXFILES),slope(MAXTEMP)
      dimension chempot(MAXCOMP,MAXFILES),energy(MAXENTRY)
      dimension chempotorig(MAXCOMP,MAXFILES)
      dimension midpt(MAXTEMP),nmol_gas(MAXCOMP),nmol_liq(MAXCOMP)
      dimension xising(MAXISING),yising(MAXISING),avgnum(MAXCOMP)
      dimension dens(MAXCOMP,0:MAXN),forder(0:MAXCRIT),weight(MAXFILES)

      logical lweight,lphase,lcritical,lpvt,lconstantp

      common /increweight/ nentry
      common /increweight1/ nmol
      common /increweight2/ num_part
      common /increweight3/ ntotentry,ntotfiles,ncomp,nising
      common /increweight4/ oldweight,beta,weight
      common /increweight5/ slope,midpt
      common /increweight6/ nmol_gas,nmol_liq,avgnum
      common /increweight7/ chempot
      common /increweight8/ energy
      common /increweight9/ xising
      common /increweight10/ yising
      common /increweight11/ dens
      common /increweight12/ forder
      common /increweight13/ lweight,lphase,lcritical,lpvt,lconstantp
c     --- local variables
      integer nfiles,ifile,noldfiles,icomp
      integer ncompi,ientry,niter,jfile,ndump,ntemp
      integer itemp,vector_dim,ntemp_crit,niter2
      integer nflip,oldsign,icount,newsign,inum
      integer ipart,jcomp,idim,jdim
  
      double precision maxdev,iter_dev,ylog,prob,ylogtmp
      double precision chempoti,chempot2min,chempot2max,chempot2incr
      double precision temperature,min_energy,max_energy
      double precision trial_vector,temp_crit,chempot_crit,fieldmix
	double precision chempot_crit_orig
      double precision volume,logz,pset,chempot_pvt,ncrit
      double precision betai,chempotincr1,chempotincrn,pressdiff
      double precision nbelow,nbelow1,pressure,sign,energy_gas
      double precision avgn,part_function,opt_in,opt_try
      double precision opt_y,xmin,xmax,xavg,energy_liq,chempotincr
      double precision part_function_gas,part_function_liq
      double precision variance,critvalue
      double precision chempot_deBrog,mass,debroglie
      
      dimension chempoti(MAXCOMP),chempot2min(MAXTEMP)
      dimension chempot_deBrog(MAXCOMP)
      dimension chempot2max(MAXTEMP),chempot2incr(MAXTEMP)
      dimension temperature(MAXTEMP),trial_vector(MAXVECTOR)
      dimension temp_crit(MAXTEMP),chempot_crit(MAXCOMP,MAXTEMP)
	dimension chempot_crit_orig(MAXCOMP,MAXTEMP)
      dimension fieldmix(MAXCOMP,MAXTEMP),ncrit(MAXTEMP)
      dimension chempot_pvt(MAXCOMP,MAXTEMP)
      dimension opt_in(MAXVECTOR+1,MAXVECTOR)
      dimension opt_try(MAXVECTOR),opt_y(MAXVECTOR+1)
      dimension chempotincr(MAXTEMP)
      dimension mass(MAXCOMP,MAXFILES)

      logical lold,lnew,lauto

      character*30 filename,temp_file,file_ising
      character*30 ftemp,ftemp2,fphase,fphase2,fname2,fname
      character*30 ftemp3,file_old
      character*30 filephaseout,filepvt2
      dimension filename(MAXFILES)

c     --- Read in input script for reweight.F
      open (1,file='input_analyse_histogram')

c     --- Initializing values
      lauto = .false.
      lold = .false.
      lweight = .false. 
      lphase = .false.
      lcritical = .false.
      lnew = .false.
      lpvt = .false.
      lconstantp = .false.

      noldfiles = 0
      nfiles = 0
      do icomp = 1,MAXCOMP
         num_part(1,icomp) = MAXN+1
         num_part(2,icomp) = -1
      enddo
      min_energy = 1.0e6
      max_energy = -1.0e6


      read(1,*)
c     --- ncomp is the number of components in the system
      read(1,*) ncomp
      read(1,*)
c     --- volume is the volume of the system
      read(1,*) volume
c     volume input in nm^3 used in Angstroms^3
      volume = volume * 1000.0d0
      read(1,*)       
c     --- lweight is the logical which tells if weights are to be determined
c     --- lweight = .true. determine weights, .false. don't
      read(1,*) lweight
      if (lweight) then
c        --- ndump is the frequency at which the iteration info 
c        --- and weights will be dumped during the iteration
         read(1,*)
         read(1,*) ndump
      endif
      if (.not. lweight .and. ncomp .gt. 2) then
         write(6,*) 'Outside of the determining of the weights,',
     &        ' this code has not been tested for systems with more',
     &        ' than two components.'
         stop
      endif
      read(1,*)
c     --- lphase is the logical which tells if phase coexistence of the 
c     --- system is to be calculated
c     --- lphase = .true. calculate phase coexistence, .false. don't
      read(1,*) lphase 

      read(1,*) 
c     --- lpvt is the logical which tells whether or not to perform 
c     --- pvt calculations
      read(1,*) lpvt     
      read(1,*)
c     --- error check: not able to do phase coexistence and PVT at the 
c     --- same time
      if (lphase .and. lpvt) then
         write(6,*) ' Can not do phase coexistence and PVT',
     &        ' calculations at the same time.'
         stop
      endif
c     --- lcritical is the logical which tells if the critical point is 
c     --- to be determined
c     --- lcritical = .true. calculate critical points, .false. don't
c      read(1,*) lcritical     
c      if (lcritical) then
c         read(1,*)
c        --- logical which tells whether or not you want it to use an automatic critical
c        --- point determination or not
c         read(1,*) lauto
c         read(1,*)
c        --- nising is the number of points in the ising distribution
c         read(1,*) nising
c        --- error check
c         if (nising .gt. MAXISING) then
c            write(6,*) 'the number of points in the ising distribution',
c     &          ' (',nising,') is greater than MAXISING (',MAXISING,').'
c            write(6,*) 'Increase MAXISING in reweight.F'
c            stop
c         endif
c         read(1,*)
c        --- file_ising is the name of the file which includes the ising
c        ---  distribution
c         read(1,*) file_ising
c         read(1,*)
c         vector_dim = (2*ncomp)+1
c        --- error check
c         if (vector_dim .gt. MAXVECTOR) then
c            write(6,*) 'vector dimension (',vector_dim,') is ',
c     &           'greater than MAXVECTOR (',MAXVECTOR,').'
c            write(6,*) 'Increase MAXVECTOR in reweight.F'
c            stop
c         endif
c        --- trial vector
c         read(1,*) (trial_vector(idim),idim=1,vector_dim)
c         read(1,*)
c        --- the number of temperatures used in determination of 
c        --- critical point
c         read(1,*) ntemp_crit
c         read(1,*)
c         do itemp = 1,ntemp_crit
c           --- read in temperatures, measured chemical potentials, 
c           --- field mixing parameters, and the number of points to 
c           --- use for the order parameter distribution
c            read(1,*) temp_crit(itemp),
c     &           (chempot_crit_orig(icomp,itemp),icomp=1,ncomp),
c     &           (fieldmix(icomp,itemp),icomp=1,ncomp),ncrit(itemp)
c           --- error check
c            if (ncrit(itemp) .gt. MAXCRIT) then
c               write(6,*) 'ncrit for temperature #',itemp,' is ',
c     &              'greater than MAXCRIT'
c               stop
c            endif
c         enddo   
c      endif
c      read(1,*) 
c     --- lold is a logical which will be true if there are old
c     --- histograms which already have determined weights in the 
c     --- data set
      read(1,*) lold
      if (lold) then
         read(1,*)
c        --- noldfiles is the number of files for which weights have 
c        --- been previously determined
         read(1,*) noldfiles
         if ((noldfiles) > MAXFILES) then
            write(6,*) 'noldfiles (',noldfiles,
     &           ') is greater than MAXFILES (',MAXFILES,')'
            write(6,*) 'Increase MAXFILES in reweight.F'
            stop
         endif
         read(1,*) 
         read(1,*) file_old
      endif
c     --- lnew is a logical which will be true if there are new
c     --- histograms which need to have the weights determined
      read(1,*)
      read(1,*) lnew
c     --- error check: if lhist == true and lnew == false, then no real
c     --- need to determine new weights when not patching any more 
c     --- information
c      if (lweight .and. (.not. lnew)) then
c         write(6,*) 'If you do not have any new histogram data, then',
c     &        ' no need to recalculate the weights.'
c         write(6,*) 'If you do have new histogram data, then change',
c     &        ' lnew in input_reweight.'
c         stop
c      endif
      if (lnew) then
         read(1,*)
c        --- read in the number of new files that have not had weights 
c        --- determined yet
         read(1,*) nfiles        
         if ((nfiles+noldfiles) > MAXFILES) then
            write(6,*) 'nfiles (',nfiles+noldfiles,
     &           ') is greater than MAXFILES (',MAXFILES,')'
            write(6,*) 'Increase MAXFILES in reweight.F'
            stop
         endif
c        --- read in the name and the number of entries in each of the 
c        --- new files
         do ifile = noldfiles+1,noldfiles+nfiles
            read(1,*)
            read(1,*) filename(ifile)
            read(1,*) 
            read(1,*) nentry(ifile)
         enddo
      endif
      close(1)
      ntotfiles = nfiles + noldfiles

c     ------------------------------------------------------------------------
c     --- Read in data files
      if (lold) then
c        --- read in data for files which have already had weights determined
         open(14,file=file_old)
         read(14,*)
         read(14,*)
         read(14,*)
         read(14,*)
         read(14,*)
         write(6,*) 'Reading in data from old files'
         do ifile = 1,noldfiles
            read(14,*)
            read(14,*) filename(ifile)
            write(6,*) 'Reading in data from ',filename(ifile)
            read(14,*) 
            read(14,*) nentry(ifile)
            read(14,*) 
            read(14,*) oldweight(ifile)
            read(14,*) 
            read(14,*) beta(ifile)
            read(14,*)
            read(14,'(10f10.2)')(chempotorig(icomp,ifile),icomp=1,ncomp)
            beta(ifile) = 1.0d0/beta(ifile)
         enddo
         close(14)
      endif
      ntotentry = 0
c     --- read in histogram from each of the files
      write(6,*) 'Reading in histogram information'
      do ifile = 1,ntotfiles
         temp_file = filename(ifile)
         write(6,*) 'Reading in data from ',temp_file
         open(10,file=temp_file)
c        --- read in the header of the histogram file which contains beta,
c        --- the number of components in the system
c        --- the chemical potential for each component and
c        --- the mass for each component
         read(10,*) beta(ifile),ncompi,
     &	(chempot(icomp,ifile),icomp=1,ncompi),
     &	(mass(icomp,ifile),icomp=1,ncompi)
         beta(ifile) = 1.0d0/beta(ifile)
         do icomp=1,ncompi
            chempotorig(icomp,ifile) = chempot(icomp,ifile)
            debroglie = 17.458/sqrt(mass(icomp,ifile)/beta(ifile))
            chempot(icomp,ifile) = chempot(icomp,ifile)
     &							-3.0*log(debroglie)/beta(ifile)
         enddo

         if (ncomp .ne. ncompi) then
            write(6,*) 'Number of components in input_reweight.dat',
     &           ' (',ncomp,') and'
            write(6,*) 'in histogram #',ifile,' (',ncompi,') do',
     &           ' not match.'
            stop
         endif
c        --- read in the histogram data from the file
         do ientry = 1,nentry(ifile)
            ntotentry = ntotentry + 1
            if (ntotentry .gt. MAXENTRY) then
               write(6,*) 'Maximum number of data entries has been',
     &              ' exceeded',ntotentry,MAXENTRY
               write(6,*) 'Increase MAXENTRY in reweight.F'
               stop
            endif
            read(10,*) (nmol(icomp,ntotentry),icomp=1,ncompi),
     &           energy(ntotentry)
            min_energy = min(min_energy,energy(ientry))
            max_energy = max(max_energy,energy(ientry))
            do icomp = 1,ncomp
               num_part(1,icomp) = min(nmol(icomp,ntotentry),
     &              num_part(1,icomp))
               num_part(2,icomp) = max(nmol(icomp,ntotentry),
     &              num_part(2,icomp))
            enddo
         enddo
         if (ifile .gt. noldfiles) then
            oldweight(ifile) = dble(nentry(ifile))
         endif
         close(10)
      enddo

      if (lphase) then
c     --- read in data for phase coexistence calculations
         open (11,file='file_phase') 
         write(6,*) 'Reading information for phase coexistence ',
     &        'calculations from file_phase'
         read(11,*)
         if (lconstantp) then
            read(11,*) chempoti(1)
            read(11,*) 
            read(11,*) logz
            read(11,*) 
            read(11,*) pset             
            pset = pset/100.0d0
c            if (ncomp .gt. 2) then
c               read(11,*) chempoti(1),chempoti(ncomp-1),logz,pset
c            else
c               read(11,*) chempoti(1),chempoti(ncomp),logz,pset
c            endif
         else
c           --- for phase coexsistence calculations if not constant pressure, 
c           --- read in initial guess
            read(11,*) chempoti(1)
         endif
         read(11,*)
         read(11,*) ntemp
c        --- error check
         if (ntemp .gt. MAXTEMP) then
            write(6,*) 'Number of temperatures sampled in phase',
     &           ' coexistence data file is greater than MAXTEMP.'
            write(6,*) 'Increase MAXTEMP in reweight.F to ',ntemp
            stop
         endif
         if (ncomp .eq. 1) then
            do itemp = 1,ntemp 
               read(11,*)
               read(11,*) temperature(itemp)
               read(11,*) 
               read(11,*) midpt(itemp)
               read(11,*) 
               read(11,*) slope(itemp)
            enddo
            debroglie = 17.458/sqrt(mass(1,1)*temperature(1))
            chempoti(1) = chempoti(1)-3.0*temperature(1)*log(debroglie)
         else
            do itemp = 1,ntemp 
               read(11,*)
               read(11,*) temperature(itemp)
               read(11,*) 
               read(11,*) midpt(itemp)
               read(11,*) 
               read(11,*) slope(itemp)
               read(11,*) 
               read(11,*) chempot2min(itemp)
               read(11,*) 
               read(11,*) chempot2max(itemp)
               read(11,*) 
               read(11,*) chempot2incr(itemp)
               debroglie = 17.458/sqrt(mass(1,1)*temperature(itemp))
               chempot2min(itemp) = chempot2min(itemp)-
     &					3.0*temperature(itemp)*log(debroglie)
               chempot2max(itemp) = chempot2max(itemp)-
     &					3.0*temperature(itemp)*log(debroglie)
            enddo
            close(11)
         endif
      endif
      
      
      
      if (lpvt) then
c        --- read in the data for the pvt calculations
         open(12,file='file_pvt')
         write(6,*) 'Reading in information for pvt calculations',
     &        ' from file_pvt.'
         read(12,*)
         read(12,*) ntemp
c        --- error check
         if (ntemp .gt. MAXTEMP) then
            write(6,*) 'Number of temperatures sampled in pvt',
     &           ' data file is greater than MAXTEMP.'
            write(6,*) 'Increase MAXTEMP in reweight.F to ',ntemp
            stop
         endif
         do itemp = 1,ntemp
            read(12,*)
            read(12,*) temperature(itemp)
            do icomp = 1,ncomp
               read(12,*)
               read(12,*) chempot_pvt(itemp,icomp)
               debroglie = 17.458/sqrt(mass(icomp,1)*temperature(itemp))
               chempot_pvt(itemp,icomp) = chempot_pvt(itemp,icomp)
     &				 - 3.0 * temperature(itemp)*log(debroglie)
            enddo 
         enddo
         close(12)
      endif

      if (lcritical) then
c        --- read in data for critical property calculation
         open(13,file=file_ising)
         write(6,*) 'Reading in information for critical property',
     &        ' calculations from ',file_ising
         do ientry = 1,nising
            read(1,*) xising(ientry),yising(ientry)
         enddo
         close(13)
      endif


c     ------------------------------------------------------------------------
c     --- iterate to find the weights using the Ferrenberg/Swendsen method
c     --- (Note: first will just put in simple substitution that Jeff had 
c     ---        then will try to make faster with some combination of math
c     ---        methods (i.e. simple substitution + Newton method)
      if (lweight) then
         write(6,*) 'Starting determination of the weights'
         niter = 0
         maxdev = 100.0d0
         do while (maxdev .gt. MAXTOL.and. niter.lt.1e4)
            if (niter .gt. 0) then
               write(6,*) 'Iteration: ',niter,' Deviation: ',maxdev
            endif
            iter_dev = 0.0d0
            maxdev = 0.0d0
            niter = niter + 1
            do ifile = 1,ntotfiles
               weight(ifile) = 0.0d0
               do ientry = 1,ntotentry
                  ylog = -1.0e9
                  do jfile = 1,ntotfiles
                     prob = -(beta(jfile)-beta(ifile))*energy(ientry)
                     do icomp = 1,ncomp
                        prob = prob + 
     &                       (((beta(jfile)*chempot(icomp,jfile))
     &                       - (beta(ifile)*chempot(icomp,ifile)))
     &                       * nmol(icomp,ientry))
                     enddo
c                    --- use logs in order to avoid overflow                   
                     ylogtmp = prob + 
     &                    log(nentry(jfile)/oldweight(jfile))
                     ylog = max(ylog,ylogtmp) + 
     &                    log(1.0d0 + exp(-1.0d0*dabs(ylog-ylogtmp)))
                  enddo
                  weight(ifile) = weight(ifile) + exp(-1.0d0*ylog)
               enddo
               iter_dev = dabs(weight(ifile)/oldweight(ifile)-1.0d0)
               if (iter_dev .gt. maxdev) maxdev = iter_dev
            enddo
            do ifile = 2,ntotfiles
               weight(ifile) = weight(ifile)/weight(1)
               oldweight(ifile) = weight(ifile)
            enddo
            weight(1) = 1.0d0
            oldweight(1) = weight(1)
            if (mod(niter,ndump) .eq. 0) then
               open(unit=20, file = 'converge.dat',status='unknown')
               write(20,'(a12,2x,i5,5x,a12,2x,f12.8)') 'iterations = '
     &              ,niter,'deviation = ',maxdev
               close(20)
               open(unit=21, file = 'towhee_weights',status='unknown')
               write(21,*) '*** Still iterating ... ***'
               write(21,'(a23,2x,i5)') 'Number of iterations = ',niter
               write(21,'(a23,2x,i5)') 'Number of components = ',ncomp
               write(21,'(a18,2x,i5)') 'Number of files = ',ntotfiles
               write(21,*) '*** Individual data for each file ***'
               write(21,*) ' Filename, # Entries, Weight, Temp.,',
     &              ' Chemical potential of each component' 
               do ifile = 1,ntotfiles
                  write(21,*) 'Name of file #',ifile
                  write(21,*) filename(ifile)
                  write(21,*) 'Number of entries in file'
                  write(21,*) nentry(ifile)
                  write(21,*) 'Determined weight for file'
                  write(21,'(g15.6)') weight(ifile)
                  write(21,*) 'Temperature'
                  write(21,*) 1.0d0/beta(ifile)
                  write(21,*) 'Chemical potentials (comp 1 to n)'
                  write(21,'(10f14.4)') 
     &                 (chempotorig(icomp,ifile),icomp = 1,ncomp)
               enddo
               close(21)
            endif
         enddo
         write(6,*) 'New weights :'
         do ifile = 1,ntotfiles
            write(6,*) ifile,weight(ifile)
         enddo
         open(unit=21, file = 'towhee_weights',status='unknown')
         write(21,*) '*** Final weights ***'
         write(21,'(a23,2x,i5)') 'Number of iterations = ',niter
         write(21,'(a23,2x,i5)') 'Number of components = ',ncomp
         write(21,'(a18,2x,i5)') 'Number of files = ',ntotfiles
         write(21,*) '*** Individual data for each file ***'
         do ifile = 1,ntotfiles
            write(21,*) 'Name of file #',ifile
            write(21,*) filename(ifile)
            write(21,*) 'Number of entries in file'
            write(21,*) nentry(ifile)
            write(21,*) 'Determined weight for file'
            write(21,'(g15.6)') weight(ifile)
            write(21,*) 'Temperature'
            write(21,*) 1.0d0/beta(ifile)
            write(21,*) 'Chemical potentials (comp 1 to n)'
            write(21,'(10f14.4)')
     &			(chempotorig(icomp,ifile),icomp = 1,ncomp)
         enddo
         close(21)
      endif

c     --- perform phase coexistence calculations
      if (lphase) then
         write(6,*) 'Starting phase coexistence calculations'
         open(30,file='phase_coex.dat')
         do itemp = 1,ntemp
            betai = 1.0d0/temperature(itemp)
            write(6,*) 'Temperature #',itemp,' out of ',ntemp
            niter = 0
            icount = 0
            if (ncomp .gt. 2) then
               if (chempotincr(itemp) .gt. 0.0d0) then
                  chempoti(ncomp) = chempot2min(itemp)
               elseif (chempotincr(itemp) .lt. 0.0d0) then
                  chempoti(ncomp) = chempot2max(itemp)
               endif
            elseif(ncomp .eq. 2) then
               chempoti(ncomp) = chempot2min(itemp)
            endif
            do while((ncomp.eq.1 .and. niter .eq. 0) .or.
     &           ((chempoti(ncomp) .lt. chempot2max(itemp)+0.00001d0) 
     &           .and.
     &           (chempoti(ncomp).gt.chempot2min(itemp)-0.00001d0)))
               chempotincr1 = chempot(1,1)/1000.0d0
               if (lconstantp) then 
                  chempotincrn = chempot(ncomp,1)/1000.0d0
               endif
               niter2 = 0
               nflip = 0
               oldsign = 1
               pressdiff = 1.0e6
               icount = icount + 1
               do while (pressdiff .gt. MAXPRESSTOL .and. 
     &              niter2 .lt. MAXITER)
                  niter = 0
                  nbelow = 0.0d0
                  do while(abs(0.5d0-nbelow) > 1.0e-6 .and.
     &                 niter .le. MAXITER)
                     if (niter .gt. 0) then
                        chempoti(1) = chempoti(1) + chempotincr1
                        call calculatez(chempoti,nbelow1,betai,
     &                       energy_liq,energy_gas,part_function,
     &                       part_function_liq,part_function_gas,
     &	                     itemp)
                        chempoti(1) = chempoti(1) - 
     &                       (log(nbelow1)-log(1.0d0-nbelow1))
     &                       *chempotincr1/(log(nbelow1)-log(nbelow)
     &                       -log(1.0d0-nbelow1)+log(1.0d0-nbelow))
                     endif
                     call calculatez(chempoti,nbelow,betai,energy_liq,
     &                    energy_gas,part_function,
     &                    part_function_liq,part_function_gas,itemp)
                     niter = niter + 1
                  enddo
                  pressure = (log(part_function_liq)-logz)*CONV2BAR/
     &                 betai/volume
                  if (.not. lconstantp) then
                     pressdiff = 1.0e-6
                  else
                     sign = pressure-pset
                     pressdiff = abs(pressure-pset)
                     if (nflip .gt. 0) then
                        chempotincrn = chempotincrn/2.0d0
                        nflip = 0
                     endif
                     if (sign .gt. 0)then
                        chempoti(ncomp) = chempoti(ncomp) + chempotincrn
                        newsign = -1
                     else
                        chempoti(ncomp) = chempoti(ncomp) - chempotincrn
                        newsign = 1
                     endif
                     if (newsign .ne. oldsign) then
                        oldsign = newsign
                        nflip = nflip + 1
                     endif
                  endif
                  niter2 = niter2 + 1
               enddo
               do icomp = 1,ncomp 
                  nmol_gas(icomp) = nmol_gas(icomp)/nbelow
                  nmol_liq(icomp) = nmol_liq(icomp)/(1.0d0-nbelow)
                  debroglie = 17.458/sqrt(mass(icomp,1)*
     &							temperature(itemp))
                  chempot_deBrog(icomp) = chempoti(icomp)+
     &							3.0*temperature(itemp)*log(debroglie)
               enddo
               energy_gas = energy_gas/(nbelow)
               energy_liq = energy_liq/(1.0d0-nbelow)
               write(30,*) 'Temperature: ',temperature(itemp)
               write(30,*) 'Chemical potential',
     &              ' (for components 1 to n): '
               write(30,'(5(f14.4,2x))')
     &				(chempot_deBrog(icomp),icomp=1,ncomp)
               write(30,*) 'Energy (liquid) = ',energy_liq
               write(30,*) 'Energy (gas) = ',energy_gas
               write(30,*) 'Num. molecules in liquid',
     &              ' (from components 1 to n):'
               write(30,'(5(f14.4,2x))')(nmol_liq(icomp),icomp=1,ncomp)
               write(30,*)'Num. molecules in gas',
     &              ' (from components 1 to n):'
               write(30,'(5(f14.4,2x))')(nmol_gas(icomp),icomp=1,ncomp)
               write(30,*) 'ln(Z(liq)) = ',log(part_function_liq)
               write(30,*) 'ln(Z(gas)) = ',log(part_function_gas)
               write(30,*) '**************************************'
               write(ftemp,*) itemp
               read(ftemp,*) fphase
               do icomp = 1,ncomp
                  write(ftemp2,*) icomp
                  read(ftemp2,*) fphase2
                  filephaseout = 'phasecomp'//
     &                 fphase2(1:len_trim(fphase2))//'temp'
     &                 //fphase(1:len_trim(fphase))//'.dat'
                  open(31,file=filephaseout)
                  write(31,*) '#Temperature : ',temperature(itemp)
                  write(31,*) '#Component : ',icomp
                  write(31,*) '#Chemical potential : ',
     &				chempot_deBrog(icomp)
                  write(31,*) '#Volume (nm^3) :',volume/1000.0d0
                  do ipart = num_part(1,icomp),num_part(2,icomp)
                     write(31,'(i6,2x,g15.4e3)') ipart,dens(icomp,ipart)
                  enddo
                  close(31)
               enddo
               chempoti(2) = chempoti(2) + chempot2incr(itemp)
            enddo
         enddo
         close(30)
      endif

c     --- perform pvt calculations
      if (lpvt) then
         write(6,*) 'Starting pvt calculations'
         inum = 1
         open(40,file='pvt.dat')
         do itemp = 1,ntemp
            write(6,*) 'Temperature #',itemp,' out of ',ntemp
            betai = 1.0d0/temperature(itemp)
            do icomp = 1,ncomp
               chempoti(icomp) = chempot_pvt(itemp,icomp)
            enddo
            midpt(itemp) = num_part(2,1)
            call calculatez(chempoti,nbelow,betai,energy_liq,energy_gas,
     &           part_function,part_function_liq,
     &           part_function_gas,itemp)
            avgn = 0.0d0
            do icomp = 1,ncomp
               avgn = avgn + avgnum(icomp)
            enddo
            write(40,*) '********************'
            write(40,*) 'Temperature #',itemp
            write(40,*) temperature(itemp)
            write(40,*) '<N>'
            write(40,*) avgn
            write(40,*) 'ln(Z)'
            write(40,*) log(part_function)
            write(ftemp,*) itemp
            read(ftemp,*) fname
            
            do icomp = 1,ncomp
               debroglie = 17.458/sqrt(mass(icomp,1)*temperature(itemp))
               chempot_deBrog(icomp) = chempoti(icomp)+
     &							3.0*temperature(itemp)*log(debroglie) 
               write(ftemp3,*) icomp
               read(ftemp3,*) fname2
               filepvt2 = 'pvtcomp'//fname2(1:len_trim(fname2))//'temp'
     &              //fname(1:len_trim(fname))//'.dat'
               open(unit=41,file=filepvt2,status='unknown')
               write(41,*) '#Component #',icomp
               write(41,*) '#Chemical potenital: ',chempot_deBrog(icomp)
               write(41,*) '#Temperature: ',temperature(itemp)
               write(41,*) '#Volume (nm^3): ',volume/1000.0d0
               do ipart = num_part(1,icomp),num_part(2,icomp)
                  write(41,'(i6,2x,g15.4e3)') ipart,dens(icomp,ipart)
               enddo
               close(41)
            enddo      
         enddo
         close(40)
      endif

c     --- perform critical point calculations
      if (lcritical) then
         do itemp = 1,ntemp_crit
            betai = 1.0d0/temp_crit(itemp)
            if (.not. lauto) then
               do icomp = 1,ncomp
                  opt_in(1,icomp) = chempot_crit(icomp,itemp)
                  opt_in(1,icomp+ncomp) = fieldmix(icomp,itemp)
                  debroglie = 17.458/sqrt(mass(icomp,ifile)/betai)
				chempot_crit(icomp,itemp) = 
     &			chempot_crit_orig(icomp,itemp)-3.*log(debroglie)/betai
			 enddo
               opt_in(1,vector_dim) = ncrit(itemp)
            else
               do icomp = 1,(2*ncomp)+2
                  do jcomp = 1,ncomp
                     opt_in(icomp,jcomp) = chempot_crit(jcomp,itemp)
                     opt_in(icomp,jcomp+ncomp) = fieldmix(jcomp,itemp)
                  enddo
                  opt_in(icomp,vector_dim) = ncrit(itemp)
               enddo
               do icomp = 2,(2*ncomp)+2
                  do jcomp = 1,(2*ncomp)+1
                     opt_in(icomp,jcomp) = opt_in(icomp,jcomp) + 
     &                    trial_vector(jcomp)
                  enddo
               enddo
               do idim = 1,vector_dim+1
                  do jdim = 1,vector_dim
                     opt_try(jdim) = opt_in(idim,jdim)
                     call calculatecrit(opt_try,betai,xavg,xmin,xmax,
     &                    variance,part_function,critvalue)
                     opt_y(idim) = critvalue
c                     opt_y(idim) = calculatecrit(opt_try,betai,xavg,
c     &                    xmin,xmax,variance,part_function)
                  enddo
               enddo
               call amoeba(opt_in,opt_y,vector_dim,betai,
     &              xavg,xmin,xmax,variance,part_function)
               do icomp = 1,ncomp
                  chempot_crit(icomp,itemp) = opt_in(1,icomp)
                  fieldmix(icomp,itemp) = opt_in(1,icomp+ncomp)
               enddo
               ncrit(itemp) = opt_in(1,vector_dim)
            endif
            do idim = 1,vector_dim
               opt_try(idim) = opt_in(1,idim)
            enddo
            call calculatecrit(opt_try,betai,xavg,xmin,xmax,
     &                    variance,part_function,critvalue)
c            opt_y(1) = calculatecrit(opt_try,betai,xavg,xmin,xmax,
c     &           variance,part_function)
            opt_y(1) = critvalue
c           don't think need another round of calculations THINK ABOUT THIS    
            open(50,file='orderparameters.dat')
            write(50,'(a6,g14.3)') '<x> = ',xavg 
     	    write(50,'(a12,g14.3)')  'variance = ',variance
            do idim = 0,int(ncrit(itemp))
               write(50,'(g12.3,2x,g12.3)') 
     &              (xmin+(idim*(xmax-xmin)/ncrit(itemp)-xavg))/
     &              sqrt(variance),forder(idim)*(ncrit(itemp)/
     &              (xmax-xmin))*sqrt(variance)
            enddo
            close(50)
            write(6,*) 'Deviation = ',real(opt_y(1)*100)
            open(51,file='criticalparameters.dat')
            write(51,'(f9.5,f12.6,20f12.5)') temp_crit(itemp),
     &           log(part_function),opt_y(1)*100,
     &           (chempot_crit_orig(icomp,itemp),icomp=1,ncomp),
     &           (avgnum(icomp),icomp=1,ncomp),
     &           (fieldmix(icomp,itemp),icomp=1,ncomp)
            close(51)
         enddo
      endif

      end


      subroutine calculatez(chempoti,nbelow,betai,energy_liq,energy_gas,
     &     part_function,part_function_liq,part_function_gas,
     &     itemp)
c     ******************************************************************
c     * this subroutine computes the partition function Z              *
c     * This code has been modified from a program originally written  *
c     * by Jeff Potoff on 04/29/1998.                                  *
c     *                                                                *
c     * originally written 04-19-2004 by C.D. Lorenz                   *
c     * last modified 05-11-2004 by C.D. Lorenz                        *
c     ******************************************************************      
      
      implicit none
c     --- parameters
      integer MAXCOMP,MAXVECTOR,MAXN,MAXFILES,MAXENTRY,MAXTEMP,MAXISING
      integer MAXCRIT,NUMBIN,MAXITER
      double precision MAXTOL,MAXPRESSTOL,CONV2BAR

c     MAXCOMP == the maximum number of components in a system
c     MAXVECTOR == the maximum number of optimizing dimensions used
c                  to determine critical properties (2*MAXCOMP)+1
c     MAXN == the maximum number of molecules of any component
c     MAXFILES == the maximum number of independent histogram files to be
c                 patched together
      parameter(MAXCOMP=10,MAXVECTOR=21,MAXN=3000,MAXFILES=30)
c     MAXENTRY == the maximum number of total entries in all files combined
c     MAXTOL == the tolerance for convergence of the weights
c     MAXPRESSTOL == the tolerance for the convergence of the pressures
      parameter(MAXENTRY=600000,MAXTOL=5.0e-5,MAXPRESSTOL=1.0e-3)
c     MAXTEMP == the maximum number of temperature entries in the phase
c                coexistence data file
c     MAXISING == the maximum number of data points in the ising.dat file
c     MAXCRIT == the maximum number of bins in the discretization of the
c                field mixing parameter
c     NUMBIN == the number of energy bins used in critical point
c               distribution calculations
      parameter(MAXTEMP=100,MAXISING=1000,MAXCRIT=1000,NUMBIN=1001)
c     MAXITER == the maximum number of iterations allowed during simulation
c     CONV2BAR == the conversion to pressure in bar
      parameter(MAXITER=5000,CONV2BAR=138.0d0)

c --- global variables for reweight.f
      integer nmol,nentry,ntotentry,ntotfiles,ncomp,nising,num_part
      dimension nmol(MAXCOMP,MAXENTRY),nentry(MAXFILES)
      dimension num_part(2,MAXCOMP)

      double precision beta,energy,chempot,slope,oldweight,weight
      double precision forder,xising,yising,dens,avgnum,nmol_liq
      double precision nmol_gas,midpt
      dimension oldweight(MAXFILES),beta(MAXFILES),slope(MAXTEMP)
      dimension chempot(MAXCOMP,MAXFILES),energy(MAXENTRY)
      dimension midpt(MAXTEMP),nmol_gas(MAXCOMP),nmol_liq(MAXCOMP)
      dimension xising(MAXISING),yising(MAXISING),avgnum(MAXCOMP)
      dimension dens(MAXCOMP,0:MAXN),forder(0:MAXCRIT),weight(MAXFILES)

      logical lweight,lphase,lcritical,lpvt,lconstantp

      common /increweight/ nentry
      common /increweight1/ nmol
      common /increweight2/ num_part
      common /increweight3/ ntotentry,ntotfiles,ncomp,nising
      common /increweight4/ oldweight,beta,weight
      common /increweight5/ slope,midpt
      common /increweight6/ nmol_gas,nmol_liq,avgnum
      common /increweight7/ chempot
      common /increweight8/ energy
      common /increweight9/ xising
      common /increweight10/ yising
      common /increweight11/ dens
      common /increweight12/ forder
      common /increweight13/ lweight,lphase,lcritical,lpvt,lconstantp

c     --- variables passed to/from subroutine
      integer itemp
      double precision nbelow,betai,energy_liq,energy_gas
      double precision part_function_liq,part_function_gas
      double precision part_function,chempoti
      dimension chempoti(MAXCOMP)      
c     --- local variables
      integer ientry,jfile,icomp,nmoli
      
      double precision ylog,prob,ylogtmp,expmylog,nsplit
      
               
c     --- initialize variables
      part_function = 0.0d0
      part_function_liq = 0.0d0
      part_function_gas = 0.0d0
      nbelow = 0.0d0
      nsplit = 0.0d0
      energy_gas = 0.0d0
      energy_liq = 0.0d0
      do icomp = 1,ncomp
         avgnum(icomp) = 0.0d0
         nmol_gas(icomp) = 0.0d0
         nmol_liq(icomp) = 0.0d0
         do ientry = 0,MAXN
            dens(icomp,ientry) = 0.0d0
         enddo
      enddo
       
      do ientry = 1,ntotentry
         ylog = -1.0e9
         do jfile = 1,ntotfiles
            prob = -1.0d0*(beta(jfile)-betai)*energy(ientry)
            do icomp = 1,ncomp
               prob = prob + (((beta(jfile)*chempot(icomp,jfile))
     &              - (betai*chempoti(icomp)))
     &              * nmol(icomp,ientry))
            enddo
c           --- use logs in order to avoid overflow                   
            ylogtmp = prob + log(nentry(jfile)/oldweight(jfile))
            ylog = max(ylog,ylogtmp) + 
     &           log(1.0d0 + exp(-1.0d0*dabs(ylog-ylogtmp)))
         enddo
         expmylog = exp(-1.0d0*ylog)
         part_function = part_function + expmylog

         nsplit = 0.0d0
         if (lphase) then
            do icomp = 2,ncomp
               nsplit = nsplit + (nmol(icomp,ientry)*slope(itemp))
            enddo
         endif
c         write(6,*) "check1: ",icomp,ientry,itemp,nmol(icomp,ientry),
c     &        midpt(itemp)
         if(nmol(1,ientry) .lt. (midpt(itemp)+nsplit)) then
            part_function_gas = part_function_gas+expmylog
            nbelow = nbelow + expmylog
c            write(6,*) "check2: ",nbelow,icomp,ientry,itemp
            do icomp = 1,ncomp
               nmol_gas(icomp) = nmol_gas(icomp) + nmol(icomp,ientry)
     &              *expmylog
            enddo
            energy_gas = energy_gas + energy(ientry)*expmylog
         else
            part_function_liq = part_function_liq + expmylog
            do icomp = 1,ncomp
               nmol_liq(icomp) = nmol_liq(icomp) + nmol(icomp,ientry)
     &              *expmylog
            enddo
            energy_liq = energy_liq + energy(ientry)*expmylog
         endif
         do icomp = 1,ncomp
            nmoli = nmol(icomp,ientry)
            avgnum(icomp) = avgnum(icomp) + nmoli*expmylog
            dens(icomp,nmoli) = dens(icomp,nmoli) + expmylog
         enddo
      enddo
      do icomp = 1,ncomp
         nmol_gas(icomp) = nmol_gas(icomp)/part_function
         nmol_liq(icomp) = nmol_liq(icomp)/part_function
         avgnum(icomp) = avgnum(icomp)/part_function
         do ientry = 0,MAXN
            dens(icomp,ientry) = dens(icomp,ientry)/part_function
         enddo
      enddo

      energy_gas = energy_gas/part_function
      energy_liq = energy_liq/part_function
      nbelow = nbelow/part_function
         
      return
      end

      subroutine  calculatecrit(opt_try,betai,xavg,xmin,xmax,variance,
     &     part_function,critvalue)
c     ******************************************************************
c     * this function calculates the difference between the ising and  *
c     * the real distribution.                                         *
c     * This code has been modified from a program originally written  *
c     * by Jeff Potoff on 04/29/1998.                                  *
c     *                                                                *
c     * originally written 05-04-2004 by C.D. Lorenz                   *
c     * last modified 05-04-2004 by C.D. Lorenz                        *
c     ******************************************************************
      
      implicit none
c     --- parameters
      integer MAXCOMP,MAXVECTOR,MAXN,MAXFILES,MAXENTRY,MAXTEMP,MAXISING
      integer MAXCRIT,NUMBIN,MAXITER
      double precision MAXTOL,MAXPRESSTOL,CONV2BAR

c     MAXCOMP == the maximum number of components in a system
c     MAXVECTOR == the maximum number of optimizing dimensions used
c                  to determine critical properties (2*MAXCOMP)+1
c     MAXN == the maximum number of molecules of any component
c     MAXFILES == the maximum number of independent histogram files to be
c                 patched together
      parameter(MAXCOMP=10,MAXVECTOR=21,MAXN=3000,MAXFILES=30)
c     MAXENTRY == the maximum number of total entries in all files combined
c     MAXTOL == the tolerance for convergence of the weights
c     MAXPRESSTOL == the tolerance for the convergence of the pressures
      parameter(MAXENTRY=600000,MAXTOL=5.0e-5,MAXPRESSTOL=1.0e-3)
c     MAXTEMP == the maximum number of temperature entries in the phase
c                coexistence data file
c     MAXISING == the maximum number of data points in the ising.dat file
c     MAXCRIT == the maximum number of bins in the discretization of the
c                field mixing parameter
c     NUMBIN == the number of energy bins used in critical point
c               distribution calculations
      parameter(MAXTEMP=100,MAXISING=1000,MAXCRIT=1000,NUMBIN=1001)
c     MAXITER == the maximum number of iterations allowed during simulation
c     CONV2BAR == the conversion to pressure in bar
      parameter(MAXITER=5000,CONV2BAR=138.0d0)

c     --- global variables for reweight.f
      integer nmol,nentry,ntotentry,ntotfiles,ncomp,nising,num_part
      dimension nmol(MAXCOMP,MAXENTRY),nentry(MAXFILES)
      dimension num_part(2,MAXCOMP)

      double precision beta,energy,chempot,slope,oldweight,weight
      double precision forder,xising,yising,dens,avgnum,nmol_liq
      double precision nmol_gas,midpt,critvalue
      dimension oldweight(MAXFILES),beta(MAXFILES),slope(MAXTEMP)
      dimension chempot(MAXCOMP,MAXFILES),energy(MAXENTRY)
      dimension midpt(MAXTEMP),nmol_gas(MAXCOMP),nmol_liq(MAXCOMP)
      dimension xising(MAXISING),yising(MAXISING),avgnum(MAXCOMP)
      dimension dens(MAXCOMP,0:MAXN),forder(0:MAXCRIT),weight(MAXFILES)

      logical lweight,lphase,lcritical,lpvt,lconstantp

      common /increweight/ nentry
      common /increweight1/ nmol
      common /increweight2/ num_part
      common /increweight3/ ntotentry,ntotfiles,ncomp,nising
      common /increweight4/ oldweight,beta,weight
      common /increweight5/ slope,midpt
      common /increweight6/ nmol_gas,nmol_liq,avgnum
      common /increweight7/ chempot
      common /increweight8/ energy
      common /increweight9/ xising
      common /increweight10/ yising
      common /increweight11/ dens
      common /increweight12/ forder
      common /increweight13/ lweight,lphase,lcritical,lpvt,lconstantp

c     --- variables that are passed to/from the subroutine
      double precision opt_try,betai,xavg,xmin,xmax,variance
      double precision part_function
      dimension opt_try(MAXVECTOR)

c     --- local variables	
      integer icomp,idim,ientry,jfile,index,ipart,ipoint,inonzero
      integer nmoli

      double precision chempoti,fieldmixi,ncrit,xcrit,ylog
      double precision prob,ylogtmp,expmylog,xpoint,ypoint,yvalue
      dimension chempoti(MAXCOMP),fieldmixi(MAXCOMP)
      
      

      do icomp = 1,ncomp
         chempoti(icomp) = opt_try(icomp)
         fieldmixi(icomp) = opt_try(icomp+ncomp)
      enddo
      ncrit = opt_try(2*ncomp+1)
      part_function = 0.0d0
      do idim = 0,int(ncrit)
         forder(idim) = 0.0d0
      enddo
      xavg = 0.0d0
         
      xmin = 1.0e6
      xmax = -1.0e6
      do ientry = 1,ntotentry         
         xcrit = nmol(1,ientry) - (fieldmixi(1)*energy(ientry))
         do icomp = 2,ncomp
            xcrit = xcrit - (fieldmixi(icomp)*nmol(icomp,ientry))
         enddo
         xmin = min(xmin,xcrit)
         xmax = max(xmax,xcrit)
      enddo
      do ientry = 1,ntotentry
         ylog = -1.0e9
         do jfile = 1,ntotfiles
            prob = -1.0d0*(beta(jfile)-betai)*energy(ientry)
            do icomp = 1,ncomp
               prob = prob + (((beta(jfile)*chempot(icomp,jfile))
     &              - (betai*chempoti(icomp)))
     &              * nmol(icomp,ientry))
            enddo
c            --- use logs in order to avoid overflow                   
            ylogtmp = prob + log(nentry(jfile)/oldweight(jfile))
            ylog = max(ylog,ylogtmp) + 
     &           log(1.0d0 + exp(-1.0d0*dabs(ylog-ylogtmp)))
         enddo
         expmylog = exp(-1.0d0*ylog)
         part_function = part_function + expmylog
         xcrit = 0.0d0
c        --- calculations for critical point calculations (iflag == 2)               
c        --- calculate the position of the critical point on the x-axis
         xcrit = nmol(1,ientry)-(fieldmixi(1)*energy(ientry))
         do icomp = 2,ncomp
            xcrit = xcrit - (fieldmixi(icomp)*
     &           nmol(icomp,ientry))
         enddo
         index = int(ncrit*(xcrit-xmin)/(xmax-xmin))
         forder(index) = forder(index) + expmylog
         xavg = xavg + (xcrit*expmylog)
         do icomp = 1,ncomp
            nmoli = nmol(icomp,ientry)
            dens(icomp,nmoli) = dens(icomp,nmoli) + expmylog
         enddo
      enddo
      do icomp = 1,ncomp
         do ipart = num_part(1,icomp),num_part(2,icomp)
            avgnum(icomp) = avgnum(icomp) + 
     &           (ipart*(dens(icomp,ipart)/part_function))
            nmol_gas(icomp) = nmol_gas(icomp) + 
     &           (ipart*(dens(icomp,ipart)/part_function))
         enddo
      enddo

      xavg = xavg/part_function
      variance = 0.0d0
      do idim = 0,MAXISING
         forder(idim) = forder(idim)/part_function
         variance = variance + forder(idim)*(xmin+
     &        (idim*(xmax-xmin)/ncrit-xavg))**2
      enddo   

      critvalue = 0.0d0
      ipoint = 1
      inonzero = 0
      do idim = 0,int(ncrit)+1
         xpoint = (xmin+(idim*(xmax-xmin)/(ncrit-xavg)))/sqrt(variance)
         ypoint = (forder(idim)*sqrt(variance))/(xmax-xmin)*ncrit
         do while (ipoint .le. nising-2 .and. 
     &        xising(ipoint+1) .lt. xpoint)
            ipoint = ipoint + 1
         enddo
         yvalue = yising(ipoint) + ((yising(ipoint+1)-yising(ipoint))/
     &        (xising(ipoint+1)-xising(ipoint)))*(xpoint-xising(ipoint))
         critvalue = critvalue + dabs(ypoint-yvalue)
         if (ypoint .gt. 1.0e-6 .and. yvalue .gt. 1.0e-6) 
     &        inonzero = inonzero + 1
      enddo
      critvalue = critvalue/inonzero
      
      return
      end

      subroutine amoeba(opt_in,opt_y,ndim,betai,
     &     xavg,xmin,xmax,variance,part_function)
c     ******************************************************************
c     * this subroutine conducts an amoeba search                      *
c     * This code has been taken from the numerical recipes book.      *
c     *                                                                *
c     * originally written 05-04-2004 by C.D. Lorenz                   *
c     * last modified 05-05-2004 by C.D. Lorenz                        *
c     ******************************************************************
      
      implicit none
c     --- parameters
      integer MAXCOMP,MAXVECTOR,MAXN,MAXFILES,MAXENTRY,MAXTEMP,MAXISING
      integer MAXCRIT,NUMBIN,MAXITER
      double precision MAXTOL,MAXPRESSTOL,CONV2BAR

c     MAXCOMP == the maximum number of components in a system
c     MAXVECTOR == the maximum number of optimizing dimensions used
c                  to determine critical properties (2*MAXCOMP)+1
c     MAXN == the maximum number of molecules of any component
c     MAXFILES == the maximum number of independent histogram files to be
c                 patched together
      parameter(MAXCOMP=10,MAXVECTOR=21,MAXN=3000,MAXFILES=30)
c     MAXENTRY == the maximum number of total entries in all files combined
c     MAXTOL == the tolerance for convergence of the weights
c     MAXPRESSTOL == the tolerance for the convergence of the pressures
      parameter(MAXENTRY=600000,MAXTOL=5.0e-5,MAXPRESSTOL=1.0e-3)
c     MAXTEMP == the maximum number of temperature entries in the phase
c                coexistence data file
c     MAXISING == the maximum number of data points in the ising.dat file
c     MAXCRIT == the maximum number of bins in the discretization of the
c                field mixing parameter
c     NUMBIN == the number of energy bins used in critical point
c               distribution calculations
      parameter(MAXTEMP=100,MAXISING=1000,MAXCRIT=1000,NUMBIN=1001)
c     MAXITER == the maximum number of iterations allowed during simulation
c     CONV2BAR == the conversion to pressure in bar
      parameter(MAXITER=5000,CONV2BAR=138.0d0)

c     --- global variables for reweight.f
      integer nmol,nentry,ntotentry,ntotfiles,ncomp,nising,num_part
      dimension nmol(MAXCOMP,MAXENTRY),nentry(MAXFILES)
      dimension num_part(2,MAXCOMP)

      double precision beta,energy,chempot,slope,oldweight,weight
      double precision forder,xising,yising,dens,avgnum,nmol_liq
      double precision nmol_gas,midpt
      dimension oldweight(MAXFILES),beta(MAXFILES),slope(MAXTEMP)
      dimension chempot(MAXCOMP,MAXFILES),energy(MAXENTRY)
      dimension midpt(MAXTEMP),nmol_gas(MAXCOMP),nmol_liq(MAXCOMP)
      dimension xising(MAXISING),yising(MAXISING),avgnum(MAXCOMP)
      dimension dens(MAXCOMP,0:MAXN),forder(0:MAXCRIT),weight(MAXFILES)

      logical lweight,lphase,lcritical,lpvt,lconstantp

      common /increweight/ nentry
      common /increweight1/ nmol
      common /increweight2/ num_part
      common /increweight3/ ntotentry,ntotfiles,ncomp,nising
      common /increweight4/ oldweight,beta,weight
      common /increweight5/ slope,midpt
      common /increweight6/ nmol_gas,nmol_liq,avgnum
      common /increweight7/ chempot
      common /increweight8/ energy
      common /increweight9/ xising
      common /increweight10/ yising
      common /increweight11/ dens
      common /increweight12/ forder
      common /increweight13/ lweight,lphase,lcritical,lpvt,lconstantp

c     --- variables passed to/from this subroutine
      integer ndim

      double precision opt_in,opt_y,betai,xavg,xmin,xmax
      double precision variance,part_function
      dimension opt_in(MAXVECTOR+1,MAXVECTOR),opt_y(MAXVECTOR+1)

c     --- local variables
      integer ndim1,niter,idim,ilo,ihi,inhi,jdim

      double precision sum,sum_opt,rtol,swap,ytry,ysave
      double precision critvalue
      dimension sum_opt(MAXVECTOR)

      ndim1 = ndim+1
      niter = 0
 1    do 12 idim = 1,ndim
         sum = 0.0d0
         do 11 jdim = 1,ndim1
            sum = sum + opt_in(jdim,idim)
 11      continue
         sum_opt(idim) = sum
 12   continue
 2    ilo = 1
      if (opt_y(1) .gt. opt_y(2))then
         ihi = 1
         inhi = 2
      else
         ihi = 2
         inhi = 1
      endif
      do 13 idim = 1,ndim1
         if (opt_y(idim) .le. opt_y(ilo)) ilo=idim
         if (opt_y(idim) .gt. opt_y(ihi)) then
            inhi = ihi
            ihi = idim
         elseif (opt_y(idim) .gt. opt_y(inhi)) then
            if (idim .ne. ihi) inhi = idim
         endif
 13   continue
      rtol = 2.0d0*abs(opt_y(ihi)-opt_y(ilo))/
     &     (abs(opt_y(ihi))+abs(opt_y(ilo)))
      if (rtol .lt. MAXTOL) then
         swap = opt_y(1) 
         opt_y(1) = opt_y(ilo)
         opt_y(ilo) = swap
         do 14 idim = 1,ndim
            swap = opt_in(1,idim)
            opt_in(1,idim) = opt_in(ilo,idim)
            opt_in(ilo,idim) = swap
 14      continue
         return
      endif
      if (niter .ge. MAXITER) then
         write(*,*) 'MAXITER exceeded in amoeba search'
         stop
      endif
      niter = niter + 2
      call amotry(opt_in,opt_y,sum_opt,ndim,ihi,-1.0d0,
     &     betai,xavg,xmin,xmax,variance,part_function,ytry)
c      ytry = amotry(opt_in,opt_y,sum_opt,ndim,ihi,-1.0,
c     &     betai,xavg,xmin,xmax,variance,part_function)
      if (ytry .le.opt_y(ilo)) then
         call amotry(opt_in,opt_y,sum_opt,ndim,ihi,2.0d0,
     &     betai,xavg,xmin,xmax,variance,part_function,ytry)
c         ytry = amotry(opt_in,opt_y,sum_opt,ndim,ihi,-1.0,
c     &     betai,xavg,xmin,xmax,variance,part_function)
      elseif (ytry .ge. opt_y(inhi)) then
         ysave =  opt_y(ihi)
         call amotry(opt_in,opt_y,sum_opt,ndim,ihi,0.5d0,
     &     betai,xavg,xmin,xmax,variance,part_function,ytry)
c         ytry = amotry(opt_in,opt_y,sum_opt,ndim,ihi,-1.0,
c     &     betai,xavg,xmin,xmax,variance,part_function)
         if (ytry .ge. ysave) then
            do 16 idim = 1,ndim1
               if (idim .ne. ilo) then
                  do 15 jdim = 1,ndim
                     sum_opt(jdim) = 0.5d0*(opt_in(idim,jdim)+
     &                    opt_in(ilo,jdim))
                     opt_in(idim,jdim) = sum_opt(jdim)
 15               continue
                  call calculatecrit(sum_opt,betai,xavg,xmin,xmax,
     &                 variance,part_function,critvalue)
c                  opt_y(idim) = calculatecrit(sum_opt,betai,xavg,xmin,
c     &                 xmax,variance,part_function)
                  opt_y(idim) = critvalue
               endif
 16         continue
            niter = niter + ndim
            goto 1
         endif
      else
         niter = niter + 1
      endif
      goto 2
      end
      

      subroutine amotry(opt_in,opt_y,sum_opt,ndim,ihi,factor,betai,xavg,
     &     xmin,xmax,variance,part_function,atry)
c     ******************************************************************
c     * this function returns a new try for the amoeba search          *
c     * This code has been taken from the numerical recipes book.      *
c     *                                                                *
c     * originally written 05-04-2004 by C.D. Lorenz                   *
c     * last modified 05-05-2004 by C.D. Lorenz                        *
c     ******************************************************************

      implicit none
c     --- parameters
      integer MAXCOMP,MAXVECTOR,MAXN,MAXFILES,MAXENTRY,MAXTEMP,MAXISING
      integer MAXCRIT,NUMBIN,MAXITER
      double precision MAXTOL,MAXPRESSTOL,CONV2BAR

c     MAXCOMP == the maximum number of components in a system
c     MAXVECTOR == the maximum number of optimizing dimensions used
c                  to determine critical properties (2*MAXCOMP)+1
c     MAXN == the maximum number of molecules of any component
c     MAXFILES == the maximum number of independent histogram files to be
c                 patched together
      parameter(MAXCOMP=10,MAXVECTOR=21,MAXN=3000,MAXFILES=30)
c     MAXENTRY == the maximum number of total entries in all files combined
c     MAXTOL == the tolerance for convergence of the weights
c     MAXPRESSTOL == the tolerance for the convergence of the pressures
      parameter(MAXENTRY=600000,MAXTOL=5.0e-5,MAXPRESSTOL=1.0e-3)
c     MAXTEMP == the maximum number of temperature entries in the phase
c                coexistence data file
c     MAXISING == the maximum number of data points in the ising.dat file
c     MAXCRIT == the maximum number of bins in the discretization of the
c                field mixing parameter
c     NUMBIN == the number of energy bins used in critical point
c               distribution calculations
      parameter(MAXTEMP=100,MAXISING=1000,MAXCRIT=1000,NUMBIN=1001)
c     MAXITER == the maximum number of iterations allowed during simulation
c     CONV2BAR == the conversion to pressure in bar
      parameter(MAXITER=5000,CONV2BAR=138.0d0)

c     --- global variables for reweight.f
      integer nmol,nentry,ntotentry,ntotfiles,ncomp,nising,num_part
      dimension nmol(MAXCOMP,MAXENTRY),nentry(MAXFILES)
      dimension num_part(2,MAXCOMP)

      double precision beta,energy,chempot,slope,oldweight,weight
      double precision forder,xising,yising,dens,avgnum,nmol_liq
      double precision nmol_gas,midpt,atry
      dimension oldweight(MAXFILES),beta(MAXFILES),slope(MAXTEMP)
      dimension chempot(MAXCOMP,MAXFILES),energy(MAXENTRY)
      dimension midpt(MAXTEMP),nmol_gas(MAXCOMP),nmol_liq(MAXCOMP)
      dimension xising(MAXISING),yising(MAXISING),avgnum(MAXCOMP)
      dimension dens(MAXCOMP,0:MAXN),forder(0:MAXCRIT),weight(MAXFILES)

      logical lweight,lphase,lcritical,lpvt,lconstantp

      common /increweight/ nentry
      common /increweight1/ nmol
      common /increweight2/ num_part
      common /increweight3/ ntotentry,ntotfiles,ncomp,nising
      common /increweight4/ oldweight,beta,weight
      common /increweight5/ slope,midpt
      common /increweight6/ nmol_gas,nmol_liq,avgnum
      common /increweight7/ chempot
      common /increweight8/ energy
      common /increweight9/ xising
      common /increweight10/ yising
      common /increweight11/ dens
      common /increweight12/ forder
      common /increweight13/ lweight,lphase,lcritical,lpvt,lconstantp

c     --- variables passed to/from this subroutine
      integer ndim,ihi

      double precision opt_in,opt_y,sum_opt,factor,opt_try
      double precision betai,xavg,xmin,xmax,variance,part_function
      dimension opt_in(MAXVECTOR+1,MAXVECTOR),opt_y(MAXVECTOR+1)
      dimension sum_opt(MAXVECTOR),opt_try(MAXVECTOR)

c     --- local variables
      integer idim
      
      double precision factor1,factor2,ytry,critvalue

      factor1 = (1.0d0-factor)/ndim
      factor2 = factor1-factor
      do 11 idim = 1,ndim
         opt_try(idim) = (sum_opt(idim)*factor1) - 
     &        (opt_in(ihi,idim)*factor2)
 11   continue
      call calculatecrit(opt_try,betai,xavg,xmin,xmax,
     &     variance,part_function,critvalue)
c      ytry = calculatecrit(opt_try,betai,xavg,xmin,xmax,variance,
c     &     part_function)
      ytry = critvalue
      if (ytry .lt. opt_y(ihi)) then
         opt_y(ihi)=ytry
         do 12 idim = 1,ndim
            sum_opt(idim) = sum_opt(idim) - opt_in(ihi,idim) + 
     &           opt_try(idim)
            opt_in(ihi,idim) = opt_try(idim)
 12      continue
      endif
      atry = ytry
      return
      end
