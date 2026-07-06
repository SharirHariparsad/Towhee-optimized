#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2020 Marcus G. Martin                       *
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
      subroutine twh_writetowhee(lfailure,lprotein,imoltyflag)
c     ******************************************************************
c     * wrapper that determines a few constants needed to dimension    *
c     * arrays and then call the main writetowhee                      *
c     *                                                                *
c     * originally written 01-16-2018 by M.G. Martin                   *
c     * last modified 04-03-2018 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,lprotein
      integer imoltyflag
c     --- local variables
      integer numboxes,nvmmax,nvmsim,nmolty

c     --- retrieve some constants
      call twh_numboxes(GLB_GET,numboxes)
      call twh_nvmmax(GLB_GET,nvmmax)
      call twh_nvmsim(GLB_GET,nvmsim)
      call twh_nmolty(GLB_GET,nmolty)

      call twh_writetowhee_main(lfailure,lprotein,imoltyflag,numboxes
     &     ,nvmmax,nvmsim,nmolty)

      return
      end

      subroutine twh_writetowhee_main(lfailure,lprotein,imoltyflag
     &     ,numboxes,nvmmax,nvmsim,nmolty)
c     ******************************************************************
c     * outputs information suitable for inpstyle 0 to towhee_altinp   *
c     *                                                                *
c     * rewritten from assemble 03-19-2002 by M.G. Martin              *
c     * split from the original writetowhee 01-16-2018 by M.G. Martin  *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_NATIVE_LAMBDA
#define FUNCTION_GET_SCALING_STYLE
#define FUNCTION_GET_SCALING_STYLE_STRING
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lprotein
      integer imoltyflag,numboxes,nvmmax,nvmsim,nmolty
c     --- local variables
      character*1 hist_suffix
      character*(ELENAMLEN) elename
      character*3 aminoshort
      character*4 pdbname
      character*8 stepstyle
      character*10 tenchar
      character*11 fchar
      character*20 routput,cvalue,coulombstyle,fieldstring,initboxtype
      character*20 initstyle,initlattice,interpolatestyle,controlstyle
      character*20 isarray,pressure_thermo_style
      dimension isarray(nmolty)
      character*30 ensemble_str
      character*50 eform
      character*(CLAPOTLEN) localcp
      character*(MAXDIRLENGTH) ffname
      character*(MIXRULLEN) classical_mixrule
      character*(FFNAMELEN) nbname
c     --- logical scalars
      logical lfound,ldoinitstyle,loutdft,louthist,lshift
      logical lwritemixrule,lwritelshift,lwriteltailc,lwritermin
      logical lwritercut,lwritercutin,lwriteinterpolatestyle,lwriterpd
      logical lwritescaling,lvalue,loutdlpoly,loutlammps,ltailc
c     --- character string arrays
      character*10 pair_list
      dimension pair_list(CROSSTYPEMAX,2)
c     --- small integer scalars
      integer*4 random_seed,luxlevel
c     --- integer scalars
      integer iatom,ivib,ibend,itor,iaa,ibox,jmolty,idim,jdim,iorder
      integer kmolty,ifile,ifld,itype,ntii,ientry,ipair,rwfile
      integer nboxpair,nunit,kmax,nfield,ivalue,ffnumber,nmolpairs
      integer aminonum,rng_code,ijaa0,ijaa1,ijaa2
      integer npair,inix,iniy,iniz,ensemble,itbang,iobang
      integer irange,pt_number,ivolmove,itbentemp,ijben2,ijben3
      integer ijimp2,ijimp3,ijimp4,ijtor2,ijtor3,ijtor4,ijvib,inaa
      integer inbend,inimproper,intorsion,invibration,itaa,itimproper
      integer ittor,itvib,ntype
c     --- integer arrays
      integer idummy
      dimension idummy(nmolty)
c     --- double precision scalars
      double precision dvalue,qqfact,temperature,ewald_prec,kalp,rcut
      double precision fieldtwo,fieldthree,fieldfour,fnbparam,pmcb,rmin
      double precision rcutin,qqatom,tavol
c     --- double precision arrays
      double precision dpboxarray
      dimension dpboxarray(numboxes+nvmsim)
      double precision dmarray
      dimension dmarray(nmolty)
      double precision dpairarray
      dimension dpairarray(1+(nmolty*(nmolty-1))/2)
      double precision dboxarray
      dimension dboxarray(numboxes)
      double precision dimtemp
      dimension dimtemp(3)
      double precision fatwo
      dimension fatwo(2)

c      write(6,*) 'starting WRITETOWHEE ',imoltyflag
c     --- retrieve constants
      call twh_ensemble(GLB_GET,ensemble)
      call twh_nmolpairs(GLB_GET,nmolpairs)

      if ( imoltyflag .eq. 0 ) then
c        --- output all of the towhee_input file except for the 
c        --- molecule information
         write(15,'(a11)') 'inputformat'
         write(15,'(a1,a6,a1)') char(39),'Towhee',char(39)
         call twh_random_code(GLB_GET,rng_code)
         write(15,'(a23)') 'random_number_generator'
         if ( rng_code .eq. RNG_RANLUX ) then
            write(15,'(a1,a6,a1)') char(39),'RANLUX',char(39)
            call twh_random_luxlevel(GLB_GET,luxlevel)
            write(15,'(a15)') 'random_luxlevel'
            write(15,*) luxlevel
         elseif ( rng_code .eq. RNG_DX_1597_2_7 ) then
            write(15,'(a1,a11,a1)') char(39),'DX-1597-2-7',char(39)
         elseif ( rng_code .eq. RNG_KISS99 ) then
            write(15,'(a1,a6,a1)') char(39),'KISS99',char(39)
         elseif ( rng_code .eq. RNG_MRG32K3A ) then
            write(15,'(a1,a8,a1)') char(39),'MRG32k3a',char(39)
         elseif ( rng_code .eq. RNG_DEBUG ) then
            write(15,'(a1,a5,a1)') char(39),'debug',char(39)
         else
            write(6,*) 'WRITETOWHEE: unknown random_code'
            write(6,*) rng_code
            lfailure = .true.
            return
         endif
         call twh_random_seed(GLB_GET,random_seed)
         write(15,'(a11)') 'random_seed'
         write(15,*) random_seed
         call twh_random_allow_restart(GLB_GET,lvalue)
         write(15,'(a20)') 'random_allow_restart'
         write(15,'(L1)') lvalue

         if (ensemble .eq. ENS_NPT) ensemble_str='npt'
         if (ensemble .eq. ENS_NVT) ensemble_str='nvt'
         if (ensemble .eq. ENS_muVT) ensemble_str='uvt'
         write(15,'(a8)') 'ensemble'
         write(15,'(a1,a30,a1)') char(39),ensemble_str,char(39)

         write(15,'(a11)') 'temperature'
         call twh_temperature(GLB_GET,temperature)
         write(15,*) temperature
         if ( ensemble .eq. ENS_NPT ) then
            write(15,'(a8)') 'pressure'
            call twh_express(GLB_GET,dvalue)
            dvalue = dvalue/CON_KPATOSIM
            write(15,*) dvalue
         endif
         write(15,'(a6)') 'nmolty'
         write(15,*) nmolty
         write(15,'(a9)') 'nmolectyp'
         do jmolty = 1,nmolty
            call twh_nmolectyp(GLB_GET,jmolty,idummy(jmolty))
         enddo
         write(15,*) (idummy(jmolty),jmolty=1,nmolty)
         if ( ensemble .eq. ENS_muVT ) then
            write(15,'(a7)') 'chempot'
            do jmolty = 1,nmolty
               call twh_chempot(GLB_GET,jmolty,dmarray(jmolty))
            enddo
            write(15,*) (dmarray(jmolty),jmolty=1,nmolty)
         endif
         write(15,'(a8)') 'numboxes'
         write(15,*) numboxes
         write(15,'(a9)') 'stepstyle'
         call twh_stepstyle(GLB_GET,stepstyle)
         write(15,'(a1,a8,a1)') char(39),stepstyle,char(39)
         if ( stepstyle .eq. 'cycles' .or. stepstyle .eq. 'moves' ) then
            write(15,'(a5)') 'nstep'
            call twh_nstep(GLB_GET,ivalue)
            write(15,*) ivalue
         else
            write(6,*) 'WRITETOWHEE: unknown stepstyle'
         endif
         write(15,'(a12)') 'controlstyle'
         call twh_controlstyle(GLB_GET,controlstyle)
         write(15,'(a1,a20,a1)') char(39),controlstyle,char(39)
         write(15,'(a9)') 'printfreq'
         call twh_printfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a9)') 'blocksize'
         call twh_blocksize(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a9)') 'moviefreq'
         call twh_moviefreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a10)') 'backupfreq'
         call twh_backupfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a11)') 'restartfreq'
         call twh_restartfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a9)') 'runoutput'
         call twh_runoutput(GLB_GET,routput)
         write(15,'(a1,a20,a1)') char(39),routput,char(39)
         write(15,'(a15)') 'pdb_output_freq'
         call twh_pdb_output_freq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a7)') 'loutdft'
         call twh_loutdft(GLB_GET,loutdft)
         write(15,*) loutdft
         write(15,'(a10)') 'loutlammps'
         call twh_loutlammps(GLB_GET,loutlammps)
         write(15,*) loutlammps
         write(15,'(a10)') 'loutdlpoly'
         call twh_loutdlpoly(GLB_GET,loutdlpoly)
         write(15,*) loutdlpoly
         if ( ensemble .eq. ENS_muVT ) then
            write(15,'(a8)') 'louthist'
            call twh_louthist(GLB_GET,louthist)
            write(15,*) louthist
            if ( louthist ) then
               write(15,'(a10)') 'hist_label'
               call twh_hist_label(GLB_GET,ivalue)
               write(15,*) ivalue
               write(15,'(a11)') 'hist_suffix'
               call twh_hist_suffix(GLB_GET,hist_suffix)
               write(15,'(a1,a1,a1)') char(39),hist_suffix,char(39)
               write(15,'(a11)') 'hist_nequil'
               call twh_hist_nequil(GLB_GET,ivalue)
               write(15,*) ivalue
               write(15,'(a12)') 'histcalcfreq'
               call twh_histcalcfreq(GLB_GET,ivalue)
               write(15,*) ivalue
               write(15,'(a12)') 'histdumpfreq'
               call twh_histdumpfreq(GLB_GET,ivalue)
               write(15,*) ivalue
            endif
         endif
         write(15,'(a20)') 'pressure_virial_freq'
         call twh_pressure_virial_freq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a20)') 'pressure_thermo_freq'
         call twh_pressure_thermo_freq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a21)') 'pressure_thermo_style'
         call twh_pressure_thermo_style(GLB_GET,pressure_thermo_style)
         write(15,'(a1,a20,a1)') char(39),pressure_thermo_style,char(39)
         call twh_pressure_thermo_number(GLB_GET,pt_number)
         if ( pt_number .eq. PTHERMO_RANGE ) then
            write(15,'(a21)') 'pressure_thermo_range'
            do irange = 1,2
               call twh_pressure_thermo_range(GLB_GET,irange
     &              ,fatwo(irange))
            enddo
            write(15,*) (fatwo(irange),irange=1,2)
         endif
         write(15,'(a13)') 'trmaxdispfreq'
         call twh_trmaxdispfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a14)') 'volmaxdispfreq'
         call twh_volmaxdispfreq(GLB_GET,ivalue)
         write(15,*) ivalue
         write(15,'(a14)') 'chempotperstep'
         do jmolty = 1,nmolty
            call twh_chempotperstep(GLB_GET,jmolty,idummy(jmolty))
         enddo
         write(15,*) (idummy(jmolty),jmolty=1,nmolty)
         write(15,'(a15)') 'loutchempotdata'
         call twh_loutchempotdata(GLB_GET,lvalue)
         write(15,*) lvalue
         write(15,'(a5)') 'linit'
         call twh_linit(GLB_GET,lvalue)
         write(15,*) lvalue
         
         write(15,'(a11)') 'initboxtype'
         call twh_initboxtype(GLB_GET,initboxtype)
         write(15,'(a1,a20,a1)') char(39),initboxtype,char(39)
         if ( initboxtype .eq. 'dimensions' ) then
            ldoinitstyle = .true.
         elseif ( initboxtype .eq. 'number density' ) then
            ldoinitstyle = .true.
         elseif ( initboxtype .eq. 'unit cell' ) then
            ldoinitstyle = .false.
         else
            write(6,*) 'WRITETOWHEE: unknown initboxtype'
            write(6,*) initboxtype
            stop
            return
         endif

         if ( ldoinitstyle ) then
            write(15,'(a9)') 'initstyle'
            do ibox = 1,numboxes
               do jmolty = 1,nmolty
                  call twh_initstyle(GLB_GET,ibox,jmolty,initstyle)
                  isarray(jmolty) = initstyle
               enddo
               write(15,*) (char(39),isarray(jmolty),char(39)
     &              ,' ',jmolty=1,nmolty)
            enddo
c           --- helix_* variables are only output if the initstyle
c           --- for that molecule type is 'COM helix' or 'helix'
            do jmolty = 1,nmolty
               lfound = .false.
               do ibox = 1,numboxes
                  call twh_initstyle(GLB_GET,ibox,jmolty,initstyle)
                  if ( initstyle .eq. 'COM helix' .or.
     &                 initstyle .eq. 'helix' ) then
                     lfound = .true.
                  endif
               enddo
               if ( lfound ) then
c                 --- output the helix_* varibles
                  write(15,'(a12)') 'helix_moltyp'
                  write(15,*) jmolty
                  write(15,'(a12)') 'helix_radius'
                  call twh_helix_radius(GLB_GET,jmolty,dvalue)
                  write(15,*) dvalue
                  write(15,'(a11)') 'helix_angle'
                  call twh_helix_angle(GLB_GET,jmolty,dvalue)
                  write(15,*) dvalue
                  write(15,'(a13)') 'helix_keytype'
                  call twh_helix_keytype(GLB_GET,jmolty,tenchar)
                  write(15,'(a1,a10,a1)') char(39),tenchar,char(39)
                  write(15,'(a13)') 'helix_keyname'
                  call twh_helix_keyname(GLB_GET,jmolty,tenchar)
                  write(15,'(a1,a10,a1)') char(39),tenchar,char(39)
                  write(15,'(a12)') 'helix_conlen'
                  call twh_helix_conlen(GLB_GET,jmolty,dvalue)
                  write(15,*) dvalue
                  write(15,'(a11)') 'helix_phase'
                  call twh_helix_phase(GLB_GET,jmolty,dvalue)
                  write(15,*) dvalue
               endif
            enddo
            write(15,'(a11)') 'initlattice'
            do ibox = 1,numboxes
               do jmolty = 1,nmolty
                  call twh_initlattice(GLB_GET,ibox,jmolty,initlattice)
                  isarray(jmolty) = initlattice
               enddo
               write(15,22) (char(39),isarray(jmolty),char(39)
     &              ,jmolty=1,nmolty)
            enddo
 22         format(NTMAX(a1,a20,a1,1x))
            write(15,'(a7)') 'initmol'
            do ibox = 1,numboxes
               do jmolty = 1,nmolty
                  call twh_initmol(GLB_GET,ibox,jmolty,idummy(jmolty))
               enddo
               write(15,*) (idummy(jmolty),jmolty=1,nmolty)
            enddo
         endif
         write(15,'(a14)') 'inix iniy iniz'
         do ibox = 1,numboxes
            call twh_inixyz(GLB_GET,ibox,inix,iniy,iniz)
            write(15,*) inix,iniy,iniz
         enddo

         if ( initboxtype .eq. 'dimensions' ) then
            ldoinitstyle = .true.
            write(15,'(a7)') 'hmatrix'
            do ibox = 1,numboxes
               do idim = 1,3
                  do jdim = 1,3
                     call twh_hmatrix(GLB_GET,BXD_REAL,ibox,idim,jdim
     &                    ,dimtemp(jdim))
                  enddo
                  write(15,*) (dimtemp(jdim),jdim=1,3)
               enddo
            enddo
         elseif ( initboxtype .eq. 'number density' ) then
            ldoinitstyle = .true.
            write(15,'(a18)') 'box_number_density'
            do ibox = 1,numboxes
               call twh_box_number_density(GLB_GET,ibox,dvalue)
               write(15,*) dvalue
            enddo
         endif

         call twh_nboxpair(GLB_GET,nboxpair)

         if ( ( ensemble .eq. ENS_NVT .and. numboxes .gt. 1 ) 
     &        .or. ensemble .eq. ENS_NPT ) then
               
            write(15,'(a5)') 'pmvol'
            call twh_pmvol(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a6)') 'pmvlpr'
            if ( nvmmax .eq. 0 ) then
               write(15,100) 1.0d0
            else
               do ivolmove = 1,nvmsim
                  call twh_pmvlpr(GLB_GET,ivolmove,dvalue)
                  dpboxarray(ivolmove) = dvalue
               enddo
               write(15,100) (dpboxarray(ivolmove),ivolmove=1,nvmsim)
            endif
            write(15,'(10x,a5)') 'rmvol'
            call twh_rmvol(GLB_GET,1,dvalue)
            write(15,100) dvalue
            write(15,'(10x,a5)') 'tavol'
            call twh_tavol(GLB_GET,tavol)
            write(15,100) tavol
            write(15,'(a6)') 'pmcell'
            call twh_pmcell(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a8)') 'pmcellpr'
            if ( nvmsim .eq. 0 ) then
               write(15,100) 1.0d0
            else
               do ivolmove=1,nvmsim
                  call twh_pmcellpr(GLB_GET,ivolmove,dvalue)
                  dpboxarray(ivolmove) = dvalue
               enddo
               write(15,100) (dpboxarray(ivolmove),ivolmove=1,nvmsim)
            endif
            write(15,'(10x,a8)') 'pmcellpt'
            if ( nboxpair .eq. 0 ) then
               write(15,100) 1.0d0
            else
               do ibox = 1,nboxpair
                  call twh_pmcellpt(GLB_GET,ibox,dvalue)
                  write(15,100) dvalue
               enddo
            endif
            write(15,'(10x,a6)') 'rmcell'
            call twh_rmcell(GLB_GET,1,1,1,dvalue)
            write(15,100) dvalue
            write(15,'(10x,a6)') 'tacell'
            call twh_tacell(GLB_GET,dvalue)
            write(15,100) dvalue
         endif
         if ( numboxes .gt. 1 ) then
c           --- 2 box rotational-bias swap move
            write(15,'(a12)') 'pm2boxrbswap'
            call twh_pm2boxrbswap(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a9)') 'pm2rbswmt'
            do jmolty = 1,nmolty
               call twh_pm2rbswmt(GLB_GET,jmolty,dvalue)
               dmarray(jmolty) = dvalue
            enddo
            write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
            write(15,'(10x,a9)') 'pm2rbswpr'
            do ibox = 1,nboxpair
               call twh_pm2rbswpr(GLB_GET,ibox,dvalue)
               dpboxarray(ibox) = dvalue
            enddo
            write(15,100) (dpboxarray(ibox),ibox=1,nboxpair)
         endif
         if ( numboxes .gt. 1 ) then
c           --- 2 box configurational-bias swap move
            write(15,'(a12)') 'pm2boxcbswap'
            call twh_pm2boxcbswap(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a9)') 'pm2cbswmt'
            do jmolty = 1,nmolty
               call twh_pm2cbswmt(GLB_GET,jmolty,dvalue)
               dmarray(jmolty) = dvalue
            enddo
            write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
            write(15,'(10x,a9)') 'pm2cbswpr'
            do ibox = 1,nboxpair
               call twh_pm2cbswpr(GLB_GET,ibox,dvalue)
               dpboxarray(ibox) = dvalue
            enddo
            write(15,100) (dpboxarray(ibox),ibox=1,nboxpair)
         endif
         if ( numboxes .gt. 1 ) then
c           ---  2 box center-of-mass switch move
            write(15,'(a15)') 'pm2boxcomswitch'
            call twh_pm2boxcomswitch(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a15)') 'pm2comswboxpair'
            do ibox = 1,nboxpair
               call twh_pm2comswboxpair(GLB_GET,ibox,dvalue)
               dpboxarray(ibox) = dvalue
            enddo
            write(15,100) (dpboxarray(ibox),ibox=1,nboxpair)
            write(15,'(10x,a15)') 'pm2comswmolpair'
            do ipair = 1,nmolpairs
               call twh_pm2comswmolpair(GLB_GET,ipair,dvalue)
               dpairarray(ipair) = dvalue
            enddo
            write(15,100) (dpairarray(ipair),ipair=1,nmolpairs)
         endif
         if ( ensemble .eq. ENS_muVT ) then
            write(15,'(a11)') 'pmuvtcbswap'
            call twh_pmuvtcbswap(GLB_GET,dvalue)
            write(15,*) dvalue
            write(15,'(a9)') 'pmuvtcbmt'
            do jmolty = 1,nmolty
               call twh_pmuvtcbmt(GLB_GET,jmolty,dmarray(jmolty))
            enddo
            write(15,*) (dmarray(jmolty),jmolty=1,nmolty)
         endif

         write(15,'(a12)') 'pm1boxcbswap'
         call twh_pm1boxcbswap(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a9)') 'pm1cbswmt'
         do jmolty = 1,nmolty
            call twh_pm1cbswmt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)

         if ( nmolpairs .ge. 1 ) then
            write(15,'(a15)') 'pm1boxcomswitch'
            call twh_pm1boxcomswitch(GLB_GET,dvalue)
            write(15,44) dvalue
            write(15,'(10x,a11)') 'pm1comswbox'
            do ibox = 1,numboxes
               call twh_pm1comswbox(GLB_GET,ibox,dvalue)
               dpboxarray(ibox) = dvalue
            enddo
            write(15,100) (dpboxarray(ibox),ibox=1,numboxes)
            write(15,'(10x,a12)') 'pm1comswpair'
            do ipair = 1,nmolpairs
               call twh_pm1comswpair(GLB_GET,ipair,dvalue)
               dpairarray(ipair) = dvalue
            enddo
            write(15,100) (dpairarray(ipair),ipair=1,nmolpairs)
         endif

         write(15,'(a6)') 'pmavb1'
         call twh_pmavb1(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a8)') 'pmavb1in'
         call twh_pmavb1in(GLB_GET,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a8)') 'pmavb1mt'
         do jmolty = 1,nmolty
            call twh_pmavb1mt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a8)') 'pmavb1ct'
         do jmolty = 1,nmolty
            do kmolty = 1,nmolty
               call twh_pmavb1ct(GLB_GET,jmolty,kmolty,dvalue)
               dmarray(kmolty) = dvalue
            enddo
            write(15,100) (dmarray(kmolty),kmolty=1,nmolty)
         enddo
         write(15,'(10x,a7)') 'avb1rad'
         call twh_avb1rad(GLB_GET,dvalue)
         write(15,101) dvalue
         
         write(15,'(a6)') 'pmavb2'
         call twh_pmavb2(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a8)') 'pmavb2in'
         call twh_pmavb2in(GLB_GET,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a8)') 'pmavb2mt'
         do jmolty = 1,nmolty
            call twh_pmavb2mt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a8)') 'pmavb2ct'
         do jmolty = 1,nmolty
            do kmolty = 1,nmolty
               call twh_pmavb2ct(GLB_GET,jmolty,kmolty,dvalue)
               dmarray(kmolty) = dvalue
            enddo
            write(15,100) (dmarray(kmolty),kmolty=1,nmolty)
         enddo
         write(15,'(10x,a7)') 'avb2rad'
         call twh_avb2rad(GLB_GET,dvalue)
         write(15,101) dvalue
         
         write(15,'(a6)') 'pmavb3'
         call twh_pmavb3(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a8)') 'pmavb3mt'
         do jmolty = 1,nmolty
            call twh_pmavb3mt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a8)') 'pmavb3ct'
         do jmolty = 1,nmolty
            do kmolty = 1,nmolty
               call twh_pmavb3ct(GLB_GET,jmolty,kmolty,dvalue)
               dmarray(kmolty) = dvalue
            enddo
            write(15,100) (dmarray(kmolty),kmolty=1,nmolty)
         enddo
         write(15,'(10x,a7)') 'avb3rad'
         call twh_avb3rad(GLB_GET,dvalue)
         write(15,101) dvalue
         
         write(15,'(a4)') 'pmcb'
         call twh_pmcb(GLB_GET,pmcb)
         write(15,44) pmcb
         write(15,'(10x,a6)') 'pmcbmt'
         do jmolty = 1,nmolty
            call twh_pmcbmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)

         write(15,'(10x,a5)') 'pmall'
         do jmolty = 1,nmolty
            call twh_pmall(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         
         write(15,'(a6)') 'pmback'
         call twh_pmback(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmbkmt'
         do jmolty = 1,nmolty
            call twh_pmbkmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)

         write(15,'(a8)') 'pmcbside'
         call twh_pmcbside(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a10)') 'pmcbsidemt'
         do jmolty = 1,nmolty
            call twh_pmcbsidemt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         
         write(15,'(a7)') 'pmpivot'
         call twh_pmpivot(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a7)') 'pmpivmt'
         do jmolty = 1,nmolty
            call twh_pmpivmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         
         write(15,'(a8)') 'pmconrot'
         call twh_pmconrot(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmcrmt'
         do jmolty = 1,nmolty
            call twh_pmcrmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)

         write(15,'(a8)') 'pmcrback'
         call twh_pmcrback(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a7)') 'pmcrbmt'
         do jmolty = 1,nmolty
            call twh_pmcrbmt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         
         write(15,'(a7)') 'pmplane'
         call twh_pmplane(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a10)') 'pmplanebox'
         do ibox = 1,numboxes
            call twh_pmplanebox(GLB_GET,ibox,dboxarray(ibox))
         enddo
         write(15,100) (dboxarray(ibox),ibox=1,numboxes)
         write(15,'(10x,a10)') 'planewidth'
         call twh_planewidth(GLB_GET,dvalue)
         write(15,100) dvalue
         
         write(15,'(a5)') 'pmrow'
         call twh_pmrow(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a8)') 'pmrowbox'
         do ibox = 1,numboxes
            call twh_pmrowbox(GLB_GET,ibox,dboxarray(ibox))
         enddo
         write(15,100) (dboxarray(ibox),ibox=1,numboxes)
         write(15,'(10x,a8)') 'rowwidth'
         call twh_rowwidth(GLB_GET,dvalue)
         write(15,100) dvalue
         
         write(15,'(a7)') 'pmtraat'
         call twh_pmtraat(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmtamt'
         do jmolty = 1,nmolty
            call twh_pmtamt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a6)') 'rmtraa'
         call twh_rmtraa(GLB_GET,1,1,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a6)') 'tatraa'
         call twh_tatraa(GLB_GET,dvalue)
         write(15,100) dvalue

         write(15,'(a11)') 'pmcomposite'
         call twh_pmcomposite(GLB_GET,dvalue)
         write(15,100) dvalue
         write(15,'(a6)') 'pmcomt'
         do jmolty = 1,nmolty
            call twh_pmcomt(GLB_GET,jmolty,dmarray(jmolty))
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(a8)') 'rmcomtra'
         call twh_rmcomtra(GLB_GET,1,1,dvalue)
         write(15,100) dvalue
         write(15,'(a8)') 'rmcomrot'
         call twh_rmcomrot(GLB_GET,1,1,dvalue)
         write(15,100) dvalue

         write(15,'(a7)') 'pmtracm'
         call twh_pmtracm(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmtcmt'
         do jmolty =1,nmolty
            call twh_pmtcmt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a6)') 'rmtrac'
         call twh_rmtrac(GLB_GET,1,1,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a6)') 'tatrac'
         call twh_tatrac(GLB_GET,dvalue)
         write(15,100) dvalue
         
         write(15,'(a8)') 'pmrotate'
         call twh_pmrotate(GLB_GET,dvalue)
         write(15,44) dvalue
         write(15,'(10x,a6)') 'pmromt'
         do jmolty = 1,nmolty
            call twh_pmromt(GLB_GET,jmolty,dvalue)
            dmarray(jmolty) = dvalue
         enddo
         write(15,100) (dmarray(jmolty),jmolty=1,nmolty)
         write(15,'(10x,a5)') 'rmrot'
         call twh_rmrot(GLB_GET,1,1,dvalue)
         write(15,100) dvalue
         write(15,'(10x,a5)') 'tarot'
         call twh_tarot(GLB_GET,dvalue)
         write(15,100) dvalue

         rwfile = 15
         call twh_rwcbmc(lfailure,WRITE_FLAG,rwfile)
         if ( lfailure ) then
            write(6,*) 'WRITETOWHEE: problem with cbmc output'
            return
         endif

c        --- start of potential output
         write(15,'(a8)') 'ffnumber'
         call twh_ffnumber(GLB_GET,ffnumber)
         write(15,*) ffnumber
         if ( ffnumber .gt. 0 ) then
            write(15,'(a11)') 'ff_filename'
            do ifile = 1,ffnumber
               call twh_ff_filename(GLB_GET,ifile,ffname)
               write(15,11) ffname
            enddo
         endif
         write(15,'(a19)') 'classical_potential'
         call twh_classical_potential(GLB_GET,localcp)
         write(15,'(a1,a30,a1)') char(39),localcp,char(39)
         
         call twh_setpotentiallogic(lfailure,lwritemixrule
     &        ,lwritelshift,lwriteltailc,lwritermin,lwritercut
     &        ,lwritercutin,lwriteinterpolatestyle,lwriterpd
     &        ,lwritescaling)
         if ( lfailure ) return
         
         if ( lwritemixrule ) then
            write(15,'(a17)') 'classical_mixrule'
            call twh_classical_mixrule(GLB_GET,classical_mixrule)
            write(15,'(a1,a30,a1)') char(39)
     &           ,classical_mixrule,char(39)
            if ( lwritescaling ) then
               write(15,'(a20)') 'cmix_rescaling_style'
               write(15,'(a1,a30,a1)') char(39),
     &              scp_get_scaling_style_string() ,char(39)
               if ( scp_get_scaling_style() .eq. SCP_SC ) then
                  write(15,'(a14)') 'cmix_lambda_lj'
                  write(15,*) scp_get_native_lambda(LAMBDA_LJ)
                  write(15,'(a13)') 'cmix_alpha_lj'
                  write(15,*) scp_get_native_lambda(ALPHA_LJ)
                  write(15,'(a13)') 'cmix_lambda_c'
                  write(15,*) scp_get_native_lambda(LAMBDA_C)
                  write(15,'(a15)') 'cmix_lprintdudl'
                  call twh_scp_eval_derivative(GLB_GET,lvalue)
                  write(15,*) lvalue
               endif
               if (scp_get_scaling_style() .eq. SCP_SC .or.
     &              scp_get_scaling_style() .eq. SCP_WCA ) then
                  call scp_get_scaled_atoms(npair, pair_list)
                  write(15,'(a10)') 'cmix_npair'
                  write(15,*) npair
                  write(15,'(a14)') 'cmix_pair_list'
                  do ipair = 1,npair
                     write(15,'(2(a1,a10,a1,1x))') (char(39)
     &                    ,pair_list(ipair,ientry),char(39)
     &                    ,ientry=1,2)
                  enddo
               elseif ( scp_get_scaling_style() .eq. SCP_NONE ) then
c                 --- nothing to output
               else
                  write(6,*) 
     &                 'WRITETOWHEE: unknown cmix_rescaling_style'
                  stop
               endif
            endif
         endif
         if ( lwriteinterpolatestyle ) then
            write(15,'(a16)') 'interpolatestyle'
            call twh_interpolatestyle(GLB_GET,interpolatestyle)
            write(15,'(a1,a20,a1)') char(39),interpolatestyle
     &           ,char(39)
         endif
         if ( lwritelshift ) then
            write(15,'(a6)') 'lshift'
            call twh_lshift(GLB_GET,lshift)
            write(15,*) lshift
         endif
         if ( lwriteltailc ) then
            write(15,'(a6)') 'ltailc'
            call twh_ltailc(GLB_GET,ltailc)
            write(15,*) ltailc
         endif
         if ( lwritermin ) then
            write(15,'(a4)') 'rmin'
            call twh_rmin(GLB_GET,rmin)
            write(15,*) rmin
         endif
         if ( lwritercut ) then
            write(15,'(a4)') 'rcut'
            call twh_rcut(GLB_GET,rcut)
            write(15,*) rcut
         endif
         if ( lwritercutin ) then
            write(15,'(a6)') 'rcutin'
            call twh_rcutin(GLB_GET,rcutin)
            write(15,*) rcutin
         endif
         if ( lwriterpd ) then
            write(15,'(a21)') 'radial_pressure_delta'
            call twh_radial_pressure_delta(GLB_GET,dvalue)
            write(15,*) dvalue
         endif
         write(15,'(a18)') 'electrostatic_form'
         call twh_electrostatic_form(GLB_GET,eform)
         write(15,'(a1,a50,a1)') char(39),eform,char(39)
         if ( eform .eq. 'none' ) then
c           --- nothing else to write here, but not an error
         elseif ( eform .eq. 'coulomb' ) then
            write(15,'(a12)') 'coulombstyle'
            call twh_coulombstyle(GLB_GET,coulombstyle)
            write(15,'(a1,a20,a1)') char(39),coulombstyle,char(39)
            if ( coulombstyle .eq. 'ewald_fixed_kmax' ) then
               write(15,'(a4)') 'kalp'
               call twh_ewald_kalp(GLB_GET,kalp)
               write(15,*) kalp
               call twh_ewald_kmax(GLB_GET,BXD_REAL,1,kmax)
               write(15,'(a4)') 'kmax'
               write(15,*) kmax
            elseif ( coulombstyle .eq. 'ewald_fixed_cutoff' ) then
               write(15,'(a10)') 'ewald_prec'
               call twh_ewald_prec(GLB_GET,ewald_prec)
               write(15,*) ewald_prec
               write(15,'(a7)') 'rcelect'
               call twh_rcelect(GLB_GET,BXD_REAL,1,dvalue)
               write(15,*) dvalue
            elseif ( coulombstyle .eq. 'minimum image' ) then
c              --- not an error, but nothing to output
            else
               write(6,*) 'WRITETOWHEE: unknown coulombstyle'
               lfailure = .true.
               return
            endif
            write(15,'(a7)') 'dielect'
            call twh_constant_qqfact(GLB_GET,qqfact)
            write(15,*) CON_EETOKCAL*CON_KCALTOK/qqfact
            write(15,'(a17)') 'lcompute_madelung'
            call twh_lcompute_madelung(GLB_GET,lvalue)
            write(15,*) lvalue
         else
            write(6,*) 'WRITETOWHEE: unknown electrostatic_form'
            write(6,*) eform
            lfailure = .true.
            return
         endif
         
         write(15,'(a6)') 'nfield'
         call twh_nfield(GLB_GET,nfield)
         write(15,*) nfield
         do ifld = 1,nfield
            write(15,'(a9)') 'fieldtype'
            call twh_fieldtypestring(GLB_GET,ifld,fieldstring)
            write(15,'(a1,a20,a1)') char(39),fieldstring,char(39)
            if ( fieldstring .eq. 'Hard Wall' ) then
               write(15,'(a6)') 'hrdbox'
               call twh_fieldint(GLB_GET,ifld,1,ivalue)
               write(15,*) ivalue
               write(15,'(a6)') 'hrdxyz'
               call twh_fieldchar(GLB_GET,ifld,1,fchar)
               write(15,'(3a1)') char(39),fchar,char(39)
               write(15,'(a6)') 'hrdcen'
               call twh_fielddp(GLB_GET,ifld,1,dvalue)
               write(15,*) dvalue
               write(15,'(a6)') 'hrdrad'
               call twh_fielddp(GLB_GET,ifld,2,dvalue)
               write(15,*) dvalue
               write(15,'(a19)') 'hrd_repulsion_style'
               call twh_fieldchar(GLB_GET,ifld,2,fchar)
               write(15,'(a1,a11,a1)') char(39),fchar,char(39)
               write(15,'(a15)') 'hrd_energy_type'
               call twh_fieldchar(GLB_GET,ifld,3,fchar)
               write(15,'(a1,a11,a1)') char(39),fchar,char(39)
               if ( fchar .eq. 'finite' ) then
                  write(15,'(a15)') 'hrd_wall_energy'
                  call twh_fielddp(GLB_GET,ifld,4,dvalue)
                  write(15,*) dvalue
               endif
            elseif ( fieldstring .eq. 'Hard Ball' ) then
               write(15,'(a6)') 'hrbbox'
               call twh_fieldint(GLB_GET,ifld,FLD_HRBBOX,ivalue)
               write(15,*) ivalue
               write(15,'(a6)') 'hrbpos'
               call twh_fielddp(GLB_GET,ifld,FLD_HRBPOSX,fieldtwo)
               call twh_fielddp(GLB_GET,ifld,FLD_HRBPOSY,fieldthree)
               call twh_fielddp(GLB_GET,ifld,FLD_HRBPOSZ,fieldfour)
               write(15,*) fieldtwo, fieldthree, fieldfour
               write(15,'(a6)') 'hrbrad'
               call twh_fielddp(GLB_GET,ifld,FLD_HRBRAD,dvalue)
               write(15,*) dvalue
               write(15,'(a6)') 'hrbfac'
               call twh_fielddp(GLB_GET,ifld,FLD_HRBFAC,dvalue)
               write(15,*) dvalue
               write(15,'(a15)') 'hrb_energy_type'
               call twh_fieldchar(GLB_GET,ifld,1,fchar)
               write(15,'(a1,a11,a1)') char(39),fchar,char(39)
               if ( fchar .eq. 'finite' ) then
                  write(15,'(a15)') 'hrb_energy_type'
                  call twh_fielddp(GLB_GET,ifld,FLD_HRBENERGY,dvalue)
                  write(15,*) dvalue
               endif
            elseif ( fieldstring .eq. 'Hooper Umbrella' ) then
               write(15,'(a6)') 'umbbox'
               call twh_fieldint(GLB_GET,ifld,1,ivalue)
               write(15,*) ivalue
               write(15,'(a6)') 'umbxyz'
               call twh_fieldchar(GLB_GET,ifld,1,fchar)
               write(15,'(3a1)') char(39),fchar,char(39)
               write(15,'(a9)') 'umbcenter'
               call twh_fielddp(GLB_GET,ifld,1,dvalue)
               write(15,*) dvalue
               write(15,'(a4)') 'umba'
               call twh_fielddp(GLB_GET,ifld,2,dvalue)
               write(15,*) dvalue
            elseif ( fieldstring .eq. 'LJ 9-3 Wall' ) then
               write(15,'(a6)') 'ljfbox'
               call twh_fieldint(GLB_GET,ifld,1,ivalue)
               write(15,*) ivalue
               write(15,'(a6)') 'ljfxyz'
               call twh_fieldchar(GLB_GET,ifld,1,fchar)
               write(15,'(3a1)') char(39),fchar,char(39)
               write(15,'(a6)') 'ljfcen'
               call twh_fielddp(GLB_GET,ifld,1,dvalue)
               write(15,*) dvalue
               write(15,'(a6)') 'ljfdir'
               call twh_fieldint(GLB_GET,ifld,2,ivalue)
               write(15,*) ivalue
               write(15,'(a6)') 'ljfcut'
               call twh_fielddp(GLB_GET,ifld,2,dvalue)
               write(15,*) dvalue
               write(15,'(a8)') 'ljfshift'
               call twh_fieldlogical(GLB_GET,ifld,lvalue)
               write(15,*) lvalue
               write(15,'(a6)') 'ljfrho'
               call twh_fielddp(GLB_GET,ifld,3,dvalue)
               write(15,*) dvalue
               write(15,'(a9)') 'ljfntypes'
               call twh_fieldint(GLB_GET,ifld,1,ivalue)
               write(15,*) ivalue
               do itype = 1,ivalue
                  call twh_fieldnbint(GLB_GET,ifld,itype,1,ntii)
                  write(15,'(a7)') 'ljfname'
                  call twh_nbname(GLB_GET,ntii,nbname)
                  write(15,'(a1,a10,a1)') char(39),nbname,char(39)
                  write(15,'(a6)') 'ljfsig'
                  call twh_fieldnbparam(GLB_GET,ifld,ntii,4,fnbparam)
                  write(15,*) fnbparam
                  write(15,'(a6)') 'ljfeps'
                  call twh_fieldnbparam(GLB_GET,ifld,ntii,5,fnbparam)
                  write(15,*) fnbparam
               enddo
            elseif ( fieldstring .eq. 'Steele Wall' ) then
               write(15,'(a10)') 'steele box'
               call twh_fieldint(GLB_GET,ifld,1,ivalue)
               write(15,*) ivalue
               write(15,'(a10)') 'steele xyz'
               call twh_fieldchar(GLB_GET,ifld,1,fchar)
               write(15,'(3a1)') char(39),fchar,char(39)
               write(15,'(a14)') 'steele surface'
               call twh_fielddp(GLB_GET,ifld,1,dvalue)
               write(15,*) dvalue
               write(15,'(a10)') 'steele dir'
               call twh_fieldint(GLB_GET,ifld,2,ivalue)
               write(15,*) ivalue
               write(15,'(a13)') 'steele cutoff'
               call twh_fielddp(GLB_GET,ifld,2,dvalue)
               write(15,*) dvalue
               write(15,'(a12)') 'steele shift'
               call twh_fieldlogical(GLB_GET,ifld,lvalue)
               write(15,*) lvalue
               write(15,'(a12)') 'steele delta'
               call twh_fielddp(GLB_GET,ifld,3,dvalue)
               write(15,*) dvalue
               write(15,'(a12)') 'steele rho_s'
               call twh_fielddp(GLB_GET,ifld,4,dvalue)
               write(15,*) dvalue
               write(15,'(a12)') 'steele ntype'
               call twh_fieldint(GLB_GET,ifld,1,ivalue)
               write(15,*) ivalue
               do itype = 1,ivalue
                  call twh_fieldnbint(GLB_GET,ifld,itype,1,ntii)
                  write(15,'(a11)') 'steele name'
                  call twh_nbname(GLB_GET,ntii,nbname)
                  write(15,'(a1,a10,a1)') char(39),nbname,char(39)
                  write(15,'(a8)') 'sigma_sf'
                  call twh_fieldnbparam(GLB_GET,ifld,ntii,2,fnbparam)
                  write(15,*) fnbparam
                  write(15,'(a10)') 'epsilon_sf'
                  call twh_fieldnbparam(GLB_GET,ifld,ntii,5,fnbparam)
                  write(15,*) fnbparam
               enddo
            elseif ( fieldstring .eq. 'Harmonic Attractor' ) then
               write(15,'(a6)') 'hafbox'
               call twh_fieldint(GLB_GET,ifld,1,ivalue)
               write(15,*) ivalue
               write(15,'(a4)') 'hafk'
               call twh_fielddp(GLB_GET,ifld,1,dvalue)
               write(15,*) dvalue
               write(15,'(a11)') 'hafentries'
               call twh_fieldint(GLB_GET,ifld,2,ivalue)
               write(15,*) ivalue
               write(15,'(a9)') 'hafrefpos'
               call twh_fieldchar(GLB_GET,ifld,1,fchar)
               write(15,'(a1,a7,a1)') char(39),fchar,char(39)
               if ( fchar .eq. 'Global') then
                  write(15,'(a10)') 'hafglobxyz'
                  call twh_fielddp(GLB_GET,ifld,2,fieldtwo)
                  call twh_fielddp(GLB_GET,ifld,3,fieldthree)
                  call twh_fielddp(GLB_GET,ifld,4,fieldfour)
                  write(15,'(3f10.6)') fieldtwo,fieldthree,fieldfour
               elseif ( fchar .ne. 'Initial') then
                  write(6,*) 'WRITETOWHEE: unknown hafrefpos'
                  stop
                  return
               endif
               write(15,'(a6)') 'hafkey'
               call twh_fieldchar(GLB_GET,ifld,2,fchar)
               write(15,'(a1,a7,a1)') char(39),fchar,char(39)
               call twh_fieldint(GLB_GET,ifld,2,ivalue)
               do ientry = 1,ivalue
                  write(15,'(a8)') 'hafmolec'
                  call twh_fieldnbint(GLB_GET,ifld,ientry,1,ivalue)
                  write(15,*) ivalue
                  call twh_fieldnbint(GLB_GET,ifld,ientry,2,ntii)
                  call twh_elementname(GLB_GET,ntii,elename)
                  if ( fchar .eq. 'Element') then
                     write(15,'(a10)') 'hafelement'
                     write(15,'(a1,a2,a1)') char(39)
     &                    ,elename,char(39)
                  elseif ( fchar .eq. 'FFtype') then
                     write(15,'(a7)') 'hafname'
                     call twh_nbname(GLB_GET,ntii,nbname)
                     write(15,'(a1,a10,a1)') char(39),nbname,char(39)
                  endif
               enddo
            else
               write(6,*) 'WRITETOWHEE: unknown fieldtype'
               write(6,*) fieldstring
               lfailure = .true.
               return
            endif
         enddo

         write(15,'(a15)') 'solvation_style'
         call twh_solvation_style(GLB_GET,cvalue)
         write(15,'(a1,a20,a1)') char(39),cvalue,char(39)
         if ( cvalue .ne. 'none' ) then
            call twh_solvation_type(GLB_GET,cvalue)
            write(15,'(a14)') 'solvation_type'
            write(15,'(a1,a20,a1)') char(39),cvalue,char(39)
         endif

      else
c        --- output the molecule information for towhee_input
         write(15,'(a11)') 'input_style'
         write(15,'(a8)') 'explicit'
         write(15,'(a5)') 'nunit'
         call twh_nunit(GLB_GET,imoltyflag,nunit)
         write(15,*) nunit
         write(15,'(a8)') 'nmaxcbmc'
         call twh_nmaxcbmc(GLB_GET,imoltyflag,ivalue)
         write(15,*) ivalue
         write(15,'(a9)') 'lpdbnames'
         write(15,*) lprotein
         
         do iatom = 1,nunit
            write(15,'(a17)') 'unit ntype qqatom'
            call twh_ntype(GLB_GET,imoltyflag,iatom,ntype)
            call twh_qqatom(GLB_GET,imoltyflag,iatom,qqatom)
            write(15,*) iatom,ntype,qqatom
            if ( lprotein ) then
c              --- include pdb stuff
               write(15,'(a27)') 'pdbname aminonum aminoshort'
               call twh_aminoshort(GLB_GET,imoltyflag,iatom,aminoshort)
               call twh_aminonum(GLB_GET,imoltyflag,iatom,aminonum)
               call twh_pdbname(GLB_GET,imoltyflag,iatom,pdbname)
               write(15,*) char(39),pdbname,char(39),' '
     &              ,aminonum,' ',char(39),aminoshort,char(39)
            endif
            write(15,'(a9)') 'vibration'
            call twh_invib(GLB_GET,imoltyflag,iatom,invibration)
            write(15,*) invibration
            do ivib = 1,invibration
               call twh_ijvib(GLB_GET,imoltyflag,iatom,ivib,ijvib)
               call twh_itvib(GLB_GET,imoltyflag,iatom,ivib,itvib)
               write(15,*) ijvib,itvib
            enddo
            write(15,'(a7)') 'bending'
            call twh_inben(GLB_GET,imoltyflag,iatom,inbend)
            write(15,*) inbend
            do ibend = 1,inbend
               call twh_orderben(GLB_GET,imoltyflag,iatom,ibend,iorder)
               call twh_itbang(GLB_GET,imoltyflag,iatom,ibend,itbang)
               call twh_orderbang(GLB_GET,imoltyflag,iatom,ibend,iobang)
               call twh_itben(GLB_GET,imoltyflag,iatom,ibend,itbentemp)
               call twh_ijben2(GLB_GET,imoltyflag,iatom,ibend,ijben2)
               call twh_ijben3(GLB_GET,imoltyflag,iatom,ibend,ijben3)
               write(15,*) ijben2,ijben3,itbentemp,iorder,itbang,iobang
            enddo
            write(15,'(a7)') 'torsion'
            call twh_intor(GLB_GET,imoltyflag,iatom,intorsion)
            write(15,*) intorsion
            do itor = 1,intorsion
               call twh_ijtor2(GLB_GET,imoltyflag,iatom,itor,ijtor2)
               call twh_ijtor3(GLB_GET,imoltyflag,iatom,itor,ijtor3)
               call twh_ijtor4(GLB_GET,imoltyflag,iatom,itor,ijtor4)
               call twh_ittor(GLB_GET,imoltyflag,iatom,itor,ittor)
               write(15,*) ijtor2,ijtor3,ijtor4,ittor
            enddo
            write(15,'(a11)') 'angle-angle'
            call twh_inaa(GLB_GET,imoltyflag,iatom,inaa)
            write(15,*) inaa
            do iaa = 1,inaa
               call twh_ijaa0(GLB_GET,imoltyflag,iatom,iaa,ijaa0)
               call twh_ijaa1(GLB_GET,imoltyflag,iatom,iaa,ijaa1)
               call twh_ijaa2(GLB_GET,imoltyflag,iatom,iaa,ijaa2)
               call twh_itaa(GLB_GET,imoltyflag,iatom,iaa,itaa)
               write(15,*) ijaa0,ijaa1,ijaa2,itaa
            enddo
            write(15,'(a16)') 'improper torsion'
            call twh_inimprop(GLB_GET,imoltyflag,iatom,inimproper)
            write(15,*) inimproper
            do itor = 1,inimproper
               call twh_ijimprop2(GLB_GET,imoltyflag,iatom,itor,ijimp2)
               call twh_ijimprop3(GLB_GET,imoltyflag,iatom,itor,ijimp3)
               call twh_ijimprop4(GLB_GET,imoltyflag,iatom,itor,ijimp4)
               call twh_itimprop(GLB_GET,imoltyflag,iatom,itor
     &              ,itimproper)
               write(15,*) ijimp2,ijimp3,ijimp4,itimproper
            enddo
         enddo
      endif
c     ---  format statements
 11   format(a MAXDIRLENGTH)
 44   format(f6.4)
 100  format(10x,20(e10.5,1x))
 101  format(10x,20(f7.4,1x))
      
c     write(6.*) 'finished WRITETOWHEE'
         
      return
      end
