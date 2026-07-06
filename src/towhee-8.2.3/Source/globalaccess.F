#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2020 Marcus G. Martin                       *
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

c     ******************************************************************
c     * this file contains all of the functions and special subroutines*
c     * for accessing global data along with subroutines that bundle   *
c     * the adjustment of multiple global data structures              *
c     *                                                                *
c     * originally written 03-24-2006 by M.G. Martin                   *
c     * last modified 03-30-2020 by M.G. Martin                        *
c     ******************************************************************

      subroutine twh_create_torsion(imolty,iunit,itor2,itor3,itor4
     &     ,itype)
c     ******************************************************************
c     * creates a new torsion by incrementing the number of torsions,  *
c     * assigning the values of ijtor2, ijtor3, ijtor4; and setting the*
c     * new torsion type                                               *
c     *                                                                *
c     * originally written 03-26-2006 by M.G. Martin                   *
c     * last modified 03-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer imolty,iunit,itor2,itor3,itor4,itype
c     --- local variables
      integer itemp

      call twh_intor(GLB_INCR,imolty,iunit,1)
      call twh_intor(GLB_GET,imolty,iunit,itemp)
      if ( itemp .gt. MAXTOR ) then
         write(6,*) 'CREATE_TORSION: new torsion exceeds MAXTOR'
         stop
      endif
      call twh_ijtor2(GLB_SET,imolty,iunit,itemp,itor2)
      call twh_ijtor3(GLB_SET,imolty,iunit,itemp,itor3)
      call twh_ijtor4(GLB_SET,imolty,iunit,itemp,itor4)
      call twh_ittor(GLB_SET,imolty,iunit,itemp,itype)

      return
      end

      subroutine twh_set_aacoeff(itype,coeffs)
c     ******************************************************************
c     * assigns the aacoeff variables for this type                    *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffs
      dimension coeffs(MINAACOEFF:MAXAACOEFF)
c     --- local variables
      integer index

      do index = MINAACOEFF,MAXAACOEFF
         call twh_aacoeff(GLB_SET,itype,index,coeffs(index))
c        --- reset the coeff back to zero for the next assignment
         coeffs(index) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_aanames(itype,nsame,namearray)
c     ******************************************************************
c     * sets the number of quartets with the same parameters and also  *
c     * set the names of those quartets                                *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMEAA,4)
c     --- local variables
      integer isame,index

      call twh_naasame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do index = 1,4
            call twh_aanames(GLB_SET,itype,isame,index
     &           ,namearray(isame,index))
c           --- reset for next assignment
            namearray(isame,index) = '-'
         enddo
      enddo
      nsame = 0

      return
      end

      subroutine twh_set_bencoeff(itype,coeffs)
c     ******************************************************************
c     * assigns the bencoeff variables for this type                   *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-05-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffs
      dimension coeffs(MINBENDCOEFF:MAXBENDCOEFF)
c     --- local variables
      integer index

      do index = MINBENDCOEFF,MAXBENDCOEFF
         call twh_bencoeff(GLB_SET,itype,index,coeffs(index))
c        --- reset the coeff back to zero for the next assignment
         coeffs(index) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_bendnames(itype,nsame,namearray)
c     ******************************************************************
c     * assigns all of the bendnames and sets nbsame for this type of  *
c     * bending angle                                                  *
c     *                                                                *
c     * originally written 04-05-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMEB,MAXBNAME)
c     --- local variables
      integer isame,iatom

      call twh_nbsame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do iatom = 1,MAXBNAME
            call twh_bendnames(GLB_SET,itype,isame,iatom
     &           ,namearray(isame,iatom))
c           --- now reset the namearray for the next allocation
            namearray(isame,iatom) = '-'
         enddo
      enddo

      return
      end

      subroutine twh_set_eamdensnum(clocal,ipair,jpair)
c     ******************************************************************
c     * translates the eam_dens_style into a defined integer to speed  *
c     * the compares                                                   *
c     *                                                                *
c     * originally written 07-20-2013 by M.G. Martin                   *
c     * last modified 07-21-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) clocal
      integer ipair,jpair
c     --- local variables
      integer eam_dens_num

c     --- these definitions must match those in preproc.h
c     --- note all tabular potentials are numbered > 100
      if ( clocal .eq. 'Ackland cubic sum' ) then
         eam_dens_num = EAM_DENS_ACKLAND
      elseif ( clocal .eq. 'exponential' ) then
         eam_dens_num = EAM_DENS_EXPON
      elseif ( clocal .eq. 'table' ) then
         eam_dens_num = EAM_DENS_TABLE
      elseif ( clocal .eq. 'none' ) then
         eam_dens_num = EAM_DENS_NONE
      else 
         write(6,*) 'SET_EAMDENSNUM: invalid eamdensnum'
         write(6,'(a)') clocal
         stop
      endif

      call twh_eam_dens_num(GLB_SET,ipair,jpair,eam_dens_num)

      return
      end

      subroutine twh_set_eamembednum(clocal,ipair)
c     ******************************************************************
c     * translates the eam_embed_style into a defined integer to speed *
c     * the compares                                                   *
c     *                                                                *
c     * originally written 07-21-2013 by M.G. Martin                   *
c     * last modified 07-21-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) clocal
      integer ipair
c     --- local variables
      integer itest

c     --- these definitions must match those in preproc.h
c     --- note all tabular potentials are numbered > 100
      if ( clocal .eq. 'Belashchenko' ) then
         itest = EAM_EMBED_BELASH
      elseif ( clocal .eq. 'Belashchenko 10 rho' ) then
         itest = EAM_EMBED_BELASHTRHO
      elseif ( clocal .eq. 'logarithmic' ) then
         itest = EAM_EMBED_LOGARITHMIC
      elseif ( clocal .eq. 'power 0.5 and 2' ) then
         itest = EAM_EMBED_POWHT
      elseif ( clocal .eq. 'power 0.5, 2, and 4' ) then
         itest = EAM_EMBED_POWHTF
      elseif ( clocal .eq. 'squareroot' ) then
         itest = EAM_EMBED_SR
      elseif ( clocal .eq. 'table' ) then
         itest = EAM_EMBED_TABLE
      else 
         write(6,*) 'SET_EAMEMBEDNUM: invalid eamembedstyle'
         write(6,'(a)') clocal
         stop
      endif

      call twh_eam_embed_num(GLB_SET,ipair,itest)

      return
      end

      subroutine twh_set_eampairnum(clocal,ipair)
c     ******************************************************************
c     * translates the eam_pair_style into a defined integer to speed  *
c     * the compares                                                   *
c     *                                                                *
c     * originally written 07-20-2013 by M.G. Martin                   *
c     * last modified 07-21-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) clocal
      integer ipair
c     --- local variables
      integer eam_pair_num

c     --- these definitions must match those in preproc.h
c     --- note all tabular potentials are numbered > 100
      if ( clocal .eq. 'Ackland 3-part' ) then
         eam_pair_num = EAM_PAIR_ACKLANDTP
      elseif ( clocal .eq. 'Ackland Power' ) then
         eam_pair_num = EAM_PAIR_ACKLANDPOW
      elseif ( clocal .eq. 'exponential' ) then
         eam_pair_num = EAM_PAIR_EXPON
      elseif ( clocal .eq. 'morse' ) then
         eam_pair_num = EAM_PAIR_MORSE
      elseif ( clocal .eq. 'table' ) then
         eam_pair_num = EAM_PAIR_TABLE
      elseif ( clocal .eq. 'Belashchenko 5-8' ) then
         eam_pair_num = EAM_PAIR_BELASHFE
      elseif ( clocal .eq. 'none' ) then
         eam_pair_num = EAM_PAIR_NONE
      else 
         write(6,*) 'SET_EAMPAIRNUM: invalid eampairnum'
         write(6,'(a)') clocal
         stop
      endif

      call twh_eam_pair_number(GLB_SET,ipair,eam_pair_num)

      return
      end

      subroutine twh_set_hbondcoeff(itype,coeffarray)
c     ******************************************************************
c     * sets up all of the coefficients for a single special hbond type*
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffarray
      dimension coeffarray(MINHBONDCOEFF:MAXHBONDCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINHBONDCOEFF,MAXHBONDCOEFF
         call twh_hbondcoeff(GLB_SET,itype,icoeff,coeffarray(icoeff))
c        --- now reset the coeffarray for the next allocation
         coeffarray(icoeff) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_hbondnames(itype,namearray)
c     ******************************************************************
c     * sets up all of the atom names for a single special hbond type  *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      character*(FFNAMELEN) namearray
      dimension namearray(MAXHBONDNAME)
c     --- local variables
      integer icoeff

      do icoeff = 1,2
         call twh_hbondnames(GLB_SET,itype,icoeff,namearray(icoeff))
c        --- now reset the namearray for the next allocation
         namearray(icoeff) = '-'
      enddo

      return
      end

      subroutine twh_set_impcoeff(itype,coeffarray)
c     ******************************************************************
c     * assigns all of the improper torsion coefficients and then      *
c     * initializes the coeffarray for the next batch                  *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffarray
      dimension coeffarray(MINIMPCOEFF:MAXIMPCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINIMPCOEFF,MAXIMPCOEFF
         call twh_impcoeff(GLB_SET,itype,icoeff,coeffarray(icoeff))
c        --- now reset the coeffarray for the next allocation
         coeffarray(icoeff) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_impnames(itype,nsame,namearray)
c     ******************************************************************
c     * assigns all of the impnames and sets nimpsame for this type of *
c     * improper torsion angle                                         *
c     *                                                                *
c     * originally written 04-13-2006 by M.G. Martin                   *
c     * last modified 04-13-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMEIMP,MAXIMPNAME)
c     --- local variables
      integer isame,iatom

      call twh_nimpsame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do iatom = 1,MAXIMPNAME
            call twh_impnames(GLB_SET,itype,isame,iatom
     &           ,namearray(isame,iatom))
c           --- now reset the namearray for the next allocation
            namearray(isame,iatom) = '-'
         enddo
      enddo

      return
      end

      subroutine twh_set_isolvtype(lfailure,solvation_style
     &     ,solvation_type)
c     ******************************************************************
c     * sets the integer code for isolvtype based on the               *
c     * solvation_style and the solvation_type                         *
c     *                                                                *
c     * originally written 01-12-2007 by M.G. Martin                   *
c     * last modified 01-12-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      character*(*) solvation_style,solvation_type
c     --- local variables
      integer isolv_type

      if ( solvation_style .eq. 'none' ) then
         isolv_type = SOLV_NONE
      else
         if ( solvation_type .eq. 'tramonto' ) then
            isolv_type = SOLV_TRAMONTO
         elseif ( solvation_type .eq. 'SASA' ) then
            isolv_type = SOLV_SASA
         elseif ( solvation_type .eq. 'EEF1' ) then
            isolv_type = SOLV_EEF1
         else
            write(6,*) 'SET_ISOLVTYPE: unknown solvation_type'
            lfailure = .true.
            return
         endif
      endif

      call twh_isolvtype(GLB_SET,isolv_type)

      return
      end

      subroutine twh_set_names(itype,nbname,bondname,anglename
     &     ,torsionname)
c     ******************************************************************
c     * assigns all of the names for each nonbonded atom type for use  *
c     * in the automatic assignment of parameters                      *
c     *                                                                *
c     * originally written 04-12-2006 by M.G. Martin                   *
c     * last modified 04-17-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      character*(*) nbname,bondname,anglename,torsionname

c     --- safety check type
      if ( itype .gt. NNTYPE ) then
         write(6,*) 'insufficient memory for this nonbond assignment'
         write(6,*) 'requested element:',itype
         write(6,*) 'NNTYPE:',NNTYPE
         stop
      endif

      call twh_nbname(GLB_SET,itype,nbname)
      call twh_cbname(GLB_SET,itype,bondname)
      call twh_caname(GLB_SET,itype,anglename)
      call twh_ctname(GLB_SET,itype,torsionname)
c     --- reset all of the names to error values
      nbname = 'init-error'
      bondname = 'init-error'
      anglename = 'init-error'
      torsionname = 'init-error'

      return
      end

      subroutine twh_set_ofnames(itype,namearray)
c     ******************************************************************
c     * assigns all of the one-five names for a single type            *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      character*(FFNAMELEN) namearray
      dimension namearray(MAXOFNAME)
c     --- local variables
      integer icoeff

      do icoeff = 1,MAXOFNAME
         call twh_ofnames(GLB_SET,itype,icoeff,namearray(icoeff))
c        --- reset all of the names to error values
         namearray(icoeff) = '-'
      enddo

      return
      end

      subroutine twh_set_qbinames(itype,namearray)
c     ******************************************************************
c     * assigns all of the qbinames for this type                      *
c     *                                                                *
c     * originally written 04-14-2006 by M.G. Martin                   *
c     * last modified 04-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      character*(FFNAMELEN) namearray
      dimension namearray(2)
c     --- local variables
      integer iatom

      do iatom = 1,2
         call twh_qbinames(GLB_SET,itype,iatom,namearray(iatom))
c        --- now reset the namearray for the next allocation
         namearray(iatom) = '-'
      enddo

      return
      end

      subroutine twh_set_nbcoeff(itype,coeffarray)
c     ******************************************************************
c     * assigns all of the nonbonded coefficients for this cross type  *
c     *                                                                *
c     * originally written 04-11-2006 by M.G. Martin                   *
c     * last modified 04-11-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffarray
      dimension coeffarray(MINNBCOEFF:MAXNBCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINNBCOEFF,MAXNBCOEFF
         call twh_nbcoeff(GLB_SET,icoeff,itype,coeffarray(icoeff))
c        --- reset for the next assignments
         coeffarray(icoeff) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_torcoeff(itype,coeffarray)
c     ******************************************************************
c     * assigns all of the torsion coefficients and then initializes   *
c     * the coeffarray for the next batch                              *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision coeffarray
      dimension coeffarray(MINTORCOEFF:MAXTORCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINTORCOEFF,MAXTORCOEFF
         call twh_torcoeff(GLB_SET,itype,icoeff,coeffarray(icoeff))
c        --- now reset the coeffarray for the next allocation
         coeffarray(icoeff) = 0.0d0
      enddo

      return
      end

      subroutine twh_set_tornames(itype,nsame,namearray)
c     ******************************************************************
c     * assigns all of the tornames and sets ntsame for this type of   *
c     * torsion angle                                                  *
c     *                                                                *
c     * originally written 04-10-2006 by M.G. Martin                   *
c     * last modified 04-10-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMET,MAXTNAME)
c     --- local variables
      integer isame,iatom

      call twh_ntsame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do iatom = 1,MAXTNAME
            call twh_tornames(GLB_SET,itype,isame,iatom
     &           ,namearray(isame,iatom))
c           --- now reset the namearray for the next allocation
            namearray(isame,iatom) = '-'
         enddo
      enddo

      return
      end

      subroutine twh_set_vibcoeffs(itype,vibarray)
c     ******************************************************************
c     * assigns all of the vibcoeffs for a single vibration type       *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype
      double precision vibarray
      dimension vibarray(MINVIBCOEFF:MAXVIBCOEFF)
c     --- local variables
      integer icoeff

      do icoeff = MINVIBCOEFF,MAXVIBCOEFF
         call twh_vibcoeff(GLB_SET,itype,icoeff,vibarray(icoeff))
      enddo

      return
      end

      subroutine twh_set_vibnames(itype,nsame,namearray)
c     ******************************************************************
c     * assigns all of the vibnames and sets nvsame for this type of   *
c     * vibration                                                      *
c     *                                                                *
c     * originally written 03-30-2006 by M.G. Martin                   *
c     * last modified 03-30-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer itype,nsame
      character*(FFNAMELEN) namearray
      dimension namearray(MAXSAMEV,MAXVNAME)
c     --- local variables
      integer isame,iatom

      call twh_nvsame(GLB_SET,itype,nsame)
      do isame = 1,nsame
         do iatom = 1,MAXVNAME
            call twh_vibnames(GLB_SET,itype,isame,iatom
     &           ,namearray(isame,iatom))
c           --- now reset the namearray for the next allocation
            namearray(isame,iatom) = '-'
         enddo
      enddo

      return
      end

c     ******************************************************************
c     * functions that primarily make it easier to access the global   *
c     * data structures in a manner that is similar to accessing an    *
c     * array                                                          *
c     * also includes functions that resemble the fetching of global   *
c     * data structures                                                *
c     ******************************************************************

      function twh_get_laahere(iaatype)
c     ******************************************************************
c     * returns true if this angle-angle type appears in a molecule    *
c     * that is being used in this simulation                          *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-28-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LAAHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer iaatype
c     --- local variables
      integer imolty,iunit,iaa,nmolty,nunit,inaa,itaa

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_inaa(GLB_GET,imolty,iunit,inaa)
            do iaa = 1,inaa
               call twh_itaa(GLB_GET,imolty,iunit,iaa,itaa)
               if ( itaa .eq. iaatype ) then
                  twh_get_laahere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_laahere = .false.

      return
      end

      function twh_get_lbendhere(bendtype)
c     ******************************************************************
c     * returns true if this angle type appears in a molecule that     *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-27-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LBENDHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer bendtype
c     --- local variables
      integer imolty,iunit,ibend,nmolty,nunit,batype,ibttest,inbend

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_inben(GLB_GET,imolty,iunit,inbend)
            do ibend = 1,inbend
               call twh_itben(GLB_GET,imolty,iunit,ibend,ibttest)
               if ( ibttest .eq. bendtype ) then
                  twh_get_lbendhere = .true.
                  return
               endif
               call twh_itbang(GLB_GET,imolty,iunit,ibend,batype)
               if ( batype .eq. bendtype ) then
                  twh_get_lbendhere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_lbendhere = .false.

      return
      end

      function twh_get_lhere(nbtype)
c     ******************************************************************
c     * returns true if this nonbonded type appears in a molecule that *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer nbtype
c     --- local variables
      integer imolty,iunit,nmolty,nunit,itype

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_ntype(GLB_GET,imolty,iunit,itype)
            if ( itype .eq. nbtype ) then
               twh_get_lhere = .true.
               return
            endif
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_lhere = .false.

      return
      end

      function twh_get_limphere(imptype)
c     ******************************************************************
c     * returns true if this improper type appears in a molecule that  *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-28-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LIMPHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer imptype
c     --- local variables
      integer imolty,iunit,iimp,nmolty,nunit,inimproper,itimproper

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_inimprop(GLB_GET,imolty,iunit,inimproper)
            do iimp = 1,inimproper
               call twh_itimprop(GLB_GET,imolty,iunit,iimp,itimproper)
               if ( itimproper .eq. imptype ) then
                  twh_get_limphere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_limphere = .false.

      return
      end

      function twh_get_lofhere(oftype)
c     ******************************************************************
c     * returns true if this one-five type appears in a molecule       *
c     * that is being used in this simulation                          *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-28-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LOFHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer oftype
c     --- local variables
      integer imolty,iunit,iof,nmolty,nunit,inonefive,itonefive

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_inof(GLB_GET,imolty,iunit,inonefive)
            do iof = 1,inonefive
               call twh_itof(GLB_GET,imolty,iunit,iof,itonefive)
               if ( itonefive .eq. oftype ) then
                  twh_get_lofhere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_lofhere = .false.

      return
      end

      function twh_get_loftor(ittype)
c     ******************************************************************
c     * gets the loftor data                                           *
c     *                                                                *
c     * originally written 03-29-2006 by M.G. Martin                   *
c     * last modified 03-29-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LOFTOR
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ittype

      call twh_loftor(GLB_GET,ittype,twh_get_loftor)
  
      return
      end

      function twh_get_ltorhere(tortype)
c     ******************************************************************
c     * returns true if this torsion type appears in a molecule that   *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-28-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LTORHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer tortype
c     --- local variables
      integer imolty,iunit,itor,nmolty,nunit,intorsion,ittorsion

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_intor(GLB_GET,imolty,iunit,intorsion)
            do itor = 1,intorsion
               call twh_ittor(GLB_GET,imolty,iunit,itor,ittorsion)
               if ( abs(ittorsion) .eq. abs(tortype) ) then
                  twh_get_ltorhere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_ltorhere = .false.

      return
      end

      function twh_get_lvibhere(vibtype)
c     ******************************************************************
c     * returns true if this vibration type appears in a molecule that *
c     * is being used in this simulation                               *
c     * returns false otherwise                                        *
c     * this used to be handled as global data, but was modified to    *
c     * replace the memory with a moderately expensive loop as this    *
c     * function/data is rarely used                                   *
c     *                                                                *
c     * originally written 03-31-2006 by M.G. Martin                   *
c     * last modified 03-20-2020 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_LVIBHERE
#include "functions.h"
c     --- variables passed to/from the function
      integer vibtype
c     --- local variables
      integer imolty,iunit,ivib,nmolty,nunit,invibration,itvibtest

      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
         do iunit = 1,nunit
            call twh_invib(GLB_GET,imolty,iunit,invibration)
            do ivib = 1,invibration
               call twh_itvib(GLB_GET,imolty,iunit,ivib,itvibtest)
               if ( itvibtest .eq. vibtype ) then
                  twh_get_lvibhere = .true.
                  return
               endif
            enddo
         enddo
      enddo
c     --- if we made it this far then there were no matches
      twh_get_lvibhere = .false.

      return
      end

      function twh_limited_double(dvalue,limit)
c     ******************************************************************
c     * check that the magnitude of a parameters has an absolute value *
c     * that is sufficiently large to reduce the noise between machines*
c     *                                                                *
c     * originally written 12-17-2007 by M.G. Martin                   *
c     * last modified 12-22-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_LIMITED_DOUBLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision dvalue,limit

      if ( dabs(dvalue) .lt. limit ) then
         twh_limited_double = 0.0d0
      else
         twh_limited_double = dvalue
      endif
  
      return
      end

      function twh_safe_double(dvalue)
c     ******************************************************************
c     * checks to see that the double precision dvalue is larger than  *
c     * SMALLEST.  If it is larger then it returns the value of dvalue *
c     * If it is smaller than SMALLEST it returns 0.0d0                *
c     * This is used when outputing information to reduce machine      *
c     * precision noise                                                *
c     *                                                                *
c     * originally written 10-30-2007 by M.G. Martin                   *
c     * last modified 10-30-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_SAFE_DOUBLE
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision dvalue

      if ( dabs(dvalue) .lt. SMALLEST ) then
         twh_safe_double = 0.0d0
      else
         twh_safe_double = dvalue
      endif
  
      return
      end

