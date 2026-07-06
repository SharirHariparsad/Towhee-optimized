c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2020 Marcus G. Martin                       *
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
      program xmg2towhee
c     ******************************************************************
c     * this program converts the embedded atom potential files from   *
c     * the xmd format into the Towhee format                          *
c     * for more information on the xmd format please see their web    *
c     * http://www.ims.uconn.edu/centers/simul/pot/potformat.htm       *
c     * in order to make this work you have to do a bit of reordering  *
c     * of their files so all the pair potentials are first, then the  *
c     * embedding and then the densities                               *
c     *                                                                *
c     * originally written 11-26-2002 by M.G. Martin                   *
c     * last modified 03-02-2020 by M.G. Martin                        *
c     ******************************************************************
#define MAXDATA 2000
#define MAXTYPES 3
      implicit none
      character*2 typenames
      dimension typenames(MAXTYPES)
      character*5 eunit,energy,junk
      character*10 forcefield
      character*17 morejunk
      character*30 filename
      integer iname,jname,index,ntypes,idata,icount,nextra,ncycle
      integer ndata_pair,ndata_dens,ndata_embed
      dimension ndata_pair(MAXTYPES,MAXTYPES)
      dimension ndata_dens(MAXTYPES),ndata_embed(MAXTYPES)
      integer icycle,nxstep,itype,jtype
      double precision eam_pair,eam_dens,eam_embed
      dimension eam_pair(2,MAXDATA,MAXTYPES,MAXTYPES)
      dimension eam_dens(2,MAXDATA,MAXTYPES)
      dimension eam_embed(2,MAXDATA,MAXTYPES)
      double precision xdiff,xincr

c     --- filename passed in when running the code
      write(6,*) 'Please enter a filename'
      read(5,*) filename
      open(unit=10,file=filename,form="formatted")
      open(unit=12,file="towhee_eam",form="formatted")

c     --- skip the first line (avaiable for comments)
      read(10,*)
      read(10,'(a5,1x,a5)') eunit,energy
      write(6,*) eunit,' ',energy
      read(10,'(a5,1x,i10)') junk,ntypes
      write(6,*) junk,ntypes
      do itype = 1,ntypes
         read(10,'(a2)') typenames(itype)
      enddo

c     --- read in all the pair potentials
      do itype = 1,ntypes
         do jtype = itype,ntypes
c           --- read in the first line with the endpoints
            read(10,1) morejunk,iname,jname,ndata_pair(itype,jtype)
     &           ,eam_pair(1,1,itype,jtype)
     &           ,eam_pair(1,ndata_pair(itype,jtype),itype,jtype)
 1          format(a16,1x,i3,i3,i5,2(e13.6,1x))

c           --- safety check ndata
            if ( ndata_pair(itype,jtype) .gt. MAXDATA ) then
               write(6,*) 'ndata ',ndata_pair(itype,jtype)
     &              ,' exceeds MAXDATA',MAXDATA
               stop
            endif
c           ---- safety check iname and jname
            if ( iname .ne. itype .or. jname .ne. jtype ) then
               write(6,*) 'mismatch of atom types',iname,itype
               write(6,*) jname,jtype
               stop
            endif
c           --- fill in the x-values
            xdiff = eam_pair(1,ndata_pair(itype,jtype),itype,jtype) 
     &           - eam_pair(1,1,itype,jtype)
            nxstep = ndata_pair(itype,jtype)-1
            xincr = xdiff/dble(nxstep)
            do idata = 2,ndata_pair(itype,jtype)
               eam_pair(1,idata,itype,jtype) 
     &              = eam_pair(1,idata-1,itype,jtype) + xincr
            enddo
         
            ncycle = (ndata_pair(itype,jtype))/4
            nextra = mod(ndata_pair(itype,jtype),4)

            do icycle = 1,ncycle
               icount = (icycle-1)*4
               read(10,2) (eam_pair(2,idata,itype,jtype)
     &              ,idata=icount+1,icount+4)
            enddo
            if ( nextra .ne. 0 ) then
               icount = (ncycle-1)*4
               read(10,2) (eam_pair(2,idata,itype,jtype)
     &              ,idata=icount+1,icount+nextra)
            endif
 2          format(4(1x,e13.6))
         enddo
      enddo

c     --- read in all the desity functions
      do itype = 1,ntypes
c        --- read in the first line with the endpoints
         read(10,11) morejunk,iname,ndata_dens(itype)
     &        ,eam_dens(1,1,itype)
     &        ,eam_dens(1,ndata_dens(itype),itype)
 11      format(a16,1x,i3,i5,2(e13.6,1x))
         
c        --- safety check ndata_dens
         if ( ndata_dens(itype) .gt. MAXDATA ) then
            write(6,*) 'ndata_dens ',ndata_dens(itype)
     &           ,' exceeds MAXDATA',MAXDATA
            stop
         endif
c        ---- safety check iname
         if ( iname .ne. itype ) then
            write(6,*) 'mismatch of atom types',iname,itype
            stop
         endif
c        --- fill in the x-values
         xdiff = eam_dens(1,ndata_dens(itype),itype) 
     &        - eam_dens(1,1,itype)
         nxstep = ndata_dens(itype)-1
         xincr = xdiff/dble(nxstep)
         do idata = 2,ndata_dens(itype)
            eam_dens(1,idata,itype) 
     &           = eam_dens(1,idata-1,itype) + xincr
         enddo
         
         ncycle = (ndata_dens(itype))/4
         nextra = mod(ndata_dens(itype),4)

         do icycle = 1,ncycle
            icount = (icycle-1)*4
            read(10,2) (eam_dens(2,idata,itype)
     &           ,idata=icount+1,icount+4)
         enddo
         if ( nextra .ne. 0 ) then
            icount = (ncycle-1)*4
            read(10,2) (eam_dens(2,idata,itype)
     &           ,idata=icount+1,icount+nextra)
         endif
      enddo

c     --- read in all the embedding functions
      do itype = 1,ntypes
c        --- read in the first line with the endpoints
         read(10,21) morejunk,iname,ndata_embed(itype)
     &        ,eam_embed(1,1,itype)
     &        ,eam_embed(1,ndata_embed(itype),itype)
 21      format(a17,1x,i3,i5,2(e13.6,1x))
c        --- safety check ndata_embed
         if ( ndata_embed(itype) .gt. MAXDATA ) then
            write(6,*) 'ndata_embed ',ndata_embed(itype)
     &           ,' exceeds MAXDATA',MAXDATA
            stop
         endif
c        ---- safety check iname
         if ( iname .ne. itype ) then
            write(6,*) 'mismatch of atom types',iname,itype
            stop
         endif
c        --- fill in the x-values
         xdiff = eam_embed(1,ndata_embed(itype),itype) 
     &        - eam_embed(1,1,itype)
         nxstep = ndata_embed(itype)-1
         xincr = xdiff/dble(nxstep)
         do idata = 2,ndata_embed(itype)
            eam_embed(1,idata,itype) 
     &           = eam_embed(1,idata-1,itype) + xincr
         enddo
         
         ncycle = (ndata_embed(itype))/4
         nextra = mod(ndata_embed(itype),4)

         do icycle = 1,ncycle
            icount = (icycle-1)*4
            read(10,2) (eam_embed(2,idata,itype)
     &           ,idata=icount+1,icount+4)
         enddo
         if ( nextra .ne. 0 ) then
            icount = (ncycle-1)*4
            read(10,2) (eam_embed(2,idata,itype)
     &           ,idata=icount+1,icount+nextra)
         endif
      enddo

c     --- output to the towhee_eam file
      write(12,'(a10)') 'EAMversion'
      write(12,'(i1)') 1
      write(12,'(a10)') 'Forcefield'
      forcefield = 'Chen_Co'
      write(12,'(a10)') forcefield

      write(12,'(a5)') 'Types'
      write(12,*) ntypes
      write(12,'(a4)') 'Atom'
      do itype = 1,ntypes
         write(12,'(a2)') typenames(itype)
      enddo
      write(12,'(a4)') 'Mass'
      do itype = 1,ntypes
         write(12,'(f7.4)') 58.9332d0
      enddo

      do itype = 1,ntypes
         do jtype = itype,ntypes
            write(12,'(a8)') 'eam_pair'
            write(12,*) itype,jtype,ndata_pair(itype,jtype)
            do idata = 1,ndata_pair(itype,jtype)
               write(12,3) eam_pair(1,idata,itype,jtype)
     &              ,eam_pair(2,idata,itype,jtype)
            enddo
         enddo
      enddo
 3    format(f13.8,1x,e13.6)

      do itype = 1,ntypes
         write(12,'(a8)') 'eam_dens'
         write(12,*) itype,ndata_dens(itype)
         do idata = 1,ndata_dens(itype)
            write(12,3) eam_dens(1,idata,itype)
     &           ,eam_dens(2,idata,itype)
         enddo
      enddo
      
      do itype = 1,ntypes
         write(12,'(a9)') 'eam_embed'
         write(12,*) itype,ndata_embed(itype)
         do idata = 1,ndata_embed(itype)
            write(12,3) eam_embed(1,idata,itype)
     &           ,eam_embed(2,idata,itype)
         enddo
      enddo
      
      end
