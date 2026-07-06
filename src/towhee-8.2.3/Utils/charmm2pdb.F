c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004 Marcus G. Martin                            *
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
      program charmm2pdb
c     ******************************************************************
c     * this conversion routine takes the Charmm style pdb file for a  *
c     * polypeptide and puts the atoms back into the "standard" order  *
c     * the input file is specified by the user                        *
c     * the output file is standard.pdb and is suitable for further    *
c     * transformation using the pdb2towhee routine                    *
c     *                                                                *
c     * not sure when I originally wrote this routine                  *
c     * last modified 06-28-2004 by M.G. Martin                        *
c     ******************************************************************

      implicit none
      integer maxperpep
      parameter ( maxperpep = 100 )

      logical lspecial

      character*20 pdbfilename
      character*1 key
      character*6 atom
      character*4 amino,ltshort,hvshort
      dimension ltshort(maxperpep),hvshort(maxperpep)
      character*4 name,mntemp
      character*4 ltname,hvname,ltmname,hvmname
      dimension ltname(maxperpep),hvname(maxperpep)
      dimension ltmname(maxperpep),hvmname(maxperpep)

      integer aminonum,numprev,iprev,idum
      integer imax,count,light,heavy,itemp

      double precision hvxcord,hvycord,hvzcord
      dimension hvxcord(maxperpep)
      dimension hvycord(maxperpep)
      dimension hvzcord(maxperpep)
      double precision ltxcord,ltycord,ltzcord
      dimension ltxcord(maxperpep)
      dimension ltycord(maxperpep)
      dimension ltzcord(maxperpep)
      double precision xtemp,ytemp,ztemp,dduma,ddumb

      write(6,*) 'input the charmm file name'
      read(5,*) pdbfilename
      write(6,*) 'input the number of atoms in the file'
      read(5,*) imax

      open(20,file=pdbfilename,form='formatted')
      open(22,file='standard.pdb',form='formatted')

      numprev = 0
      light = 0
      heavy = 0
      iprev = 0
      do count = 1,imax+1
         if ( count .lt. imax+1) then
            read(20,*) atom,idum,name,amino,aminonum,xtemp,ytemp,ztemp
     &           ,dduma,ddumb,mntemp
         else
            aminonum = 0
         endif
         if ( numprev .ne. aminonum ) then
c           --- new peptide, output previous peptide info
c           --- if this is the first one then nothing is output
            lspecial = .false.
            if ( heavy .gt. 4 ) then
               if ( hvname(heavy-1) .eq. 'C' 
     &              .and. hvname(heavy) .eq. 'O' ) then
                  lspecial = .true.
               endif
            endif

            if ( lspecial ) then
c              --- peptide format of the C-O backbone
               do itemp = 1,2
                  iprev = iprev + 1
                  write (22,40)
     &                 'ATOM  '
     &                 ,iprev,hvname(itemp)
     &                 ,hvshort(itemp),numprev
     &                 ,hvxcord(itemp),hvycord(itemp),hvzcord(itemp)
     &                 ,1.0,2.2,hvmname(itemp)
               enddo
               do itemp = heavy-1,heavy
                  iprev = iprev + 1
                  write (22,40)
     &                 'ATOM  '
     &                 ,iprev,hvname(itemp)
     &                 ,hvshort(itemp),numprev
     &                 ,hvxcord(itemp),hvycord(itemp),hvzcord(itemp)
     &                 ,1.0,2.2,hvmname(itemp)
               enddo
               do itemp = 3,heavy-2
                  iprev = iprev + 1
                  write (22,40)
     &                 'ATOM  '
     &                 ,iprev,hvname(itemp)
     &                 ,hvshort(itemp),numprev
     &                 ,hvxcord(itemp),hvycord(itemp),hvzcord(itemp)
     &                 ,1.0,2.2,hvmname(itemp)
               enddo
            else
c              --- normal non-peptide format
               do itemp = 1,heavy
                  iprev = iprev + 1
                  write (22,40)
     &                 'ATOM  '
     &                 ,iprev,hvname(itemp)
     &                 ,hvshort(itemp),numprev
     &                 ,hvxcord(itemp),hvycord(itemp),hvzcord(itemp)
     &                 ,1.0,2.2,hvmname(itemp)
               enddo
            endif

            do itemp = 1,light
               iprev = iprev + 1
               write (22,40)
     &              'ATOM  '
     &              ,iprev,ltname(itemp)
     &              ,ltshort(itemp),numprev
     &              ,ltxcord(itemp),ltycord(itemp),ltzcord(itemp)
     &              ,1.0,2.2,ltmname(itemp)
            enddo
            heavy = 0
            light = 0
            numprev = aminonum
         endif

 40      format(A6,I5,1x,A4,1x,A3,2x,i4,4x,3f8.3
     &        ,f6.2,f6.2,4x,5x,A4) 
c        --- the first letter of the name is the key
         key = name(1:1)
         if ( key .eq. 'H' ) then
c          --- this is a light atom
            light = light + 1
            if ( light .gt. maxperpep ) then
               write(6,*) ' light exceeded maxperpep ',amino
               stop
            endif
            ltname(light) = name
            ltshort(light) = amino
            ltxcord(light) = xtemp
            ltycord(light) = ytemp
            ltzcord(light) = ztemp
            ltmname(light) = mntemp
         else
c           --- this is a heavy atom
            heavy = heavy + 1
            if ( heavy .gt. maxperpep ) then
               write(6,*) ' heavy exceeded maxperpep ',amino
               stop
            endif
            hvname(heavy) = name
            hvshort(heavy) = amino
            hvxcord(heavy) = xtemp
            hvycord(heavy) = ytemp
            hvzcord(heavy) = ztemp
            hvmname(heavy) = mntemp
         endif
      enddo

      close(20)
      close(22)

      end

