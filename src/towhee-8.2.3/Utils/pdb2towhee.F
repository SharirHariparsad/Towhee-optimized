c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2003 Marcus G. Martin                       *
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
      program pdb2towhee
c     ******************************************************************
c     * this routine converts a pdb format file into the input files   *
c     * needed to start a towhee simulation                            *
c     * last modified 05-18-2010 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- local variables
      logical lshift
      character*20 pdbfilename
      character*5 c_aminonum,c_numprev
      character*4 name
      character*3 amino
      integer icdum,aminocount,imax,count
      double precision xmin,ymin,zmin,xcord,ycord,zcord
      double precision xmax,ymax,zmax

      write(6,*) 'Please enter the type of input file (integer)'
      write(6,*) '1 - protein pdb file'
      write(6,*) '2 - generic pdb file'
      write(6,*) '3 - generic charmm output file'
      write(6,*) '4 - charmm output file (DNA)'
      write(6,*) '5 - charmm output file (RNA)'
      read(5,*) icdum
      write(6,*) 'Please enter the pdb file name (character string)'
      read(5,*) pdbfilename
      write(6,*) 'Please enter the number of atoms (integer)'
      read(5,*) imax
      write(6,*) 'Do you wish to shift the positions of the atoms?'
      read(5,*) lshift

      open(20,file=pdbfilename,form='formatted')
      open(21,file='towhee_altinp',form='formatted')
      open(22,file='towhee_coords',form='formatted')

     
      aminocount = 0
      c_numprev = '0'
      do count = 1,imax

         call readinfo(icdum,xcord,ycord,zcord,c_aminonum,amino,name)

         if ( count .eq. 1 .or. xcord .lt. xmin ) xmin = xcord
         if ( count .eq. 1 .or. xcord .gt. xmax ) xmax = xcord
         if ( count .eq. 1 .or. ycord .lt. ymin ) ymin = ycord
         if ( count .eq. 1 .or. ycord .gt. ymax ) ymax = ycord
         if ( count .eq. 1 .or. zcord .lt. zmin ) zmin = zcord
         if ( count .eq. 1 .or. zcord .gt. zmax ) zmax = zcord
         
         if ( c_aminonum .ne. c_numprev ) then
            aminocount = aminocount + 1

            if ( icdum .eq. 1 ) then
c              --- this is a new amino acid, output towhee name
 40            format(a1,a2,a1,1x,a1,a1,a1,1x,a1,1x,a1,a1,a1)
               if ( amino .eq. 'ALA' ) then
                  write(21,40) char(39),'a0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'ASN' ) then
                  write(21,40) char(39),'n0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'CYS' ) then
                  write(21,40) char(39),'c0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'CYD' ) then
                  write(21,40) char(39),'cs',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'?'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'CYP' ) then
                  write(21,40) char(39),'cp',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'GLN' ) then
                  write(21,40) char(39),'q0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'GLY' ) then
                  write(21,40) char(39),'g0',char(39)
     &                 ,char(39),'r',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'ILE' ) then
                  write(21,40) char(39),'i0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'LEU' ) then
                  write(21,40) char(39),'l0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'DLE' ) then
                  write(21,40) char(39),'l0',char(39)
     &                 ,char(39),'d',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'MET' ) then
                  write(21,40) char(39),'m0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'PHE' ) then
                  write(21,40) char(39),'f0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'PRO' ) then
                  write(21,40) char(39),'p0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'SER' ) then
                  write(21,40) char(39),'s0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'THR' ) then
                  write(21,40) char(39),'t0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'TRP' ) then
                  write(21,40) char(39),'w0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'TYR' ) then
                  write(21,40) char(39),'y0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'VAL' ) then
                  write(21,40) char(39),'v0',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'DVA' ) then
                  write(21,40) char(39),'v0',char(39)
     &                 ,char(39),'d',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'ASP' ) then
                  write(21,40) char(39),'d-',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'GLU' ) then
                  write(21,40) char(39),'e-',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'ARG' ) then
                  write(21,40) char(39),'r+',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'HIS' ) then
                  write(21,40) char(39),'h?',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'HSP' ) then
                  write(21,40) char(39),'h+',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'LYS' ) then
                  write(21,40) char(39),'k+',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'LYR' ) then
                  write(21,40) char(39),'kr',char(39)
     &                 ,char(39),'l',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'FOR' ) then
                  write(21,40) char(39),'zf',char(39)
     &                 ,char(39),'r',char(39),
     &                 '0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'ETA' ) then
                  write(21,40) char(39),'ze',char(39)
     &                 ,char(39),'r',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'ACE' ) then
                  write(21,40) char(39),'za',char(39)
     &                 ,char(39),'r',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'NMA' ) then
                  write(21,40) char(39),'zn',char(39)
     &                 ,char(39),'r',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               elseif ( amino .eq. 'NH2' ) then
                  write(21,40) char(39),'zm',char(39)
     &                 ,char(39),'r',char(39)
     &                 ,'0'
     &                 ,char(39),'-',char(39)
               else
                  write(6,*) 'unknown amino acid',amino
                  stop
               endif
            elseif (icdum .eq. 4) then
               if (amino .eq. 'CYT') then
                  write(21,*) 'dC'
               else if (amino .eq. 'GUA') then
                  write(21,*) 'dG'
               else if (amino .eq. 'ADE') then
                  write(21,*) 'dA'
               else if (amino .eq. 'THY') then
                  write(21,*) 'dT'
               else 
                  write(6,*) 'unknown deoxyribonuclease',amino
                  write(6,*) 'if meant to be a RNA strand'
                  write(6,*) 'then run again and choose input type 5'
               endif
            elseif (icdum .eq. 5) then
               if (amino .eq. 'CYT') then
                  write(21,*) 'C'
               else if (amino .eq. 'GUA') then
                  write(21,*) 'G'
               else if (amino .eq. 'ADE') then
                  write(21,*) 'A'
               else if (amino .eq. 'URA') then
                  write(21,*) 'U'
               else 
                  write(6,*) 'unknown ribonuclease',amino
                  write(6,*) 'if meant to be a DNA strand'
                  write(6,*) 'then run again and choose input type 4'
               endif
            endif                
            c_numprev = c_aminonum
         endif
      enddo
      
      rewind(20)
      write(6,*) 'original xmin,ymin,zmin',xmin,ymin,zmin
      write(6,*) 'original xmax,ymax,zmax',xmax,ymax,zmax
      if ( .not. lshift ) then
c        --- reset the min to zero so we don't shift
         xmin = 0.0d0
         ymin = 0.0d0
         zmin = 0.0d0
      endif
      write(6,*) 'number of groups for builder',aminocount
      if ( icdum .lt. 4 ) then
         write(6,*) 'Note: you must check all of the histidines'
         write(6,*) 'to replace the h? with the proper protonation'
         write(6,*) 'And, you need to change the last variable - to N'
         write(6,*) 'for the n-terminus and C for the c-terminus'
      endif

      do count = 1,imax
         call readinfo(icdum,xcord,ycord,zcord,c_aminonum,amino,name)
         write(22,100)
     &        xcord-xmin,ycord-ymin,zcord-zmin
     &        ,char(39),amino,char(39),char(39),name,char(39)
      enddo
c     --- write a dummy line
      amino = 'nul'
      name = 'nul'
      write(22,100) 0.0,0.0,0.0
     &        ,char(39),amino,char(39),char(39),name,char(39)

 100  format(3(f10.5,1x),a1,a3,a1,1x,a1,a4,a1) 

      close(20)
      close(21)
      close(22)

      write(6,*) 'new xmax,ymax,zmax',xmax-xmin,ymax-ymin,zmax-zmin

      end

      subroutine readinfo(icdum,xcord,ycord,zcord,c_aminonum,amino,name)
      implicit none
c     --- variables passed to/from the subroutine
      integer icdum
      double precision xcord,ycord,zcord
      character*5 c_aminonum,c_numprev
      character*4 name
      character*3 amino
c     --- local variables
      character*6 atom
      integer iunit

      if ( icdum  .eq. 1 .or. icdum .eq. 2 ) then
         read(20,12) atom,iunit,name,amino,c_aminonum,xcord,ycord
     &        ,zcord
 12      format(A6,I5,1x,A4,1x,A3,2x,A4,4x,3f8.3)
      elseif ( icdum .eq. 3 .or. icdum .eq. 4 .or. icdum .eq. 5 ) then
         read(20,*) iunit,c_aminonum,amino,name,xcord,ycord,zcord
      else
         write(6,*) 'invalid value for icdum'
         stop
      endif

      return
      end
