      program jre_to_towhee
c     ******************************************************************
c     * converts from the spead format to the towhee format            *
c     *                                                                *
c     * originally written 12-01-2005 by M.G. Martin                   *
c     * last modified 12-01-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
      character*5 atomname,firstpart,secondpart
      character*20 filename
      integer outfile,infile,ival
      integer iline,imain,isubclass,totline
      double precision wellvalues
      dimension wellvalues(4)
      double precision sigma,mweight,bondrad
      
      infile = 10
      write(6,*) 'Enter the file name to transform'
      read(5,*) filename
      open(infile,file=filename,form='formatted')

      outfile = 11
      filename = 'output'
      open(outfile,file=filename,form='formatted')

c     --- first line has the total number of lines
      read(infile,*) totline
      write(6,*) 'totline',totline
      do iline = 1,totline
         read(infile,*) imain,isubclass,(wellvalues(ival),ival=1,4)
     &        ,sigma,mweight,bondrad,atomname
         write(outfile,*)
         write(outfile,'(6x,a15)') 'type = type + 1'
         write(outfile,'(6x,a13,f10.5)') 'mass(type) = ',mweight
         write(outfile,'(6x,a16,a1,a1)') 'element(type) = '
     &        ,char(39),char(39)
         write(outfile,'(6x,a17,a1,a4,a1)') 'bondpatt(type) = '
     &        ,char(39),'null',char(39)
c        --- names
         write(outfile,'(6x,a19)') 'nbff(type) = ffname'
         write(outfile,'(6x,a15,a1,a5,a1)') 'nbname(type) = '
     &        ,char(39),atomname,char(39)
c        --- nonbonded names are a combination of the main and subclass
c        --- integers
         write(firstpart,'(i5)') imain
         write(secondpart,'(i5)') isubclass
         write(outfile,'(6x,a15,a1,a5,a5,a1)') 'cbname(type) = '
     &        ,char(39),firstpart,secondpart,char(39)
         write(outfile,'(6x,a15,a1,a5,a5,a1)') 'caname(type) = '
     &        ,char(39),firstpart,secondpart,char(39)
         write(outfile,'(6x,a15,a1,a5,a5,a1)') 'ctname(type) = '
     &        ,char(39),firstpart,secondpart,char(39)
         write(outfile,'(6x,a29)') 'ii = twh_getnbtype(type,type)'
c        --- sigma
c        --- convert from nm to A
         write(outfile,'(6x,a16,e15.7)') 'nbcoeff(1,ii) = ',sigma*10.0d0
         write(outfile,'(6x,a16,e15.7)') 'nbcoeff(2,ii) = '
     &        ,-wellvalues(1)
         write(outfile,'(6x,a16,e15.7)') 'nbcoeff(3,ii) = '
     &        ,-wellvalues(2)
         write(outfile,'(6x,a16,e15.7)') 'nbcoeff(4,ii) = '
     &        ,-wellvalues(3)
         write(outfile,'(6x,a16,e15.7)') 'nbcoeff(5,ii) = '
     &        ,-wellvalues(4)
c        --- stash the bondrad in a nonbond coefficient until we figure
c        --- out what to do with it
         write(outfile,'(6x,a16,e15.7)') 'nbcoeff(6,ii) = ',bondrad

      enddo

      close(infile)
      close(outfile)
      end
