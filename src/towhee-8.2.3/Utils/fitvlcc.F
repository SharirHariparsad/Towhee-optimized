c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2013 Marcus G. Martin                       *
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
c     ---- pound defines for the entire file
#define MAXDATA 500
c     --- begin actual code
      program fitvlcc
c     ******************************************************************
c     * this utility routine is used to process groups of towhee_vlcc  *
c     * files combined in a specific format and output the critical    *
c     * temperatures and densities.                                    *
c     *                                                                *
c     * Here is an example file demonstrating the file format          *
c     ******************************************************************
c     * An example input file and the script for turning it into VLCC  *
c     * fitting and plotting information for TraPPE-UA ethane can be   *
c     * found in the                                                   *
c     * Examples/VLCC_Fit                                              *
c     * directory of the code distribution                             *
c     *                                                                *
c     * This code has been passed around for over a decade and I'm not *
c     * sure who originally wrote it.  It has been extensively         *
c     * modified from its original form                                *
c     * last modified 07-26-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
      character*30 fileroot,name,fitvars
      common /filestuff/ fileroot
      logical lfitnonlinear,lfitrectilinear,lfound
      logical ltransformexp,lfit_a0_a1_b0_b1
      integer inp,ns,i,pressoption
      double precision beta,tc,dtc,dumma,dummb
      double precision slope,sm,intercept,sb
c     --- dp arrays
      double precision temperature,rhog,errg,rhol,errl
      dimension  temperature(MAXDATA),rhog(MAXDATA),errg(MAXDATA)
      dimension rhol(MAXDATA),errl(MAXDATA)
      double precision rhoa,erra
      dimension rhoa(MAXDATA),ERRA(MAXDATA)
      double precision press,errp
      dimension press(MAXDATA),errp(MAXDATA)
      double precision tmin
      double precision xxx,yyy,www
      dimension xxx(MAXDATA),yyy(MAXDATA),www(MAXDATA)

      write(6,*) 'Enter the root file name'
      read(5,*) fileroot

      inp = 11
      name = '.sim'
      call twh_openappendfile(inp,name)

      write(6,*) 'Vapor pressure analysis (0 for none, 1 for virial '
     &     ,'2 for thermodynamic pressure)'
      read(5,*) pressoption
      write(6,*) 'Input minimum temperature for plots'
      read(5,*) tmin
      write(6,*) 
      read(inp,*)                                         
      read(inp,*)                                         
      read(inp,*) beta                                    
      read(inp,*)                                         

      lfound = .false.
      i = 0
      do while ( i .lt. MAXDATA .and. .not. lfound )
         i = i + 1
         if ( pressoption .eq. 1 ) then
c           --- also read in the virial pressure information
            read(inp,*) temperature(I),rhog(I),errg(I)
     &           ,rhol(I),errl(I),press(i),errp(i)
         elseif ( pressoption .eq. 2 ) then
c           --- read in the thermodynamic pressure information
            read(inp,*) temperature(I),rhog(I),errg(I)
     &           ,rhol(I),errl(I),dumma,dummb,press(i),errp(i)
         else
c           --- just read in the density data
            read(inp,*) temperature(I),rhog(I),errg(I),rhol(I),errl(I)
         endif
         if (temperature(I) .EQ. 0.0d0 ) then
            NS=I-1                                                
            lfound = .true.
         endif                                                    
      enddo
c     --- compute the average/diameter properties
      do i = 1,ns
         rhoa(i) = ( rhog(i) + rhol(i) ) / 2.0
         erra(i) = dsqrt( errg(i)**2 + errl(i)**2 ) / 2.0d0
      enddo

c     --- open the files for the output of simulation data
      name = '.vsat'
      call twh_openappendfile(21,name)
      name = '.lsat'
      call twh_openappendfile(22,name)
      name = '.asat'
      call twh_openappendfile(23,name)
c     --- plot simulation results
      do I=1,NS
         write(21,*) rhog(i),temperature(i),errg(i)
         write(22,*) rhol(i),temperature(i),errl(i)
         write(23,*) rhoa(i),temperature(i),erra(i)
      enddo
      close(21)
      close(22)
      close(23)

c     --- vapor pressure transformation
      if ( pressoption .eq. 1 .or. pressoption .eq. 2 ) then
c        --- open the vapor pressure files
         name = '.tvp'
         call twh_openappendfile(35,name)
         name = '.cc'
         call twh_openappendfile(36,name)
         do i = 1,ns
            write(35,*) temperature(i),press(i),errp(i)
c           --- transform data for claus-clap fit
            xxx(i) = 1.0d0/temperature(i)
            yyy(i) = dlog(press(i))
            www(i) = 1.0d0/((errp(i)/press(i))**2 )
            write(36,*) xxx(i),yyy(i)
         enddo
         close(35)
         close(36)
c        --- call the least squares fitting for CC
         call ccfit(ns,xxx,yyy,www,slope,sm,intercept,sb)
      endif
      write(6,*)
      write(6,*) 'Do you wish to transform experimental data into'
     &     ,' a plotable format?'
      read(5,*) ltransformexp
      if ( ltransformexp ) call transformexpdata(beta)

c     --- see if we want to do a rectilinear fit
      write(6,*) 'Do you wish to fit to the following equations?'
      write(6,*) 'tau = 1 - T/Tc'
      write(6,*) 'rho_diff = a_0 tau^beta'
      write(6,*) 'rho_diam = rho_c + b_1 tau'
      read(5,*) lfitrectilinear
      if ( lfitrectilinear ) then
         call rectilinearfit(ns,rhog,rhol,rhoa,errg,errl,erra
     &        ,temperature,beta,tmin,tc,dtc)
c        --- output critical pressure
         if ( pressoption .ne. 0 ) then
            call critpress(tc,dtc,slope,sm,intercept,sb)
         endif
      endif

c     --- see if we want to try a nonlinear fit
      write(6,*) 'Do you wish to fit to the following equations?'
      write(6,*) 'tau = 1 - T/Tc'
      write(6,*) 'rho_diff = a_0 tau^beta'
      write(6,*) 'rho_diam = rho_c + b_0 tau^(1-alpha) + b_1 tau'
      read(5,*) lfitnonlinear
      if ( lfitnonlinear ) then
         fitvars = 'a_0 b_0 b_1'
         call nonlinearfit(ns,rhog,rhol,temperature,tmin,fitvars
     &        ,tc)
c        --- output critical pressure (use the rectilinear dtc)
         if ( pressoption .ne. 0 ) then
            call critpress(tc,dtc,slope,sm,intercept,sb)
         endif
      endif

c     --- see if we want to try an expanded fit
      write(6,*) 'Do you wish to fit to the following equations?'
      write(6,*) 'tau = 1 - T/Tc'
      write(6,*) 'rho_diff = a_0 tau^beta + a_1 tau^(beta + delta)'
      write(6,*) 'rho_diam = rho_c + b_0 tau^(1-alpha) + b_1 tau'
      read(5,*) lfit_a0_a1_b0_b1
      if ( lfit_a0_a1_b0_b1 ) then
         fitvars = 'a_0 a_1 b_0 b_1'
         call nonlinearfit(ns,rhog,rhol,temperature,tmin,fitvars
     &        ,tc)
c        --- output critical pressure (use the rectilinear dtc)
         if ( pressoption .ne. 0 ) then
            call critpress(tc,dtc,slope,sm,intercept,sb)
         endif
      endif


c     --- close the input and output files
      close(11)

      end                                                              

      subroutine fit(X,Y,ndata,SIG,A,B,SIGA,SIGB,CHI2,Q,COV)
c     ******************************************************************
c     * performs a weighted fit to the data                            *
c     *                                                                *
c     * original author unknown                                        *
c     * last modified 07-21-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ndata
      double precision x,y,sig,a,b,siga,sigb,chi2,q,cov
      DIMENSION X(ndata),Y(ndata),SIG(ndata)
c     --- local variables
      integer i
      double precision sx,sy,st2,ss,wt,sxoss,t,gammq

      SX=0.
      SY=0.
      ST2=0.
      B=0.
      SS=0.
      do I = 1,ndata
         WT = 1./(SIG(I)**2)
         SS = SS + WT
         SX = SX + X(I)*WT                
         SY = SY + Y(I)*WT                
      enddo
      SXOSS = SX/SS                      
      do I = 1,ndata                
         T = (X(I)-SXOSS)/SIG(I)        
         ST2 = ST2+T*T                  
         B = B + T*Y(I)/SIG(I)            
      enddo
      B = B/ST2                          
      A = (SY-SX*B)/SS                   
      SIGA = sqrt((1.+SX*SX/(SS*ST2))/SS)
      SIGB = sqrt(1./ST2)                
      CHI2 = 0.0d0
      COV = -SX/(SS*ST2)                 
      do I = 1,ndata                
         CHI2 = CHI2 + ((Y(I)-A-B*X(I))/SIG(I))**2
      enddo
      Q = GAMMQ(0.5*dble(ndata-2),0.5*CHI2)
      return
      end

      FUNCTION GAMMQ(A,X)              
      implicit none
c     --- variables passed to/from the function
      double precision gammq,a,x
c     --- local variables
      double precision gammcf,gln,gamser
c     --- initialize variables
      gln = 0.0d0
      gammcf = 0.0d0

      if(X.LT.0..OR.A.LE.0.) then 
         write(6,*) 'values unexpectedly small'
      endif
      if(X.LT.A+1.)then                
        call GSER(GAMSER,A,X,GLN)      
        GAMMQ=1.-GAMSER                
      else                             
        call GCF(GAMMCF,A,X,GLN)       
        GAMMQ=GAMMCF                   
      endif                            
      return                           
      end                              

      subroutine GSER(GAMSER,A,X,GLN)  
      implicit none
c     --- variables passed to/from the subroutine
      double precision gamser,a,x,gln
c     --- parameters
      integer itmax
      double precision eps
      PARAMETER (ITMAX=100,EPS=3.E-7)  
c     --- local variables
      integer n
      double precision gammln,ap,sum,del

      GLN=GAMMLN(A)                    
      if(X.LE.0.)then                  
        if(X.LT.0.) then
           write(6,*) 'values unexpectedly small'
        endif
        GAMSER=0.                      
        return                         
      endif                            
      AP=A                             
      SUM=1./A                         
      DEL=SUM                          
      do 11 N=1,ITMAX                  
        AP=AP+1.                       
        DEL=DEL*X/AP                   
        SUM=SUM+DEL                    
        if(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    continue                         
      write(6,*) 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)  
      return                           
      end   
                           
      subroutine GCF(GAMMCF,A,X,GLN)   
      implicit none
c     --- variables passed to/from the subroutine
      double precision gammcf,a,x,gln
c     --- parameters
      integer itmax
      double precision eps
      PARAMETER (ITMAX=100,EPS=3.E-7)  
c     --- local variables
      integer n
      double precision gammln,gold,a0,a1,b0,b1,fac,an,ana,anf,g

      GLN=GAMMLN(A)                    
      GOLD=0.                          
      A0=1.                            
      A1=X                             
      B0=0.                            
      B1=1.                            
      FAC=1.                           
      do 11 N=1,ITMAX                  
        AN=dble(N)                    
        ANA=AN-A                       
        A0=(A1+A0*ANA)*FAC             
        B0=(B1+B0*ANA)*FAC             
        ANF=AN*FAC                     
        A1=X*A0+ANF*A1                 
        B1=X*B0+ANF*B1                 
        if(A1.NE.0.)then               
          FAC=1./A1                    
          G=B1*FAC                     
          if(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G                       
        endif                          
11    continue       
      stop 'problem in GCF'
1     GAMMCF=EXP(-X+A*dLOG(X)-GLN)*G   
      return                           
      end                              

      FUNCTION GAMMLN(XX)              
      implicit none
c     --- variables passed to/from the function
      double precision gammln,xx
c     --- local variables
      integer j
      double precision COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE                         
      TMP=X+FPF                        
      TMP=(X+HALF)*LOG(TMP)-TMP        
      SER=ONE                          
      do 11 J=1,6
        X=X+ONE                        
        SER=SER+COF(J)/X               
11    continue                         
      GAMMLN=TMP+LOG(STP*SER)          
      return                           
      end                              

      subroutine twh_openappendfile(number,name)
c     ******************************************************************
c     * appends the name variable to the filefoot common and then opens*
c     * a file with that name                                          *
c     *                                                                *
c     * originally written some time in 2005 by M.G. Martin            *
c     * last modified 08-24-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer number
      character*30 name
c     --- common blocks
      character*30 fileroot
      common /filestuff/ fileroot
c     --- local variables
      integer flength
      character*30 filename
      flength = index(fileroot,' ') - 1
c     --- open the file
      filename = fileroot(1:flength)//name
      open(number,file=filename,form='formatted')
      return
      end

      subroutine ccfit(npoint,xxx,yyy,www,slope,sm,intercept,sb)
c     ******************************************************************
c     * fits the pressures to a linear equation                        *
c     *                                                                *
c     * last modified 07-26-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      integer npoint
      double precision xxx(npoint),yyy(npoint),www(npoint)
      double precision slope,sm,intercept,sb
c     --- local variables
      integer i
      double precision sumx,sumy,sumw,sumdw,di
      dimension di(npoint)
      double precision press,delpress,tboil,deltboil,pboil
      double precision hvap,dhvap
c     Weighted linear least squares fit
      
      write(6,*) 
      write(6,*) 'Weighted Linear Least Squares Fit for Pressures'
      sumx = 0.0d0
      sumy = 0.0d0
      sumw  = 0.0d0
      do i = 1, npoint
         sumx = sumx + xxx(i) * www(i)
         sumy = sumy + yyy(i) * www(i)
         sumw = sumw + www(i)
      enddo
      sumx = sumx / sumw
      sumy = sumy / sumw

      sumdw = 0.0d0
      do i = 1, npoint
         sumdw = sumdw + www(i) * ( xxx(i) - sumx )**2
      enddo
         
      slope = 0.0d0
      do i = 1, npoint
         slope = slope + www(i) * ( xxx(i) - sumx ) * yyy(i)
      enddo
      slope = slope / sumdw
      
      intercept = sumy - slope * sumx
         
      do i = 1, npoint
         di(i) = yyy(i) - slope * xxx(i) - intercept
      enddo
         
      sm = 0.0d0
      do i = 1, npoint
         sm = sm + www(i) * di(i)**2
      enddo
      sm = sm / ( dble(npoint-2) * sumdw )
      sm = dsqrt(sm)
         
      sb = 0.0d0
      do i = 1, npoint
         sb = sb + www(i) * di(i)**2
      enddo
      sb = ( (1.0d0/sumw) + (sumx**2/sumdw) ) * sb / dble(npoint-2)
      sb = dsqrt(sb)
         
      write(6,*)
      write(6,*) 'SLOPE =     ', slope, '  +/-  ', sm
      write(6,*) 'INTERCEPT = ', intercept,  '  +/-  ', sb
      write(6,*)

c     --- transform the slope into a heat of vaporization
      hvap = -slope*8.314d0/1000.0d0
      dhvap = dabs((sm/slope)*hvap)
      write(6,*) 'Delta Hvap [kJ/mol] ',hvap,' error',dhvap

c     --- estimate boiling point from this relationship
      pboil = 101.325d0
      tboil = slope / ( dlog(pboil) - intercept )
      deltboil = tboil*dsqrt((sm/slope)**2 
     &     + (sb/(dlog(pboil)-intercept))**2)
      write(6,*) 'Boiling Point [K] ',tboil,' error',deltboil


      return
      end

      subroutine critpress(temp,deltemp,slope,sm,intercept,sb)
c     ******************************************************************
c     * computes the critical pressure using the simple linear fit     *
c     *                                                                *
c     * split out of old routines 07-26-2013 by M.G. Martin            *
c     * last modified 07-26-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed in/out of the subroutine
      double precision temp,deltemp,slope,sm,intercept,sb
c     --- local variables
      double precision press,delpress

      press = dexp(slope/temp + intercept)
      delpress = dsqrt( ((slope/temp)**2)*
     &     (  (sm/slope)**2 + (deltemp/temp)**2 ) + sb**2 ) 
      delpress = delpress*press
      write(6,*) 'Critical pressure [kPa]:',press,' error',delpress 

      return
      end

      subroutine transformexpdata(beta)
c     ******************************************************************
c     * transforms some standard experimental data formats into the    *
c     * format used for the Towhee simulation results                  *
c     *                                                                *
c     * split out of main routines 08-23-2005 by M.G. Martin           *
c     * last modified 08-23-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision beta
c     --- local variables
      integer itest
      character*30 name
      double precision tt,um,dum,rliq,rgas,psat

      write(6,*) 'Experimental Data Format'
      write(6,*) '1: Smith and Srivastava data format'
      write(6,*) '2: Hensel and Warren data format'
      read(5,*) itest
      if ( itest .eq. 1 .or. itest .eq. 2 ) then
c        --- open the files for processing experimental data
         name = '.exp'
         call twh_openappendfile(13,name)
         name = '.resultsex'
         call twh_openappendfile(66,name)
         name = '.xexpg'
         call twh_openappendfile(28,name)
         name = '.xexpl'
         call twh_openappendfile(29,name)

         tt=1.
         read(13,*) 
         read(13,*) um
         read(13,*)
         read(13,*)
     
         write(66,*)
         write(66,*)
         write(66,*) beta
         write(66,*)
         dum=1.0d0
         write(6,*) 'molecular mass = ',um
 50      if (tt.gt.0.0000001) then
            if ( itest .eq. 1 ) then
               read(13,*) tt,rliq,rgas
               rgas=um/rgas
               rliq=um/rliq
            elseif ( itest .eq. 2 ) then
               read(13,*) tt,psat,rliq,rgas
c              --- shift from celsius to kelvin
               tt = tt + 273.15d0
            endif
            if (tt.gt.0.00001) then
               if ( rgas .ne. um ) write(28,*) rgas,tt
               if ( rliq .ne. um ) write(29,*) rliq,tt
            endif
            write(66,*) tt,rgas,dum,rliq,dum,(rliq+rgas)/2.,dum,
     +           (rliq-rgas),dum
            goto 50
         endif
c        --- close these files
         close(13)
         close(28)
         close(29)
         close(66)
      else
         write(6,*) 'unknown format'
      endif
      
      return
      end

      subroutine rectilinearfit(ns,rhog,rhol,rhoa,errg,errl,erra
     &     ,temperature,beta,tmin,tc,dtc)
c     ******************************************************************
c     * fits the data to the scaling law and the law of rectilinear    *
c     * diameters                                                      *
c     *                                                                *
c     * originally written prior to 1995 by persons unknown            *
c     * last modified 07-26-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ns
      double precision tc,beta,tmin,dtc
      double precision rhog,rhol,temperature,errg,errl
      dimension rhog(MAXDATA),rhol(MAXDATA),temperature(MAXDATA)
      dimension errl(MAXDATA),errg(MAXDATA)
      double precision rhoa,erra
      dimension rhoa(MAXDATA),ERRA(MAXDATA)
c     --- local variables
      integer npoint,i
      character*30 name
      double precision ar,br,sigar,sigbr,chi2r,qr,covr,b_1,a_0,rhoc
      double precision chi2s,qs,covs,as,bs,sigas,sigbs,tt
      double precision da,db,drhoc,sumrho,difrho
      double precision rmsd_scale,rmsd_diameter,tau
      double precision rho1,rho2,rho3,delt

      double precision delrho,rhos,errs
      dimension RHOS(MAXDATA),ERRS(MAXDATA)

c     --- open the output files for rectilinear fit
      name = '.xcrit'
      call twh_openappendfile(24,name)
      name = '.xfitg'
      call twh_openappendfile(25,name)
      name = '.xfitl'
      call twh_openappendfile(26,name)
      name = '.xfita'
      call twh_openappendfile(27,name)
      name = '.errortc'
      call twh_openappendfile(30,name)


c     --- compute and store some quantities
      do i = 1,ns
         delrho = ( rhol(I)-rhog(I) )/2.0d0
         rhos(i) = delrho
         errs(i) = 0.5d0 * dsqrt( errl(i)**2 + errg(i)**2 )
         rhos(I) = rhos(I)**(1./beta)                               
         errs(I) = rhos(i) * (1.0d0/beta) * errs(i) / delrho
      enddo

c     --- weighted fit to the scaling law where we have rewritten the
c     --- equation as
c     --- [(rhol - rhog)/2]**(1/beta) = a_0^(1/beta) - T a_0^(1/beta)/Tc
      call fit(temperature,rhos,NS,errs,AS,BS,SIGAS,SIGBS,CHI2S,QS,COVS)
c     --- extract critial temperature as the ratio of the scaling 
c     --- constant divided by the scaling slope
      tc = -as/bs
      dtc = tc * dsqrt( (sigas/as)**2 + (sigbs/bs)**2 )
c     --- extract a_0 from the scaling constant
      a_0 = (AS)**beta                                     
      da = ((AS**beta) * beta * sigas / as)

c     --- weighted fit to the rectilinear law where we rewrite the 
c     --- equation so we don't have to know Tc ahead of time
c     --- despite the fact we already know Tc at this point
c     --- rhoa = (rhoc + b_1) - T b_1/Tc
      call fit(temperature,rhoa,NS,ERRA,AR,BR,SIGAR,SIGBR,CHI2R,QR,COVR)
c     --- extract b_1 from the rectilinear slope and the critical temp
      b_1 = -br * tc
      db = b_1 * dsqrt( (sigbr/br)**2 + (dtc/tc)**2 )
c     --- extract rhoc from the rectilinear constant and b_1
      rhoc = ar - b_1
      drhoc = dsqrt( sigar**2 + db**2 )

      write(30,*) rhoc,tc,drhoc,dtc
      write(6,1003) NS,TC,DTC,RHOC,DRHOC,a_0,da,b_1,db
 1003 format('  =========== RESULTS:',/,
     +       '     NUMBER OF SAMPLES : ',I9  ,/,
     +       '  Critical Temperature : ',F12.6,' +/-',F12.6,/,
     +       '      Critical Density : ',F12.6,' +/-',F12.6,/,
     +       '                   a_0 : ',F12.6,' +/-',F12.6,/,
     +       '                   b_1 : ',F12.6,' +/-',F12.6)

      rmsd_scale = 0.0d0
      rmsd_diameter = 0.0d0
      do I = 1,NS
         tau = (tc - temperature(i))/tc
         sumrho = rhoc + b_1*tau
         rmsd_diameter = rmsd_diameter 
     &        + ( sumrho - 0.5d0*(rhol(i)+rhog(i)) )**2
         difrho = a_0 * tau**beta
         rmsd_scale = rmsd_scale
     &        + ( difrho - 0.5d0*(rhol(i)-rhog(i)) )**2
      enddo
      rmsd_scale = dsqrt(rmsd_scale/dble(ns))
      rmsd_diameter = dsqrt(rmsd_diameter/dble(ns))
      write(6,'(a25,F12.6,/)') '   Scale*Diameter RMSD : ',
     &     rmsd_scale*rmsd_diameter

c     --- start of plotres section
C     === CALCULATE SCALING CURVE AND RECTILINEAR CURVE
      NPOINT=200
      DELT=ABS((TMIN-TC)/NPOINT)
      do I=1,NPOINT-1
         TT=TMIN+DELT*dble(I-1)
         tau = (tc - tt)/tc
         sumrho = rhoc + b_1*tau
         difrho = a_0 * tau**beta
         RHO2 = sumrho+difrho
         RHO1 = sumrho-difrho
         RHO3 = sumrho
         write(25,*)RHO1,tt
         write(26,*)RHO2,tt
         write(27,*)RHO3,tt
      enddo

C     === PLOT TC AND RHOC:
      write( 6,200) RHOC,tc                    
 200  format(2X,4(2X,F12.4))                   
      write(24,*) RHOC,TC                     

c     --- close rectilinear fit output files
      close(24)
      close(25)
      close(26)
      close(27)
      close(30)

      return
      end

      subroutine nonlinearfit(ndata,rhog,rhol,temperature,tmin,fitvars
     &     ,besttemp)
c     ******************************************************************
c     * uses a combination of linear least squares plus a searching    *
c     * algorithm to find the best fit to a more complicated expression*
c     * for the trend in diameter density                              *
c     *                                                                *
c     * if fitvars is 'a_0 b_0 b_1' then                               *
c     * rho_diff = a_0 (1 - T/Tc)^beta                                 *
c     * rho_diam = rhoc + b_0 (1 - T/Tc)^(1-alpha) + b_1 (1 - T/Tc)    *
c     *                                                                *
c     * if fitvars is 'a_0 b_0 b_1' then                               *
c     * rho_diff = a_0 (1 - T/Tc)^beta + a_1 (1 - T/Tc)^(beta + delta) *
c     * rho_diam = rhoc + b_0 (1 - T/Tc)^(1-alpha) + b_1 (1 - T/Tc)    *
c     *                                                                *
c     * originally written 06-10-2005 by M.G. Martin                   *
c     * last modified 01-31-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*30 fitvars
      integer ndata
      double precision tmin
      double precision rhog,rhol,temperature
      dimension rhog(ndata),rhol(ndata),temperature(ndata)
c     --- local variables
      character*30 name
      integer idata,tempsearch,itemp,maxpoint,ivar,nvars
      double precision templow,temphigh,temptest,fitbeta,delta
      double precision besttemp,besterror,testerror,value
      double precision testone,testtwo,rhoaverage,tau,omalpha
      double precision aaa,zzz,uuu,testbeta,bestbeta
      integer maxfitpoint
      parameter (maxfitpoint = 3)
      dimension aaa(maxfitpoint)
      dimension uuu(ndata)
      dimension zzz(maxfitpoint,ndata)
      double precision avalues
      dimension avalues(maxfitpoint)
      logical lscalefit,lscale,lscanbeta
      dimension lscale(ndata)
      double precision bestavalue,bestbvalue
      dimension bestavalue(maxfitpoint),bestbvalue(maxfitpoint)

c     --- open nonlinear fit files
      if ( fitvars .eq. 'a_0 b_0 b_1' ) then
         name = '.nlfitg'
         call twh_openappendfile(40,name)
         name = '.nlfitl'
         call twh_openappendfile(41,name)
         name = '.nlfita'
         call twh_openappendfile(42,name)
         name = '.nlfitc'
         call twh_openappendfile(43,name)
         name = '.nltcvserr'
         call twh_openappendfile(44,name)
      elseif ( fitvars .eq. 'a_0 a_1 b_0 b_1' ) then
         name = '.a0a1b0b1fitg'
         call twh_openappendfile(40,name)
         name = '.a0a1b0b1fitl'
         call twh_openappendfile(41,name)
         name = '.a0a1b0b1fita'
         call twh_openappendfile(42,name)
         name = '.a0a1b0b1fitc'
         call twh_openappendfile(43,name)
         name = '.a0a1b0b1tcvserr'
         call twh_openappendfile(44,name)
      else
         write(6,*) 'Error: unknown fitvars'
         return
      endif
      
      if ( fitvars .eq. 'a_0 b_0 b_1' ) then
         write(6,*) 'Do you wish to fit for the value of beta? (T/F)'
         read(5,*) lscanbeta
      else
c        --- always utilize a single value for beta
         lscanbeta = .false.
      endif

      if ( .not. lscanbeta ) then
         write(6,*) 'Input a value for beta (typical 0.325)'
         read(5,*) fitbeta
         if ( fitbeta .le. 0.0d0 ) then
            write(6,*) 'Error: beta must be positive'
            return
         endif
      endif

      write(6,*) 'Please enter a value for 1-alpha (typical 0.89)'
      read(5,*) omalpha
      if ( fitvars .eq. 'a_0 a_1 b_0 b_1' ) then
         write(6,*) 'Please enter a value for delta (typical 0.5)'
         read(5,*) delta
         if ( delta .le. 0.0d0 ) then
            write(6,*) 'Error: delta must be positive'
            return
         endif
      endif

c     --- scan a temperature space from the highest temperature to 
c     --- 10% higher than the linear fit Tc
      templow = temperature(1)
      do idata = 2,ndata
         if ( temperature(idata) .gt. templow ) then
            templow = temperature(idata)
         endif
      enddo
      templow = templow*1.01d0
      temphigh = 2.0d0*templow
      write(6,*) 'Searching critical temp between ',templow,' and '
     &     ,temphigh
      write(6,*) 'Please enter the number of temperature search points'
      read(5,*) tempsearch
      if ( tempsearch .lt. 1 ) then
         write(6,*) 'you must choose at least 1 search point'
         tempsearch = 1
      endif

      do itemp = 0,tempsearch
c        --- initialize leastsquares variables for the first equation
         do idata = 1,ndata
            uuu(idata) = 0.0d0
            lscale(idata) = .true.
            do ivar = 1,maxfitpoint
               zzz(ivar,idata) = 0.0d0
            enddo
         enddo
c        --- set the critical temperature we are testing
         temptest = templow + dble(itemp)*(temphigh-templow)
     &        /dble(tempsearch)
         if ( lscanbeta ) then
c           --- take the log of both sides of the equation to separate
c           --- beta so it is fitable via least squares
c           --- use the least squares routine
c           --- ln[(rhol - rhog)/2] = ln(avalue) + beta ln[1 - T/T_c]
            do idata = 1,ndata
               uuu(idata) = dlog(0.5d0*(rhol(idata)-rhog(idata)))
c              --- set for the avalue: zzz(1,dp) = 1
               zzz(1,idata) = 1.0d0
c              --- set for testbeta: zzz(2,db) = ln(1 - T/T_c)
               zzz(2,idata) = dlog( 1.0d0 - temperature(idata)
     &              /temptest)
            enddo
c           --- determine the aaa values via leastsquares
            call twh_leastsquares(ndata,2,ndata,maxfitpoint,zzz,aaa
     &           ,uuu,lscalefit,lscale)
c           --- we fit two variables, but only one is an avalue
            nvars = 1
            do ivar = 1,nvars
               avalues(ivar) = dexp(aaa(ivar))
            enddo
c           --- set testbeta
            testbeta = aaa(2)
         else
            if ( fitvars .eq. 'a_0 b_0 b_1' ) then
               nvars = 1
            elseif ( fitvars .eq. 'a_0 a_1 b_0 b_1' ) then
               nvars = 2
            else
               write(6,*) 'Error: unknown fitvars in scaling law'
               return
            endif
c           --- compute data for the first, relatively simple, equation
c           --- fix beta, use linear least squares
            testbeta = fitbeta
c           --- (rhol - rhov)/2 = A (1 - T/T_c)^beta
            do idata = 1,ndata
c              --- store the difference in densities as the target
               uuu(idata) = 0.5d0*(rhol(idata)-rhog(idata))
c              --- set for the avalue: zzz(1,dp) = tau^beta
               tau = 1.0d0 - temperature(idata)/temptest
               zzz(1,idata) = tau**testbeta
               if ( fitvars .eq. 'a_0 a_1 b_0 b_1' ) then
c                  --- zzz(2,dp) = (1-T/Tc)^(beta+delta)
                  zzz(2,idata) = tau**(testbeta + delta)
               endif
            enddo
c           --- determine the aaa values via leastsquares
            call twh_leastsquares(ndata,nvars,MAXDATA,maxfitpoint
     &           ,zzz,aaa,uuu,lscalefit,lscale)
c           --- set avalues
            do ivar = 1,nvars
               avalues(ivar) = aaa(ivar)
            enddo
         endif

c        --- reinitialize leastsquares variables for the second equation
         do idata = 1,ndata
            uuu(idata) = 0.0d0
            lscale(idata) = .true.
            do ivar = 1,maxfitpoint
               zzz(ivar,idata) = 0.0d0
            enddo
         enddo
c        --- use the leastsquares routine to compute the best fit for
c        --- the nonlinear diameter law
         do idata = 1,ndata
c           --- diameter we are trying to fit to
            uuu(idata) = 0.5d0*(rhol(idata)+rhog(idata))
c           --- critical density is just a constant
            zzz(1,idata) = 1.0d0
            tau = 1.0d0 - temperature(idata)/temptest
c           --- (1-alpha) prefactor
            zzz(2,idata) = tau**omalpha
c           --- linear prefactor
            zzz(3,idata) = tau
         enddo
c        --- determine the best values via leastsquares
         call twh_leastsquares(ndata,3,MAXDATA,maxfitpoint,zzz,aaa,uuu
     &        ,lscalefit,lscale)
c        --- compute the error from these fit values
         testone = 0.0d0
         testtwo = 0.0d0
         do idata = 1,ndata
            tau = 1.0d0 - temperature(idata)/temptest
            if ( fitvars .eq. 'a_0 b_0 b_1' ) then
               testone = testone + ( 
     &              avalues(1)*tau**testbeta 
     &              - 0.5d0*(rhol(idata)-rhog(idata)) )**2
            elseif ( fitvars .eq. 'a_0 a_1 b_0 b_1' ) then
               testone = testone + ( 
     &              avalues(1)*tau**testbeta 
     &              + avalues(2)*tau**(testbeta + delta)
     &              - 0.5d0*(rhol(idata)-rhog(idata)) )**2
            else
               write(6,*) 'Error: unknown fitvars in error calc'
               return
            endif
            testtwo = testtwo + ( aaa(1)*zzz(1,idata)
     &           + aaa(2)*zzz(2,idata) + aaa(3)*zzz(3,idata) 
     &           - uuu(idata) )**2
         enddo
         testerror = dsqrt(testone/dble(ndata))
     &        *dsqrt(testtwo/dble(ndata))
         write(44,*) temptest,testerror
         if ( itemp .eq. 0 .or. testerror .lt. besterror ) then
            besterror = testerror
            besttemp = temptest
            do ivar = 1,nvars
               bestavalue(ivar) = avalues(ivar)
            enddo
            do ivar = 1,3
               bestbvalue(ivar) = aaa(ivar)
            enddo
            bestbeta = testbeta
         endif

      enddo

c     ---- output the liquid and gas coexistence densities from the fit
      maxpoint = 200
      do idata = 0,200
         temptest = tmin + (besttemp-tmin)*dble(idata)/dble(maxpoint)
c        --- compute the average from the nonlinear diameter law
         rhoaverage = bestbvalue(1) 
     &        + bestbvalue(2)*(1.0d0 - temptest/besttemp)**omalpha
     &        + bestbvalue(3)*(1.0d0 - temptest/besttemp)
         write(42,*) rhoaverage,temptest
c        --- compute half the difference between the densities using
c        --- the scaling law
         tau = 1.0d0 - temptest/besttemp
         value = bestavalue(1)*tau**bestbeta
         if ( fitvars .eq. 'a_0 a_1 b_0 b_1' ) then
c          --- add in the delta term
            value = value + bestavalue(2)*tau**(bestbeta + delta)
         endif
c        --- vapor is average minus the difference
         write(40,*) rhoaverage-value,temptest
c        --- liquid is average plus the difference
         write(41,*) rhoaverage+value,temptest
      enddo
c     --- output the critical point in plottable form
      write(43,*) bestbvalue(1),besttemp
c     --- output all of the fitted values
      write(6,*) 'Results from the nonlinear fit to the data'
      write(6,*) '   Critical Temperature: ',besttemp
      write(6,*) '      Critical Density : ',bestbvalue(1)
      write(6,*) '                   a_0 : ',bestavalue(1)
      if ( fitvars .eq. 'a_0 a_1 b_0 b_1' ) then
         write(6,*) '                   a_1 : ',bestavalue(2)
      endif
      write(6,*) '                   b_0 : ',bestbvalue(2)
      write(6,*) '                   b_1 : ',bestbvalue(3)
      if ( lscanbeta ) then
         write(6,*) '      Scaling Law beta : ',bestbeta
      endif
      write(6,*) '   Scale*Diameter RMSD : ',besterror

c     --- close nonlinear output files
      close(40)
      close(41)
      close(42)
      close(43)
      close(44)

      return
      end
