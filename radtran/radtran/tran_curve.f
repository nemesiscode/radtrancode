      subroutine tran_curve(fknu,delad,y1,T,lcalch,k_g,del_g,ng,
     1ndata,x,y,sig)
C     ************************************************************************
C     Subroutine to calculate a transmission curve with transmission varying
C     from 0.01 to 0.99 and absorber amounts equally spaced in log space
C
C     Input variables
C	fknu		real 	knu at temperature of layer
C       delad		real 	mean line spacing/(ad0*sqrt(T)
C	y1		real 	aL/aD at pressure and temperature of layer
C 	lcalch		integer	Switch to determine how transmission function
C				is calculated.
C				1  = Goody-Voigt
C				2  = Malkmus-Lorentz
C                               3  = k-distribution
C	ndata		integer	Number of points in transmission curve
C	k_g(maxg)		real	k-distribution (if required)
C	del_g(maxg)	real	g-ordinates (if required)
C	ng		integer	number of g-ordinates
C
C     Output variables
C
C	x(ndata)	real	Absorber amounts
C	y(ndata)	real	Transmission
C	sig(ndata)	real	Transmission weights
C
C     Pat Irwin		1/3/95
C
C     ***********************************************************************
      include '../includes/arrdef.f'
      real pi
      parameter (pi=3.1415927)
      real fknu,delad,y1,T,tau_goody_voigt2
      integer lcalch
      double precision SL,BL
      integer ndata
      real x(ndata),y(ndata),sig(ndata)

      real k_g(maxg),del_g(maxg)
      integer ng

C     ***********************************************************************
C     CODE
C     ***********************************************************************

C     Limits are set so that ytau(U1) = 0.99 and ytau(U2)=0.01
C     Firstly find where this is the case by linear interpolation
C     in log10(U) space


C      print*,'fknu,delad,y1,T,lcalch,k_g,del_g,ng,
C     1ndata,x,y,sig'

C      print*,fknu,delad,y1,T,lcalch,k_g,del_g,ng,
C     1ndata,x,y,sig



C     Set approximate limits for U
      U1=-4.0-log10(fknu)
      U2=U1+10.0+log10(delad)

C     Now find where ytau(U)=0.99 and 0.01
      imin=2
      imax=ndata
      lmin=0
      lmax=0

      yv=0.25*y1/delad		!Band parameter used by Malkmus-Lorentz model
      call ml_lac1(fknu,yv,SL,BL)

      do 15 i=1,ndata
        UL=U1 + (U2-U1)*real(i-1)/real(ndata-1)
        U=10.0**UL
        x(i)=U
        y(i)=0.
        if(lcalch.eq.1)then
C         Goody-Voigt
          y(i)=exp(-tau_goody_voigt2(fknu,delad,y1,T,U))
        else if (lcalch.eq.2)then
C         Malkmus-Lorentz
          y(i)=sngl(dexp(-0.5*PI*BL*(DSQRT(1. + 4.*SL*U/(PI*BL)) - 1.)))
        else if(lcalch.eq.3)then
C         K-distribution
          y(i)=tran_ck(k_g,del_g,ng,U)
        end if

C        print*,i,x(i),y(i)

        if(y(i).lt.0.99.and.lmin.eq.0)then
          imin=i
          lmin=1
        endif
        if(y(i).lt.0.01.and.lmax.eq.0)then
          imax=i
          lmax=1
        endif    
15    continue

     

      dxmin=log10(x(imin))-log10(x(imin-1))
      dymin=y(imin)-y(imin-1)
      dxmax=log10(x(imax))-log10(x(imax-1))
      dymax=y(imax)-y(imax-1)


C      print*,imin,imax
C      print*,dxmin,dxmax
C      print*,dymin,dymax

C     Now set limits in absorber space using simple linear interpolation 
      U1=log10(x(imin-1))+dxmin*(0.99-y(imin-1))/dymin
      U2=log10(x(imax-1))+dxmax*(0.01-y(imax-1))/dymax

C      print*,U1,U2

      do 25 i=1,ndata
        UL=U1 + (U2-U1)*real(i-1)/real(ndata-1)
        U=10.0**UL
        x(i)=U
        y(i)=0.
        if(lcalch.eq.1)then
C         Goody-Voigt
          y(i)=exp(-tau_goody_voigt2(fknu,delad,y1,T,U))
C          print*,x(i),y(i)
        else if (lcalch.eq.2)then
C         Malkmus-Lorentz
          y(i)=sngl(dexp(-0.5*PI*BL*(DSQRT(1. + 4.*SL*U/(PI*BL)) - 1.)))
        else if(lcalch.eq.3)then
C         K-distribution
          y(i)=tran_ck(k_g,del_g,ng,U)
        end if
        sig(i)=1e-2
25    continue

      return

      end

      real function tran_ck(k_g,del_g,ng,U)
C     *******************************************************************
C     Function to calculate the transmission of a path given the absorber
C     amount U and the K-coefficients and weights.
C
C     Pat Irwin		1/3/95
C
C     *******************************************************************
      include '../includes/arrdef.f'
      real k_g(maxg),del_g(maxg),U
      integer ng


      sum=0.
      do 10 i=1,ng
       sum=sum+exp(-k_g(i)*U)*del_g(i)
10    continue

      tran_ck=sum
      return
      end

