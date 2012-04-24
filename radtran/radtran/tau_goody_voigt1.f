      real function tau_goody_voigt1(knu0,delad,y0,EL,SFB,qrot,P,T,U,q)
C     $Id: tau_goody_voigt1.f,v 1.3 2011-06-17 15:40:27 irwin Exp $
C     *********************************************************************
C
C     Version 1.1
C     Calculate the transmission of a path using the the formulation of the
C     Goody-Voigt random band model given in Strong (1992) D.Phil.Thesis
C
C     Integration of the function over x is performed in two stages using 
C     Simpsons rule with 101 points over the ranges:
C	1)	x = 0 to 30*y
C 	2)	x = 30*y to 3030*y
C     This is the same method used in tran_eks.f by K.Strong.
C
C     ----------------------------------------------------------------------
C     
C     Version 1.2
C     Version 1.1 gives spurious results if the coefficient of V(x,y) in the
C     denominator of the integral is large. We really want to set the first
C     limit of x to be where the integrand is 1% of the integrand when x=0.
C
C     The integrand may be written:
C	B(x,y) =         V(x,y)
C 		    ---------------
C	  	      1 + A*V(x,y)
C
C     Put C = 0.01*Bmax = 0.01*B(0,y). When the function B is 1% of max
C     then:
C		V1 = V(x,y) = C/(1-AC)
C
C     In the wings, the Voigt line shape looks like a Lorentz line. Thus:
C  	 	(nu - nu0) is proportional to 1/sqrt(k)
C     Now when x=30y, V(x,y) is always less than 1% of Vmax = V(0,y). Thus
C     Putting D = V(0,y), when the factor A is large, the integration limit
C     must be scaled to:
C		x_int = 30*y*sqrt(0.01*D/V1) 
C
C     and the integration is then more appropriate. Further checks to ensure
C     numerical accuracy are in progress.
C 
C     ------------------------------------------------------------------------
C     Input Variables
C
C	knu0	real	Absorption coefficient at 296K (*1e20 (molec/cm2)**-1)
C	delad	real	line spacing / alpha_D0 (K**-0.5)
C	y0	real	alpha_L0 / alpha_D0 (K**-0.5)
C	EL	real	Lower state energy (cm-1)
C	SFB	real	Self-to-foreign broadening parameter
C	qrot	real	Rotational Partition Function Coefficient
C	P	real	Mean pressure of path (atm)
C	T	real	Mean temperature of path (K)
C	U	real	Absorber amount in path (*1e20 (molec/cm2))
C	q	real	Fractional abundance of active gas 
C
C     Output variable
C
C 	tau_goody_voigt1	real	Transmission of path
C     ------------------------------------------------------------------------
C
C	Pat Irwin	26/9/94
C
C     *********************************************************************
      implicit none
      real knu,U,delad,y,T,EL,knu0,C1,T0,y0,SFB,P,q,P0,qrot
      integer i,idim,ndim,ico
      parameter (P0 = 1.,T0=296.0, C1=1.439,idim=500)
      real SQT0
      parameter (SQT0 = 17.20465053)		! Sqrt(296.0)
      real x,dx,s(idim),A1,A2,humlic,simp_int,Aconst,lx,Cconst,Vconst
      real Bconst,Dconst,Ratio,Tconst

      if(knu0.eq.0.or.U.eq.0.)then
       tau_goody_voigt1=0.
      else

C     calculate y=alpha_L/alpha_D at path conditions

      y=y0*(P/P0)*(SQT0/T)*(q + (1-q)/SFB)

C     Calculate absorption coefficient at path temperature
      knu=(knu0*(T0/T)**qrot)*exp(C1*EL*(1/T0 - 1/T))

      if(knu.eq.0.0)then
       tau_goody_voigt1=0.
       return
      endif
      
      Dconst = humlic(0.,y)
      Aconst = knu*U*delad/sqrt(T)
      Bconst = Dconst/(1. + Aconst*Dconst)
      Cconst = 0.01*Bconst      

      Vconst = Cconst / (1. - Cconst*Aconst)

      Ratio = sqrt(0.01*Dconst/Vconst)
      
      lx = 30.*y*Ratio


      if(lx.eq.0)then
       print*,'tau_goody_voigt1: Error: lx = 0.'
       print*,'knu0,delad,y0,EL,SFB,qrot,P,T,U,q',knu0,delad,
     1		y0,EL,SFB,qrot,P,T,U,q
       print*,'setting tau to zero'
       tau_goody_voigt1 = 0.0
       return
      endif

      ico = 0
777   Tconst = humlic(lx,y)/(1. + Aconst*humlic(lx,y))
      if(Tconst.gt.Cconst)then
       lx = lx*1.5
       ico = ico+1
C       if (ico.gt.10)then
C         print*,'tau_goody_voigt1 : ico = ',ico
C       endif
       goto 777
      endif


      do 10 i=1,101
       x=lx*(i-1)/100.
       s(i)=humlic(x,y)/(1. + Aconst*humlic(x,y))
10    continue

      ndim=101
      dx = lx/100.

      A1=simp_int(s,idim,ndim,dx)
 
      do 20 i=1,101
       x=lx + (i-1)*lx
       s(i)=humlic(x,y)/(1. + Aconst*humlic(x,y))
20    continue      

      ndim=101
      dx = lx

      A2=simp_int(s,idim,ndim,dx)

      tau_goody_voigt1=2.*knu*U*(A1+A2)

      end if

      return
      end


