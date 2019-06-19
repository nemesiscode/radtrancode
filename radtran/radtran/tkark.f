      real function tkark(kap100,kap198,kap296,dline,P,T,U,q)
C     $Id:
C     *********************************************************************
C
C     Calculate the transmission of a path using the the formulation of the
C     Karkoschka and Tomasko (2009)
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
C	kap100	real	Absorption coefficient at 100 ((km-am)**-1)
C	kap198	real	Absorption coefficient at 198 ((km-am)**-1)
C	kap296	real	Absorption coefficient at 296 ((km-am)**-1)
C	dline	real	line spacing / alpha_D0 (K**-0.5)
C	P	real	Mean pressure of path (atm)
C	T	real	Mean temperature of path (K)
C	U	real	Absorber amount in path (*1e20 (molec/cm2))
C	q	real	Fractional abundance of active gas 
C
C     Output variable
C
C 	tkark	real	Transmission of path
C     ------------------------------------------------------------------------
C
C	Pat Irwin	26/9/94
C
C     *********************************************************************
      implicit none
      real knu,U,dline,y,T,kap100,kap198,kap296,T0,y0,SFB,P,q,P0,qrot,a
      integer i,idim,ndim,ico,id,iso,imod
      parameter (P0 = 1.,T0=296.0, idim=500)
      real SQT0,kmamagat
      parameter (SQT0 = 17.20465053)          
      parameter (kmamagat = 2.687e4)   
      real x,dx,s(idim),A1,A2,humlic,simp_int,Aconst,lx,Cconst,Vconst
      real Bconst,Dconst,Ratio,Tconst
      real z,tmp
      character*1 ans
 
C      print*,'kark'
C      print*,kap100,kap198,kap296,dline
C      print*,P,T,U,q

C      read(5,1)ans
1     format(a)

      if(kap296.eq.0.or.U.eq.0.)then
       tkark=0.
 
      else

C      interpolate k-coeffiecient to path temperature.
       z = (T-198.0)/98.0
       tmp = 0.5*z*(z-1.0)*alog(kap100) + (1.0-z**2)*alog(kap198)+
     &  0.5*z*(z+1.0)*alog(kap296)
       knu = exp(tmp)
C       print*,'knu = ',knu

C      calculate y=alpha_L/alpha_D at path conditions
       y0 = 125.0
       SFB = 1.4
       y=y0*(P/P0)*(SQT0/T)*(q + (1-q)/SFB)
      
C       print*,'y = ',y

       Dconst = humlic(0.,y)
C       Aconst = knu*U*delad/sqrt(T)
       Aconst = knu*(U/kmamagat)*dline/sqrt(T)
       Bconst = Dconst/(1. + Aconst*Dconst)
       Cconst = 0.01*Bconst      

       Vconst = Cconst / (1. - Cconst*Aconst)

       Ratio = sqrt(0.01*Dconst/Vconst)
      
       lx = 30.*y*Ratio


       if(lx.eq.0)then
        print*,'tkark: Error: lx = 0.'
        print*,'kap296,dline,SFB,P,T,U,q',kap296,dline,
     1		SFB,P,T,U,q
        print*,'setting tau to zero'
        tkark = 0.0
        return
       endif

       ico = 0
777    Tconst = humlic(lx,y)/(1. + Aconst*humlic(lx,y))
       if(Tconst.gt.Cconst)then
        lx = lx*1.5
        ico = ico+1
C        if (ico.gt.10)then
C         print*,'tkark : ico = ',ico
C        endif
        goto 777
       endif


       do 10 i=1,101
        x=lx*(i-1)/100.
        s(i)=humlic(x,y)/(1. + Aconst*humlic(x,y))
10     continue

       ndim=101
       dx = lx/100.

       A1=simp_int(s,idim,ndim,dx)
 
       do 20 i=1,101
        x=lx + (i-1)*lx
        s(i)=humlic(x,y)/(1. + Aconst*humlic(x,y))
20     continue      

       ndim=101
       dx = lx

       A2=simp_int(s,idim,ndim,dx)

C       print*,'output',knu,U,U/kmamagat,A1,A2
       tkark=2.*knu*(U/kmamagat)*(A1+A2)
C       print*,'tkark = ',tkark  

      end if

      return
      end


