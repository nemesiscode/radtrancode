      subroutine ml_lacis(knu0in,yvin,ELin,SFBin,qrotin,Pin,
     1 Tin,qin,SL,BL)
C     $Id: ml_lacis.f,v 1.3 2011-06-17 15:40:27 irwin Exp $
C     ***********************************************************************
C     
C     Convert EKS Malkmus-Lorentz Parameters to those used by Lacis and
C     Oinas (1991)
C
C     Pat Irwin	3/6/94
C
C     ***********************************************************************
      implicit double precision (a-h,o-z)
      integer i
      real knu0in,yvin,ELin,SFBin,qrotin,Pin,Tin,qin
      double precision knu0,yv,EL,SFB,qrot,P,T,q,pi,C1,T0,SL,BL
      double precision knu
      parameter (T0=296.0, pi=3.1415927, C1=1.439)

      knu0 = dble(knu0in)
      yv = dble(yvin)
      EL = dble(ELin)
      SFB = dble(SFBin)
      qrot = dble(qrotin)
      P = dble(Pin)
      T = dble(Tin)
      q = dble(qin)

      IF(T.LT.140.0)T=140.0 ! Code becomes unstable for small T

C      print*,'ml-lacis',knu0,yv,EL,SFB,qrot,P,T,q

      knu=(knu0*(T0/T)**qrot)*dexp(C1*EL*(1/T0 - 1/T))

C      print*,'ml-lacis',knu

      A = knu
      B = pi*A*yv*P*(q+(1-q)/SFB)*dsqrt(T0/T)

C      print*,A,B
      SL = A
      BL = 4.*B/(pi*A)

      return
      end
