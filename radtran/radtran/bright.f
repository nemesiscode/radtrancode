      real function bright(v,R)
C     $Id: bright.f,v 1.2 2011-06-17 15:40:25 irwin Exp $
C     ****************************************************************
C     Utility routine to calculate the brightness temperature
C      R is the radiance in units of W cm-2 sr-1 cm
C      v is the wavenumber in cm-1
C
C      bright is the brightness temperature in K
C
C     Pat Irwin   20/1/94
C
C     ****************************************************************
      implicit none
      real c1,c2,v,T,a,b,R
      parameter(c1=1.1911e-12,c2=1.439)

C
      a = c1*v*v*v/R
      bright = c2*v/log(a+1.0)

      return

      end   
