      real function planck(v,T)
C     ****************************************************************
C     Utility routine to calculate the Planck function.
C      T is the temperature in K
C      v is the wavenumber in cm-1
C
C      Calculated function is in units of W cm-2 sr-1 cm
C
C     Pat Irwin   20/1/94
C
C     ****************************************************************
      implicit none
      real v,T
      double precision c1,c2,a,b,tmp
      parameter(c1=1.1911e-12,c2=1.439)

      if(v.eq.0.0) then
       planck = 0.0
      else
       a = c1*v*v*v
       tmp = c2*v/T
       b = dexp(tmp) - 1.
       planck = sngl(a/b)
      endif

      return

      end   
