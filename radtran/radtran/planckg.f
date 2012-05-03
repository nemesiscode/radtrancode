      real function planckg(v,T)
C     ****************************************************************
C     Utility routine to calculate the dB/dT of Planck function.
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
      double precision c1,c2,bot,tmp,top
      parameter(c1=1.1911e-12,c2=1.439)

C
      tmp = c2*v/T
      bot = (dexp(tmp) - 1.0)**2
      top = dexp(tmp)*c1*c2*(v**4)/(T**2)
      planckg = sngl(top/bot)

      return

      end   
