      real function planckg_wave(iwave,v,T)
C     ****************************************************************
C     Utility routine to calculate the dB/dT of Planck function.
C
C      iwave indicates wavenumbers(0) or wavelength(1)
C      T is the temperature in K
C      v is the wavenumber (in cm-1) or wavelength(um)
C
C      Calculated function is in units of W cm-2 sr-1 cm
C      Calculated function is in units of W cm-2 sr-1 um-1 for iwave=1
C
C     Pat Irwin   2003		Original(!)
C     Pat irwin   25/8/04       Updated to deal with either wavenumber of
C                               wavelength
C
C     ****************************************************************
      implicit none
      real v,T,y
      integer iwave
      double precision c1,c2,bot,tmp,top,A
      parameter(c1=1.1911e-12,c2=1.439)

C
      if(iwave.eq.0)then
       y = v
       A = c1*c2*(y**4)/T**2
      else
       y=1e4/v  
       A = c1*c2*(y**6)*1e-4/T**2
      endif

      tmp = c2*y/T
      bot = (dexp(tmp) - 1.0)**2
      top = dexp(tmp)*A

      planckg_wave = sngl(top/bot)

      return

      end   
