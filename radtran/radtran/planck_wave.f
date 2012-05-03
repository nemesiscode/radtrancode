      real function planck_wave(iwave,v,T)
C     ****************************************************************
C     Utility routine to calculate the Planck function.
C
C      iwave indicates wavenumbers(0) or wavelength(1)
C      T is the temperature in K
C      v is the wavenumber (in cm-1) or wavelength(um)
C
C      Calculated function is in units of W cm-2 sr-1 cm for iwave=0
C      Calculated function is in units of W cm-2 sr-1 um-1 for iwave=1
C
C     Pat Irwin   20/1/94	Original
C     Pat irwin	  25/8/04	Updated to deal with either wavenumber of
C				wavelength
C
C     ****************************************************************
      implicit none
      real v,T,y
      integer iwave
      double precision c1,c2,a,b,tmp
      parameter(c1=1.1911e-12,c2=1.439)

 
      if(iwave.eq.0)then
       y = v
       a = c1*(y**3)
      else
       y=1e4/v
       a = c1*(y**5)/1e4
      endif

      if(y.eq.0.0) then
        planck_wave = 0.0
      else
        tmp = c2*y/T
        b = dexp(tmp) - 1.
        planck_wave = sngl(a/b)
      endif

      return

      end   
