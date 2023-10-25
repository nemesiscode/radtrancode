C-----------------------------------------------------------------------
      real function invplanck_wave(iwave,v,B)
C-----------------------------------------------------------------------
C     Function to calculate the inverse Planck function.
C
C     iwave indicates wavenumbers(0) or wavelength(1)
C      
C     B is the plank function at T and v
C         units of W cm-2 sr-1 cm for iwave=0
C         units of W cm-2 sr-1 um-1 for iwave=1
C
C     v is the wavenumber (in cm-1) or wavelength(um)
C
C	Output: Temperature in K
C
C	NB. based on inverting planck_wave.f
C
C-----------------------------------------------------------------------
C     20/10/23	Nick Teanby		Original
C-----------------------------------------------------------------------
      implicit none
      real v,B
      integer iwave
      double precision c1,c2,aa,bb,t
      parameter(c1=1.1911e-12,c2=1.439)

      if (B.eq.0.0) then
         t = 0.0
      else if (v.eq.0.0) then
         t = 0.0
      else if(iwave.eq.0)then
         aa = c2*dble(v)
         bb = dlog( (c1*dble(v)**3)/dble(B) + 1.0 ) 
         t = aa/bb
      else
         aa = c2 * 1.0e4 / dble(v)
         bb = dlog( (c1*1.0e16)/(dble(B)*dble(v)**5) + 1.0 )
         t = aa/bb
      endif
      invplanck_wave = sngl(t)

      return

      end   
