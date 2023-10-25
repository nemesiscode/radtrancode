C-----------------------------------------------------------------------
      subroutine calc_nlte_t(T,p,v,iwave,inlte_flag,Tnlte)
C-----------------------------------------------------------------------
C
C	Apply some basic non-LTE functions to the emission temperatures.
C	If the function is to be applied to the Planck function instead
C	then the invplanck_wave routine can be used to implement this.
C
C	IN:
C	  T		 real		T(K)
C	  p		 real		Pressure(atm)
C	  v		 real		wavenumber(cm-1) or wavelength(microns)
C	  iwave	 int		=0 v is in wavenumbers
C	           			=1 v is in microns
C	  inlte_flag int		flag to specify nLTE expression to use
C					=0 assume perfect LTE (ie do nothing)
C					=1 Orton form as in Radtran manual P27 approx for Saturn/Titan
C					=2 Manuel Lopez Puertas CH3 calcs for k1=1e-11 for Titan
C					=3 Manuel Lopez Puertas CH3 calcs for k1=1e-12 for Titan
C
C	OUT:
C	  Tnlte	real		equivalent non-LTE emission temperature (K)
C
C	REQUIRES:
C	  function planck_wave(iwave,v,T)		T->B
C	  function invplanck_wave(iwave,v,B)	B->T
C
C	CAUTION:
C	  This is a pretty basic way to implement non-LTE. 
C	  These functions are applied to _all_ sources of opacity 
C	  (ie all gases + continuum CIA + etc). This is probably OK for
C	  a distinct emission line superimposed on an LTE continuum if the
C	  continuum is generated at higher pressures where these functions
C	  should tend to LTE anyway.
C
C	  Enture that T and p are correct. In cirsradg_wave the press and
C	  emtemp arrays are in opposite senses and layinc must be used. See
C	  cirsradg_wave and you'll see what I mean.
C
C-----------------------------------------------------------------------
C     20/10/23	Nick Teanby		Original
C-----------------------------------------------------------------------
      implicit none
      real T,p,v,Tnlte
      integer iwave,inlte_flag
      real f,B
      real planck_wave, invplanck_wave
      external planck_wave, invplanck_wave
      
      if (inlte_flag.eq.0) then
c     ** leave T unchanged, ie assume perfect LTE **        
         Tnlte = T
      else if (inlte_flag.eq.1) then
c     ** Orton formula from Radtran manual P27, estimated for saturn and titan **
c     planck function is scaled by F=P(atm)*1e6/(1 + p(atm)*1e6) 
C	ie. some smooth function that =0.5 at 1ubar
	   f = p*1.0e6 / (1.0 + p*1.0e6)
	   B = planck_wave(iwave,v,T)
	   Tnlte = invplanck_wave(iwave,v,B*f)
      else if (inlte_flag.eq.2) then
c     ** Approx to Manuel Lopez Puertas CH3 calcs for k1=1e-11
	   f = p*1.3e8 / (1.0 + p*1.3e8)
	   B = planck_wave(iwave,v,T)
	   Tnlte = invplanck_wave(iwave,v,B*f)
      else if (inlte_flag.eq.3) then
c     ** Approx to Manuel Lopez Puertas CH3 calcs for k1=1e-12
	   f = p*1.3e7 / (1.0 + p*1.3e7)
	   B = planck_wave(iwave,v,T)
	   Tnlte = invplanck_wave(iwave,v,B*f)
      else
c     ** default to doing nothing (always a good plan, lols!) **
         Tnlte = T
      endif
      
      return
      
      end
      