      SUBROUTINE CALCTAUGRAD(NPRO,NCONT,H,P,T,DUST,MOLWT,XSEC,IRAY,VV,
     1 K_G,HEIGHT,DTAUDS,DTAUDC,DTAUR)
C     **************************************************************
C     Routine to estimate the rate of change of optical depth per km
C     at a point in the atmosphere
C
C     Input parameters
C  	NPRO	INTEGER	Number of vertical levels in profile
C	NCONT	INTEGER	Number of aerosol types
C	H(NPRO)	REAL	Profile heights
C	P(NPRO)	REAL	Profile pressures
C	T(NPRO)	REAL	Profile temperatures
C	DUST(MAXPRO/MAXCON) REAL Dust abundances
C	MOLWT	REAL	Atmospheric molecular weight
C	XSEC(NCONT) REAL	Dust particle x-sections
C	IRAY	INTEGER	Switch for Rayleigh scattering
C	VV	REAL	Current wavenumber
C	K_G(NPRO) REAL	Gas Absorption coefficient at each level
C	HEIGHT	REAL	Height in atmosphere for calculation
C
C     Output variables
C	DTAUDS	REAL	Total optical depth per km
C	DTAUDC	REAL	Total optical depth of scatterers per km
C	DTAUR	REAL	Optical depth of Rayleigh scatteing per km
C
C     Pat Irwin	22/3/13
C
C     ************************************************************** 	
      IMPLICIT NONE
      include '../includes/arrdef.f'
      include '../includes/constdef.f'
      INTEGER NPRO,IFL,K,NCONT,IRAY,J
      REAL DUST(MAXPRO,MAXCON),XSEC(MAXCON),MOLWT
      REAL HEIGHT,H(NPRO),PNOW,P(NPRO),TNOW,T(NPRO),F
      REAL K_G(NPRO),K1,RAYLEIGHJ,VV,C,TOTAM,DKDS,DNOW
      REAL DKDC,DTAUR,DTAUDC(MAXCON),DTAUDS
      REAL RAYLEIGHV,RAYLEIGHA

C     Calculate position in height array and determine local
C       temperature and pressure 
      CALL INTERP_PT(NPRO,H,P,T,HEIGHT,PNOW,TNOW,F,IFL)

      K1 = (1-F)*K_G(IFL) + F*K_G(IFL+1)
      TOTAM = MODBOLTZ*(PNOW/TNOW)
C     TOTAM is molecules/cm2 for a 1km-long path

      DKDS = TOTAM*K1*1E-20	! optical depth of gas/km

      DKDC = 0.0		! optical depth of aerosol/km

C     DUST() is in units of particles/gram. Need to multiply by g/cm2 
C     of a 1km-long path to get optical depth/km.
      DO J=1,NCONT
        DNOW = (1-F)*DUST(IFL,J) + F*DUST(IFL+1,J)
        DTAUDC(J)=XSEC(J)*DNOW*TOTAM*MOLWT/AVOGAD
        DKDC = DKDC+DTAUDC(J)
      ENDDO

      DTAUR = 0.0
      IF(IRAY.GT.0)THEN
C       print*,IRAY,VV,PNOW,TNOW,TOTAM
       IF(IRAY.EQ.1)THEN
        DTAUR = TOTAM*RAYLEIGHJ(VV,PNOW,TNOW)
       ELSEIF(IRAY.EQ.2)THEN
        DTAUR = TOTAM*RAYLEIGHV(VV,PNOW,TNOW)
       ELSE
        DTAUR = TOTAM*RAYLEIGHA(VV,PNOW,TNOW)
       ENDIF
C       print*,TOTAM,RAYLEIGHJ(VV,PNOW,TNOW),
C     1  RAYLEIGHV(VV,PNOW,TNOW),
C     2  RAYLEIGHA(VV,PNOW,TNOW)

C       print*,TOTAM*RAYLEIGHJ(VV,PNOW,TNOW),
C     1  TOTAM*RAYLEIGHV(VV,PNOW,TNOW),
C     2  TOTAM*RAYLEIGHA(VV,PNOW,TNOW)
   
C       print*,DTAUR

C       stop
       
      ENDIF

C     Calculate total optical depth/km at current conditions
      DTAUDS = DKDS + DKDC + DTAUR


      RETURN

      END

