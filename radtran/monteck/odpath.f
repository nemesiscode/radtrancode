      SUBROUTINE ODPATH(VV,IRAY,PVEC,AVEC,NPRO,NGAS,NCONT,MOLWT,
     1 RADIUS,P,T,H,DUST,TABK,IGDIST,XSEC)
C     $Id:
C     ***************************************************************
C     Subroutine to calculate equivalent CG path between two points
C     in a  spherically symmetric atmosphere. Code is based on the
C     RADTRAN layer.f routine and uses Simpson's rule integration.
C
C     Input variables:
C	VV		REAL	Calculated wavenumber (require to estimate
C				Rayleigh scattering)
C	IRAY		INTEGER 0=Rayleigh scatt off, 1 = on
C	PVEC(3)		REAL	Starting position. Convention is that
C				z-azis is the zenith at the tangent point
C				of a limb path.
C	AVEC(3)		REAL	Vector of new path.
C	NPRO		INTEGER	Number of points in atm profile.
C	NGAS		INTEGER	Number of gases
C	NCONT		INTEGER	Number of dust types 
C	MOLWT		REAL	Molecular weight of atmosphere
C	RADIUS		REAL	Radius of planet
C	P(MAXPRO)		REAL	Pressure (prf)
C	T(MAXPRO)		REAL	Temperature (prf)
C	H(MAXPRO)		REAL	Heights (prf)
C	DUST(MAXPRO,MAXCON)REAL	Dust profiles
C	TABK(MAXG,MAXPRO)	REAL	Overlapped gas k-distributions at each 
C				level in the atmosphere
C	IGDIST		INTEGER	Required g-ordinate for calculation
C	XSEC(MAXCON)	REAL	Aerosol x-sections
C
C     Output variables:
C	TAUTOT		REAL	Total optical thickness of path
C
C     Pat Irwin		1/4/99	Original
C     Pat Irwin		22/3/12	Massively revised
C 
C     ***************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'

      REAL PVEC(3),AVEC(3),TAUREC,TAUTOT,VV,TVEC(3)
      INTEGER NPRO,NGAS,J,NCONT,I,K,IFL,IRAY,IDIST,MINT
      PARAMETER (MINT=101)
      REAL MOLWT,F,DELS,W(MINT),S
      REAL TNOW,PNOW,CALCALT,RADIUS,HEIGHT
      REAL P(MAXPRO),T(MAXPRO),H(MAXPRO)
      REAL TABK(MAXG,MAXPRO),K_G(MAXPRO)
      REAL T1,T2,FSCAT,DTAU,DTAUSC,TAUREQ
      REAL DTAUDS,DTAUR,DTAUDC,XSEC(MAXCON)
      INTEGER IGDIST
      REAL DUST(MAXPRO,MAXCON)
C     ------------------------------------------------------------------
C     Check some numbers
C     ------------------------------------------------------------------
      IF(NCONT.GT.MAXCON)THEN
       PRINT*,'ODPATH: NCONT > MAXCON'
       PRINT*,NCONT,MAXCON
       STOP
      ENDIF

      IF(NPRO.GT.MAXPRO)THEN
       PRINT*,'ODPATH: NPRO > MAXPRO'
       PRINT*,NPRO,MAXPRO
       STOP
      ENDIF


C     Load up relevant k-ordinate
C      print*,'igdist',igdist
      DO I=1,NPRO
       K_G(I)=TABK(IGDIST,I)
      ENDDO

C     Initialise integration array
      DO 20 I=1,MINT
       W(I)=2.
       IF(I.EQ.2*(I/2))W(I)=4.
20    CONTINUE
      W(1)=1.
      W(MINT)=1.

C     Compute path length
      S = 0.0
      DO 10 I=1,3
       S = S + AVEC(I)**2
10    CONTINUE
      S = SQRT(S) 
      DELS=S/FLOAT(MINT-1)

      TAUTOT=0.
  
      DO 30 I=1,MINT

       DO 25 J=1,3
        TVEC(J) = PVEC(J) + FLOAT(I-1)*AVEC(J)/FLOAT(MINT-1)
25     CONTINUE

       HEIGHT = CALCALT(TVEC,RADIUS)
C      Find place in height array and local pressure, temperature
       CALL INTERP_PT(NPRO,H,P,T,HEIGHT,PNOW,TNOW,F,IFL)

C      Determine optical depths/km
       CALL CALCTAUGRAD(NPRO,NCONT,H,P,T,DUST,MOLWT,XSEC,IRAY,VV,
     1 K_G,HEIGHT,DTAUDS,DTAUDC,DTAUR)

C      Add optical depth element to integration array
       TAUTOT = TAUTOT + DTAUDS*W(I)

30    CONTINUE

      TAUTOT = TAUTOT*DELS/3.

      RETURN

      END
