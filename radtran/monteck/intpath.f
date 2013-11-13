      SUBROUTINE INTPATH(VV,IRAY,TAUREQ,PVEC,DVEC,NPRO,NGAS,NCONT,MOLWT,
     1 RADIUS,P,T,H,DUST,TABK,IGDIST,XSEC,XOMEGA,TEND,FSCAT,TAUSCAT)
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
C	TAUREQ		REAL	Required opacity
C	PVEC(3)		REAL	Starting position. Convention is that
C				z-azis is the zenith at the point where
C				the photon enters the atmosphere
C	DVEC(3)		REAL	Direction vector of new path.
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
C	XOMEGA(MAXCON)	REAL	Particle single scattering albedos
C
C     Output variables:
C	TEND		REAL	Temperature at end of path
C	PVEC(3)		REAL	Final position vector
C	FSCAT		REAL	Probability of scattering
C	TAUSCAT(MAXCON)	REAL	Scattering opacity of path for each
C				of NCONT particles. If IRAY=1, then 
C				TAU(NCONT+1) contains the rayleigh
C				scattering opacity.
C
C     Pat Irwin		1/4/99	Original
C     Pat Irwin		22/3/12	Massively revised
C 
C     ***************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'

      REAL PVEC(3),DVEC(3),TAUREC,VV,TVEC(3),PVEC1(3),AVEC(3)
      INTEGER NPRO,NGAS,J,NCONT,I,K,IFL,IRAY,IDIST,MINT
      PARAMETER (MINT=201)
      REAL MOLWT,F,DELS
      REAL TNOW,PNOW,CALCALT,RADIUS,HEIGHT
      REAL P(MAXPRO),T(MAXPRO),H(MAXPRO),TAUTOT(MINT)
      REAL TAUTOTC(MINT)
      REAL TABK(MAXG,MAXPRO),K_G(MAXPRO),TAUNOW(MINT),TAUC(MINT)
      REAL T1,T2,DTAU,DTAUSC,TAUREQ,DIST(MINT),PEND
      REAL DTAUDS,DTAUR,XSEC(MAXCON),XOMEGA(MAXCON)
      REAL TEND,FSCAT,TAUSCAT(MAXCON),RADCLOSE,S,DTAUDC(MAXCON)
      INTEGER IGDIST
      REAL DUST(MAXPRO,MAXCON),RADIUS1
C     ------------------------------------------------------------------
C     Check some numbers
C     ------------------------------------------------------------------
      IF(NCONT.GT.MAXCON)THEN
       PRINT*,'INTPATH: NCONT > MAXCON'
       PRINT*,NCONT,MAXCON
       STOP
      ENDIF

      IF(NPRO.GT.MAXPRO)THEN
       PRINT*,'INTPATH: NPRO > MAXPRO'
       PRINT*,NPRO,MAXPRO
       STOP
      ENDIF

 
1000  CONTINUE

C     Load up relevant k-ordinate
C      print*,'igdist',igdist
      DO I=1,NPRO
       K_G(I)=TABK(IGDIST,I)
C       print*,igdist,k_g(i)
      ENDDO

C     First extend trajectory to see if photon is heading for ground or out to space
      RADIUS1 = RADCLOSE(PVEC,DVEC)
      IF(RADIUS1.LT.0.0.OR.(RADIUS1.GT.RADIUS))THEN
C       Photon either moving upwards or is moving downwards but doesn't hit surface
C       Find point at which it leaves the atmosphere
        CALL HITSPHERE(PVEC,DVEC,RADIUS+H(NPRO),PVEC1)
C         print*,'photon trajectory leaves atmosphere'
C         print*,PVEC1
      ELSE
C       Photon will strike surface, find point at which is strikes
        CALL HITSPHERE(PVEC,DVEC,RADIUS,PVEC1)
C         print*,'photon trajectory strikes surface'
C         print*,PVEC1
      ENDIF
     
      HEIGHT=CALCALT(PVEC1,RADIUS)
      S=0.
      DO I=1,3
       AVEC(I)=PVEC1(I)-PVEC(I)
       S=S+(PVEC1(I)-PVEC(I))**2
      ENDDO
      S=SQRT(S)
C      PRINT*,'Length of path = ',S
      DELS=S/FLOAT(MINT-1)

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
       TAUNOW(I)=DTAUDS
       TAUC(I)=DTAUDC(1)
30    CONTINUE

      TAUTOT(1)=0.
      TAUTOTC(1)=0.
      DIST(1)=0.
      DO I=2,MINT
       TAUTOT(I)=TAUTOT(I-1)+0.5*(TAUNOW(I-1)+TAUNOW(I))*DELS
       TAUTOTC(I)=TAUTOTC(I-1)+0.5*(TAUC(I-1)+TAUC(I))*DELS
       DIST(I)=DIST(I-1)+DELS
      ENDDO
 
      IF(TAUTOT(MINT).LE.TAUREQ)THEN
C       print*,'Atmosphere too thin'
C      There is insufficient opacity in atmosphere. Photon either leaves atmosphere
C      entirely or strikes surface. Add on 10km to finishing position to make sure
C      subsequent codes realises that the photon has run out of atmosphere
       DO I=1,3
        PVEC(I)=PVEC1(I)+DVEC(I)*10.
       ENDDO
    
      ELSE
c      Interpolate the optical depth array to find where photon ends up.
       J=-1
       I=1
       DO I=2,MINT
        IF(TAUREQ.GE.TAUTOT(I-1).AND.TAUREQ.LT.TAUTOT(I))THEN
         J=I
         F=(TAUREQ-TAUTOT(I-1))/(TAUTOT(I)-TAUTOT(I-1))  
        ENDIF
       ENDDO
       IF(J.LT.0)THEN
C       Error in intpath.f. Cannot find location in array
        PRINT*,'Error in INTPATH.f. Cannot find location in array'
        PRINT*,TAUREQ
        PRINT*,(TAUTOT(I),I=1,MINT)
        STOP    
       ELSE
        DO I=1,3
         TVEC(I)=PVEC(I)+(DIST(J)+F*DELS)*DVEC(I)
        ENDDO
        HEIGHT = CALCALT(TVEC,RADIUS)
C       Find place in height array and local pressure, temperature
        CALL INTERP_PT(NPRO,H,P,T,HEIGHT,PEND,TEND,F,IFL)

C       Determine optical depths/km at this point in the atmosphere
        CALL CALCTAUGRAD(NPRO,NCONT,H,P,T,DUST,MOLWT,XSEC,IRAY,VV,
     1   K_G,HEIGHT,DTAUDS,DTAUDC,DTAUR)

        DTAUSC=0.0
        DO J=1,NCONT
         TAUSCAT(J)=DTAUDC(J)*XOMEGA(J)
         DTAUSC = DTAUSC + TAUSCAT(J)
        ENDDO
        IF(IRAY.EQ.1)THEN
         TAUSCAT(NCONT+1)=DTAUR
         DTAUSC = DTAUSC + DTAUR
        ENDIF

        IF(DTAUDS.GT.0)THEN
         FSCAT = DTAUSC/DTAUDS
        ELSE
         FSCAT = 0.0
        ENDIF

C       Set PVEC to the new position vector before return
        DO I=1,3
         PVEC(I)=TVEC(I)
        ENDDO

       ENDIF

      ENDIF

      RETURN

      END
