      SUBROUTINE NEWPATH(VV,IRAY,TAUREQ,PVEC,DVEC,NPRO,NGAS,NCONT,MOLWT,
     1 RADIUS,P,T,H,DUST,TABK,IGDIST,XSEC,XOMEGA,TEND,FSCAT,TAUSCAT)
C     $Id:
C     ***************************************************************
C     Subroutine to calculate position in a spherically symmetric 
C     atmosphere which is the required optical thickness away, in a 
C     direction specified, from a given starting position.
C    
C     Code uses Trapezium rule integration.
C
C     Input variables:
C	VV		REAL	Calculated wavenumber (require to estimate
C				Rayleigh scattering)
C	IRAY		INTEGER 0=Rayleigh scatt off, 1 = on
C	TAUREQ		REAL	Required optical path
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
C     Pat Irwin		1/8/05
C 
C     ***************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'
      REAL PVEC(3),DVEC(3),TAUREC,TAUTOT,TAUR1,TAUR2,TAUR,VV,RAYLEIGHJ
      INTEGER NPRO,NGAS,J,NCONT,I,K,IFL,IRAY
      REAL DNOW,MOLWT,F,DELS,DELS0
      REAL TNOW,PNOW,CALCALT,RADIUS,HEIGHT,TEND,DELTAU,PVEC1(3)
      REAL P(MAXPRO),T(MAXPRO),H(MAXPRO)
      REAL TABK(MAXG,MAXPRO),TAUG,TAUA,TAUS,K_G
      REAL P1,P2,T1,T2,K1,K2,TOTAM1,TOTAM2,DKDS1,DKDS2,FSCAT,F1
      REAL TAUREQ,TAU1,DKDC1,DKDC2,DTAUDS1,DTAUDS2
      INTEGER IGDIST
      REAL DUST(MAXPRO,MAXCON),CONT(MAXCON),XSEC(MAXCON),C1(MAXCON)
      REAL C2(MAXCON),TAUD,XOMEGA(MAXCON),TAUSCAT(MAXCON)

C     ------------------------------------------------------------------

C     Check some numbers
C     ------------------------------------------------------------------
      IF(NGAS.GT.MAXGAS)THEN
       PRINT*,'NEWPATH: NGAS > MAXGAS'
       PRINT*,NGAS,MAXGAS
       STOP
      ENDIF

      IF(NCONT.GT.MAXCON)THEN
       PRINT*,'NEWPATH: NCONT > MAXCON'
       PRINT*,NCONT,MAXCON
       STOP
      ENDIF

      IF(NPRO.GT.MAXPRO)THEN
       PRINT*,'NEWPATH: NPRO > MAXPRO'
       PRINT*,NPRO,MAXPRO
       STOP
      ENDIF

      DELS0=50.0		! default max. path length (km)
      TAU1=0.0			! optical depth of new path

      HEIGHT = CALCALT(PVEC,RADIUS)
      CALL INTERP_PT(NPRO,H,P,T,HEIGHT,P1,T1,F,IFL)

      K1 = (1-F)*TABK(IGDIST,IFL) + F*TABK(IGDIST,IFL+1)
      TOTAM1 = MODBOLTZ*(P1/T1)
      DKDS1 = TOTAM1*K1*1E-20	! optical depth of gas/km

      DKDC1 = 0.0		! optical depth of aerosol/km
      DO J=1,NCONT
        DNOW = (1-F)*DUST(IFL,J) + F*DUST(IFL+1,J)
        C1(J)=DNOW*TOTAM1*MOLWT/AVOGAD
        DKDC1 = DKDC1+C1(J)*XSEC(J)
      ENDDO


      IF(IRAY.EQ.1)THEN
       TAUR1 = TOTAM1*RAYLEIGHJ(VV,P1,T1)
      ELSE
       TAUR1 = 0.0
      ENDIF

C     Calculate optical depth/km at current conditions
      DTAUDS1 = DKDS1 + DKDC1 + TAUR1

132   CONTINUE

      DELS=DELS0
      IF(DTAUDS1.GT.0.0)DELS = 0.1/DTAUDS1

      IF(DELS.GT.DELS0)DELS=DELS0

      IF(DELS.EQ.0.0)THEN 
       PRINT*,'NEWPATH. DELS GONE to ZERO'
       DELS=1E-25
      ENDIF

133   CONTINUE

C     Try moving down to new position and calculate OD of path traversed
      DO J=1,3
       PVEC1(J) = PVEC(J) + DELS*DVEC(J)
      ENDDO

      HEIGHT = CALCALT(PVEC1,RADIUS)


C      IF(HEIGHT.LT.H(1).OR.HEIGHT.GT.H(NPRO))THEN
C       DO J=1,3
C        PVEC(J)=PVEC1(J)
C       ENDDO
C       GOTO 101
C      ENDIF

      CALL INTERP_PT(NPRO,H,P,T,HEIGHT,P2,T2,F,IFL)

      K2 = (1-F)*TABK(IGDIST,IFL) + F*TABK(IGDIST,IFL+1)
    
      TOTAM2 = MODBOLTZ*(P2/T2)
      DKDS2 = TOTAM2*K2*1E-20	! new optical depth of gas/km

      DKDC2=0.0	! new optical depth of aerosol/km
      DO J=1,NCONT
        DNOW = (1-F)*DUST(IFL,J) + F*DUST(IFL+1,J)
        C2(J)=DNOW*TOTAM1*MOLWT/AVOGAD
        DKDC2=DKDC2 + C2(J)*XSEC(J)

        CONT(J)=0.5*(C1(J)+C2(J))

      ENDDO

      IF(IRAY.EQ.1)THEN
       TAUR2 = TOTAM2*RAYLEIGHJ(VV,P2,T2)
      ELSE
       TAUR2 = 0.0
      ENDIF

C     Calculate optical depth/km at new position
      DTAUDS2 = DKDS2 + DKDC2 + TAUR2

      print*,height,p2,t2,dkds2,dkdc2,taur2

      TAUR = 0.5*(TAUR1+TAUR2)*DELS	! optical thickness of rayleigh 
C					  scattering

      TAUG = 0.5*(DKDS1+DKDS2)*DELS   	! optical thickness of gas

      TAUD = 0.0	! optical thickness of aerosol
      TAUS = 0.0	! scattering opacity
      DO J=1,NCONT
       TAUD = TAUD + CONT(J)*XSEC(J)*DELS
       TAUSCAT(J) = CONT(J)*XSEC(J)*XOMEGA(J)*DELS
       TAUS = TAUS + TAUSCAT(J)
      ENDDO


C     Add on Rayleigh scattering
      IF(IRAY.EQ.1)THEN
       TAUD = TAUD + TAUR
       TAUS = TAUS + TAUR
       TAUSCAT(NCONT+1)=TAUR
      ENDIF
    
      TAUA = TAUD-TAUS	! Absorption opacity of aerosol

      IF((TAUG+TAUD).NE.0.0)THEN
        FSCAT = TAUS/(TAUG+TAUD)
      ELSE
        FSCAT = 0.0
      ENDIF

      TEND = 0.5*(T1+T2)

      DELTAU = TAUG+TAUD

C     If the opacity of the path is greater than 0.1, then we have come
C     too far. We need to reduce the step length and try again
      IF(DELTAU.GT.0.1)THEN
       DELS = DELS*0.095/DELTAU
C       PRINT*,'DELTAU > 0.1',DELTAU
       GOTO 133
      ENDIF

C     Calculate revised OD of total path to new point in atmosphere
      TAUTOT = TAU1 + DELTAU

      IF(TAUTOT.LT.TAUREQ) THEN
        K1=K2
        P1=P2
        T1=T2
        TOTAM1 = TOTAM2
        DKDS1 = DKDS2
        TAU1 = TAUTOT
        DO J=1,NCONT 
         C1(J)=C2(J)
        ENDDO
        TAUR1 = TAUR2
        DTAUDS1 = DTAUDS2
        DKDC1 = DKDC2
C	Start photon from new position
        DO I=1,3
         PVEC(I)=PVEC1(I)
        ENDDO

        GOTO 132

      ELSE

        F1 = (TAUREQ-TAU1)/(TAUTOT-TAU1)

        DO 21 I=1,3
         PVEC(I) = PVEC(I) + F1*DVEC(I)*DELS
21      CONTINUE

      ENDIF

101   CONTINUE

      RETURN

      END

