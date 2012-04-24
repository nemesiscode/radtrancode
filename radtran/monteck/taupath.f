      SUBROUTINE TAUPATH(VV,IRAY,TAUREQ,PVEC,DVEC,NPRO,NGAS,NCONT,MOLWT,
     1 RADIUS,P,T,H,DUST,TABK,IGDIST,XSEC,XOMEGA,TEND,FSCAT,TAUSCAT)
C     $Id:
C     ***************************************************************
C     Subroutine to point in a spherically symmetric atmosphere which is an 
C     optical distance of TAUREQ away from the starting position PVEC(3) in 
C     the direction DVEC(3).
C     Code uses Trapezium rule integration.
C
C     Input variables:
C	VV		REAL	Calculated wavenumber (require to estimate
C				Rayleigh scattering)
C	IRAY		INTEGER 0=Rayleigh scatt off, 1 = on
C	TAUREQ		REAL	Required transmission
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
C     Pat Irwin		29/7/05
C 
C     ***************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'

      REAL PVEC(3),DVEC(3),TAUREC,TAUTOT,VV
      INTEGER NPRO,NGAS,J,NCONT,I,K,IFL,IRAY,IDIST
      REAL DNOW,MOLWT,F,DELS,DELS0,DELS1
      REAL TNOW,PNOW,CALCALT,RADIUS,HEIGHT,TEND,PVEC1(3)
      REAL P(MAXPRO),T(MAXPRO),H(MAXPRO),PEND
      REAL TABK(MAXG,MAXPRO),K_G(MAXPRO)
      REAL T1,T2,FSCAT,DTAU,DTAUSC,TAUREQ
      REAL DTAUDS,DTAUR
      INTEGER IGDIST
      REAL DUST(MAXPRO,MAXCON),CONT(MAXCON),XSEC(MAXCON)
      REAL TAUD,XOMEGA(MAXCON),TAUSCAT(MAXCON),DTAUDC(MAXCON)
      REAL DIST,OLDTAU,NEWTAU,OLDTAUGRAD,OLDTAU1,OLDDIST,OLDDIST1     
C     ------------------------------------------------------------------
C     Check some numbers
C     ------------------------------------------------------------------
      IF(NCONT.GT.MAXCON)THEN
       PRINT*,'TAUPATH: NCONT > MAXCON'
       PRINT*,NCONT,MAXCON
       STOP
      ENDIF

      IF(NPRO.GT.MAXPRO)THEN
       PRINT*,'TAUPATH: NPRO > MAXPRO'
       PRINT*,NPRO,MAXPRO
       STOP
      ENDIF


C     Load up relevant k-ordinate
C      print*,'igdist',igdist
      DO I=1,NPRO
       K_G(I)=TABK(IGDIST,I)
      ENDDO

      DELS0=(H(NPRO)-H(1))/10.0	! default max. path length (km)

      DO I=1,3
       PVEC1(I)=PVEC(I)
      ENDDO

C      print*,'PVEC',PVEC
      IDIST=1
      HEIGHT = CALCALT(PVEC1,RADIUS)      
      CALL CALCTAUGRAD(NPRO,NCONT,H,P,T,DUST,MOLWT,XSEC,IRAY,VV,
     1 K_G,HEIGHT,DTAUDS,DTAUDC,DTAUR)

      OLDTAUGRAD = DTAUDS
      DIST=0.0
      OLDTAU=0.0
      OLDTAU1=0.0
      OLDDIST=0.0
      OLDDIST1=0.0

999   DELS = 0.1/OLDTAUGRAD
      IF(DELS.GT.DELS0)DELS=DELS0


      IDIST=IDIST+1      
      DO I=1,3
       PVEC1(I)=PVEC1(I)+DELS*DVEC(I)
      ENDDO
 
      HEIGHT = CALCALT(PVEC1,RADIUS)

      IF(HEIGHT.LT.H(1).OR.HEIGHT.GT.H(NPRO))THEN
         DO I=1,3
          PVEC(I)=PVEC1(I)
         ENDDO
         RETURN
      ENDIF

      CALL CALCTAUGRAD(NPRO,NCONT,H,P,T,DUST,MOLWT,XSEC,IRAY,VV,
     1 K_G,HEIGHT,DTAUDS,DTAUDC,DTAUR)
      
      DIST = DIST + DELS
      NEWTAU = OLDTAU + 0.5*(DTAUDS+OLDTAUGRAD)*DELS

C      print*,IDIST,HEIGHT,DIST,DELS,NEWTAU,TAUREQ

      IF(NEWTAU.LT.TAUREQ)THEN
       OLDTAU1 = OLDTAU
       OLDTAU = NEWTAU
       OLDTAUGRAD = DTAUDS
       OLDDIST1 = OLDDIST
       OLDDIST = DIST
       GOTO 999
      ENDIF

C     Calculate final interpolated position
      DELS1 = DELS*(TAUREQ - OLDTAU)/(NEWTAU-OLDTAU) 
      DO I=1,3
       PVEC(I) = PVEC1(I)+DELS1*DVEC(I)
      ENDDO

C     Calculate conditions at end of path
      HEIGHT = CALCALT(PVEC,RADIUS)
      CALL INTERP_PT(NPRO,H,P,T,HEIGHT,PEND,TEND,F,IFL)

C      print*,'HEIGHT,TEND',HEIGHT,TEND
      CALL CALCTAUGRAD(NPRO,NCONT,H,P,T,DUST,MOLWT,XSEC,IRAY,VV,
     1 K_G,HEIGHT,DTAUDS,DTAUDC,DTAUR)
      
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


      RETURN

      END


      SUBROUTINE CALCTAUGRAD(NPRO,NCONT,H,P,T,DUST,MOLWT,XSEC,IRAY,VV,
     1 K_G,HEIGHT,DTAUDS,DTAUDC,DTAUR)

      IMPLICIT NONE
      include '../includes/arrdef.f'
      include '../includes/constdef.f'
      INTEGER NPRO,IFL,K,NCONT,IRAY,J
      REAL DUST(MAXPRO,MAXCON),XSEC(MAXCON),MOLWT
      REAL HEIGHT,H(NPRO),PNOW,P(NPRO),TNOW,T(NPRO),F
      REAL K_G(NPRO),K1,RAYLEIGHJ,VV,C,TOTAM,DKDS,DNOW
      REAL DKDC,DTAUR,DTAUDC(MAXCON),DTAUDS


      CALL INTERP_PT(NPRO,H,P,T,HEIGHT,PNOW,TNOW,F,IFL)

      K1 = (1-F)*K_G(IFL) + F*K_G(IFL+1)
      TOTAM = MODBOLTZ*(PNOW/TNOW)
      DKDS = TOTAM*K1*1E-20	! optical depth of gas/km

      DKDC = 0.0		! optical depth of aerosol/km
      DO J=1,NCONT
        DNOW = (1-F)*DUST(IFL,J) + F*DUST(IFL+1,J)
        DTAUDC(J)=XSEC(J)*DNOW*TOTAM*MOLWT/AVOGAD
        DKDC = DKDC+DTAUDC(J)
      ENDDO

      IF(IRAY.EQ.1)THEN
       DTAUR = TOTAM*RAYLEIGHJ(VV,PNOW,TNOW)
      ELSE
       DTAUR = 0.0
      ENDIF

C     Calculate optical depth/km at current conditions
      DTAUDS = DKDS + DKDC + DTAUR


      RETURN

      END

      SUBROUTINE INTERP_PT(NPRO,H,P,T,HEIGHT,PNOW,TNOW,F,IFL)
      IMPLICIT NONE
      INTEGER NPRO,IFL,K
      REAL HEIGHT,H(NPRO),PNOW,P(NPRO),TNOW,T(NPRO),F
      F=-1.0
      IFL = 0
      DO K=1,NPRO-1
        IF(HEIGHT.GE.H(K).AND.HEIGHT.LT.H(K+1))THEN
         F = (HEIGHT - H(K))/(H(K+1)-H(K))
         IFL=K
        ENDIF
      ENDDO
      IF(IFL.EQ.0)THEN
        IF(HEIGHT.LT.H(1))THEN
         IFL=1
         F=0.0
        ENDIF
        IF(HEIGHT.GE.H(NPRO))THEN
         IFL=NPRO-1
         F=1.0
        ENDIF
      ENDIF


      PNOW = (1-F)*P(IFL) + F*P(IFL+1)
      TNOW = (1-F)*T(IFL) + F*T(IFL+1)


      RETURN

      END

