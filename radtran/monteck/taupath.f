      SUBROUTINE TAUPATH(VV,IRAY,TAUREQ,PVEC,DVEC,NPRO,NGAS,NCONT,MOLWT,
     1 RADIUS,P,T,H,DUST,TABK,IGDIST,XSEC,XOMEGA,TEND,FSCAT,TAUSCAT)
C     $Id:
C     ***************************************************************
C     Subroutine to point in a spherically symmetric atmosphere which is an 
C     optical distance of TAUREQ away from the starting position PVEC(3) in 
C     the direction DVEC(3). Code has same function as NEWPATH but is tidier.
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
C      CHARACTER*1 ANS
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
C      print*,'AA',HEIGHT,OLDTAUGRAD
      DIST=0.0
      OLDTAU=0.0
      OLDTAU1=0.0
      OLDDIST=0.0
      OLDDIST1=0.0

999   DELS = 0.1/OLDTAUGRAD
      IF(DELS.GT.DELS0)DELS=DELS0
C      print*,'BB',DELS

      IDIST=IDIST+1      
      DO I=1,3
       PVEC1(I)=PVEC1(I)+DELS*DVEC(I)
      ENDDO
 
      HEIGHT = CALCALT(PVEC1,RADIUS)
C      print*,'CC',HEIGHT
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

C      print*,'DD',IDIST,HEIGHT,DIST,DELS,DTAUDS,NEWTAU,TAUREQ

C      READ(5,1)ANS
C1     FORMAT(A)
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


