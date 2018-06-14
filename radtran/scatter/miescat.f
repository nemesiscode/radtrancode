************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C				SUBROUTINE MIESCAT
C
C	Calculates the phase function, scattering, and extinction
C	coefficients at a a given wavelength for a specified size 
C	distribution using Mie theory.
C
C	The distribution can be "modified gamma", "log-normal", or
C	single particle size depending upon the input parameters.   
C
C	Integration is carried out over a a size distribution, calling 
C	DMIE for each size.
C
C	88-11-11  -- Adapted from program RPHASE and subroutine DMIE
C		     obtained from PNF
C	88-12-06  -- Added multiple size distribution modes.
C	89-02-15  -- Removed RR,ANR arrays, no more max. number of 
C	             points: terminate integration when N < 10**-9*Nmax
C	89-02-20  -- Added option for Log-Normal distribution;  terminate 
C		     integration when N*QS < 10**-6*(N*QS)max
C
C	Parameters:
C
C	(All R*4, will be converted to/from R*8 for computations)
C
C	Input:
C
C	XLAM	Wavelength in microns.
C	DSIZE	Particle size distribution parameters: 4 parameters
C		per mode, 10 modes allowed. 
C		For each mode:
C			If DSIZE(3) > 0
C				either standard gamma distribution, 
C				DSIZE = (A, B, ALPHA)
C				or MCS gamma distribution, 
C				DSIZE = (A, B, C)
C			ELSE
C				Log normal distribution, 
C				DSIZE = (R0,SIG,0)
C			ENDIF
C	NMODE	Number of particle modes
C	RS	Size integration limits and stepsize. If RS(2) < RS(1),
C		integration is terminated when N * QS < 10**-6 * Max, and
C		RS(2) is returned. 
C	REFINDX Real and complex indices of refraction.
C	THETA	Angles at which to return PHASE (must be < 90) and 
C		should (but need not) be ordered in ascending sequence.
C	NTHETA 	Number of points in THETA array.
C
C	Output:
C
C	XSCAT 	Mean scattering X-section (cm**2)
C	XEXT	Mean extinction X-section (cm**2)
C	PHAS	Phase function, contains 2*NTHETA or 2*NTHETA-1
C		points, depending upon whether 90 degrees was 
C		included (this arises because a value for 180-THETA 
C		as well as for THETA is returned).
C	NPHAS	Number of values in PHAS array. This should be 2*NTHETA
C		(or 2*NTHETA-1 if THETA = 90 exists).
C
C-----------------------------------------------------------------------

      SUBROUTINE MIESCAT (XLAM, ISCAT, DSIZE, NMODE, RS, REFINDX, THETA, 
     1		NTHETA, XSCAT, XEXT, PHAS, NPHAS,csratio,
     2  	RE2,TMAG2,fixtoggle)

      IMPLICIT NONE
      INTEGER		MAXTH, MAXMODE
      PARAMETER 	(MAXTH=100, MAXMODE=10)

      INTEGER	NPHAS, NTHETA, INR, J, N, M, I, K,  NMODE,  MO
      INTEGER   ISCAT(NMODE)
      REAL 	XLAM, DSIZE(3,NMODE), RS(3), REFINDX(2), PHAS(NPHAS), 
     1 		THETA(NTHETA), XSCAT, XEXT
	REAL	R2, RE2,TMAG2, csratio
	INTEGER	fixtoggle
      DOUBLE PRECISION NQMAX(MAXMODE), RMAX(MAXMODE)
      DOUBLE PRECISION RFR, RFI, RR, ANR, DELR, R1, XX, AA, BB, 
     1          ALPHA, GAMMA, CC,
     2		KSCAT, KEXT, ANORM, FUNC(4,2*MAXTH), VV, PHAS0(2*MAXTH), 
     3		CQ, X1, DX, QEXT, QSCAT, ELTRMX(4, MAXTH, 2), ANR1,
     4		THETD(MAXTH), PI,PI2,R2DOB, RE2DOB,TMAG2DOB
      LOGICAL CONT, CONT0

      PARAMETER (PI = 3.1415926535897932D0)

      PI2 = 1.D0/DSQRT(2.D0*PI)


      
C       PRINT*,'-------------------'
C       PRINT*,XLAM, (ISCAT(I),I=1,NMODE)
C       PRINT*,(DSIZE(I,1),I=1,3)
C       PRINT*,NMODE, (RS(I),I=1,3) 
C       PRINT*,(REFINDX(I),I=1,2) 
C       PRINT*,NTHETA,(THETA(I),I=1,NTHETA)
C       PRINT*,NPHAS
C       PRINT*,XSCAT,XEXT
C       PRINT*,NPHAS,(PHAS(I),I=1,NPHAS)
C       PRINT*,'-------------------'

      IF (NTHETA.GT.MAXTH) THEN
		WRITE (*,*) '  TOO MANY ANGLE POINTS: NTHETA = ',NTHETA,
     1			' MAXTH = ', MAXTH
		STOP
      ENDIF
      IF (NMODE.GT.MAXMODE) THEN
		WRITE (*,*)'  TOO MANY PARTICLE MODES; NMODE = ', NMODE,
     1			' MAXMODE = ', MAXMODE
        STOP
      ENDIF

      DO J = 1, NTHETA
       IF ((THETA(J).LT.0).OR.(THETA(J).GT.90.)) THEN
        WRITE (*,*) ' ANGLE <0 OR > 90'
        STOP
       ENDIF
       THETD(J) = THETA(J)
      ENDDO
		
C-----------------------------------------------------------------------
C
C	Size integration parameters.
C
C-----------------------------------------------------------------------

C       print*,'Miescat RS = ',RS
       R1 = RS(1)
       DELR = RS(3)
       IF (RS(2).LT.RS(1)) THEN
 	INR = 1000000001
		CONT0 = .FALSE.
       ELSE
		INR = 1+INT((RS(2)-RS(1))/RS(3))
		IF (INR.GT.1 .AND. INR/2*2.NE.INR) INR = INR+1
		CONT0 = .TRUE. 
       ENDIF

C       print*,'Miescat : ',RS
C       print*,'Miescat : ',INR
C-----------------------------------------------------------------------
C
C	Compute the peaks of the size distribution.
C
C-----------------------------------------------------------------------

	IF (.NOT.CONT0) THEN
		DO M = 1,NMODE
	  		NQMAX(M) = 0.
			IF(DSIZE(2,M).NE.0.)THEN
				AA = DSIZE(1,M)	! IF LOG-NORM = R0
				BB = DSIZE(2,M)	! IF LOG-NORM = SIG
				ALPHA = 0.0
				CC = 0.0
C    RMAX(M) is radius where distribution is maxiumum.
                                IF(ISCAT(M).EQ.1) THEN
C				 Standard Gamma
 				 ALPHA = DSIZE(3,M)	! = (1-3B)/B
				 RMAX(M) = ALPHA*AA*BB 
   				ELSEIF (ISCAT(M).EQ.2) THEN
C				 Log-normal
				 RMAX(M) = DEXP(DLOG(AA)-BB*BB)
                                ELSE
C				 MCS modified Gamma
				 CC = DSIZE(3,M)
                                 RMAX(M)=(AA/(BB*CC))**(1.0/CC)
                                ENDIF
		        ENDIF
		ENDDO
	ENDIF

C-----------------------------------------------------------------------
C
C	Initialise variables for use in DMIE.
C
C-----------------------------------------------------------------------

	DO J = 1,NPHAS
		PHAS0(J) = 0.0
	ENDDO
	KSCAT = 0.0
	KEXT = 0.0
	ANORM = 0.0

	RFR = REFINDX(1)
	RFI = REFINDX(2)

C-----------------------------------------------------------------------
C
C	Main loop: Calculation of phase function for each dropsize.
C
C-----------------------------------------------------------------------

	DO M = 1,INR

		RR = R1 + (M-1)* DELR

C		print*,M,INR,RR,XX
C		Use homogeneous sphere scattering model unless otherwise specified
		if(csratio.eq.-1.0)then
			XX = 2.0*PI*RR/XLAM
			CALL DMIE (XX, RFR, RFI, THETD, NTHETA, QEXT,
     1				QSCAT, CQ, ELTRMX)

C		Use coated sphere scattering model if Maltmieser explicitly specified
		else
			R2DOB = R2		!convert to double precision
			RE2DOB = RE2
			TMAG2DOB = TMAG2
			XX = 2.0*PI/XLAM
			R2 = SNGL(RR*(1.0-csratio)**(1.0/3.0))!find core rad
			if(fixtoggle.eq.1)then!fix core ref index and retrieve shell ref index
C			print*, 'MM fixed core: Calling DMIESS'
			CALL DMIESS(RR, RFR, RFI, THETD, NTHETA, QEXT, 
     1				 QSCAT, CQ, ELTRMX, 
     2				R2DOB, RE2DOB, TMAG2DOB, XX)	
			elseif(fixtoggle.eq.0)then!fix shell ref index and retrieve core ref index
C			print*, 'MM fixed shell: Calling DMIESS'
			CALL DMIESS(RR, RE2DOB, TMAG2DOB,
     1				 THETD, NTHETA, QEXT, QSCAT,
     2				 CQ, ELTRMX, R2DOB, RFR, RFI, XX)
			else
				print*,'Error miescat:'
				print*,'fixtoggle must equal 0 or 1'
				print*,csratio,fixtoggle
				stop
			endif 
		endif

		DO J = 1, NTHETA		! THETA <,= 90
			DO I = 1,4
			   IF(QEXT.GE.0.0)THEN
				FUNC(I,J) = ELTRMX(I,J,1)
			   ELSE
				FUNC(I,J) = -999.9
			   END IF
			ENDDO
		ENDDO
		DO J = NTHETA + 1, NPHAS	! THETA > 90
			DO I = 1,4
			   IF(QEXT.GE.0.0)THEN
				FUNC(I,J) = ELTRMX(I,NPHAS+1-J,2)
			   ELSE
				FUNC(I,J) = -999.9
			   END IF
			ENDDO
		ENDDO

C-----------------------------------------------------------------------
C
C	Compute particle distribution functions and check for size 
C	cut-off. 
C
C-----------------------------------------------------------------------

		ANR = 0.
		CONT = CONT0
		DO MO = 1, NMODE
			IF(DSIZE(2,MO).NE.0.)THEN
				AA = DSIZE(1,MO)	!R0 IF LOG-NORMAL
				BB = DSIZE(2,MO)	!SIG "  "    " 
				ALPHA = 0.0
				CC = 0.0
                                IF(ISCAT(MO).EQ.1) THEN ! STANDARD-GAMMA
				 ALPHA = DSIZE(3,MO)	!(1-3B)/B IF GAMMA
				 ANR1 = (RR**ALPHA) * 
     1						EXP(-RR/(AA*BB))
				ELSEIF(ISCAT(MO).EQ.3) THEN ! MCS GAMMA
				 CC=DSIZE(3,MO)
				 ANR1 = (RR**AA) * 
     1						DEXP(-BB*RR**CC)
				ELSE			! LOG-NORMAL
				 ANR1 = PI2/(RR*BB) * EXP(-0.5D0
     1						* ((DLOG(RR)- DLOG(AA))
     2						/BB)**2)
				ENDIF
			ELSE
				ANR1=1.			!SINGLE PARTICLE
			END IF
			ANR = ANR + ANR1
			NQMAX(MO) = DMAX1(NQMAX(MO),(ANR1*QSCAT))
                        IF(.NOT.CONT)THEN
 			 IF ((RR.LT.RMAX(MO)).OR.(ANR1*QSCAT.GT.
     1				1.E-6*NQMAX(MO))) CONT = .TRUE.
			 ENDIF
		ENDDO
		IF ((.NOT.CONT0).AND.(.NOT.CONT)) THEN
C		 WRITE (*,*) ' Miescat: size integration terminated'
C		 WRITE (*,*) ' Wavelength: ', XLAM, ' r =', RR

C		 RS(2) = RR
		 GOTO 10
		ENDIF

C-----------------------------------------------------------------------
C
C	Initialise for integration.
C
C-----------------------------------------------------------------------

		IF (M/2*2.NE.M) THEN
			VV = 2.0*DELR/3.0	! M IS ODD
		ELSE
			VV = 4.0*DELR/3.0	! M IS EVEN
		ENDIF
		IF ((M.EQ.1).OR.(M.EQ.INR)) VV = DELR/3.0

C-----------------------------------------------------------------------
C 
C	Integration over drop size by Simpson's rule.
C
C-----------------------------------------------------------------------
       
		IF(QEXT.GE.0)THEN
 		 DO J = 1, NPHAS
			PHAS0(J) = PHAS0(J) + 0.5 * ANR * VV *
     1				(FUNC(1,J) + FUNC(2,J))
		 ENDDO
		 KSCAT = KSCAT + PI * RR * RR * QSCAT * ANR * VV
		 KEXT = KEXT + PI * RR * RR * QEXT  * ANR * VV
        	 ANORM = ANORM + ANR * VV
		END IF
	ENDDO

C-----------------------------------------------------------------------
C
C	Normalise integrations. Cross sections are returned as cm2.
C
C-----------------------------------------------------------------------

10	CONTINUE
        IF(ANORM.GT.0.0)THEN
         XSCAT = SNGL(KSCAT/ANORM * 1.e-8)
 	 XEXT = SNGL(KEXT/ANORM * 1.e-8)
        ELSE
         XSCAT = 0.0
         XEXT = 0.0
         KSCAT = 1.0	!Dummy to stop PHAS going silly
        END IF
	DO J = 1, NPHAS
		PHAS(J) = XLAM * XLAM * sngl(PHAS0(J)/(PI*KSCAT))
	ENDDO

	RETURN

	END

************************************************************************
************************************************************************
