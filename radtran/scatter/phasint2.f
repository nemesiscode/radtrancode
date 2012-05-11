      SUBROUTINE PHASINT2( NF, MU, NMU, NPHI, ISCAT,
     1 CONS, NCONS, ICONT, NCONT, VWAVE)
C     $Id: phasint2.f,v 1.3 2011-06-17 15:57:54 irwin Exp $
C     ******************************************************************
C
C     Subroutine to calculate the phase function integrated over azimuth angle
C     for a prespecified set of angles.
C
C     Input variables:
C	NF		INTEGER Number  of fourier components
C	MU(MAXMU)	DOUBLE 	Ordinates of cos(zenith)
C	NMU		INTEGER	Number of zenith ordinates
C	NPHI		INTEGER	Number of azimuth ordinates
C	ISCAT		INTEGER	Scattering ID
C	CONS(MAXSCATPAR)	DOUBLE	Constants array for phase function calculation
C	NCONS		INTEGER	Number of values in array used
C	ICONT		INTEGER	Aerosol ID
C	NCONT		INTEGER	Number of Aerosols in calculation
C	VWAVE		REAL	Wavenumber of calculation	
C
C     Output variables (to PSHARE common block):
C	PTPL(MAXCON,MAXF,MAXMU,MAXMU)	DOUBLE	Integrated phase function coefficients 
C					in 'plus' direction (i.e. downwards)
C					for up to 40 fourier components
C	PTMI(MAXCON,MAXF,MAXMU,MAXMU)	DOUBLE	Integrated phase function coefficients
C					in 'minus' direction (i.e. upwards)
C					for up to 40 fourier components
C
C      26DEC96 .. Adapted from PHASINT1 to calculate all fourier componenents
C		  in one go.

C     ******************************************************************

      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NMU, NPHI, ISCAT, I, J, K,
     1 NCONS,KL,ICONT,NCONT,NF,IDUMP
      DOUBLE PRECISION CONST
      DOUBLE PRECISION PTPL(MAXCON,MAXF,MAXMU,MAXMU)
      DOUBLE PRECISION PTMI(MAXCON,MAXF,MAXMU,MAXMU)
      DOUBLE PRECISION CONS(MAXSCATPAR), MU(MAXMU),
     1 DPHI, PI, PLX, PMX, STHI, STHJ
      REAL VWAVE,SUM
      DOUBLE PRECISION PHI,WPHI,PHASL(200,3),PHASM(200,3)
      DOUBLE PRECISION CPL,CMI,PL,PM
      CHARACTER*30 FOURIER
     
      COMMON /SCATDUMP/ IDUMP
      COMMON /PSHARE/ PTPL,PTMI
C
C     ******************************************************************


C       PRINT*,'PHASINT2'
       PI = 4.0D0*DATAN(1.0D0)
       CONST = 1.0

       IF (NPHI.LE.NF) THEN 
        PRINT*,"WARNING: NPHI has been set less than NF, which will"
        PRINT*,"result in an unphysical oscillitory solution."
        STOP
       ENDIF

       DPHI = 2.0*PI/NPHI

       DO J = 1,NMU
	DO I = 1,NMU
          DO KL=1,NF+1 
  	   PTPL(ICONT,KL,I,J) = 0.0D0
 	   PTMI(ICONT,KL,I,J) = 0.0D0
          END DO
          STHI = DSQRT(1.0D0 - MU(I)*MU(I))	! sin(theta(i))
          STHJ = DSQRT(1.0D0 - MU(J)*MU(J))	! sin(theta(j))
	  DO K = 0,NPHI
            PHI = K*DPHI
	    CPL = STHI*STHJ*DCOS(PHI) + MU(I)*MU(J) 
	    CALL PHASE1(CPL,PL,ISCAT,CONS,NCONS,ICONT,NCONT,VWAVE)
	    CMI = STHI*STHJ*DCOS(PHI) - MU(I)*MU(J) 
	    CALL PHASE1(CMI,PM,ISCAT,CONS,NCONS,ICONT,NCONT,VWAVE)

            PHASL(K+1,1)=CPL
            PHASL(K+1,2)=PL
            PHASM(K+1,1)=CMI
            PHASM(K+1,2)=PM
            PHASL(K+1,3)=0
            PHASM(K+1,3)=0

            DO KL = 0,NF

             PLX = PL*DCOS(KL*PHI)
             PMX = PM*DCOS(KL*PHI)
            
             WPHI = 1.0*DPHI
             IF(K.EQ.0.OR.K.EQ.NPHI)WPHI=0.5*DPHI

             IF(KL.EQ.0)THEN
              WPHI = WPHI/(2.0*PI)
             ELSE
              WPHI = WPHI/PI
             END IF

             PTPL(ICONT,KL+1,I,J)=PTPL(ICONT,KL+1,I,J)+WPHI*PLX
             PTMI(ICONT,KL+1,I,J)=PTMI(ICONT,KL+1,I,J)+WPHI*PMX
            END DO
          END DO


C         Now check that the calculated Fourier coefficients provide an
C	  adequate fit to the real phase function
          DO K = 0,NPHI

           PHI = K*DPHI
           DO KL=0,NF
            PHASL(K+1,3)=PHASL(K+1,3)+
     &                        PTPL(ICONT,KL+1,I,J)*COS(KL*PHI)
            PHASM(K+1,3)=PHASM(K+1,3)+
     &                        PTMI(ICONT,KL+1,I,J)*COS(KL*PHI)
           END DO

          END DO

         
          IF(I.EQ.NMU.AND.J.EQ.NMU)THEN
            IF(IDUMP.EQ.1)THEN
C	     Output complete phase function and fourier fit for worst
C	     possible case of maximum illumination and viewing angles
             FOURIER='fourier*.dat'
             FOURIER(8:8) = char(icont+48)
             OPEN(48,file=FOURIER,status='unknown')          
              write(48,*)nphi
              do k=0,nphi
               PHI = K*DPHI
               write(48,*)phi,phasl(k+1,1),phasl(k+1,2),phasl(k+1,3),
     &                       phasm(k+1,1),phasm(k+1,2),phasm(k+1,3)
              end do
              close(48)
             ELSEIF(IDUMP.EQ.2) THEN
C  	      Quality of fourier fit for worst possible case of maximum
C 	      illumination and viewing angles    
              sum = 0.0
	      do k=0,nphi
		sum=sum + sngl(phasl(k+1,2) - phasl(k+1,3))**2
              end do
              sum=sum/(1.0*(nphi+1))
              k = nphi/2
              write(*,1000)icont,phasl(1,2),phasl(k,2),sum
              write(*,1000)icont,phasl(1,3),phasl(k,3),sum
1000  format(' Phasint2: Scat = ',i2,' Extreme phase = ',f7.3,f7.3,
     &' Fourier rms : ',e9.3)
             END IF
          ENDIF

	ENDDO
      END DO


      RETURN
      END

  
