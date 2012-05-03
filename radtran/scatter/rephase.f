      SUBROUTINE REPHASE(MU, NMU, NPHI, ISCAT,
     1 CONS, NCONS, ICONT, NCONT, VWAVE, PTPL, PTMI)
C     $Id: rephase.f,v 1.2 2011-06-17 15:57:54 irwin Exp $
C     ******************************************************************
C
C     Subroutine to calculate the phase function integrated over azimuth angle
C     for a prespecified set of angles.
C
C     Input variables:
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
C     Output variables:
C	PTPL(5,41,MAXMU,MAXMU)	DOUBLE	Integrated phase function coefficients 
C					in 'plus' direction (i.e. downwards)
C					for up to 40 fourier components
C	PTMI(5,41,MAXMU,MAXMU)	DOUBLE	Integrated phase function coefficients
C					in 'minus' direction (i.e. upwards)
C					for up to 40 fourier components
C
C      26DEC96 .. Adapted from PHASINT1 to calculate all fourier componenents
C		  in one go.
C      05MAY96    Adapted to calculate fourier components using matrices
C     ******************************************************************

      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NMU, NPHI, ISCAT, I, J, K, IP, NITER
      INTEGER NCONS,KL,J1,ICONT,NCONT,NF,IDUMP,IDIM
      INTEGER K1,K2,ITAB,LF
      PARAMETER (IDIM=101)
      DOUBLE PRECISION TABLEP(MAXMU,MAXMU,MAXF,MAXF)
      DOUBLE PRECISION TABLEM(MAXMU,MAXMU,MAXF,MAXF)
      DOUBLE PRECISION PTPL(MAXCON,MAXF,MAXMU,MAXMU),
     1 PTMI(MAXCON,MAXF,MAXMU,MAXMU),CONS(MAXSCATPAR), CONS1(MAXSCATPAR), MU(MAXMU),
     2 XI, XJ, P, XFAC,DPHI, CALPHA, PI, X0, X1, Y0, Y1, C0, C1, 
     3 PLI0, STEST, XM,PLI1, PLJ0, PLJ1, XF, YF, PLX, PMX, SSUM, 
     4 STHI, STHJ
      REAL VWAVE,SUM,SUM1
      INTEGER NY1,NY2,NA1,NA2,NC1,NC2
      DOUBLE PRECISION PHI,WPHI,PHASL(200,3),PHASM(200,3)
      DOUBLE PRECISION CPL(IDIM,IDIM),CPM(IDIM,IDIM),PL,COEF(MAXF),
     1  CPHI,TPA(IDIM)
      DOUBLE PRECISION APL(IDIM,IDIM),APM(IDIM,IDIM),
     1  YY(IDIM,IDIM),XPA(IDIM)
      CHARACTER*30 FOURIER 

      DOUBLE PRECISION CPL1,CMI1,PL1,PM1
 
      COMMON /FTABLES/TABLEP,TABLEM
      COMMON /FTAB1/ITAB
      COMMON /SCATDUMP/ IDUMP

      PARAMETER (STEST = 0.5)
     
C     The variable IDUMP has been added but has not yet been included in
C     the calling parameter list. This will be done when things are quiet
C     again!
C
C     ******************************************************************


       PI = 4.0D0*DATAN(1.0D0)

       DPHI = PI/NPHI

       NF=MAXF-1

C      Set up conversion table if ITAB <> -1
       IF(ITAB.NE.-1)THEN
        PRINT*,'Rephase: Calculating Fourier table'
        CALL FTABLE( MU, NMU, NPHI, NF)
        PRINT*,'Table calculated'
       ENDIF

C      Fourier analyse phase function

       DO KL=1,MAXF
        COEF(KL)=0.0
       END DO

       XM = 0.0
       DO 200 K = 0,NPHI
        PHI = K*DPHI
        CPHI = DCOS(PHI)
        CALL PHASE1(CPHI,PL,ISCAT,CONS,NCONS,ICONT,NCONT,VWAVE)
        TPA(K+1)=PL
        XM = XM + 1.0*PL/(NPHI+1)
	XPA(K+1)=0.0
200    CONTINUE

C      Keep calculating the fourier coefficients until the sum of
C      ((phase - calc_phase)/mean_phase)^2 is less than a prespecified
C      value
       XFAC = 0.5
       DO 100 KL=0,MAXF
         LF = KL
         COEF(KL+1) = 0.0
         IF(KL.GT.0)XFAC=1.0
         DO 150 K=0,NPHI
           PHI = K*DPHI
           PLX = TPA(K+1)*DCOS(KL*PHI)
           WPHI = DPHI*XFAC/PI
           IF(K.EQ.0.OR.K.EQ.NPHI)WPHI=WPHI*0.5          
           COEF(KL+1)=COEF(KL+1)+2*WPHI*PLX
150      CONTINUE

         SSUM=0.0
         DO 160 K=0,NPHI
          PHI = K*DPHI
          XPA(K+1)=XPA(K+1) + COEF(KL+1)*DCOS(KL*PHI)
          SSUM = SSUM+((TPA(K+1)-XPA(K+1))/XM)**2
160      CONTINUE
         SSUM = 100.0*SQRT(SSUM/(NPHI+1))
         IF(IDUMP.GT.0)PRINT*,KL,SSUM
         IF(SSUM.LT.STEST)GOTO 400
100    CONTINUE


400    DO K2=1,LF+1
       DO K1 = 1,LF+1
        YY(K2,K1)=0
       END DO
       YY(K2,1)=COEF(K2)
      END DO
      ny1 = nf+1
      ny2 = 1

       DO J = 1,NMU
	DO I = 1,NMU
         DO K1=1,LF+1
          DO K2=1,LF+1
           APL(K1,K2) = TABLEP(I,J,K1,K2)
           APM(K1,K2) = TABLEM(I,J,K1,K2)
          END DO
         END DO
         na1 = nf+1
         na2 = nf+1

         CALL DMULT_MAT(IDIM,APL,NA1,NA2,YY,NY1,NY2,CPL,NC1,NC2)
         CALL DMULT_MAT(IDIM,APM,NA1,NA2,YY,NY1,NY2,CPM,NC1,NC2)


         DO K2=1,LF+1
          PTPL(ICONT,K2,I,J) = CPL(K2,1)
          PTMI(ICONT,K2,I,J) = CPM(K2,1)
         END DO

         DO K2=LF+2,MAXF
          PTPL(ICONT,K2,I,J) = 0.0
          PTMI(ICONT,K2,I,J) = 0.0
         END DO

        END DO
       END DO

       IF(IDUMP.GT.0)THEN

C      Now check that the calculated Fourier coefficients provide an
C      adequate fit to the real phase function by analysing the extreme
C      case where I=NMU, J=NMU
       
       I=NMU
       J=NMU
       STHI = DSQRT(1.0D0 - MU(I)*MU(I))     ! sin(theta(i))
       STHJ = DSQRT(1.0D0 - MU(J)*MU(J))     ! sin(theta(j))
       DO K = 0,NPHI
            PHI = K*DPHI
            CPL1 = STHI*STHJ*DCOS(PHI) + MU(I)*MU(J)
            CALL PHASE1(CPL1,PL1,ISCAT,CONS,NCONS,ICONT,NCONT,VWAVE)
            CMI1 = STHI*STHJ*DCOS(PHI) - MU(I)*MU(J)
            CALL PHASE1(CMI1,PM1,ISCAT,CONS,NCONS,ICONT,NCONT,VWAVE)

            PHASL(K+1,1)=CPL1
            PHASL(K+1,2)=PL1
            PHASM(K+1,1)=CMI1
            PHASM(K+1,2)=PM1
            PHASL(K+1,3)=0
            PHASM(K+1,3)=0
       END DO


       DO K = 0,NPHI
         PHI = K*DPHI
         DO KL=0,LF
         PHASL(K+1,3)=PHASL(K+1,3)+
     &                        PTPL(ICONT,KL+1,NMU,NMU)*COS(KL*PHI)
         PHASM(K+1,3)=PHASM(K+1,3)+
     &                        PTMI(ICONT,KL+1,NMU,NMU)*COS(KL*PHI)
         END DO
       END DO

         
       IF(IDUMP.EQ.2)THEN
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
       ELSEIF(IDUMP.EQ.1) THEN
C  	      Quality of fourier fit for worst possible case of maximum
C 	      illumination and viewing angles    
              sum = 0.0
              sum1 = 0.0
	      do k=0,nphi
		sum=sum + (phasl(k+1,2) - phasl(k+1,3))**2
              end do
              sum=sqrt(sum/(1.0*(nphi+1)))
	      do k=0,nphi
		sum1=sum1 + (phasm(k+1,2) - phasm(k+1,3))**2
              end do
              sum=sqrt(sum/(1.0*(nphi+1)))
              sum1=sqrt(sum1/(1.0*(nphi+1)))
              k = 0.5*nphi
              write(*,1000)icont,phasl(1,2),phasl(k,2),sum
              write(*,1001)icont,phasl(1,3),phasl(k,3)
              write(*,1000)icont,phasm(1,2),phasm(k,2),sum1
              write(*,1001)icont,phasm(1,3),phasm(k,3)
              write(*,*)'Number of coeffs used : ',LF
1000  format(' Rephase: Scat = ',i2,' Real Extreme phase = ',f7.3,
     &f7.3,' Fourier rms : ',e9.3)
1001  format(' Rephase: Scat = ',i2,' Fit Extreme phase = ',f7.3,
     &f7.3)
      ENDIF

      ENDIF

      RETURN
      END

  
