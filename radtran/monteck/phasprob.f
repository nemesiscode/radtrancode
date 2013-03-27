      SUBROUTINE PHASPROB(NCONT,XHG,NPHASE,THETA)
C     ******************************************************************
C     Subroutine calculates probability particle with H-G phase 
C     function defined by XHG, scattering into different angles.
C     Method from Hansen and Travis (1974)
C     Routine first calculates at high resolution and then interpolates 
C     results on to equal probability grid THETA.
C
C     Input parameters
C	NCONT		INTEGER		Number of particle types
C	XHG(MAXCON,3)	REAL		Combined H-G function
C	NPHASE		INTEGER		Number of points in THETA
C
C     Output parameter
C	THETA(MAXCON,100) REAL		Scattering angles (radians) spaced
C					with equal probability
C
C     Pat Irwin		Original 	11/3/99
C			Tidied		28/7/05
C
C     ******************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NPHASE,MPHASE,I,J,NCONT,ICONT
      REAL PI
      PARAMETER(PI=3.1415927,MPHASE=1000)
      REAL THETA(MAXCON,100),XHG(MAXCON,3)
      REAL TEMP(MPHASE),HENYEYMC,XTEMP(MPHASE),XTH,DTH
      REAL XR,R,XF,THG(3),CALPHA

      DTH = PI/(1.0*MPHASE-1.0)


      DO 999 ICONT=1,NCONT

       DO I=1,3
        THG(I)=XHG(ICONT,I)
       ENDDO

       IF(THG(1).LT.0.0)THEN
        PRINT*,'PHASPROB: Set exact dipole scattering'
       ENDIF

       
       DO 10 I=1,MPHASE
        XTH = (I-1)*DTH
C       If fraction in Combined H-G is +ve then use, otherwise assume
C	dipole phase function.
        IF(THG(1).GT.0.0)THEN
         TEMP(I) = HENYEYMC(XTH,THG)*SIN(XTH)
        ELSE
         CALPHA = COS(XTH)
         TEMP(I) = 0.75*(1.0 + CALPHA**2)*SIN(XTH)
        ENDIF        
10     CONTINUE

       XTEMP(1)=0
       DO 20 I=2,MPHASE
        XTEMP(I) = XTEMP(I-1)+0.5*(TEMP(I)+TEMP(I-1))*DTH
20     CONTINUE
      
       DO 30 I=1,MPHASE
        XTEMP(I)=XTEMP(I)/XTEMP(MPHASE)
30     CONTINUE


       XR = 1.0*(NPHASE-1)
       THETA(ICONT,1)=0.0
       THETA(ICONT,NPHASE)=PI

       J=1
       DO 40 I=2,NPHASE-1
        R = FLOAT(I-1)/XR
45      IF(R.GE.XTEMP(J).AND.R.LT.XTEMP(J+1))THEN
         XF = (R-XTEMP(J))/(XTEMP(J+1)-XTEMP(J))
         THETA(ICONT,I) = DTH*(J-1+XF)
        ELSE
         J=J+1
         GOTO 45
        ENDIF
40     CONTINUE

999   CONTINUE

C      open(12,file='phasprob.dat',status='unknown')
C      write(12,*)ncont,nphase
C      do i=1,ncont
C       write(12,*),'Type: ',i
C       write(12,*),(xhg(i,j),j=1,3)
C       do j=1,nphase
C        write(12,*)j,theta(i,j)
C       enddo
C      enddo 
C      close(12)

      RETURN

      END
   

