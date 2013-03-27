      SUBROUTINE DEFLECT(IDUM,THETA,NPHASE,ALPHA,PHI)
C     *******************************************************************
C     Subroutine to calculate new vector of photon after being scattered
C     by a particle.
C
C     Input parameters:
C	IDUM		INTEGER	Initialization switch for random number
C				generator
C	THETA(NPHASE)	REAL	Tabulated scattering angles (radians)
C	NPHASE		INTEGER THETA array length
C
C     Output parameters
C	ALPHA		REAL	Forward scattering angle (radians)
C	PHI		REAL	Associated azimuth deflection (radians)
C
C     Pat Irwin		5/3/99
C			1/4/99
C
C     *******************************************************************

      IMPLICIT NONE
      INTEGER IDUM,I,NPHASE
      REAL RAN11,PI,PHI,THETA(NPHASE)
      REAL XR,XF,ALPHA,TMP
      PARAMETER (PI=3.1415927)
      
C     Calculate aziumuth part of scattering
      PHI=2*PI*RAN11(IDUM)

C     Calculate the zenith part of scattering
      XR = 1.0 + (NPHASE-1)*RAN11(IDUM)
      I=INT(XR)
      XF = XR-I
      ALPHA = (1.0-XF)*THETA(I) + XF*THETA(I+1)


      RETURN

      END

      
