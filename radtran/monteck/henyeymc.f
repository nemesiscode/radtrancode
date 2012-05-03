      REAL FUNCTION HENYEYMC(THETA,XHG)
C     *****************************************************************
C     Function to calculate the phase function at alpha = cos(theta) using
C     the combined H-G function.
C
C     Input variables:
C	theta		real	scattering angle
C	xhg(1):f	real	Fraction of forward scattering peak
C	xhg(2):g1	real	assymetry factor of forward peak (-1 < g1 < 1)
C	xhg(3):g2	real	assymetry factor of backward peak (-1 < g2 < 1)
C
C     Output variable
C	henyeymc	real	Phase function
C
C     Pat Irwin		5/3/99
C
C     *****************************************************************

      IMPLICIT NONE
      REAL THETA,XHG(3),X1,X2,Y,PI
      PARAMETER (PI=3.1415927)
      REAL F,G1,G2,ALPHA


      F =XHG(1)
      G1=XHG(2)
      G2=XHG(3)

      ALPHA = COS(THETA)

      X1 = (1.0-G1**2)/((1.0+G1**2 - 2*G1*ALPHA)**1.5)
      X2 = (1.0-G2**2)/((1.0+G2**2 - 2*G2*ALPHA)**1.5)

      Y = F*X1 + (1.0-F)*X2

      HENYEYMC=Y

      RETURN
      END
