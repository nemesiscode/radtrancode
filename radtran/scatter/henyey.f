      real function henyey(alpha,f,g1,g2)
C     *****************************************************************
C     Function to calculate the phase function at alpha = cos(theta) using
C     the combined H-G function.
C
C     Input variables:
C	alpha	real	cos(scattering angle)
C	f	real	Fraction of forward scattering peak
C	g1	real	assymetry factor of forward peak (-1 < g1 < 1)
C	g2	real	assymetry factor of backward peak (-1 < g2 < 1)
C
C     Output variable
C	henyey	real	Phase function
C
C     Pat Irwin		12/5/97
C
C     *****************************************************************

      real alpha,f,g1,g2,x1,x2,pi

      pi=3.1415927

      x1 = (1.0-g1**2)/((1.0+g1**2 - 2*g1*alpha)**1.5)
      x2 = (1.0-g2**2)/((1.0+g2**2 - 2*g2*alpha)**1.5)

      y = f*x1 + (1.0-f)*x2

C      henyey = y/(4*pi)
      henyey = y
      return
      end
