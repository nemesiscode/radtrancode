      REAL FUNCTION CALCALT(PVEC,RADIUS)
C     ****************************************************************
C     Function to find the altitude of a point (x,y,z) in the atmosphere 
C     of a spherical planet.
C
C     Input variables
C	PVEC(3)	REAL	Position vector (x,y,z) (km)
C	RADIUS  REAL    Plantetary radius (km)
C
C     Output variable
C	CALCALT REAL	Altitude (km)
C
C     Documented	Pat Irwin	1/1/05
C
C     ****************************************************************

      REAL PVEC(3),RADIUS,ALTITUDE
      INTEGER I

      ALTITUDE = 0.0
      DO I=1,3 
        ALTITUDE = ALTITUDE + PVEC(I)**2
      ENDDO

      CALCALT = SQRT(ALTITUDE) - RADIUS

      RETURN

      END
