      SUBROUTINE HITSPHERE(PVEC,DVEC,RADIUS,PVEC1)
C     ***************************************************************
C     Find point of intersection between a line with R=PVEC+LAMBDA*DVEC
C     and a sphere of radius = RADIUS
C     Input variables
C	PVEC(3)	REAL	Starting position of line
C	DVEC(3) REAL	Directional unit vector of line
C	RADIUS	REAL	Radius of sphere
C
C     Output variables
C	PVEC1(3) REAL	Position vector of point of intersection
C
C     Pat Irwin	22/3/13
C
C     ***************************************************************

      REAL PVEC(3),DVEC(3),PVEC1(3),RADIUS,A,B,C
      REAL PVEC2(3),LAMBDA,LAMBDA1,LAMBDA2
      INTEGER I

      A=DVEC(1)**2 + DVEC(2)**2 + DVEC(3)**2
      B=2.0*(PVEC(1)*DVEC(1)+PVEC(2)*DVEC(2)+PVEC(3)*DVEC(3))
      C=PVEC(1)**2+PVEC(2)**2+PVEC(3)**2-RADIUS**2

C     There are two possible intercepts between a line and a sphere. Need to find the right one
      LAMBDA1 = (-B+SQRT(B**2-4*A*C))/(2*A)
      LAMBDA2 = (-B-SQRT(B**2-4*A*C))/(2*A)

      IF(LAMBDA1.LT.0.0.OR.LAMBDA2.LT.0.0)THEN
C      If either LAMBDA1 or LAMBDA2 is negative then this is in the opposite direction to the 
C      trajectory of the photon. Hence we need to take the one on the right direction.
       IF(LAMBDA1.GT.0) THEN
        LAMBDA = LAMBDA1
       ELSE
        LAMBDA = LAMBDA2
       ENDIF
      ELSE
C      If both LAMBDA1 and LAMBDA2 are positive, then we need to take the closest one.
       IF(LAMBDA1.LT.LAMBDA2)THEN
        LAMBDA=LAMBDA1
       ELSE
        LAMBDA=LAMBDA2
       ENDIF
      ENDIF

      DO I=1,3
       PVEC1(I)=PVEC(I)+DVEC(I)*LAMBDA
      ENDDO

      RETURN

      END

