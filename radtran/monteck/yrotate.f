      SUBROUTINE YROTATE(DVEC,THETA)
C     ***************************************************************
C     Routine to rotate a vector DVEC and angle THETA about the y-axis.
C
C     Input variables
C	DVEC(3)	REAL 	Input vector
C	THETA	REAL	Required angle of rotation
C
C     Outut variable
C	DVEC(3)	REAL	Rotated vector
C
C     Pat Irwin		Original	1/8/05
C
C     ***************************************************************

      REAL DVEC(3),THETA,X,Y,Z,PI
      PARAMETER (PI=3.1415927)

      DTR = PI/180.0

      X=DVEC(1)
      Y=DVEC(2)
      Z=DVEC(3)

      DVEC(1) = X*COS(THETA*DTR)+Z*SIN(THETA*DTR)
      DVEC(2) = Y
      DVEC(3) = Z*COS(THETA*DTR)-X*SIN(THETA*DTR)

      RETURN

      END
