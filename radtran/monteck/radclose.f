      REAL FUNCTION RADCLOSE(PVEC,DVEC)
C     ********************************************************
C     Simple routine to find the closest approach to the origin
C     of a line with initial position PVEC and direction DVEC
C     ********************************************************
      IMPLICIT NONE
      REAL PVEC(3),DVEC(3),TVEC(3),SUM1,SUM2
      REAL RMIN,LMIN,A,B,C
      INTEGER I

C     Need to see if photon is moving towards or away from the closest approach
C     to the origin

      SUM1=0.
      SUM2=0.
      DO I=1,3
        SUM1=SUM1+PVEC(I)**2
        SUM2=SUM2+(PVEC(I)+DVEC(I)*10.)**2
      ENDDO

      IF(SUM2.GT.SUM1)THEN
C      Photon is moving upwards and so cannot hit the surface

       RADCLOSE=-1.

      ELSE
C      Photon is moving downwards and so could potentially strike the surface
C      Find closest approach

C      R=PVEC+L*DVEC. Distance from origin = R^2. Differentiate
C      with respect to L and find value of L when distance is minumum.
C      Then calculate the minimum distance. 
       A=DVEC(1)**2 + DVEC(2)**2 + DVEC(3)**2
       B=2.0*(PVEC(1)*DVEC(1)+PVEC(2)*DVEC(2)+PVEC(3)*DVEC(3))
       C=PVEC(1)**2+PVEC(2)**2+PVEC(3)**2

       LMIN = -B/(2.*A)
       RMIN = A*LMIN**2 + B*LMIN+C

        RADCLOSE=SQRT(RMIN)

      ENDIF
 
      RETURN

      END
      


