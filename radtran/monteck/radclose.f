      REAL FUNCTION RADCLOSE(PVEC,DVEC)
C     ********************************************************
C     Simple routine to find the closest approach to the origin
C     of a line with initial position PVEC and direction DVEC
C     ********************************************************
      IMPLICIT NONE
      REAL PVEC(3),DVEC(3),TVEC(3),SUM1,SUM2,SUM3,SUM4
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
       TVEC(1)=PVEC(2)*DVEC(3)-PVEC(3)*DVEC(2)
       TVEC(2)=PVEC(3)*DVEC(1)-PVEC(1)*DVEC(3)
       TVEC(3)=PVEC(1)*DVEC(2)-PVEC(2)*DVEC(1)

       SUM3=0.
       SUM4=0.
       DO I=1,3
         SUM3=SUM3+TVEC(I)**2
         SUM4=SUM4+DVEC(I)**2
       ENDDO

       RADCLOSE=SQRT(SUM3/SUM4)

      ENDIF
 


      RETURN

      END
      


