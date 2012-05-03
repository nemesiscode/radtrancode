      SUBROUTINE MODVEC(VEC,ALPHA,PHI)
C     *******************************************************************
C     Subroutine to calculate new vector of photon after being scattered
C     by a particle.
C
C     Input parameters:
C	VEC(3)		REAL	Initial vector
C	ALPHA		REAL	Forward scattering angle
C	PHI		REAL	Associated azimuth deflection
C     Output parameters
C	VEC(3)		REAL	Modified vector
C
C     Pat Irwin		5/3/99
C
C     *******************************************************************

      IMPLICIT NONE
      INTEGER I
      REAL VEC(3),OUTVEC(3),PHI,ALPHA,XT,PI,YT,ARCTAN
      REAL TVEC(3),R,THIN,PHIN,DTR
      PARAMETER (DTR=1.7453293E-2,PI=3.1415927)
      

      R=0.0
      DO 10 I=1,3
       R = R+VEC(I)**2
10    CONTINUE
      R = SQRT(R)

C     For input vector, calculate vector in spherical polar coordinates
C     R, THIN, PHIN
      THIN = ACOS(VEC(3)/R)
      IF(THIN.EQ.0.0)THEN
       PHIN = 0.0
      ELSE
       XT = VEC(1)
       YT = VEC(2)
       PHIN = ARCTAN(YT,XT)
      ENDIF

C      print*,'R,THETA,PHIN (deg) = ',R,THIN/DTR,PHIN/DTR
C      print*,'alpha,phi (deg) : ',ALPHA/DTR,PHI/DTR

      OUTVEC(1)=R*SIN(ALPHA)*COS(PHI)
      OUTVEC(2)=R*SIN(ALPHA)*SIN(PHI)
      OUTVEC(3)=R*COS(ALPHA)

      TVEC(1)=COS(THIN)*OUTVEC(1)+SIN(THIN)*OUTVEC(3)
      TVEC(2)=OUTVEC(2)
      TVEC(3)=-SIN(THIN)*OUTVEC(1)+COS(THIN)*OUTVEC(3)

      VEC(1)=COS(PHIN)*TVEC(1)-SIN(PHIN)*TVEC(2)
      VEC(2)=SIN(PHIN)*TVEC(1)+COS(PHIN)*TVEC(2)
      VEC(3)=TVEC(3)

      RETURN
      END
