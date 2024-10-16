      SUBROUTINE CONRAY(IRAY,V0,DV,PRESS,TEMP,fheh2,fch4h2,fnh3,
     & TOTAM,POLY)
C     *****************************************************************
C     Subroutine to calculate Rayleigh scattering optical depth in
C     Jovian-type atmosphere
C
C     Pat Irwin 	27/7/99
C
C     *****************************************************************
      INCLUDE '../includes/arrdef.f'
      INTEGER IP,NX,IORDP,IRAY
      REAL V,DV,PRESS,TEMP,TOTAM,POLY(IORDP1)
      REAL XX(IORDP1),YY(IORDP1),XMIN,XMAX
      REAL fheh2,fch4h2,fnh3

      DO 10 IP=1,IORDP1

          FF = V0 + FLOAT(IP-1)*0.5*DV
          IF(IRAY.EQ.1)THEN
           POLY(IP) = TOTAM*RAYLEIGHJ(FF,PRESS,TEMP)*1E20
          ELSEIF(IRAY.EQ.2)THEN
           POLY(IP) = TOTAM*RAYLEIGHV_IGNATIEV(FF,PRESS,TEMP)*1E20
          ELSEIF(IRAY.EQ.3)THEN
           POLY(IP) = TOTAM*RAYLEIGHA(FF,PRESS,TEMP)*1E20
          ELSE
           POLY(IP) = TOTAM*RAYLEIGHLS(FF,fheh2,fch4h2,fnh3)*1E20
          ENDIF
10    CONTINUE

C     Convert polynomial calculations to a polynomial fit.
  
      DO IP=1,IORDP1
         XX(IP)=0.5*DV*(IP-1)
         YY(IP)=POLY(IP) 
      END DO
      XMIN=XX(1)
      XMAX=XX(IORDP1)
      NX=IORDP1
              
      CALL CALC_PCOEFF(NX,YY,XX,XMIN,XMAX,POLY)

      RETURN

      END
