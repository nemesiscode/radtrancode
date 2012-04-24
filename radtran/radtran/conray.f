      SUBROUTINE CONRAY(IRAY,V0,DV,IORDP1,PRESS,TEMP,TOTAM,POLY)
C     *****************************************************************
C     Subroutine to calculate Rayleigh scattering optical depth in
C     Jovian-type atmosphere
C
C     Pat Irwin 	27/7/99
C
C     *****************************************************************
      INTEGER IORDP1,IP,NX,IORDP,IRAY
      PARAMETER (IORDP = 3)
      REAL V,DV,PRESS,TEMP,TOTAM,POLY(IORDP1)
      REAL XX(IORDP),YY(IORDP),XMIN,XMAX

      IF(IORDP.NE.IORDP1)THEN
       PRINT*,'CONRAY: Error. IORDP <> IORDP1 ',IORDP,IORDP1
      ENDIF

      DO 10 IP=1,IORDP1

          FF = V0 + FLOAT(IP-1)*0.5*DV
          IF(IRAY.EQ.1)THEN
           POLY(IP) = TOTAM*RAYLEIGHJ(FF,PRESS,TEMP)*1E20
          ELSEIF(IRAY.EQ.2)THEN
           POLY(IP) = TOTAM*RAYLEIGHV(FF,PRESS,TEMP)*1E20
          ELSE
           POLY(IP) = TOTAM*RAYLEIGHA(FF,PRESS,TEMP)*1E20
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
              
      CALL CALC_PCOEFF(NX,YY,XX,XMIN,XMAX,IORDP1,POLY)

      RETURN

      END
