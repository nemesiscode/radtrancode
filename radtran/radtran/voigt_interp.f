      SUBROUTINE VOIGT_INTERP(ABSCO,AL,AD,EQW,RODGERS)
C     $Id: voigt_interp.f,v 1.2 2011-06-17 15:40:28 irwin Exp $
C     *********************************************************************
C
C     Interpolates the VOIGT Equivalent width look-up table 
C
C     Pat Irwin		16/2/94
C
C     *********************************************************************
      IMPLICIT NONE
      REAL ABSCO,AL,AD,EQW,VWIDTH(117,161),CONV,X,Y,RX,RY
      REAL y1,y2,y3,y4,t,u,VW,CALC_LOR_WIDTH,CALC_DOP_WIDTH
      REAL WL,WD,eqw1
      DOUBLE PRECISION RODGERS(4,0:7)
      INTEGER I,J

      PARAMETER(CONV=2.3025851)
      COMMON /VOIGTDATA/VWIDTH

C     CALCULATE ORDINATES

      X=LOG(AL/AD)/CONV
      Y=LOG(ABSCO/AD)/CONV

      RX=1.+(X+20.0)/0.25
      RY=1.+(Y+20.0)/0.25

      IF(RX.GE.0)THEN
       I=INT(RX)
      ELSE
       I=INT(RX-1.)
      END IF
      IF(RY.GE.0)THEN
       J=INT(RY)
      ELSE
       J=INT(RY-1.)
      END IF

      IF(RX.GE.1.AND.RX.LT.117.AND.RY.GE.1.AND.RY.LT.161)THEN

C      Ordinates are within the table and may be interpolated
       y1=VWIDTH(I,J)
       y2=VWIDTH(I+1,J)
       y3=VWIDTH(I+1,J+1)  
       y4=VWIDTH(I,J+1)
       t=RX-I
       u=RY-J
       VW=(1.0-t)*(1.0-u)*y1 + t*(1.0-u)*y2 +t*u*y3 + (1.0-t)*u*y4    
       IF(VW.NE.-999.)THEN
         EQW=AD*EXP(VW)
       ELSE
         WL = CALC_LOR_WIDTH(ABSCO,AL,2,RODGERS)
         WD = CALC_DOP_WIDTH(ABSCO,AD,2,RODGERS)
         EQW = SQRT(WL*WL + WD*WD - (WL*WD/ABSCO)**2)
       END IF

       IF(EQW.EQ.0.0)THEN
         WL = CALC_LOR_WIDTH(ABSCO,AL,2,RODGERS)
         WD = CALC_DOP_WIDTH(ABSCO,AD,2,RODGERS)
         EQW = SQRT(WL*WL + WD*WD - (WL*WD/ABSCO)**2)
       END IF

      ELSE

       PRINT*,'VOIGT_INTERP : PARAMETERS OUT OF BOUNDS'
       PRINT*,'AL,AD,ABSCO,RX,RY'
       PRINT*,AL,AD,ABSCO,RX,RY

C      CALCULATE MIXED WIDTH USING RODGERS APPROXIMATIONS

       WL = CALC_LOR_WIDTH(ABSCO,AL,2,RODGERS)
       WD = CALC_DOP_WIDTH(ABSCO,AD,2,RODGERS)
       EQW = SQRT(WL*WL + WD*WD - (WL*WD/ABSCO)**2)
      END IF
      
      RETURN
      END
