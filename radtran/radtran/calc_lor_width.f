      REAL FUNCTION CALC_LOR_WIDTH(ABSCO,AL,LOR_MOD,RODGERS)
C     $Id:
C     ***********************************************************
C     Calculates the Lorentz width of a line assuming different
C     assumptions.
C
C     Pat Irwin	?/?/??	Original
C     Pat Irwin	26/4/12	Commented.
C
C     ***********************************************************
      REAL ABSCO,AL,PI
      DOUBLE PRECISION RODGERS(4,0:7),SUM,Z
      PARAMETER (PI=3.1415927)
      INTEGER LOR_MOD

      IF(LOR_MOD.EQ.0)THEN
C     Weak equivalent width
       CALC_LOR_WIDTH = ABSCO

      ELSEIF(LOR_MOD.EQ.1)THEN
C     Strong Lorentz equivalent width
       CALC_LOR_WIDTH = 2*SQRT(ABS(ABSCO*AL))
   
      ELSE
C     Combined Lorentz equivalent width using Rodgers+Williams data
       Z = DBLE(ABSCO/(2*PI*AL))
       IF(Z.LE.0.)THEN
        CALC_LOR_WIDTH = 0.
        print*,'Error. Calc_Lor_Width z.le.0'
       ELSE
        SUM=0.
        IF(Z.LE.1.98)THEN
         DO 10 N=0,6
          SUM=SUM+RODGERS(1,N)*Z**N
10       CONTINUE
         SUM=SUM*Z
        ELSE
         DO 20 N=0,5
          SUM=SUM+RODGERS(2,N)*Z**(-N)
20       CONTINUE
         SUM=SUM*DSQRT(Z)
        ENDIF
        CALC_LOR_WIDTH = 2*PI*AL*SNGL(SUM)
       END IF
      END IF

      RETURN

      END


C      ELSEIF(LOR_MOD.EQ.2)THEN
C     Combined Lorentz equivalent width
C       CALC_LOR_WIDTH = ABSCO * ( 1. + (0.25*ABSCO/AL)**1.25 )**(-0.4)

