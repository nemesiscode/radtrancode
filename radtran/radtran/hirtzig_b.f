C***********************************************************************
C
C  PROGRAM        HIRTZIG_B    SUBROUTINE          
C
C  PURPOSE        COMPUTE CH4 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS
C
C		  Alternative modification to CH4 lineshape recommended by HIRTZIG (2013)
C
C  VERSION        1.0	M. T. Roman	27-JULY-2023
C
C***********************************************************************
       
      REAL FUNCTION HIRTZIG_B(X,Y,DV1)
      IMPLICIT NONE

      REAL X,Y,DV1,CHI,HUMLIC0,HUMLIC

      HUMLIC0=HUMLIC(X,Y)

      IF (DV1.LT.25.6) THEN
           CHI = 1.0
      ELSE
       IF(DV1.LT.60.)THEN
           CHI = 1.2378*EXP(-DV1/120.)
       ELSE
           CHI = 0.757046/0.367879*EXP(-DV1/60.)
       ENDIF
      ENDIF

      HIRTZIG_B=CHI*HUMLIC0

      RETURN                                              
      END
      
