C***********************************************************************
C
C  PROGRAM        HARTMANNCH4B    SUBROUTINE          
C
C  PURPOSE        COMPUTE CH4 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS
C
C		  Modification to CH4 lineshape recommended by Sromovsky (2012)
C
C  VERSION        1.0	P.G.J. Irwin	07/06/2012
C
C***********************************************************************
       
      REAL FUNCTION HARTMANNCH4B(X,Y,DV1)
      IMPLICIT NONE

      REAL X,Y,DV1,CHI,HUMLIC0,HUMLIC

      HUMLIC0=HUMLIC(X,Y)

      IF (DV1.LT.26.) THEN
           CHI = 1.0
      ELSE 
       IF(DV1.LT.60.)THEN
           CHI = EXP(1.0-DV1/26.)
       ELSE
           CHI = 0.335*EXP(-DV1/280.)
       ENDIF
      ENDIF

      HARTMANNCH4B=CHI*HUMLIC0

      RETURN                                              
      END
      
