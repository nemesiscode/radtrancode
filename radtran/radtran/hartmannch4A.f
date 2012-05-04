C***********************************************************************
C
C  PROGRAM        HARTMANNCH4A    SUBROUTINE          
C
C  PURPOSE        COMPUTE CH4 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS
C
C		  Original coeffiecients recommended by Hartmann (2002)
C
C  VERSION        1.0	P.G.J. Irwin	21/06/2011
C
C***********************************************************************
       
      REAL FUNCTION HARTMANNCH4A(X,Y,DV1)
      IMPLICIT NONE

      REAL X,Y,DV1,CHI,HUMLIC0,HUMLIC

      HUMLIC0=HUMLIC(X,Y)

      IF (DV1.LT.26.) THEN
           CHI = 1.0
      ELSE 
       IF(DV1.LT.60.)THEN
           CHI = 8.72*EXP(-DV1/12.)
       ELSE
           CHI = 0.0684*EXP(-DV1/393.)
       ENDIF
      ENDIF

      HARTMANNCH4A=CHI*HUMLIC0

      RETURN                                              
      END
      
