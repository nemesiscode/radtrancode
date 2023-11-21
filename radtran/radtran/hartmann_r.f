C***********************************************************************
C
C  PROGRAM        HARTMANN_R   SUBROUTINE          
C
C  PURPOSE        COMPUTE CH4 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS
C
C		  Modification to CH4 lineshape recommended by Hartmann, 
C                 somewhat arbitrarily modified for a sharper drop at
C                 lower distance from the line center
C                 
C
C  VERSION        1.0	M. Roman	09-NOV-2023
C
C***********************************************************************
       
      REAL FUNCTION HARTMANN_R(X,Y,DV1)
      IMPLICIT NONE

      REAL X,Y,DV1,CHI,HUMLIC0,HUMLIC

      HUMLIC0=HUMLIC(X,Y)

      IF (DV1.LT.26.) THEN
           CHI = 1.0
      ELSE 
       IF(DV1.LT.60.)THEN
           CHI = EXP(1.0-(DV1/26.)**2.)
       ELSE
           CHI = 0.019543294*EXP(-DV1/200.)
       ENDIF
      ENDIF

      HARTMANN_R=CHI*HUMLIC0

      RETURN                                              
      END
      
