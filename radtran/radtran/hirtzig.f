C***********************************************************************
C
C  PROGRAM        HIRTZIG    SUBROUTINE          
C
C  PURPOSE        COMPUTE CH4 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS
C
C		  Modification to CH4 lineshape recommended by HIRTZIG (2013)
C
C  VERSION        1.0	M. T. Roman	27-JULY-2023
C
C***********************************************************************
       
      REAL FUNCTION HIRTZIG(X,Y,DV1)
      IMPLICIT NONE

      REAL X,Y,DV1,CHI,HUMLIC0,HUMLIC

      HUMLIC0=HUMLIC(X,Y)

      IF (DV1.LT.25.6) THEN
           CHI = 1.0
      ELSE 
           CHI = 1.2378*EXP(-DV1/120.)
      ENDIF

      HIRTZIG=CHI*HUMLIC0

      RETURN                                              
      END
      
