C***********************************************************************
C
C  PROGRAM        BAILLYNH3    SUBROUTINE          
C
C  PURPOSE        COMPUTE NH3 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS
C
C		  Modified from Bailly (2004) paper
C
C  VERSION        1.0	R.S. GILES	06/06/2014
C
C***********************************************************************
      
      REAL FUNCTION BAILLYNH3(X,Y,DV)
      IMPLICIT NONE

      REAL X,Y,DV,DV1,CHI,HUMLIC0,HUMLIC

      HUMLIC0=HUMLIC(X,Y)
      DV1=ABS(DV)

      IF (DV1.LT.10.0) THEN
           CHI = 1.0*EXP(-DV*0.00243)
      ENDIF
      IF (DV1.GE.10.0.AND.DV1.LT.65.0) THEN
           CHI = 1.42*EXP(-DV1*0.0354)*EXP(-DV*0.00243)
      ENDIF
      IF (DV1.GE.65.0) THEN
	   CHI = 0.172*EXP(-DV1*0.00296)*EXP(-DV*0.00243)
      ENDIF

      BAILLYNH3=CHI*HUMLIC0

      RETURN                                              
      END
      
