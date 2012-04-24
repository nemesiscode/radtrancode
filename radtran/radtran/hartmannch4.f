C***********************************************************************
C
C  PROGRAM        HARTMANNCH4    SUBROUTINE          
C
C  PURPOSE        COMPUTE CH4 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS
C
C		  Modified from raw Hartmann (2002) paper recommendations
C		  by C. de Bergh for Titan.
C
C  VERSION        1.0	P.G.J. Irwin	21/06/2011
C
C***********************************************************************
       
      REAL FUNCTION HARTMANNCH4(X,Y,DV1)

      REAL X,Y,DV1,CHI,HUMLIC0,HUMLIC

      HUMLIC0=HUMLIC(X,Y)

      IF (DV1.LT.25.6) THEN
           CHI = 1.0
      ELSE 
           CHI = 1.1757*EXP(-DV1/158.4)
      ENDIF

      HARTMANNCH4=CHI*HUMLIC0

      RETURN                                              
      END
      
