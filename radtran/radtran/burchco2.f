C***********************************************************************
C
C  PROGRAM        BURCHCO2    SUBROUTINE          
C
C  PURPOSE        COMPUTE CO2 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS AND LINE MIXING FOR VENUS
C		  4.3 micron self-broadening re: Burch et al.
C
C  VERSION        1.0	P.G.J. Irwin	17/07/2005
C                 1.1	N Teanby	26/9/06	added polynomial for
C					3 to 6 cm-1 to remove discontinuity.
C
C		  Ref:
C		  a)	Burch, D. E. et al., 1969, Absorption of 
C                       infrared radiant energy by CO2 and H2O. IV. 
C                       Shapes of collision-broadened CO2 lines. 
C                       J. Opt. Soc. America, 59, 267-280
C
C***********************************************************************
C
       
      FUNCTION BURCHCO2(X,Y,DV1)
      IMPLICIT NONE

      REAL X,Y,DV1,CHI,HUMLIC0,BURCHCO2,HUMLIC

      HUMLIC0=HUMLIC(X,Y)

      IF (DV1 .LT. 3.0) THEN
           CHI = 1.0
      ELSE IF (DV1 .GE. 3.0 .AND. DV1 .LT. 6.0) THEN
         CHI = -1.40253 +1.92162*DV1 -0.479585*DV1**2 +0.0353706*DV1**3
	ELSE IF (DV1 .GE. 6.0 .AND. DV1 .LT. 46.0) THEN
	   CHI = 10**(-0.233-1.087*DV1/100.0)
	ELSE IF (DV1 .GE. 46.0 .AND. DV1 .LE. 136.0) THEN
	   CHI = 10**(-0.133-1.3*DV1/100.0)
	ELSE IF (DV1 .GT. 136.0) THEN
	   CHI = 10**(-0.7-0.88*DV1/100.0)
	ENDIF

      BURCHCO2=CHI*HUMLIC0

      RETURN                                              
      END
      
