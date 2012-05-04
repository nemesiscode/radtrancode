C***********************************************************************
C
C  PROGRAM        CHICO2    SUBROUTINE          
C
C  PURPOSE        COMPUTE CO2 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS AND LINE MIXING FOR VENUS
C		  NEAR-INFRARED WINDOWSC
C
C  VERSION        1.0	C.C.C. Tsang	09/03/2005
C
C  DESCRIPTION    THIS ROUTINE CALCULATES THE CHI FACTOR MODIFIED
C                 VOIGT FUNCTION OVER A FREQUENCY MESH DUE TO
C                 A LINE CENTRED AT SWN CM-1 USING TEMPERTAURE INDEPENDENT
C		  CHI FORM FACTOR FOR CO2 USING ONE OF A NUMBER OF DIFFERNET
C		  CORRECTION FACTORS BASED ON THE 2.7UM BAND HEAD.
C 		
C		  Ref:
C		  a)	Collard A., 1993, Thermal Emission from the Nightside
C			of Venus, DPhil. Thesis
C		  
C		  b)	Pollack et. al., 1993, Near-Infrared Light from
C			Venus' Nightside: A Spectroscopic Analysis, Icarus,
C			103, 1-42
C		
C		  c)	Meadows V.S and Crisp D., 1996, Ground-based near-
C			infrared observations of the Venus nightside: The
C			thermal structure and water abundance near the surface,
C			Journal of Geophysical Research, 101, 4595-4622
C			!TEMPERATURE DEPENDENT!
C
C		  d)	Tonkov M.V. et al., 1996, Measurements and empirical
C			modeling of pure CO2 absorption in the 2.3um region at
C			room temperature: far wings, allowed and collision-
C			induced bands, Applied Optics, 35, 24, 4863-4870,
C			CHI VALUES FROM TONKOV HAS BEEN USED HERE(09/03/05)
C
C  ARGUMENTS      NUM    I*4 I/P NUMBER OF FREQUENCY INTERVALS
C                 FF     R*4 I/P FINE WAVENUMBER GRID [cm-1]
C                 SWN    R*4 I/P WAVENUMBER OF LINE CENTRE [cm-1] 
C                 STRPAR R*4 I/P PATH ADJUSTED LINE STRENGTH [cm-1]
C                 WIDPAR R*4 I/P PATH ADJUSTED LORENTZ HALFWIDTH [cm-1] 
C                 DOPWID R*4 I/P DOPPLER LINE WIDTH [cm-1]
C                 YMIX   R*4 I/P LINE MIXING Y COEFFICIENT
C                 T      R*4 I/P PATH TEMPERATURE [K]
C                 XABS   R*4 O/P FREQUENCY MESH POINT ABSORPTIONS
C  
C SUBROUTINES    VOIVEC
C***********************************************************************
C
       
       FUNCTION CHICO2(X,Y,DV1)
       IMPLICIT NONE

       REAL X,Y,DV1,CHI,HUMLIC0,CHICO2,HUMLIC

        HUMLIC0=HUMLIC(X,Y)

        IF (DV1 .LT. 3.0) THEN
           CHI = 1.0
	   CHICO2=CHI*HUMLIC0
        ELSE  
	
        IF (DV1 .GE. 3.0 .AND. DV1 .LT. 150.0) THEN
	   CHI = 1.084*EXP(-0.027*DV1)
	   CHICO2=CHI*HUMLIC0
	ELSE
	   
	IF (DV1 .GE. 150.0 .AND. DV1 .LT. 300.0) THEN
	   CHI = 0.208*EXP(-0.016*DV1)
	   CHICO2=CHI*HUMLIC0
	ELSE
	 
	IF (DV1 .GE. 300.0) THEN
	   CHI = 0.025*EXP(-0.009*DV1)
	   CHICO2=CHI*HUMLIC0
	   
	   ENDIF
	  ENDIF
         ENDIF
	ENDIF

       RETURN                                              
       END
       
