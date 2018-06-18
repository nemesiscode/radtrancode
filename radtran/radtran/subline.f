      REAL FUNCTION SUBLINE(IDGAS,PRESS,TEMP,IPROC,VV,VLIN,ABSCO,X,Y,
     1 AD,FNH3,FH2,LLQ,DOUBV)
C     ****************************************************************
C     Function to calculate the line contribution from multiple line
C     broadening types.
C
C     Input variables:
C	IDGAS	INTEGER	Gas ID
C	PRESS	REAL	Pressure
C	TEMP	REAL	Temperature
C	IPROC	INTEGER	Line processing parameter
C	VV	REAL*8	Calculation wavenumber
C	VLIN	REAL*8	Line wavenumber
C	ABSCO	REAL	Computed line strength 
C	X	REAL	DV/AD
C	Y	REAL	AL/AD
C	AD	REAL 	AD i.e. Doppler width
C	FNH3	REAL	Fraction of NH3 (needed of models 3 and 6)
C	FH2	REAL	Fraction of H2 (needed of models 3 and 6)
C	LLQ	CHARACTER*15 Lower state quantum numbers
C	DOUBV	REAL	Distance to other doublet member (cm-1) for ammonia
C			lines only (models 3 and 6)


C	
C     ****************************************************************
      IMPLICIT NONE
      INTEGER IPROC,IDGAS
      
      REAL HUMLIC,GVOICO2,CHICO2,BURCHCO2,BURCHCO2_NORM
      REAL BAILLYNH3,HARTMANNCH4,HARTMANNCH4A,HARTMANNCH4B
      REAL ABSCO,X,Y,AD,TRATIO,DOUBV,GAMMA,VVWEISS
      REAL ALIN,SBLIN,ELIN,TDW,TDWS,TCORDW,TCORS1,TCORS2
      REAL NH2,NHE,YPH3_H2,YPH3_HE,YPH3,PI,DV,TSTIM,TS1,TS2
      REAL PRESS,TEMP,GETNH3,FNH3,FH2,WY,DPEXP,LNABSCO
      DOUBLE PRECISION VV,VLIN
      CHARACTER*15 LLQ

C     AD is the Doppler-broadened line width
C     Y is the Lorentz/collision-broadened line width divided by
C      the Doppler-broadened line width.
C     HUMLIC, GVOICO2 and PARTF are all programs that return a REAL value for
C      use in calculating the lineshapes.
C     NH2, NHe are the temperature exponents for H2 and He broadening from
C      Levy, A. 1994, J.Mol.Spec., 166, 20-31
C     YPH3_H2, YPH3_HE are the temperature-dependent halfwidths for PH3 due
C      to H2 and He collision-broadening.
C     YPH3 is the combined (added) temperature-dependent halfwidths for PH3
C      due to H2 and He collision-broadening.

      PARAMETER (BURCHCO2_NORM = 1.006956255,PI=3.1415927)
C     BURCHCO2_NORM = ratio of area under voigt / area under burch lineshape
C     used to renomalise the linestrengths

      TRATIO=296.0/TEMP
      DV=VV-VLIN


      IF(IPROC.EQ.14)THEN
C      IPROC=14 :: van Vleck-Weisskopf-Ben-Reuven NH3 lineshape
       IF(IDGAS.NE.11)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' van Vleck-Weisskopf-Ben-Reuven broadening'
        WRITE(*,*)' for gases other than NH3.'
        WRITE(*,*)' IDGAS = ',IDGAS
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF
       IF(FNH3.LT.0.OR.FH2.LT.0)THEN
C       Hard-wiring FNH3 and FH2 to typical Saturn conditions 
        FNH3=1e-4		
        FH2 = 0.881
       ENDIF
       WY = GETNH3(PRESS,TEMP,FNH3,FH2,LLQ,DOUBV)
       IF(WY.EQ.0)THEN
C       Default to Vleck-Weisskopf lineshape
        GAMMA=Y*AD
        SUBLINE = ABSCO*VVWEISS(GAMMA,VV,VLIN)
       ELSE 
        SUBLINE = (Y*AD + WY*DV)/((VV-VLIN)**2 + (Y*AD)**2)
        SUBLINE = SUBLINE + (Y*AD + WY*DV)/((VV+VLIN)**2 + (Y*AD)**2)
        SUBLINE = SUBLINE*ABSCO*(VV/VLIN)/PI
       ENDIF
      

      ELSE IF(IPROC.EQ.13)THEN
C      IPROC=13::Sub-lorentzian lineshape from Bailly et al. (2004)
       IF(IDGAS.NE.11)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' Bailly sublorentzian shape '
        WRITE(*,*)' for gases other than NH3.'
        WRITE(*,*)' IDGAS = ',IDGAS
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF
	SUBLINE = ABSCO*BAILLYNH3(ABS(X),Y,DV)/AD

      ELSE IF(IPROC.EQ.12)THEN
C      IPROC=12::Doppler broadening only
       SUBLINE = ABSCO*EXP(-X**2)/(AD*SQRT(PI))

      ELSE IF(IPROC.EQ.11)THEN
C      IPROC=11:: Modification to Hartmann sub-lorentzian lineshape for methane
C 		 in hydrogen atmosphere recommended by Sromovsky et al. (2012)

       IF(IDGAS.NE.6)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' Sromovsky CH4 Hartmann sublorentzian shape '
        WRITE(*,*)' for gases other than CH4.'
        WRITE(*,*)' IDGAS = ',IDGAS
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF

       SUBLINE = ABSCO*
     1                   HARTMANNCH4B(ABS(X),Y,ABS(DV))/AD


      ELSE IF(IPROC.EQ.10)THEN
C      IPROC=10:: Original Hartmann sub-lorentzian lineshape for methane
C 		 in hydrogen atmosphere
       IF(IDGAS.NE.6)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' Original CH4 Hartmann sublorentzian shape '
        WRITE(*,*)' for gases other than CH4.'
        WRITE(*,*)' IDGAS = ',IDGAS
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF

       SUBLINE = ABSCO*
     1                   HARTMANNCH4A(ABS(X),Y,ABS(DV))/AD


      ELSE IF(IPROC.EQ.9)THEN
C      IPROC=9:: Hartmann sub-lorentzian lineshape for methane modified
C      for Titan by C. de Bergh
       IF(IDGAS.NE.6)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' C. de Bergh CH4 Hartmann sublorentzian shape '
        WRITE(*,*)' for gases other than CH4.'
        WRITE(*,*)' IDGAS = ',IDGAS
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF

       SUBLINE = ABSCO*
     1                   HARTMANNCH4(ABS(X),Y,ABS(DV))/AD
       

      ELSEIF(IPROC.EQ.81)THEN
C      IPROC=81:: Burch CO2 sub-lorentzian line shape (renormalised)
       IF(IDGAS.NE.2)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' CO2 Burch sublorentzian shape '
        WRITE(*,*)' for gases other than CO2.'
        WRITE(*,*)' IDGAS = ',IDGAS
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF
       SUBLINE = BURCHCO2_NORM*ABSCO*
     1                     BURCHCO2(ABS(X),Y,ABS(DV))/AD

      ELSE IF(IPROC.EQ.8)THEN
C      IPROC=8:: Burch CO2 sub-lorentzian line shape
       IF(IDGAS.NE.2)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' CO2 Burch sublorentzian shape '
        WRITE(*,*)' for gases other than CO2.'
        WRITE(*,*)' IDGAS = ',IDGAS
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF

       SUBLINE = ABSCO*BURCHCO2(ABS(X),Y,ABS(DV))/AD

      ELSE IF(IPROC.EQ.7)THEN
C      IPROC=7 :: CO2 Sublorentzian line shape for Venus
       IF(IDGAS.NE.2)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' CO2 sublorentzian shape '
        WRITE(*,*)' for gases other than CO2.'
        WRITE(*,*)' IDGAS = ',IDGAS
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF
       SUBLINE = ABSCO*CHICO2(ABS(X),Y,ABS(DV))/AD

      ELSE IF(IPROC.EQ.6)THEN
C      IPROC=6 :: Rosenkrantz-Ben-Reuven NH3 lineshape
C      same as IPROC=3 except default is voigt lineshape
       IF(IDGAS.NE.11)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' Rosenkrantz-Ben-Reuven broadening'
        WRITE(*,*)' for gases other than NH3.'
        WRITE(*,*)' IDGAS = ',IDGAS
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF
       IF(FNH3.LT.0.OR.FH2.LT.0)THEN
C       Hard-wiring FNH3 and FH2 to typical Saturn conditions 
        FNH3=1e-4		
        FH2 = 0.881
       ENDIF
       WY = GETNH3(PRESS,TEMP,FNH3,FH2,LLQ,DOUBV)
       IF(WY.EQ.0)THEN
C       IPROC=0 :: Default to Voigt lineshape
        SUBLINE = ABSCO*HUMLIC(ABS(X),Y)/AD
       ELSE 
        SUBLINE = ABSCO*(Y*AD + WY*DV)/(PI*(DV**2 + (Y*AD)**2))
       ENDIF

      ELSE IF(IPROC.EQ.5)THEN
C      IPROC=5 :: Levy etal 1994 PH3 broadening
       IF(IDGAS.NE.28)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the Levy, A.'
        WRITE(*,*)' etal (1994) broadening for gases'
        WRITE(*,*)' other than PH3.'
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF
       NH2 = 0.732  ! H2, He broadening temperature exponents
       NHE = 0.303  ! (Levy, A. 1994,J.Mol.Spec.,166,20-31).
       YPH3_H2 = ALIN*PRESS*(TRATIO**NH2)
       YPH3_HE = ALIN*PRESS*(TRATIO**NHe)
C      Multiply the broadening linewidths by Jovian H2 and He abundances
C      (0.863, 0.134) respectively
       YPH3 = (0.863*YPH3_H2) + (0.134*YPH3_HE)
       SUBLINE = ABSCO*YPH3/(PI*(DV**2 + (YPH3)**2))

      ELSE IF(IPROC.EQ.4)THEN
C      IPROC=4 :: Lorentz broadening only
       SUBLINE = ABSCO*Y*AD/(PI*(DV**2 + (Y*AD)**2))
C       print*,'L'
      ELSE IF(IPROC.EQ.3)THEN
C      IPROC=3 :: Rosenkrantz-Ben-Reuven NH3 lineshape
       IF(IDGAS.NE.11)THEN
        WRITE(*,*)' SUBLINE.f :: Can_t use the'
        WRITE(*,*)' Rosenkrantz-Ben-Reuven broadening'
        WRITE(*,*)' for gases other than NH3.'
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF
       IF(FNH3.LT.0.OR.FH2.LT.0)THEN
C       Hard-wiring FNH3 and FH2 to typical Saturn conditions 
        FNH3=1e-4		
        FH2 = 0.881
       ENDIF
       WY = GETNH3(PRESS,TEMP,FNH3,FH2,LLQ,DOUBV)
       IF(WY.EQ.0)THEN
        SUBLINE = ABSCO*GVOICO2(ABS(X),Y,TEMP,AD)/AD
       ELSE 
        SUBLINE = ABSCO*(Y*AD + WY*DV)/(PI*(DV**2 + (Y*AD)**2))
       ENDIF

      ELSE IF(IPROC.EQ.2)THEN
C      IPROC=2 :: VanVleck-Weisskopf lineshape
       GAMMA = Y*AD
       SUBLINE = ABSCO*VVWEISS(GAMMA,VV,VLIN)

      ELSE IF(IPROC.EQ.1)THEN
C      IPROC=1 :: Sub-Lorentzian Lineshape
       IF(IDGAS.NE.2)THEN
        WRITE(*,*)' SUBLINE.f :: Can_T use this sub-'
        WRITE(*,*)' Lorentzian broadening for gases'
        WRITE(*,*)' other than CO2.'
        WRITE(*,*)' '
        WRITE(*,*)' Stopping program.'
        STOP
       ENDIF
       SUBLINE = ABSCO*GVOICO2(ABS(X),Y,TEMP,AD)/AD

      ELSE IF(IPROC.EQ.0)THEN  
C      IPROC=0 :: Default to Voigt lineshape
C       SUBLINE = ABSCO*HUMLIC(ABS(X),Y)/AD
       SUBLINE = ABSCO*HUMLIC(X,Y)/AD

      ELSE
       PRINT*,'Line processing code not recognised in subline.f'
       PRINT*,'Requested IPROC = ',IPROC
       STOP
      ENDIF


      RETURN

      END
