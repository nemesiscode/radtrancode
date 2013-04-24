      REAL FUNCTION CALCMOLWT(NGAS,VMR,ID,ISO)
C     ****************************************************************
C     Subroutine to calculate molecular weight at a level in the atmosphere
C     (assuming the vmrs) add up to 1.0.
C     Code assumes the GASFIL data has already been read in elsewhere by
C     RDGAS subroutine
C
C     Input variables
C	NGAS		INTEGER	Number of gases
C	VMR(NGAS)	REAL	VMR of each gas
C	ID(NGAS)	INTEGER	Gas IDs
C	ISO(NGAS)	INTEGER	Gas Isotope Numbers
C
C     Output variable
C	CALCMOLWT	REAL	Molecular weight.	
C
C     Pat Irwin		Original	9/5/12
C
C     ****************************************************************
      IMPLICIT NONE
      INTEGER NGAS,ID(NGAS),ISO(NGAS),IGAS
      REAL VMR(NGAS),GETMASS,SUM,XTEST,TEST,XXMASS
      PARAMETER(XTEST=0.001)

      XXMASS=0.0
      SUM=0.0
   
      DO 100 IGAS=1,NGAS
        XXMASS=XXMASS+GETMASS(ID(IGAS),ISO(IGAS))*VMR(IGAS)
        SUM=SUM+VMR(IGAS)
100   CONTINUE

      TEST=ABS(SUM-1.0)

      IF(TEST.GT.XTEST)THEN
       PRINT*,'Warning in CALCMOLWT - Sum of vmrs <> 1.0',SUM
      ENDIF
      
      CALCMOLWT = XXMASS

      RETURN

      END 
