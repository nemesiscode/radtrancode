      REAL FUNCTION LINECONTRIB(IPROC,IDGAS,VV,TCORDW,TCORS1,TCORS2,
     1 PRESS,TEMP,FRAC,VLIN,SLIN,ELIN,ALIN,SBLIN,TDW,TDWS,LLQ,DOUBV,
     2 FNH3,FH2)
C     ****************************************************************
C     Function to calculate the line contribution from multiple line
C     broadening types.
C
C     Input variables:
C	IPROC	INTEGER	Line processing parameter
C	IDGAS	INTEGER	Gas ID
C	VV	REAL	Calculation wavenumber
C	TCORDW	REAL	Doppler line width coefficient
C	TCORS1 	REAL	Temperature coefficient 1
C	TCORS2 	REAL	Temperature coefficient 2
C	PRESS	REAL	Pressure
C	TEMP	REAL	Temperature
C	FRAC	REAL	Fraction
C	VLIN	REAL	Line wavenumber
C	SLIN	REAL	Line strength at STP
C	ELIN	REAL 	Lower state energy
C	ALIN	REAL	Air-broadened width
C	SBLIN	REAL	Self-broadened width
C	TDW	REAL	Temperature dependance of air-broadening width
C	TDWS	REAL	Temperature dependance of self-broadening width
C	LLQ	CHARACTER*15 Lower state quantum numbers
C	FNH3	REAL	Fraction of NH3 (needed of models 3 and 6)
C	FH2	REAL	Fraction of H2 (needed of models 3 and 6)
C	
C     ****************************************************************
      IMPLICIT NONE
      INTEGER IPROC,IDGAS
      
      REAL SUBLINE,VV,VLIN,ABSCO,X,Y,AD,TRATIO,DOUBV
      REAL SLIN,ALIN,SBLIN,ELIN,TDW,TDWS,TCORDW,TCORS1,TCORS2
      REAL DV,TSTIM,TS1,TS2
      REAL PRESS,TEMP,FRAC,FNH3,FH2,WY,DPEXP,LNABSCO
      CHARACTER*15 LLQ

C     AD is the Doppler-broadened line width
C     Y is the Lorentz/collision-broadened line width divided by
C      the Doppler-broadened line width.

C     Doppler line width
      AD=TCORDW*VLIN

C     Stimulated emission coefficient.
      TS1 = (1.0 - DPEXP(-1.439*VLIN/TEMP))
      TS2 = (1.0 - DPEXP(-1.439*VLIN/296.0))
      TSTIM=1.0
      IF(TS2.NE.0.)TSTIM=TS1/TS2

      LNABSCO=LOG(SLIN)+LOG(TCORS1)+TCORS2*ELIN+LOG(TSTIM)
      ABSCO=EXP(LNABSCO)

      TRATIO = 296./TEMP
      Y = (ALIN*(1.-FRAC)*TRATIO**TDW+(ALIN-SBLIN)*FRAC*
     1 TRATIO**TDWS)*PRESS/AD

      DV = VV-VLIN
      X = DV/AD

C      print*,IDGAS,PRESS,TEMP,IPROC,VV,VLIN,ABSCO,X,Y,
C     1 FNH3,FH2,LLQ,DOUBV

      LINECONTRIB=SUBLINE(IDGAS,PRESS,TEMP,IPROC,VV,VLIN,ABSCO,X,Y,
     1 AD,FNH3,FH2,LLQ,DOUBV)

C      print*,'LINE = ',LINECONTRIB
      
      RETURN

      END
