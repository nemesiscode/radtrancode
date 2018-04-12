C     $Id: lincom.f,v 1.6 2010-11-01 11:22:44 irwin Exp $
C***********************************************************************
C_TITL: LINCOM
C
C_DESC:
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************
C-----------------------------------------------------------------------
C
C	Line data variables.
C
C-----------------------------------------------------------------------
      INTEGER MAXLIN
      PARAMETER (MAXLIN=2000000)
C MAXLIN: maximum number of lines which can be stored.

      INTEGER IDLIN(MAXLIN),LINMAX,NXTLIN
C IDLIN: Air Force Geospace Lab. identifier.

      DOUBLE PRECISION SLIN(MAXLIN),VLIN(MAXLIN)
      REAL ALIN(MAXLIN),PSHIFT(MAXLIN)
C VLIN: Line position [cm^-1].
C SLIN: Line strength [cm^-1 molecule^-1 cm^-2] at STP.
C ALIN: Lorentz half width [cm^-1] at STP.
C PSHIFT: line pressure-shift parameter.

      REAL ABSCO_ARR(MAXLIN),AD_ARR(MAXLIN),Y_ARR(MAXLIN)
      REAL TSTIM_ARR(MAXLIN)
C ABSCO_ARR: Calculated line strength at calculation T,P 
C AD_ARR:    Calculated doppler width at calculation T,P 
C Y_ARR:     Calculated AL/AD at calculation T,P 
C TSTIM_ARR: Calculated TSTIM coeff at calculation T,P

      REAL ELIN(MAXLIN),SBLIN(MAXLIN),TDW(MAXLIN),TDWS(MAXLIN)
C ELIN: Lower-state energy line-position [cm^-1].
C SBLIN: Line self-broadening coefficients. NOTE: the self broadening
C coefficient used in this program is the mixing gas 'air'-broadened
C half width minus the self-broadened half width.

      REAL DOUBV(MAXLIN)
C DOUBV: Inversion doublet separation (only applies to longwave NH3 lines
C at present.

      CHARACTER*15 LLQ(MAXLIN)
C LLQ: lower state local quanta. (15 elements for Hitran04, 9 otherwise)

      COMMON /LINCOM/ VLIN,SLIN,ALIN,ELIN,SBLIN,PSHIFT,DOUBV,TDW,TDWS,
     1 LLQ,IDLIN,LINMAX,NXTLIN,ABSCO_ARR,AD_ARR,Y_ARR,TSTIM_ARR
