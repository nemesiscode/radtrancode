C     $Id: lincomc.f,v 1.1 2011-06-17 14:48:20 irwin Exp $
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
      PARAMETER (MAXLIN=10000)
C MAXLIN: maximum number of lines which can be stored.

      INTEGER IDLIN(2,MAXLIN),NLINE(2),NXTREC
      REAL VLIN(2,MAXLIN),SLIN(2,MAXLIN),ALIN(2,MAXLIN)
      REAL PSHIFT(2,MAXLIN)
C VLIN: Line position [cm^-1].
C SLIN: Line strength [cm^-1 molecule^-1 cm^-2] at STP.
C ALIN: Lorentz half width [cm^-1] at STP.
C PSHIFT: line pressure-shift parameter.

      REAL ELIN(2,MAXLIN),SBLIN(2,MAXLIN),TDW(2,MAXLIN)
      REAL TDWS(2,MAXLIN)
C ELIN: Lower-state energy line-position [cm^-1].
C SBLIN: Line self-broadening coefficients. NOTE: the self broadening
C coefficient used in this program is the mixing gas 'air'-broadened
C half width minus the self-broadened half width.

      REAL DOUBV(2,MAXLIN)
C DOUBV: Inversion doublet separation (only applies to longwave NH3 lines
C at present.

      CHARACTER*15 LLQ(2,MAXLIN)
C LLQ: lower state local quanta. (15 elements for Hitran04, 9 otherwise)

      COMMON /LINCOMC/ VLIN,SLIN,ALIN,ELIN,SBLIN,PSHIFT,DOUBV,TDW,TDWS,
     1 LLQ,IDLIN,NLINE,NXTREC
