C     $Id:
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
      INTEGER MAXLINSEQ
      PARAMETER (MAXLINSEQ=20000000)
      CHARACTER*256 SEQBUF(MAXLINSEQ)
C MAXLINSEQ: maximum number of lines which can be stored.

      COMMON /LINCOMSEQ/SEQBUF
