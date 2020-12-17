C***********************************************************************
C_TITL: ELECTRON
C
C_DESC:	For passing electron partial pressure to ngascon
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************

      REAL PELEC, AMOUNTHMIN, PPRESSHMIN
      INTEGER EFLAG

      COMMON /ELECTRON/ PELEC, AMOUNTHMIN, PPRESSHMIN
      COMMON /ELECTRON/ EFLAG
