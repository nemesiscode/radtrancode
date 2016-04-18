C***********************************************************************
C_TITL: EMCEE
C
C_DESC:	For passing emcee inputs to a priori.
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************

      INTEGER VMRflag, GRAVflag
      COMMON /MCMCflag/ VMRflag, GRAVflag

      REAL MCMCmass, MCMCrad, MCMCvmr(50), MCMCtemp(200)
      COMMON /MCMCmr/ MCMCmass, MCMCrad, MCMCvmr, MCMCtemp

