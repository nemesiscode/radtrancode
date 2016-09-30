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

      INTEGER VMRflag, GRAVflag, CLOUDflag
      COMMON /MCMCflag/ VMRflag, GRAVflag, CLOUDflag

      REAL MCMCmass, MCMCrad, MCMCvmr(50), MCMCtemp(200)
      COMMON /MCMCmr/ MCMCmass, MCMCrad, MCMCvmr, MCMCtemp

      REAL MCMCheight(200), MCMCcont(200)
      COMMON /MCMCmr2/ MCMCheight, MCMCcont

      REAL MCMCpr, MCMCpvar, MCMCimag, MCMCreal
      COMMON /MCMCmr3/ MCMCpr, MCMCpvar, MCMCimag, MCMCreal

      REAL MCMChknee, MCMCdeep, MCMCfsh  
      COMMON /MCMCmr4/ MCMChknee, MCMCdeep, MCMCfsh

      INTEGER MCtemplen
      COMMON /MCMCint/ MCtemplen
