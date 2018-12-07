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

      INTEGER MCMCflag, PHASEflag, Cyflag
      COMMON /MCMCflag/ MCMCflag, PHASEflag, Cyflag

      REAL MCMCmass, MCMCrad, MCMCvmr(200,50), MCMCtemp(200)
      COMMON /MCMCmr/ MCMCmass, MCMCrad, MCMCvmr, MCMCtemp

      REAL MCMCheight(200), MCMCcont(200), MCMCpres(200)
      COMMON /MCMCmr2/ MCMCheight, MCMCcont, MCMCpres

      REAL MCMCpr, MCMCpvar, MCMCimag(50), MCMCreal
      COMMON /MCMCmr3/ MCMCpr, MCMCpvar, MCMCimag, MCMCreal

      REAL MCMChknee, MCMCdeep, MCMCfsh  
      COMMON /MCMCmr4/ MCMChknee, MCMCdeep, MCMCfsh

      REAL MCMChknee2, MCMCdeep2, MCMCfsh2
      COMMON /MCMCmr5/ MCMChknee2, MCMCdeep2, MCMCfsh2

      REAL MCMCflat(50), MCMCflon(50), MCMCsolzen(50)
      COMMON /MCMCmr6/ MCMCflat, MCMCflon, MCMCsolzen

      REAL MCMCemzen(50),MCMCazi(50),MCMCwt(50)
      COMMON /MCMCmr7/ MCMCemzen,MCMCazi,MCMCwt

      REAL MCMCtop, MCMCscat, MCMCshape
      COMMON /MCMCmr9/ MCMCtop, MCMCscat, MCMCshape

      INTEGER MCtemplen, MCMCnav
      COMMON /MCMCint/ MCtemplen, MCMCnav

      REAL kl_g_Cy(300,50,200)
      COMMON /MCMCmr8/ kl_g_Cy

