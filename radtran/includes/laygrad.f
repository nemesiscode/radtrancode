C***********************************************************************
C_TITL: LAYGRAD
C
C_DESC:	Parameters for holding rate of change of path-calculated layers
C	with vmr,T,Fp and dust.
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************

      REAL DTE(MAXLAY,MAXPRO),DFP(MAXLAY,MAXPRO),DAM(MAXLAY,MAXPRO)
      REAL DCO(MAXLAY,MAXPRO),DFC(MAXLAY,MAXPRO)

      COMMON /LGRAD/ DTE,DFP,DFC,DAM,DCO
