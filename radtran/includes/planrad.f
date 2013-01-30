C     $Id: molwtcom.f,v 1.5 2011-06-17 14:47:50 irwin Exp $
C***********************************************************************
C_TITL: PLANRAD
C
C_DESC:	Way of passing adjusted radius to forward model, also 
C	    for passing radius to forwarddisc for disc-integrated calcs
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************
C-----------------------------------------------------------------------
C
C	
C	
C
C-----------------------------------------------------------------------
	  REAL RADIUS2
	  INTEGER jradf
      COMMON /PLANRAD/RADIUS2,jradf