C     $Id: parcom.f,v 1.3 2002-11-23 18:35:52 parrish Exp $
C***********************************************************************
C_TITL: PARCOM
C
C_DESC:	Parameter inclusion file.
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************

C-----------------------------------------------------------------------
C
C	Adaptive integration variables.
C
C-----------------------------------------------------------------------
      INTEGER MAXINT
      PARAMETER (MAXINT=2000)
C MAXINT: each integration interval is subdivided into up to MAXINT 
C intervals and the integration adaptively performed for each sub
C interval.

C-----------------------------------------------------------------------
C
C	Continuum variables.
C
C-----------------------------------------------------------------------
      INTEGER NWAV
      PARAMETER(NWAV=100)
C NWAV is number of wavenumber points which can be used to compute 
C polynomials for continuum.
