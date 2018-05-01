C     $Id: dbcom.f,v 1.6 2011-06-17 14:49:05 irwin Exp $
C***********************************************************************
C_TITL: LCOCOM
C
C_DESC:	Common variables used by line continuum routines.
C***********************************************************************

      INTEGER MLCO,IJLCO,NBINLCO,IDLCO,ISOLCO
      PARAMETER (MLCO=3000)
      DOUBLE PRECISION VLCO(MLCO),SLCO(MLCO),LCOBINSIZE,VMINLCO,
     &  VMAXLCO
      REAL LCLSE(MLCO),LCWIDA(MLCO),LCWIDS(MLCO),LCTDEP(MLCO)

      COMMON /LCODP/VMINLCO,VMAXLCO,LCOBINSIZE,VLCO,SLCO
      COMMON /LCOREAL/LCLSE,LCWIDA,LCWIDS,LCTDEP
      COMMON /LCINT/IJLCO,NBINLCO,IDLCO,ISOLCO


