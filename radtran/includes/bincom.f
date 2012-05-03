C     $Id: bincoms.f,v 1.7 2005-06-13 10:46:26 irwin Exp $
C***********************************************************************
C_TITL: BINCOMS
C
C_DESC:
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************

      INTEGER CURBIN,FSTBIN,LSTBIN,NBIN
C CURBIN: number of the bin currently in use.

      INTEGER NLINES(MAXBIN),FSTLIN(MAXBIN),LSTLIN(MAXBIN)
C NLINES: number of lines in each bin.
C FSTLIN: location (in the line arrays) of the first line for each stored
C (`local') bin.
C LSTLIN: location of the last line (if no lines then LSTLIN= FSTLIN-1).

      REAL VBIN(MAXBIN),VBINB(MAXBIN)
C VBIN: minimum wavenumber for each bin.
C VBINB: minimum wavenumber for each broadband bin.

      COMMON /BINCOM/ VBIN,VBINB,CURBIN,FSTBIN,LSTBIN,
     1 NLINES,FSTLIN,LSTLIN,NBIN,BANDPAR,BANDTYP

C-----------------------------------------------------------------------
C
C	Band parameters
C
C-----------------------------------------------------------------------
      INTEGER MAXBGAS
      PARAMETER(MAXBGAS=15)
      REAL BANDPAR(MAXBIN,MAXBGAS,7)
      INTEGER BANDTYP(MAXBGAS)
