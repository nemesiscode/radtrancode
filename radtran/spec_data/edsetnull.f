      SUBROUTINE EDSETNULL
C     $Id: edsetnull.f,v 1.1 2010-06-09 16:33:16 irwin Exp $
C--------------------------------------------------------------
C
C_HIST:	
C		NT	20.jan.2004 Added option to broaden PH3 for Saturn
C				h2 and he abundances.
C--------------------------------------------------------------
C     common variables used by all linedata routines
C     mixed data types are avoided
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
      LOGICAL WSBTAB,CO2TAB,H2OTAB,SCOFLD,WEDAD,BEZARD,PMIRR
      LOGICAL H2HePH3J,H2HePH3S
      COMMON /EDLOG/WSBTAB,CO2TAB,H2OTAB,SCOFLD,WEDAD,BEZARD,PMIRR,
     $      H2HePH3J,H2HePH3S
C--------------------------------------------------------------
      LOGICAL ASKYN
      CHARACTER*256 BUFFER

      WSBTAB=.FALSE.
      CO2TAB=.FALSE.
      H2OTAB=.FALSE.
      SCOFLD=.FALSE.
      WEDAD=.FALSE.
      BEZARD=.FALSE.
      PMIRR=.FALSE.
      H2HePH3J=.FALSE.
      H2HePH3S=.FALSE.

      RETURN

      END
