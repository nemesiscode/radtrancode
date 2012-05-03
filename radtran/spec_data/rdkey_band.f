      SUBROUTINE RDKEY_BAND(LUN)
C     $Id: rdkey_band.f,v 1.3 2011-06-17 15:53:02 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  RDKEY: reads in line data key file
C
C_KEYS:   SUBR,LINEDATA
C
C_DESCR:  opens key file, prompting for name if needed
C         and reads in all information in file
C
C_ARGS:   LUN:INTEGER  logical unit for keyfile
C
C_FILES :
C         unit LUN - keyfile
C
C_CALLS:  FILE    sets file extension in text string
C
C_BUGS:
C
C_HIST:
C        22feb91 SBC Original version
C	 25oct94 PGJI modified for band files
C
C_END:
C--------------------------------------------------------------
      INCLUDE '../includes/dbcom.f' 
C--------------------------------------------------------------
      INTEGER LUN

      DBLUN=LUN
20    FORMAT(A)
      CALL FILE(KEYFIL,KEYFIL,'key')

      OPEN(UNIT=LUN,FILE=KEYFIL,STATUS='OLD')

      READ(LUN,20)DBFILE
      CALL REMSP(DBFILE)
      WRITE(*,20)DBFILE

      READ(LUN,20)GASFIL
      CALL REMSP(GASFIL)
      WRITE(*,20)GASFIL

      CLOSE(LUN)

      RETURN
      END
