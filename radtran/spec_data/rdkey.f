      SUBROUTINE RDKEY(LUN)
C     $Id: rdkey.f,v 1.5 2011-06-17 15:53:02 irwin Exp $
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
C
C_END:
C--------------------------------------------------------------
      INCLUDE '../includes/dbcom.f' 
C--------------------------------------------------------------
      INTEGER LUN
      CHARACTER*56 TMPNAM
      DBLUN=LUN
      CALL FILE(KEYFIL,KEYFIL,'key')
      OPEN(UNIT=LUN,FILE=KEYFIL,STATUS='OLD')
      READ(LUN,20)DBNAME
C     note not overwriting KEYFIL with original name
      READ(LUN,*)
      READ(LUN,20)TMPNAM
      IPFILE(1:56)=TMPNAM(1:56)
      CALL REMSP(IPFILE)
      IPFILE(57:100)='                                            '
      WRITE(*,20)IPFILE
      READ(LUN,20) TMPNAM
      DBFILE(1:56)=TMPNAM(1:56)
      DBFILE(57:100)='                                            '
      CALL REMSP(DBFILE)
      WRITE(*,20)DBFILE
      READ(LUN,20) TMPNAM
      INDFIL(1:56)=TMPNAM(1:56)
      INDFIL(57:100)='                                            '
      CALL REMSP(INDFIL)
      WRITE(*,20)INDFIL
      READ(LUN,20) TMPNAM
      GASFIL(1:56)=TMPNAM(1:56)
      GASFIL(57:100)='                                            '
      CALL REMSP(GASFIL)
      WRITE(*,20)GASFIL
      READ(LUN,20) TMPNAM
      ISOFIL(1:56)=TMPNAM(1:56)
      ISOFIL(57:100)='                                            '
      CALL REMSP(ISOFIL)
      WRITE(*,20)ISOFIL

C     relies on size of variables to avoid comments at end of line
20    FORMAT(A56)

      READ(LUN,*)DBSIZ
      READ(LUN,*)DBFORM
      READ(LUN,*)INDRL
      INDRL = ISYS()
      READ(LUN,*)DBRECL
      CLOSE(LUN)

      RETURN
      END
