      SUBROUTINE RDKEY_CORRK(LUN,NGAS,IDGAS,ISOGAS,CORNAME)
C     $Id: rdkey_corrk.f,v 1.2 2011-06-17 15:40:27 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  RDKEY_CORRK: reads in line data key file
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
      INCLUDE '../includes/arrdef.f' 
      INCLUDE '../includes/dbcom.f' 
C--------------------------------------------------------------
      INTEGER LUN,NGAS,ID,ISO
      INTEGER IDGAS(NGAS),ISOGAS(NGAS)
      CHARACTER*100 CORNAME(10),TMPNAM,OUTNAM
      DBLUN=LUN

      CALL FILE(KEYFIL,KEYFIL,'key')
      OPEN(UNIT=LUN,FILE=KEYFIL,STATUS='OLD')
      DO 34 I=1,NGAS
       READ(LUN,*)ID,ISO
       IF(ID.NE.IDGAS(I))THEN
        PRINT*,'rdkey_corrk: ID <> IDGAS'
        PRINT*,'ID = ',ID
        PRINT*,'I, IDGAS(I) = ',I,IDGAS(I)
        STOP
       END IF
       IF(ISO.NE.ISOGAS(I))THEN
        PRINT*,'rdkey_corrk: ISO <> ISOGAS'
        PRINT*,'ISO = ',ISO
        PRINT*,'I, ISOGAS(I) = ',I,ISOGAS(I)
        STOP
       END IF

       READ(LUN,20)TMPNAM
       CALL REMSP(TMPNAM)
       CALL TwiddleExpand(TMPNAM, OUTNAM)
       CORNAME(I) = OUTNAM

34    CONTINUE
      READ(LUN,20)TMPNAM
      CALL REMSP(TMPNAM)
      CALL TwiddleExpand(TMPNAM, GASFIL)

20    FORMAT(A)
      CLOSE(LUN)

      RETURN
      END
