      SUBROUTINE FNDWAV(WAVE)
C     $Id: fndwav.f,v 1.2 2011-06-17 15:53:01 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  FNDWAV: searches data base file for wavenumber
C
C_KEYS:   SUBR,LINEDATA
C
C_DESCR:  performs a binary chop search of an already open line
C         data base file. On exit DBREC points to the record containing the
C         next transition at a wavenumber greater or equal to the one specified.
C         The line parameters (in common) are undefined on exit.
C         If a wavenumeber higher than any in the file is requested record
C         DBSIZ+1 is returned.
C
C_FILES :
C         unit DBLUN - already open on entry
C
C_ARGS:   WAVE:DOUBLE PRECISION     wavenumber to search for
C         DBLUN:COMMON  units number of data base
C         DBREC:COMMON  record pointer (OUTPUT)
C
C_CALLS:  RDLINE        performs internal read on line data record
C
C_BUGS:
C
C_HIST:
C         22feb91 SBC Original version
C
C_END:
C--------------------------------------------------------------
      INCLUDE '../includes/dbcom.f' 
C--------------------------------------------------------------
      CHARACTER*256 BUFFER
      DOUBLE PRECISION WAVE
      REAL LOG2
      PARAMETER (LOG2=0.30103)
      INTEGER DELREC,N
      N=INT(LOG10(FLOAT(DBSIZ))/LOG2)-2
      DBREC=0
      DELREC=2**N
10    DBREC=DBREC+DELREC
      IF(DBREC.GT.DBSIZ)GOTO 20
      READ(DBLUN,30,REC=DBREC)BUFFER(1:DBRECL)
30    FORMAT(A)
C     # denotes a HEADER comment so searching next
      IF(BUFFER(2:2).EQ.'#'.OR.BUFFER(1:1).EQ.'#')GOTO 10
      CALL RDLINE(BUFFER)
      IF(LNWAVE.LT.WAVE)GOTO 10
20    CONTINUE
      DBREC=DBREC-DELREC
      IF(DELREC.GT.1)THEN
        DELREC=DELREC/2
        GOTO 10
        END IF
      DBREC=DBREC+1
      RETURN
      END
