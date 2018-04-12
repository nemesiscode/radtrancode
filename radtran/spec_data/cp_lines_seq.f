      PROGRAM CP_LINES_SEQ
C     $Id:
C--------------------------------------------------------------
C_TITLE:  COPY: copies a subset of sequential-access line data base to a new data base
C
C_KEYS:   PROG,LINEDATA
C
C_DESCR:  copies line data bases (eg HITRAN or GEISA) from an existing data
C         base into a sequential ascii file which can be turned into a new data
C         base using makedb. Unlike SELECT.FOR it copies ALL lines
C         between given wavenumber limits rather than selecting them by
C         gas, isotope and strength.
C
C_FILES :
C         unit 2 - data base files
C         unit 3 - output file
C
C_BUGS:   
C
C_HIST:
C         26apr92 SBC Original version
C	  10apr18 PGJI Updated for sequential access
C_END:
C--------------------------------------------------------------
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
C
      DOUBLE PRECISION VMIN,VMAX
      CHARACTER*256 BUFFER
      CHARACTER*100 OPNAME
C
C     open existing sequential access database
      CALL PROMPT('Enter name of sequential-access data base:')
      READ(*,1)KEYFIL
1     FORMAT(A)
      OPEN(DBLUN,FILE=KEYFIL,STATUS='OLD',READONLY)

C     read in parameters
2     CALL PROMPT('wavenumber limits to extract?')
      READ(*,*)VMIN,VMAX
      IF(VMAX.LE.VMIN)GOTO 2

      CALL PROMPT('Enter name of output file:')
      READ(*,1)OPNAME
      OPEN(UNIT=3,FILE=OPNAME,STATUS='UNKNOWN')

      CALL PROMPT('Enter record length : ')
      READ*,DBRECL

      IF(DBRECL.NE.160)THEN
       PRINT*,'Code not set up for DBRECL = ',DBRECL
       GOTO 998
      ENDIF

100   READ(DBLUN,1,END=999)BUFFER
      READ(BUFFER,101,ERR=998)LNID,LNISO,LNWAVE
101   FORMAT(I2,I1,F12.6)
      IF(LNWAVE.GE.VMIN.AND.LNWAVE.LE.VMAX)WRITE(3,1)BUFFER(1:DBRECL)
      IF(LNWAVE.GT.VMAX)GOTO 999
      GOTO 100

998   PRINT*,'Error in reading. Aborting'
      CLOSE(3)
      CLOSE(DBLUN)
      STOP

999   PRINT*,'Transfer complete'
      CLOSE(3)
      CLOSE(DBLUN)


      END
