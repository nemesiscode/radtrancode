      PROGRAM COPY
C     $Id: cp_lines.f,v 1.2 2011-06-17 15:53:01 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  COPY: copies a subset of a line data base to a new data base
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
C_CALLS:  RDGAS    reads in the gas file information
C         RDISO    reads in the isotope file information
C         WTEXT    writes test to screen
C         PROMPT   prompts user for input
C         RDKEY    reads in details from the key file
C         FNDWAV   searches data base for wavenumber entry
C
C_BUGS:   
C
C_HIST:
C         26apr92 SBC Original version
C_END:
C--------------------------------------------------------------
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
C
      INTEGER FSTLIN,LSTLIN,LINE
      DOUBLE PRECISION VMIN,VMAX
      CHARACTER*256 BUFFER
      CHARACTER*100 OPNAME
C
C     open data base
      CALL PROMPT('name of data base key')
      READ(*,1)KEYFIL
1     FORMAT(A)
      CALL RDKEY(2)
      CALL RDGAS
      CALL RDISO
C     read in parameters
C     parameters are:_
C      VMIN            minimum wavenumber
C      VMAX            maximum wavenumber
2     CALL PROMPT('wavenumber limits?')
      READ(*,*)VMIN,VMAX
      IF(VMAX.LE.VMIN)GOTO 2
C
      CALL PROMPT('new file?')
      READ(*,21)OPNAME
21    FORMAT(A)
      OPEN(UNIT=3,FILE=OPNAME,STATUS='UNKNOWN')
C
      OPEN(DBLUN,FILE=DBFILE,STATUS='OLD',FORM='FORMATTED',
     1ACCESS='DIRECT',RECL=DBRECL)
C
      WRITE(BUFFER,114)
114   FORMAT(' # data records written by routine COPY')
      WRITE(3,111)BUFFER(1:DBRECL)
      WRITE(BUFFER,201)
201   FORMAT(' # original data base files:')
      WRITE(3,111)BUFFER(1:DBRECL)
      WRITE(BUFFER,202)DBFILE
202   FORMAT(' #  ',A)
      WRITE(3,111)BUFFER(1:DBRECL)
      WRITE(BUFFER,202)KEYFIL
      WRITE(3,111)BUFFER(1:DBRECL)
      WRITE(BUFFER,203)
203   FORMAT(' # original data base header if any:')
      WRITE(3,111)BUFFER(1:DBRECL)
      LINE=1
204   READ(DBLUN,111,REC=LINE)BUFFER(1:DBRECL)
      IF(BUFFER(2:2).EQ.'#'.OR.BUFFER(1:1).EQ.'#')THEN
        WRITE(3,111)BUFFER(1:DBRECL)
        LINE=LINE+1
        GOTO 204
        END IF
      WRITE(BUFFER,205)
205   FORMAT(' # selection criteria:')
      WRITE(3,111)BUFFER(1:DBRECL)
      WRITE(BUFFER,206)VMIN,VMAX
206   FORMAT(' # wavenumber range =',F10.3,' -',F10.3)
      WRITE(3,111)BUFFER(1:DBRECL)
C
      CALL FNDWAV(VMIN)
      FSTLIN=DBREC-1
      IF(FSTLIN.LT.1)FSTLIN=1
      CALL FNDWAV(VMAX)
      LSTLIN=DBREC-1
C
C     copy required lines
      DO 150 LINE=FSTLIN,LSTLIN
      READ(DBLUN,111,REC=LINE)BUFFER(1:DBRECL)
111   FORMAT(A)
      WRITE(3,111)BUFFER(1:DBRECL)
150   CONTINUE
C
      LINE=LSTLIN-FSTLIN+1
      WRITE(*,112)LINE
112   FORMAT(1X,I8,' lines copied')
      STOP
      END
