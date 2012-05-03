      PROGRAM SCAN
C     $Id: scan.f,v 1.2 2011-06-17 15:53:03 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  SCAN: lists line data to screen - ascii as stored in file
C
C_KEYS:   PROG,LINEDATA
C
C_DESCR:  searches for a wavelength entry in data base files then lists
C         transitions to the screen in format as stored in file
C
C_FILES :
C         unit 1 - misc
C         unit 2 - data base
C
C_CALLS:  RDKEY     reads key file
C         FNDWAV    searches data base for wavenumber entry
C         RDLINE    uses internal read to translate line data record
C         PROMPT    prompts user for input
C         REMSP     removes leading spaces from text string
C         WTEXT     writes text to screen
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
      CHARACTER*100 TEXT
      REAL WAVE
      INTEGER IDIS
      LOGICAL ABO,ASKYN
      CALL PROMPT('data base key?')
      READ(*,102)KEYFIL
102   FORMAT(A)
      CALL RDKEY(2)
      print*,'DBRECL = ',DBRECL
      OPEN(UNIT=DBLUN,FILE=DBFILE,STATUS='OLD',RECL=DBRECL,
     1ACCESS='DIRECT',FORM='FORMATTED')
      CALL PROMPT('from wavenumber?')
      READ(*,*)WAVE
      CALL FNDWAV(WAVE)
      WRITE(*,23)DBREC
23    FORMAT(' begining at record:',I10)
      IDIS=0
10    READ(DBLUN,12,REC=DBREC)BUFFER(1:DBRECL)
12      FORMAT(A)
      IF(DBFORM.EQ.0)THEN
        WRITE(TEXT,20)DBREC,BUFFER(1:67)
        CALL WTEXT(TEXT)
        WRITE(TEXT,21)BUFFER(68:DBRECL)
        CALL WTEXT(TEXT)
20      FORMAT(I10,':',A)
21      FORMAT(12X,A)
        IDIS=IDIS+2
       ELSE IF(DBFORM.EQ.1)THEN
C        print*,'buffer = ',
C        print*,buffer
C        print*,DBRECL
        IF(DBRECL.EQ.80)THEN
 	 WRITE(TEXT,22)DBREC
         CALL WTEXT(TEXT)
	 WRITE(TEXT,24)BUFFER(1:DBRECL)
	 CALL WTEXT(TEXT)
        ELSE IF(DBRECL.EQ.120)THEN
 	 WRITE(TEXT,20)DBREC,BUFFER(1:67)
         CALL WTEXT(TEXT)
	 WRITE(TEXT,21)BUFFER(68:DBRECL)
	 CALL WTEXT(TEXT)
        END IF
22      FORMAT(I10,': -> ')
24      FORMAT(A)
         IDIS=IDIS+2
       ELSE
        WRITE(*,11)
11      FORMAT(' invalid data format')
        STOP
        END IF
      DBREC=DBREC+1
      IF(IDIS.EQ.22)THEN
        IDIS=0
        ABO=ASKYN('abort?')
        IF(ABO)STOP
        END IF
      IF(DBREC.LE.DBSIZ)GOTO 10
      STOP
      END
