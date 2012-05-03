      PROGRAM INV_SCAN
C     $Id: inv_scan.f,v 1.2 2011-06-17 15:53:01 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  
C
C_KEYS:   PROG,LINEDATA
C
C_DESCR:  
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
C
C_END:
      IMPLICIT NONE
C--------------------------------------------------------------
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
C
      INTEGER FSTLIN,LSTLIN,LINE,ID,LINE1,LSCAN,J,LTOT
      INTEGER FSTREC
      REAL VMIN,VMAX,LNWAVE1,DWAVE
      CHARACTER*256 BUFFER,BUFFER1,XNAME,YNAME,XNAME1,YNAME1
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
2     CALL PROMPT('Enter wavenumber limits?')
      READ(*,*)VMIN,VMAX
      IF(VMAX.LE.VMIN)GOTO 2
C
      CALL PROMPT('Enter name of new file?')
      READ(*,21)OPNAME
21    FORMAT(A)

      CALL PROMPT('Enter gas id : ')
      READ(*,*)ID

      OPEN(UNIT=3,FILE=OPNAME,STATUS='UNKNOWN')
C
      OPEN(DBLUN,FILE=DBFILE,STATUS='OLD',FORM='FORMATTED',
     1ACCESS='DIRECT',RECL=DBRECL)

      OPEN(UNIT=4,FILE=INDFIL,STATUS='OLD',ACCESS='DIRECT',
     1FORM='UNFORMATTED',RECL=INDRL)

C
C
      CALL FNDWAV(VMIN)
      FSTLIN=DBREC-1
      IF(FSTLIN.LT.1)FSTLIN=1
      CALL FNDWAV(VMAX)
      LSTLIN=DBREC-1
      LSCAN = 0
      LTOT = 0

      PRINT*,'FSTLIN,LSTLIN = ',FSTLIN,LSTLIN
C
C     copy required lines

C     find first line of gas required
      DO 140 LINE=FSTLIN,LSTLIN
        READ(DBLUN,111,REC=LINE)BUFFER(1:DBRECL)
        READ(BUFFER,98,ERR=99)LNID
        IF(LNID.EQ.ID) THEN
         FSTREC = LINE
         GOTO 145
        ENDIF
140   CONTINUE


145   READ(DBLUN,111,REC=LINE)BUFFER(1:DBRECL)
C      WRITE(*,*)'LINE : ',LINE
      READ(BUFFER,113,ERR=99)LNWAVE
        
      XNAME = BUFFER(4:55)
      YNAME = BUFFER(1:3)//BUFFER(56:100)
          
        LINE1 = FSTREC
        DWAVE = 0.0

150     READ(DBLUN,111,REC=LINE1)BUFFER1(1:DBRECL)
C        WRITE(*,*)'        LINE1 : ',LINE1

        READ(BUFFER1,113,ERR=99)LNWAVE1
        XNAME1 = BUFFER1(4:55)
        YNAME1 = BUFFER1(1:3)//BUFFER1(56:100)
         
        IF((XNAME.NE.XNAME1).AND.(YNAME.EQ.YNAME1))THEN
             WRITE(*,*)'Inversion doublet found'
             DWAVE = 0.5*(LNWAVE-LNWAVE1)
             LSCAN=LSCAN+1
	     GOTO 160
        ENDIF

C       Find record number of next gas entry
        READ(4,REC=LINE1)J
        LINE1 = J

        IF(LINE1.LE.LSTLIN)GOTO 150

160   LTOT=LTOT+1
    
      LDOUBV = DWAVE

      WRITE(3,114)BUFFER(1:DBRECL),LDOUBV



      READ(4,REC=LINE)J
      LINE = J
      IF(LINE.LE.LSTLIN)GOTO 145

      WRITE(*,112)LTOT,LSCAN

98    FORMAT(I2)

111   FORMAT(A)

112   FORMAT(1X,I8,' lines copied',I8,' inversion lines')

113   FORMAT(3X,F12.6)
114   FORMAT(A,F10.6)

      CLOSE(DBLUN)
      CLOSE(3)

115   FORMAT(A,F10.6)

      STOP

99    CONTINUE
      WRITE(*,*)'I/O error'
      CLOSE(DBLUN)
      CLOSE(3)


      END
