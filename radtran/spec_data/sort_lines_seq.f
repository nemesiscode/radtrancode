      PROGRAM SORT_LINES_SEQ
C     $Id:
C--------------------------------------------------------------
C_TITLE:  COPY: copies a subset of sequential-access line data base to a new data base
C
C_KEYS:   PROG,LINEDATA
C
C_DESCR:  
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
      IMPLICIT NONE
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
C
      DOUBLE PRECISION VMIN,VMAX,STRMIN,V1,V2
      INTEGER I,NLINE,ISWAP
      CHARACTER*256 BUFFER
      CHARACTER*100 OPNAME
      CHARACTER*128 LDATA(1000000),LSWAP
C
C     open existing sequential access database
      CALL PROMPT('Enter name of sequential-access data base:')
      READ(*,1)KEYFIL
1     FORMAT(A)
      OPEN(DBLUN,FILE=KEYFIL,STATUS='OLD',READONLY)

      CALL PROMPT('Enter name of output file:')
      READ(*,1)OPNAME
      OPEN(UNIT=3,FILE=OPNAME,STATUS='UNKNOWN')

      CALL PROMPT('Enter record length : ')
      READ*,DBRECL

      I=0
      PRINT*,'Reading in data'
100   READ(DBLUN,1,END=999)BUFFER
      I=I+1
      LDATA(I)=BUFFER(1:DBRECL)
      GOTO 100

999   CONTINUE
      NLINE=I
      CLOSE(DBLUN)

      PRINT*,'Sorting...'
300   CONTINUE
      ISWAP=0

      DO 201 I=1,NLINE-1 
       BUFFER=LDATA(I)
       CALL RDLINE(BUFFER)
       V1=LNWAVE
       BUFFER=LDATA(I+1)
       CALL RDLINE(BUFFER)
       V2=LNWAVE
       IF(V1.GT.V2)THEN
        LSWAP=LDATA(I)
        LDATA(I)=LDATA(I+1)
        LDATA(I+1)=LSWAP
        ISWAP=1
       ENDIF
201   CONTINUE
      IF(ISWAP.EQ.1) GOTO 300

      DO 400 I=1,NLINE
       BUFFER=LDATA(I)
       WRITE(3,1)BUFFER(1:DBRECL)
400   CONTINUE 

      PRINT*,'Transfer complete'
      CLOSE(3)

      END

