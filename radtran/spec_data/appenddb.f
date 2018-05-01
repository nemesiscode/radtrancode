      PROGRAM APPENDDB
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INCLUDE '../includes/dbcom.f' 
C ../includes/dbcom.f stores the linedata base variables.

      INTEGER IREC
      CHARACTER*100 OPNAME
      CHARACTER*256 BUFFER

C******************************** CODE *********************************

C Open database ...
      CALL PROMPT('name of existing data base key')
      READ(*,1)KEYFIL
      CALL REMSP(KEYFIL)
1     FORMAT(A)
      CALL RDKEY(2)

      OPEN(DBLUN,FILE=DBFILE,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='DIRECT',RECL=DBRECL)


      CALL PROMPT('Enter name of sequential access database : ')
      READ(*,1)OPNAME

      OPEN(12,FILE=OPNAME,STATUS='OLD')
      IREC=DBSIZ

201   READ(12,1,END=999)BUFFER
      IREC=IREC+1
      WRITE(DBLUN,1,REC=IREC)BUFFER(1:DBRECL)
      GOTO 201
   
999   CLOSE(12)
      CLOSE(DBLUN)
      print*,'Old number of records = ',DBSIZ
      print*,'New number of records = ',IREC
      print*,''
      print*,'Update the key file accordingly'

      END
 
