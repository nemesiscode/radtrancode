      SUBROUTINE RDGAS
C     $Id: rdgas.f,v 1.4 2011-06-17 15:53:02 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  RDGAS: reads in details of gases in data base
C
C_KEYS:   SUBR, LINEDATA
C
C_DESCR:  opens gas information file (file name already read into common by
C         RDKEY) and reads in gas information, i.e. identifier, name,
C         and number of isotopes; then
C         molecular mass, isotope identifier and relative abundance of each
C         isotope.
C	  Then routine reads in the 4th order polynomial coefficients of
C	  the total Partition function (rot+vib) for each isotope, where:
C		Q(T) = A +B*T + C*T^2 + D*T^3
C
C_ARGS:   
C
C_FILES : unit DBLUN - gas file
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:   4feb91 SBC Original version
C
C_END:
C--------------------------------------------------------------
      INCLUDE '../includes/dbcom.f' 
C--------------------------------------------------------------
      INTEGER I,J,K,L
      CHARACTER*10 BUFFER
      OPEN(UNIT=DBLUN,FILE=GASFIL,STATUS='OLD')

C     Skip past header : 
233   READ(DBLUN,10)BUFFER
      IF(BUFFER(1:2).EQ.'# ')THEN
C       WRITE(*,10)BUFFER
       GOTO 233
      ENDIF
      L=0
30    READ(DBLUN,*,END=100)
      READ(DBLUN,*,END=100)I
      READ(DBLUN,10)GASNAM(I)
C      WRITE(*,10)GASNAM(I)
10    FORMAT(A)
      READ(DBLUN,*)DBNISO(I)
C      WRITE(*,*)DBNISO(I)
      DO 20 J=1,DBNISO(I)
      READ(DBLUN,*)K,RELABU(K,I),MASSNO(K,I)
      READ(DBLUN,*)DBQTA(K,I),DBQTB(K,I),DBQTC(K,I),DBQTD(K,I)
C      WRITE(*,*)K,RELABU(K,I),MASSNO(K,I)
C      WRITE(*,*)DBQTA(K,I),DBQTB(K,I),DBQTC(K,I),DBQTD(K,I)
20    CONTINUE
      L=L+1
      GOTO 30
100   CONTINUE
C      WRITE(*,40)L
40    FORMAT(I3,' gases recognised by software')
      CLOSE(DBLUN)
      RETURN
      END
