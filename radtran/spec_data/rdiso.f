      SUBROUTINE RDISO
C     $Id: rdiso.f,v 1.4 2011-06-17 14:57:07 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  RDISO: Reads in isotope details for line data base.
C
C_KEYS:   SUBR,LINEDATA
C
C_DESCR:  Opens isotope data file (file name already read into common by
C         RDKEY) and reads in isotope mapping. i.e. the N isotopes of a gas
C         are refered to by the integers 1 to N. The same integers are used
C         irrespective of data base. The isotope identifiers in the data base
C         file will generally be different. The isotope data file contains the
C         mapping between the two.
C        
C_ARGS:   
C
C_FILES : unit DBLUN - the isotope data file
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
      INTEGER I, J, K, ISO(MAXISO), ST, EN
      CHARACTER*100 ISOSTR,BUFFER

      WRITE(*,*)'Reading Isotope file :' 
      OPEN(UNIT = DBLUN, FILE = ISOFIL, STATUS = 'OLD')
      K = 0
 
5     READ(DBLUN,15)BUFFER
      IF(BUFFER(1:2).eq.'# ')GOTO 5   

     
10    READ (DBLUN, 13, END = 100, ERR = 100) I1,I2, ISOSTR

13    FORMAT (I6, I5, A)
15    FORMAT(A)

      IF(I1.EQ.0) GOTO 100

      
      I = I1
      LOCID(I)=I2

      J = 0
      ST = 1
      EN = 1
      DO WHILE (ST .LT. LEN(ISOSTR))
         DO WHILE (ST .LT. LEN(ISOSTR) .AND. ISOSTR(ST:ST) .EQ. ' ')
            ST = ST + 1
         END DO
         EN = ST
         DO WHILE (EN .LT. LEN(ISOSTR) .AND. ISOSTR(EN:EN) .NE. ' ')
            EN = EN + 1
         END DO
         EN = EN - 1
         IF (ST .LT. LEN(ISOSTR) .AND. EN .LE. LEN(ISOSTR)) THEN
            J = J + 1
            READ (ISOSTR(ST:EN), 14) ISO(J)
 14         FORMAT (I5)
         END IF
         ST = EN + 1
      END DO

      DO 30 J = 1, DBNISO(LOCID(I))
C     note- some isotopes (CH3D) have separate IDs in the data base
C     these are mapped onto a single local ID but only the non zero
C     isotope codes are mapped to avoid overwriting the main entries
      IF (ISO(J) .GT. 0) DBISO(J, LOCID(I)) = ISO(J)

 30   CONTINUE
      K = K + 1
      PRINT*,K,I,LOCID(I),DBNISO(LOCID(I)),(DBISO(J,LOCID(I)),
     1   J=1,DBNISO(LOCID(I)))
      GOTO 10



 100  CONTINUE
      WRITE (*, 20) K
 20   FORMAT(I3, ' gases in data base')

      CLOSE (DBLUN)

C      DO 200 I=1,K
C       PRINT*,I,LOCID(I),DBNISO(LOCID(I)),(DBISO(J,LOCID(I)),
C     1   J=1,DBNISO(LOCID(I)))
C200   CONTINUE
      
      RETURN
      END
