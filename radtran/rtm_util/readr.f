      REAL FUNCTION READR(TEXT,IERROR)
C     $Id: readr.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C     ***********************************************************************
C     Function to read a real number from a text string of indeterminate length
C
C       1/1/90  SBC     Original version
C       3/10/94 PGJI    Added $Id: readr.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C
C     ***********************************************************************
      CHARACTER TEXT*(*) 
      CHARACTER*8 FORM
      INTEGER I,J,IERROR
      I=LEN(TEXT)
      J=INDEX(TEXT,'.')
      IF(J.EQ.0)THEN
        WRITE(FORM,10)I
10      FORMAT('(I',I2.2,')')
        READ(TEXT,FMT=FORM,ERR=99)J
        READR=FLOAT(J)
       ELSE
        J=I-J
        WRITE(FORM,20)I,J
20      FORMAT('(F',I2.2,'.',I1,')')
        READ(TEXT,FMT=FORM,ERR=99)READR
        END IF
      TEXT=' '
      IERROR=0
      RETURN
99    READR=0
      IERROR=-1
      RETURN
      END
