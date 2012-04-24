      INTEGER FUNCTION READI(TEXT,IERROR)
C     $Id: readi.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C     ***********************************************************************
C     Function to read an integer from a text string of indeterminate length
C
C	1/1/90	SBC	Original version
C	3/10/94	PGJI	Added $Id: readi.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C
C     ***********************************************************************
      CHARACTER TEXT*(*) 
      CHARACTER*8 FORM
      INTEGER J,IERROR
      J=LEN(TEXT)
      WRITE(FORM,10)J
10    FORMAT('(I',I2.2,')')
      READ(TEXT,FMT=FORM,ERR=99)READI
      IERROR=0
      RETURN
99    IERROR=1
      RETURN
      END




