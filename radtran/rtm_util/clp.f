      SUBROUTINE CLP(TEXT,NPAR,FIRST,LAST)
C     $Id: clp.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C     ********************************************************************
C     A simple command line processor
C     Routine simply finds all parameters and returns their location
C
C     On entry NPAR is the maximum number of parameters allowed
C     On exit NPAR is the number of parameters found
C     FIRST and LAST are the positions of the first and last character
C     for each parameter
C
C     Only spaces are used to delimit parameters
C
C     1/1/90    Original Version:       SBC
C     3/10/94   Updated Header          PGJI
C
C     ********************************************************************
      CHARACTER TEXT*(*)
      INTEGER NPAR,FIRST(NPAR),LAST(NPAR)

      INTEGER I,J,MAXPAR

      J=LEN(TEXT)
      I=0
      MAXPAR=NPAR
      NPAR=0
100   I=I+1
      IF(I.GT.J)RETURN
C     removing control characters
      IF(ICHAR(TEXT(I:I)).LE.26)TEXT(I:I)=' '
      IF(TEXT(I:I).NE.' ')THEN
C       here if found begining of a parameter
        NPAR=NPAR+1
        IF(NPAR.LE.MAXPAR)FIRST(NPAR)=I
200     I=I+1
        IF(I.GT.J)THEN
          IF(NPAR.LE.MAXPAR)LAST(NPAR)=J
          RETURN
         ELSE IF(TEXT(I:I).EQ.' ')THEN 
          IF(NPAR.LE.MAXPAR)LAST(NPAR)=I-1
          GOTO 100
          END IF
        GOTO 200
        END IF
      GOTO 100
      END
