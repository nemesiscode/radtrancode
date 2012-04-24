      SUBROUTINE FCAT(IFILE1,IFILE2,OFILE)
C     $Id: fcat.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C     ********************************************************************
C     Subroutine to concatenate two files together
C
C     ********************************************************************
      INTEGER L
      CHARACTER*(*) IFILE1,IFILE2,OFILE
C
C     first copy the 1st file name to the output array and remove any
C     leading spaces
      CALL REMSP(IFILE1)
      CALL REMSP(IFILE2)
      OFILE=IFILE1
C     now stepping back from the end of the name through the last four
C     characters
      LE1=LEN(IFILE1)
      DO 40 L=LE1,1,-1
      IF(IFILE1(L:L).NE.' ')THEN
        LE1=L
        GOTO 50
        END IF
40    CONTINUE
      RETURN
50    CONTINUE


      LE2=LEN(IFILE2)
      DO 60 L=LE2,1,-1
      IF(IFILE2(L:L).NE.' ')THEN
        LE2=L
        GOTO 70
        END IF
60    CONTINUE
      RETURN
70    CONTINUE

      DO 80 L=1,LE2
       L1 = LE1+L
       OFILE(L1:L1)=IFILE2(L:L)
80    CONTINUE

      RETURN

      END 
