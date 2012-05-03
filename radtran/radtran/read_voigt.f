      SUBROUTINE READ_VOIGT
C     $Id: read_voigt.f,v 1.2 2011-06-17 15:40:27 irwin Exp $
C     **********************************************************************
C
C     Reads in and takes the log of the Voigt EW look up table
C
C     Pat Irwin		16/2/94
C 
C     **********************************************************************
      IMPLICIT NONE
      REAL VWIDTH(117,161)
      REAL TWIDTH(117,161)
      COMMON /VOIGTDATA/VWIDTH
      INTEGER I,J
      CHARACTER*100 ANAME

      ANAME = 'voigt_ew.tab'
      CALL DATARCHIVE(ANAME)
 
      OPEN(15,FILE=ANAME,STATUS='OLD',FORM='UNFORMATTED')
       READ(15)TWIDTH
      CLOSE(15)

      DO 10 I=1,117

       DO 20 J=1,161

        IF(TWIDTH(I,J).LE.0.)THEN
C         PRINT*,I,J,TWIDTH(I,J)
         VWIDTH(I,J)=-999.
        ELSE
         VWIDTH(I,J)=LOG(TWIDTH(I,J))
        END IF

20     CONTINUE

10    CONTINUE

      RETURN

      END
      
