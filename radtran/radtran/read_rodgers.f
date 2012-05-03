      SUBROUTINE READ_RODGERS(RODGERS)
C     $Id: read_rodgers.f,v 1.2 2011-06-17 15:40:27 irwin Exp $
C     ***********************************************************************
C     
C     Reads in the Rodgers and Williams (1974) data file containing numerical
C     coefficients to the Lorentz and Doppler widths
C
C     Pat Irwin		16/2/94
*
C     ***********************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION RODGERS(4,0:7)
      CHARACTER*10 HEAD
      CHARACTER*100 ANAME
      INTEGER I,N
1     FORMAT(A)

      ANAME = 'rodgers_williams.dat'
      CALL DATARCHIVE(ANAME)
      OPEN(31,FILE=ANAME,STATUS='OLD')

      DO 10 I=1,8
       READ(31,1)HEAD
10    CONTINUE
      DO 20 I=0,6
       READ(31,*)N,RODGERS(1,I),RODGERS(2,I)
20    CONTINUE
      DO 30 I=1,6
       READ(31,1)HEAD
30    CONTINUE
      DO 40 I=0,7
       READ(31,*)N,RODGERS(3,I),RODGERS(4,I)
40    CONTINUE
      CLOSE(31)

      RETURN

      END
