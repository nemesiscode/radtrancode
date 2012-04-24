      SUBROUTINE READ_YAM(ACO2,NCO2,AN2,NN2,ECO2,EN2)
C     $Id: read_yam.f,v 1.2 2011-06-17 15:53:03 irwin Exp $
C     *************************************************************
C     SUBROUTINE TO READ IN YAMAMOTO CO2 DATA FROM FILE
C     DATA PROCESSED BY J.T.SCHOFIELD TO GET N VALUES
C
C     PAT IRWIN    10/2/93
C     *************************************************************
      REAL ACO2(200),NCO2(200),AN2(200),NN2(200),ECO2(200),EN2(200)
      INTEGER M,I,K
      CHARACTER*10 DUMMY
      CHARACTER*100 ANAME

      ANAME = 'yamamoto.dat'
      CALL DATARCHIVE(ANAME)
      OPEN(12,FILE=ANAME,STATUS='OLD')

      DO 10 I=1,15
       READ(12,100)DUMMY
10    CONTINUE
      DO 20 I=1,200
       READ(12,105) K,ACO2(I),NCO2(I),AN2(I),NN2(I),ECO2(I),EN2(I)
C       WRITE(6,105) K,ACO2(I),NCO2(I),AN2(I),NN2(I),ECO2(I),EN2(I)
20    CONTINUE
      CLOSE(12)
100   FORMAT(A10)
105   FORMAT(1X,I4,2(F12.4,F12.3),2F10.6)
      RETURN
      END
