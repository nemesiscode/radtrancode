      SUBROUTINE INTERP(X,Y,N,YOUT,XIN)
C     $Id: interp.f,v 1.1.1.1 2000-08-17 09:26:54 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  INTERP: subroutine to perform interpolation
C
C_KEYS:   RADTRAN,PROG,VMS
C
C_DESCR:  general interpolation routine for GENLBL support programs
C
C_ARGS:   
C
C_FILES : none
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:   25aug92 SBC copied from VERINT
C
C_END:
C--------------------------------------------------------------
C     simple linear interpolation for preliminary use
      INTEGER I,N
      REAL X(N),Y(N),XIN,YOUT
C
      IF(X(1).LT.X(N))THEN
        DO 10 I=1,N
        IF(X(I).GT.XIN)GOTO 30
10      CONTINUE
       ELSE
        DO 20 I=1,N
        IF(X(I).LT.XIN)GOTO 30
20      CONTINUE
        END IF
      I=N
30    CONTINUE
      IF(I.EQ.1)I=2
      YOUT=Y(I-1)+(Y(I)-Y(I-1))*(XIN-X(I-1))/(X(I)-X(I-1))
      RETURN
      END 
