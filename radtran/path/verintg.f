C     $Id:
C--------------------------------------------------------------
C_TITLE:  VERINT: subroutines to perform vertical interpolation
C
C_KEYS:   RADTRAN,PROG,VMS
C
C_DESCR:  two subroutines to perform vertical interpolation and integration
C         of atmospheric profiles. They are stored together because the
C        integration must implicitly assume the same interpolation scheme
C
C_ARGS:   described in each routine
C
C_FILES : none
C
C_CALLS:  FILE,PROMPT,ASKYN
C
C_BUGS:
C
C_HIST:   2jul92 SBC Original version
C
C_END:
C--------------------------------------------------------------
C
      SUBROUTINE VERINTG(X,Y,N,YOUT,XIN,J,F)
C     simple linear interpolation for preliminary use
C     note that it is important that interpolation used is
C     consistent with that implicit in the integration
C     scheme
      INTEGER I,N,J
      REAL X(N),Y(N),XIN,YOUT,F
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
      J=I-1
      F = (XIN-X(I-1))/(X(I)-X(I-1))
      YOUT=(1.0-F)*Y(I-1)+F*Y(I)
      RETURN
      END
