C     $Id: verint.f,v 1.1.1.1 2000-08-17 09:26:56 irwin Exp $
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
      SUBROUTINE VERINT(X,Y,N,YOUT,XIN)
C     simple linear interpolation for preliminary use
C     note that it is important that interpolation used is
C     consistent with that implicit in the integration
C     scheme
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
      if(X(I-1).eq.X(I))then YOUT=Y(I)
      if (isnan(YOUT)) then
       YOUT=Y(I)
C       print*,"XIN,XI,XI-1,YI-1,YI",XIN,X(I),X(I-1),Y(I-1),Y(I)
C       print*,'YOUT is not a number!'
C       print*,'I = ',I
C       do I=1,N
C        print*,I,X(I),Y(I)
C       enddo
C       stop
      end if

      RETURN
      END
