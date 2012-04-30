      SUBROUTINE VERINT(X,Y,N,YOUT,XIN)
    !   simple linear interpolation for preliminary use
    !   note that it is important that interpolation used is
    !   consistent with that implicit in the integration
    !   scheme
      INTEGER I,N
      REAL X(N),Y(N),XIN,YOUT

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
