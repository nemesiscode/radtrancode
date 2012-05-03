      INTEGER FUNCTION PINDEX(X,NX,X1)
C     **********************************************************
C     Function for looking up closest point beliw X1 in array X(NX)
C
C     Pat Irwin	2/3/12	Updated for Radtrans2.0
C
C     **********************************************************
      INTEGER NX
      REAL X(NX),X1
      IF(X1.GT.X(NX))THEN
C       PRINT*,'Warning in INDEX. X1>XMAX',X1,X(NX)
       X1=X(NX)
      END IF
      IF(X1.LT.X(1))THEN
C       PRINT*,'Warning in INDEX. X1>XMIN',X1,X(1)
       X1=X(1)
      END IF

      J=0

      DO 10 I=1,NX-1
       IF(X1.GE.X(I))THEN
        J=I
       ELSE
        GOTO 20
       END IF
10    CONTINUE
20    CONTINUE

      PINDEX=J
      RETURN
      END
