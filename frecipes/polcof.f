      SUBROUTINE POLCOF(XA,YA,N,COF)
      PARAMETER (NMAX=15)
      DIMENSION XA(N),YA(N),COF(N),X(NMAX),Y(NMAX)
      DO 11 J=1,N
        X(J)=XA(J)
        Y(J)=YA(J)
11    CONTINUE
      DO 14 J=1,N
        CALL POLINT(X,Y,N+1-J,0.,COF(J),DY)
        XMIN=1.E38
        K=0
        DO 12 I=1,N+1-J
          IF (ABS(X(I)).LT.XMIN)THEN
            XMIN=ABS(X(I))
            K=I
          ENDIF
          IF(X(I).NE.0.)Y(I)=(Y(I)-COF(J))/X(I)
12      CONTINUE
        IF (K.LT.N+1-J) THEN
          DO 13 I=K+1,N+1-J
            Y(I-1)=Y(I)
            X(I-1)=X(I)
13        CONTINUE
        ENDIF
14    CONTINUE
      RETURN
      END
