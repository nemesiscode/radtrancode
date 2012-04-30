      SUBROUTINE VANDER(X,W,Q,N)
      PARAMETER (NMAX=100,ZERO=0.0,ONE=1.0)
      DIMENSION X(N),W(N),Q(N),C(NMAX)
      IF(N.EQ.1)THEN
        W(1)=Q(1)
      ELSE
        DO 11 I=1,N
          C(I)=ZERO
11      CONTINUE
        C(N)=-X(1)
        DO 13 I=2,N
          XX=-X(I)
          DO 12 J=N+1-I,N-1
            C(J)=C(J)+XX*C(J+1)
12        CONTINUE
          C(N)=C(N)+XX
13      CONTINUE
        DO 15 I=1,N
          XX=X(I)
          T=ONE
          B=ONE
          S=Q(N)
          K=N
          DO 14 J=2,N
            K1=K-1
            B=C(K)+XX*B
            S=S+Q(K1)*B
            T=XX*T+B
            K=K1
14        CONTINUE
          W(I)=S/T
15      CONTINUE
      ENDIF
      RETURN
      END
