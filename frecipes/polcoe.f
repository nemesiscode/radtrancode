      SUBROUTINE POLCOE(X,Y,N,COF)
      PARAMETER (NMAX=15)
      DIMENSION X(N),Y(N),COF(N),S(NMAX)
      DO 11 I=1,N
        S(I)=0.
        COF(I)=0.
11    CONTINUE
      S(N)=-X(1)
      DO 13 I=2,N
        DO 12 J=N+1-I,N-1
          S(J)=S(J)-X(I)*S(J+1)
12      CONTINUE
        S(N)=S(N)-X(I)
13    CONTINUE
      DO 16 J=1,N
        PHI=N
        DO 14 K=N-1,1,-1
          PHI=K*S(K+1)+X(J)*PHI
14      CONTINUE
        FF=Y(J)/PHI
        B=1.
        DO 15 K=N,1,-1
          COF(K)=COF(K)+B*FF
          B=S(K)+X(J)*B
15      CONTINUE
16    CONTINUE
      RETURN
      END
