      SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)
      PARAMETER (NMAX=100)
      DIMENSION A(NP,NP),ALUD(NP,NP),INDX(N),B(N),X(N),R(NMAX)
      REAL*8 SDP
      DO 12 I=1,N
        SDP=-B(I)
        DO 11 J=1,N
          SDP=SDP+DBLE(A(I,J))*DBLE(X(J))
11      CONTINUE
        R(I)=SDP
12    CONTINUE
      CALL LUBKSB(ALUD,N,NP,INDX,R)
      DO 13 I=1,N
        X(I)=X(I)-R(I)
13    CONTINUE
      RETURN
      END
