      SUBROUTINE POLDIV(U,N,V,NV,Q,R)
      DIMENSION U(N),V(NV),Q(N),R(N)
      DO 11 J=1,N
        R(J)=U(J)
        Q(J)=0.
11    CONTINUE
      DO 13 K=N-NV,0,-1
        Q(K+1)=R(NV+K)/V(NV)
        DO 12 J=NV+K-1,K+1,-1
          R(J)=R(J)-Q(K+1)*V(J-K)
12      CONTINUE
13    CONTINUE
      R(NV)=0.
      RETURN
      END
