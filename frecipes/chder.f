      SUBROUTINE CHDER(A,B,C,CDER,N)
      DIMENSION C(N),CDER(N)
      CDER(N)=0.
      CDER(N-1)=2*(N-1)*C(N)
      IF(N.GE.3)THEN
        DO 11 J=N-2,1,-1
          CDER(J)=CDER(J+2)+2*J*C(J+1)
11      CONTINUE
      ENDIF
      CON=2./(B-A)
      DO 12 J=1,N
        CDER(J)=CDER(J)*CON
12    CONTINUE
      RETURN
      END
