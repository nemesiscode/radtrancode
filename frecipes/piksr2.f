      SUBROUTINE PIKSR2(N,ARR,BRR)
      DIMENSION ARR(N),BRR(N)
      DO 12 J=2,N
        A=ARR(J)
        B=BRR(J)
        DO 11 I=J-1,1,-1
          IF(ARR(I).LE.A)GO TO 10
          ARR(I+1)=ARR(I)
          BRR(I+1)=BRR(I)
11      CONTINUE
        I=0
10      ARR(I+1)=A
        BRR(I+1)=B
12    CONTINUE
      RETURN
      END
