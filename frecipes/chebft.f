      SUBROUTINE CHEBFT(A,B,C,N,FUNC)
      REAL*8 PI
      PARAMETER (NMAX=50, PI=3.141592653589793D0)
      REAL*8 SUM
      DIMENSION C(N),F(NMAX)
      BMA=0.5*(B-A)
      BPA=0.5*(B+A)
      DO 11 K=1,N
        Y=COS(PI*(K-0.5)/N)
        F(K)=FUNC(Y*BMA+BPA)
11    CONTINUE
      FAC=2./N
      DO 13 J=1,N
        SUM=0.D0
        DO 12 K=1,N
          SUM=SUM+F(K)*COS((PI*(J-1))*((K-0.5D0)/N))
12      CONTINUE
        C(J)=FAC*SUM
13    CONTINUE
      RETURN
      END
