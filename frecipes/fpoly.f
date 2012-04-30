      SUBROUTINE FPOLY(X,P,NP)
      DIMENSION P(NP)
      P(1)=1.
      DO 11 J=2,NP
        P(J)=P(J-1)*X
11    CONTINUE
      RETURN
      END
