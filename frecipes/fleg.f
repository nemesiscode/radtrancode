      SUBROUTINE FLEG(X,PL,NL)
      DIMENSION PL(NL)
      PL(1)=1.
      PL(2)=X
      IF(NL.GT.2) THEN
        TWOX=2.*X
        F2=X
        D=1.
        DO 11 J=3,NL
          F1=D
          F2=F2+TWOX
          D=D+1.
          PL(J)=(F2*PL(J-1)-F1*PL(J-2))/D
11      CONTINUE
      ENDIF
      RETURN
      END
