      SUBROUTINE CHINT(A,B,C,CINT,N)
      DIMENSION C(N),CINT(N)
      CON=0.25*(B-A)
      SUM=0.
      FAC=1.
      DO 11 J=2,N-1
        CINT(J)=CON*(C(J-1)-C(J+1))/(J-1)
        SUM=SUM+FAC*CINT(J)
        FAC=-FAC
11    CONTINUE
      CINT(N)=CON*C(N-1)/(N-1)
      SUM=SUM+FAC*CINT(N)
      CINT(1)=2.*SUM
      RETURN
      END
