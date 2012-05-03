      SUBROUTINE AVEVAR(DATA,N,AVE,VAR)
      DIMENSION DATA(N)
      AVE=0.0
      VAR=0.0
      DO 11 J=1,N
        AVE=AVE+DATA(J)
11    CONTINUE
      AVE=AVE/N
      DO 12 J=1,N
        S=DATA(J)-AVE
        VAR=VAR+S*S
12    CONTINUE
      VAR=VAR/(N-1)
      RETURN
      END
