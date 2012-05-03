      SUBROUTINE ECLAZZ(NF,N,EQUIV)
      LOGICAL EQUIV
      DIMENSION NF(N)
      NF(1)=1
      DO 12 JJ=2,N
        NF(JJ)=JJ
        DO 11 KK=1,JJ-1
          NF(KK)=NF(NF(KK))
          IF (EQUIV(JJ,KK)) NF(NF(NF(KK)))=JJ
11      CONTINUE
12    CONTINUE
      DO 13 JJ=1,N
        NF(JJ)=NF(NF(JJ))
13    CONTINUE
      RETURN
      END
