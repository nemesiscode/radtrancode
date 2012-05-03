      SUBROUTINE BKSUB(NE,NB,JF,K1,K2,C,NCI,NCJ,NCK)
      DIMENSION C(NCI,NCJ,NCK)
      NBF=NE-NB
      DO 13 K=K2,K1,-1
        KP=K+1
        DO 12 J=1,NBF
          XX=C(J,JF,KP)
          DO 11 I=1,NE
            C(I,JF,K)=C(I,JF,K)-C(I,J,K)*XX
11        CONTINUE
12      CONTINUE
13    CONTINUE
      DO 16 K=K1,K2
        KP=K+1
        DO 14 I=1,NB
          C(I,1,K)=C(I+NBF,JF,K)
14      CONTINUE
        DO 15 I=1,NBF
          C(I+NB,1,K)=C(I,JF,KP)
15      CONTINUE
16    CONTINUE
      RETURN
      END
