      SUBROUTINE KENDL2(TAB,I,J,IP,JP,TAU,Z,PROB)
      DIMENSION TAB(IP,JP)
      EN1=0.
      EN2=0.
      S=0.
      NN=I*J
      POINTS=TAB(I,J)
      DO 12 K=0,NN-2
        KI=K/J
        KJ=K-J*KI
        POINTS=POINTS+TAB(KI+1,KJ+1)
        DO 11 L=K+1,NN-1
          LI=L/J
          LJ=L-J*LI
          M1=LI-KI
          M2=LJ-KJ
          MM=M1*M2
          PAIRS=TAB(KI+1,KJ+1)*TAB(LI+1,LJ+1)
          IF(MM.NE.0)THEN
            EN1=EN1+PAIRS
            EN2=EN2+PAIRS
            IF(MM.GT.0)THEN
              S=S+PAIRS
            ELSE
              S=S-PAIRS
            ENDIF
          ELSE
            IF(M1.NE.0)EN1=EN1+PAIRS
            IF(M2.NE.0)EN2=EN2+PAIRS
          ENDIF
11      CONTINUE
12    CONTINUE
      TAU=S/SQRT(EN1*EN2)
      VAR=(4.*POINTS+10.)/(9.*POINTS*(POINTS-1.))
      Z=TAU/SQRT(VAR)
      PROB=ERFCC(ABS(Z)/1.4142136)
      RETURN
      END
