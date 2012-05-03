      SUBROUTINE KSONE(DATA,N,FUNC,D,PROB)
      DIMENSION DATA(N)
      CALL SORT(N,DATA)
      EN=N
      D=0.
      FO=0.
      DO 11 J=1,N
        FN=J/EN
        FF=FUNC(DATA(J))
        DT=AMAX1(ABS(FO-FF),ABS(FN-FF))
        IF(DT.GT.D)D=DT
        FO=FN
11    CONTINUE
      PROB=PROBKS(SQRT(EN)*D)
      RETURN
      END
