      SUBROUTINE KSTWO(DATA1,N1,DATA2,N2,D,PROB)
      DIMENSION DATA1(N1),DATA2(N2)
      CALL SORT(N1,DATA1)
      CALL SORT(N2,DATA2)
      EN1=N1
      EN2=N2
      J1=1
      J2=1
      FO1=0.
      FO2=0.
      D=0.
1     IF(J1.LE.N1.AND.J2.LE.N2)THEN
        IF(DATA1(J1).LT.DATA2(J2))THEN
          FN1=J1/EN1
          DT=AMAX1(ABS(FN1-FO2),ABS(FO1-FO2))
          IF(DT.GT.D)D=DT
          FO1=FN1
          J1=J1+1
        ELSE
          FN2=J2/EN2
          DT=AMAX1(ABS(FN2-FO1),ABS(FO2-FO1))
          IF(DT.GT.D)D=DT
          FO2=FN2
          J2=J2+1
        ENDIF
      GO TO 1
      ENDIF
      PROB=PROBKS(SQRT(EN1*EN2/(EN1+EN2))*D)
      RETURN
      END
