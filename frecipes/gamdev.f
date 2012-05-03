      FUNCTION GAMDEV(IA,IDUM)
      IF(IA.LT.1)PAUSE
      IF(IA.LT.6)THEN
        X=1.
        DO 11 J=1,IA
          X=X*RAN1(IDUM)
11      CONTINUE
        X=-LOG(X)
      ELSE
1         V1=2.*RAN1(IDUM)-1.
          V2=2.*RAN1(IDUM)-1.
        IF(V1**2+V2**2.GT.1.)GO TO 1
          Y=V2/V1
          AM=IA-1
          S=SQRT(2.*AM+1.)
          X=S*Y+AM
        IF(X.LE.0.)GO TO 1
          E=(1.+Y**2)*EXP(AM*LOG(X/AM)-S*Y)
        IF(RAN1(IDUM).GT.E)GO TO 1
      ENDIF
      GAMDEV=X
      RETURN
      END
