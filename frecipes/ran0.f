      FUNCTION RAN0(IDUM)
      DIMENSION V(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        ISEED=ABS(IDUM)
        IDUM=1
        DO 11 J=1,97
          DUM=RAN(ISEED)
11      CONTINUE
        DO 12 J=1,97
          V(J)=RAN(ISEED)
12      CONTINUE
        Y=RAN(ISEED)
      ENDIF
      J=1+INT(97.*Y)
      IF(J.GT.97.OR.J.LT.1)PAUSE
      Y=V(J)
      RAN0=Y
      V(J)=RAN(ISEED)
      RETURN
      END
