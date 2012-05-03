      SUBROUTINE QTRAP(FUNC,A,B,S)
      PARAMETER (EPS=1.E-6, JMAX=20)
      OLDS=-1.E30
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S,J)
        IF (ABS(S-OLDS).LT.EPS*ABS(OLDS)) RETURN
        OLDS=S
11    CONTINUE
      PAUSE 'Too many steps.'
      END
