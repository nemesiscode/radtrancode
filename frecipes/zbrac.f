      SUBROUTINE ZBRAC(FUNC,X1,X2,SUCCES)
      PARAMETER (FACTOR=1.6,NTRY=50)
      LOGICAL SUCCES
      IF(X1.EQ.X2)PAUSE 'You have to guess an initial range'
      F1=FUNC(X1)
      F2=FUNC(X2)
      SUCCES=.TRUE.
      DO 11 J=1,NTRY
        IF(F1*F2.LT.0.)RETURN
        IF(ABS(F1).LT.ABS(F2))THEN
          X1=X1+FACTOR*(X1-X2)
          F1=FUNC(X1)
        ELSE
          X2=X2+FACTOR*(X2-X1)
          F2=FUNC(X2)
        ENDIF
11    CONTINUE
      SUCCES=.FALSE.
      RETURN
      END
