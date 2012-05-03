      FUNCTION ERFC(X)
      IF(X.LT.0.)THEN
        ERFC=1.+GAMMP(.5,X**2)
      ELSE
        ERFC=GAMMQ(.5,X**2)
      ENDIF
      RETURN
      END
