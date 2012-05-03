      FUNCTION FACTLN(N)
      DIMENSION A(100)
      DATA A/100*-1./
      IF (N.LT.0) PAUSE 'negative factorial'
      IF (N.LE.99) THEN
        IF (A(N+1).LT.0.) A(N+1)=GAMMLN(N+1.)
        FACTLN=A(N+1)
      ELSE
        FACTLN=GAMMLN(N+1.)
      ENDIF
      RETURN
      END
