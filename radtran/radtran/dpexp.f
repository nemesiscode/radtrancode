      FUNCTION DPEXP(X)
C     $Id: dpexp.f,v 1.2 2011-06-17 15:40:26 irwin Exp $
C     simple function to perform exponentiation in double precision
C     needed on Decstations since exp(-large number) gives floating point
C     overflow instead of underflow so can't be just ignored.
C     Has no useful effect under VMS since exponent range is the same for
C     double and single precision and the error is an underflow anyway.
      Y = MAX(X,-500.0)
      IF(Y.GE.88.7)THEN
       Y=88.7
      ENDIF
      DPEXP=SNGL(DEXP(DBLE(Y)))
      RETURN
      END
