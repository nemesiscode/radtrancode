      SUBROUTINE DGSER(GAMSER,A,X,GLN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ITMAX=100,EPS=3.D-7)
      GLN=DGAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*DEXP(-X+A*DLOG(X)-GLN)
      RETURN
      END
