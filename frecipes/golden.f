      FUNCTION GOLDEN(AX,BX,CX,F,TOL,XMIN)
      PARAMETER (R=.61803399,C=.38196602)
      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
        X1=BX
        X2=BX+C*(CX-BX)
      ELSE
        X2=BX
        X1=BX-C*(BX-AX)
      ENDIF
      F1=F(X1)
      F2=F(X2)
1     IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2)))THEN
        IF(F2.LT.F1)THEN
          X0=X1
          X1=X2
          X2=R*X1+C*X3
          F0=F1
          F1=F2
          F2=F(X2)
        ELSE
          X3=X2
          X2=X1
          X1=R*X2+C*X0
          F3=F2
          F2=F1
          F1=F(X1)
        ENDIF
      GOTO 1
      ENDIF
      IF(F1.LT.F2)THEN
        GOLDEN=F1
        XMIN=X1
      ELSE
        GOLDEN=F2
        XMIN=X2
      ENDIF
      RETURN
      END
