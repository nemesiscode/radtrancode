      FUNCTION RTSEC(FUNC,X1,X2,XACC)
      PARAMETER (MAXIT=30)
      FL=FUNC(X1)
      F=FUNC(X2)
      IF(ABS(FL).LT.ABS(F))THEN
        RTSEC=X1
        XL=X2
        SWAP=FL
        FL=F
        F=SWAP
      ELSE
        XL=X1
        RTSEC=X2
      ENDIF
      DO 11 J=1,MAXIT
        DX=(XL-RTSEC)*F/(F-FL)
        XL=RTSEC
        FL=F
        RTSEC=RTSEC+DX
        F=FUNC(RTSEC)
        IF(ABS(DX).LT.XACC.OR.F.EQ.0.)RETURN
11    CONTINUE
      PAUSE 'RTSEC exceed maximum iterations'
      END
