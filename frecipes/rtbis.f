      FUNCTION RTBIS(FUNC,X1,X2,XACC)
      PARAMETER (JMAX=40)
      FMID=FUNC(X2)
      F=FUNC(X1)
      IF(F*FMID.GE.0.) PAUSE 'Root must be bracketed for bisection.'
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        FMID=FUNC(XMID)
        IF(FMID.LE.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
      PAUSE 'too many bisections'
      END
