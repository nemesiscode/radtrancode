      FUNCTION EL2(X,QQC,AA,BB)
      PARAMETER(PI=3.14159265, CA=.0003, CB=1.E-9)
      IF(X.EQ.0.)THEN
        EL2=0.
      ELSE IF(QQC.NE.0.)THEN
        QC=QQC
        A=AA
        B=BB
        C=X**2
        D=1.+C
        P=SQRT((1.+QC**2*C)/D)
        D=X/D
        C=D/(2.*P)
        Z=A-B
        EYE=A
        A=0.5*(B+A)
        Y=ABS(1./X)
        F=0.
        L=0
        EM=1.
        QC=ABS(QC)
1       B=EYE*QC+B
        E=EM*QC
        G=E/P
        D=F*G+D
        F=C
        EYE=A
        P=G+P
        C=0.5*(D/P+C)
        G=EM
        EM=QC+EM
        A=0.5*(B/EM+A)
        Y=-E/Y+Y
        IF(Y.EQ.0.)Y=SQRT(E)*CB
        IF(ABS(G-QC).GT.CA*G)THEN
          QC=SQRT(E)*2.
          L=L+L
          IF(Y.LT.0.)L=L+1
          GO TO 1
        ENDIF
        IF(Y.LT.0.)L=L+1
        E=(ATAN(EM/Y)+PI*L)*A/EM
        IF(X.LT.0.)E=-E
        EL2=E+C*Z
      ELSE
        PAUSE 'failure in EL2'
      ENDIF
      RETURN
      END
