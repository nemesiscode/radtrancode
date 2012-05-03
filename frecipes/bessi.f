      FUNCTION BESSI(N,X)
      PARAMETER(IACC=40,BIGNO=1.0E10,BIGNI=1.0E-10)
      IF (N.LT.2) PAUSE 'bad argument N in BESSI'
      TOX=2.0/X
      BIP=0.0
      BI=1.0
      BESSI=0.
      M=2*((N+INT(SQRT(FLOAT(IACC*N)))))
      DO 11 J=M,1,-1
        BIM=BIP+FLOAT(J)*TOX*BI
        BIP=BI
        BI=BIM
        IF (ABS(BI).GT.BIGNO) THEN
          BESSI=BESSI*BIGNI
          BI=BI*BIGNI
          BIP=BIP*BIGNI
        ENDIF
        IF (J.EQ.N) BESSI=BIP
11    CONTINUE
      BESSI=BESSI*BESSI0(X)/BI
      RETURN
      END
