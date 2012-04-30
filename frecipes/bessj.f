      FUNCTION BESSJ(N,X)
      PARAMETER (IACC=40,BIGNO=1.E10,BIGNI=1.E-10)
      IF(N.LT.2)PAUSE 'bad argument N in BESSJ'
      TOX=2./X
      IF(X.GT.FLOAT(N))THEN
        BJM=BESSJ0(X)
        BJ=BESSJ1(X)
        DO 11 J=1,N-1
          BJP=J*TOX*BJ-BJM
          BJM=BJ
          BJ=BJP
11      CONTINUE
        BESSJ=BJ
      ELSE
        M=2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
        BESSJ=0.
        JSUM=0
        SUM=0.
        BJP=0.
        BJ=1.
        DO 12 J=M,1,-1
          BJM=J*TOX*BJ-BJP
          BJP=BJ
          BJ=BJM
          IF(ABS(BJ).GT.BIGNO)THEN
            BJ=BJ*BIGNI
            BJP=BJP*BIGNI
            BESSJ=BESSJ*BIGNI
            SUM=SUM*BIGNI
          ENDIF
          IF(JSUM.NE.0)SUM=SUM+BJ
          JSUM=1-JSUM
          IF(J.EQ.N)BESSJ=BJP
12      CONTINUE
        SUM=2.*SUM-BJ
        BESSJ=BESSJ/SUM
      ENDIF
      RETURN
      END
