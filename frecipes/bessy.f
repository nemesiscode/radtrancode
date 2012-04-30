      FUNCTION BESSY(N,X)
      IF(N.LT.2)PAUSE 'bad argument N in BESSY'
      TOX=2./X
      BY=BESSY1(X)
      BYM=BESSY0(X)
      DO 11 J=1,N-1
        BYP=J*TOX*BY-BYM
        BYM=BY
        BY=BYP
11    CONTINUE
      BESSY=BY
      RETURN
      END
