      SUBROUTINE ZBRAK(FX,X1,X2,N,XB1,XB2,NB)
      DIMENSION XB1(1),XB2(1)
      NBB=NB
      NB=0
      X=X1
      DX=(X2-X1)/N
      FP=FX(X)
      DO 11 I=1,N
        X=X+DX
        FC=FX(X)
        IF(FC*FP.LT.0.) THEN
          NB=NB+1
          XB1(NB)=X-DX
          XB2(NB)=X
        ENDIF
        FP=FC
        IF(NBB.EQ.NB)RETURN
11    CONTINUE
      RETURN
      END
