      SUBROUTINE RKDUMB(VSTART,NVAR,X1,X2,NSTEP,DERIVS)
      PARAMETER (NMAX=10)
      COMMON /PATH/ XX(200),Y(10,200)
      DIMENSION VSTART(NVAR),V(NMAX),DV(NMAX)
      external DERIVS
      DO 11 I=1,NVAR
        V(I)=VSTART(I)
        Y(I,1)=V(I)
11    CONTINUE
      XX(1)=X1
      X=X1
      H=(X2-X1)/NSTEP
      DO 13 K=1,NSTEP
        CALL DERIVS(X,V,DV)
        CALL RK4(V,DV,NVAR,X,H,V,DERIVS)
        IF(X+H.EQ.X)PAUSE 'Stepsize not significant in RKDUMB.'
        X=X+H
        XX(K+1)=X
        DO 12 I=1,NVAR
          Y(I,K+1)=V(I)
12      CONTINUE
13    CONTINUE
      RETURN
      END
