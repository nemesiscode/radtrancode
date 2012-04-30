      SUBROUTINE LAGUER(A,M,X,EPS,POLISH)
      COMPLEX A(*),X,DX,X1,B,D,F,G,H,SQ,GP,GM,G2,ZERO
      LOGICAL POLISH
      PARAMETER (ZERO=(0.,0.),EPSS=6.E-8,MAXIT=100)
      DXOLD=CABS(X)
      DO 12 ITER=1,MAXIT
        B=A(M+1)
        ERR=CABS(B)
        D=ZERO
        F=ZERO
        ABX=CABS(X)
        DO 11 J=M,1,-1
          F=X*F+D
          D=X*D+B
          B=X*B+A(J)
          ERR=CABS(B)+ABX*ERR
11      CONTINUE
        ERR=EPSS*ERR
        IF(CABS(B).LE.ERR) THEN
          DX=ZERO
          RETURN
        ELSE
          G=D/B
          G2=G*G
          H=G2-2.*F/B
          SQ=CSQRT((M-1)*(M*H-G2))
          GP=G+SQ
          GM=G-SQ
          IF(CABS(GP).LT.CABS(GM)) GP=GM
          DX=M/GP
        ENDIF
        X1=X-DX
        IF(X.EQ.X1)RETURN
        X=X1
        CDX=CABS(DX)
        IF(ITER.GT.6.AND.CDX.GE.DXOLD)RETURN
        DXOLD=CDX
        IF(.NOT.POLISH)THEN
          IF(CABS(DX).LE.EPS*CABS(X))RETURN
        ENDIF
12    CONTINUE
      PAUSE 'too many iterations'
      RETURN
      END
