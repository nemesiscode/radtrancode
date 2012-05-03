      SUBROUTINE QROOT(P,N,B,C,EPS)
      PARAMETER (NMAX=20,ITMAX=20,TINY=1.0E-6)
      DIMENSION P(N),Q(NMAX),D(3),REM(NMAX),QQ(NMAX)
      D(3)=1.
      DO 12 ITER=1,ITMAX
        D(2)=B
        D(1)=C
        CALL POLDIV(P,N,D,3,Q,REM)
        S=REM(1)
        R=REM(2)
        CALL POLDIV(Q,N-1,D,3,QQ,REM)
        SC=-REM(1)
        RC=-REM(2)
        DO 11 I=N-1,1,-1
          Q(I+1)=Q(I)
11      CONTINUE
        Q(1)=0.
        CALL POLDIV(Q,N,D,3,QQ,REM)
        SB=-REM(1)
        RB=-REM(2)
        DIV=1./(SB*RC-SC*RB)
        DELB=(R*SC-S*RC)*DIV
        DELC=(-R*SB+S*RB)*DIV
        B=B+DELB
        C=C+DELC
        IF((ABS(DELB).LE.EPS*ABS(B).OR.ABS(B).LT.TINY)
     *      .AND.(ABS(DELC).LE.EPS*ABS(C)
     *      .OR.ABS(C).LT.TINY)) RETURN
12    CONTINUE
      PAUSE 'too many iterations in QROOT'
      END
