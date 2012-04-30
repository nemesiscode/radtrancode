      SUBROUTINE SOR(A,B,C,D,E,F,U,JMAX,RJAC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(JMAX,JMAX),B(JMAX,JMAX),C(JMAX,JMAX),
     *    D(JMAX,JMAX),E(JMAX,JMAX),F(JMAX,JMAX),U(JMAX,JMAX)
      PARAMETER(MAXITS=1000,EPS=1.D-5,ZERO=0.D0,HALF=.5D0,QTR=.25D0,ONE=
     *1.D0)
      ANORMF=ZERO
      DO 12 J=2,JMAX-1
        DO 11 L=2,JMAX-1
          ANORMF=ANORMF+ABS(F(J,L))
11      CONTINUE
12    CONTINUE
      OMEGA=ONE
      DO 15 N=1,MAXITS
        ANORM=ZERO
        DO 14 J=2,JMAX-1
          DO 13 L=2,JMAX-1
            IF(MOD(J+L,2).EQ.MOD(N,2))THEN
              RESID=A(J,L)*U(J+1,L)+B(J,L)*U(J-1,L)+
     *            C(J,L)*U(J,L+1)+D(J,L)*U(J,L-1)+
     *            E(J,L)*U(J,L)-F(J,L)
              ANORM=ANORM+ABS(RESID)
              U(J,L)=U(J,L)-OMEGA*RESID/E(J,L)
            ENDIF
13        CONTINUE
14      CONTINUE
        IF(N.EQ.1) THEN
          OMEGA=ONE/(ONE-HALF*RJAC**2)
        ELSE
          OMEGA=ONE/(ONE-QTR*RJAC**2*OMEGA)
        ENDIF
        IF((N.GT.1).AND.(ANORM.LT.EPS*ANORMF))RETURN
15    CONTINUE
      PAUSE 'MAXITS exceeded'
      END
