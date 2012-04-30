      SUBROUTINE ADI(A,B,C,D,E,F,G,U,JMAX,K,ALPHA,BETA,EPS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(JJ=100,KK=6,NRR=32,MAXITS=100,ZERO=0.D0,TWO=2.D0,HALF=.5
     *D0)
      DIMENSION A(JMAX,JMAX),B(JMAX,JMAX),C(JMAX,JMAX),D(JMAX,JMAX),
     *    E(JMAX,JMAX),F(JMAX,JMAX),G(JMAX,JMAX),U(JMAX,JMAX),
     *    AA(JJ),BB(JJ),CC(JJ),RR(JJ),UU(JJ),PSI(JJ,JJ),
     *    ALPH(KK),BET(KK),R(NRR),S(NRR,KK)
      IF(JMAX.GT.JJ)PAUSE 'Increase JJ'
      IF(K.GT.KK-1)PAUSE 'Increase KK'
      K1=K+1
      NR=2**K
      ALPH(1)=ALPHA
      BET(1)=BETA
      DO 11 J=1,K
        ALPH(J+1)=SQRT(ALPH(J)*BET(J))
        BET(J+1)=HALF*(ALPH(J)+BET(J))
11    CONTINUE
      S(1,1)=SQRT(ALPH(K1)*BET(K1))
      DO 13 J=1,K
        AB=ALPH(K1-J)*BET(K1-J)
        DO 12 N=1,2**(J-1)
          DISC=SQRT(S(N,J)**2-AB)
          S(2*N,J+1)=S(N,J)+DISC
          S(2*N-1,J+1)=AB/S(2*N,J+1)
12      CONTINUE
13    CONTINUE
      DO 14 N=1,NR
        R(N)=S(N,K1)
14    CONTINUE
      ANORMG=ZERO
      DO 16 J=2,JMAX-1
        DO 15 L=2,JMAX-1
          ANORMG=ANORMG+ABS(G(J,L))
          PSI(J,L)=-D(J,L)*U(J,L-1)+(R(1)-E(J,L))*U(J,L)
     *        -F(J,L)*U(J,L+1)
15      CONTINUE
16    CONTINUE
      NITS=MAXITS/NR
      DO 27 KITS=1,NITS
        DO 24 N=1,NR
          IF(N.EQ.NR)THEN
            NEXT=1
          ELSE
            NEXT=N+1
          ENDIF
          RFACT=R(N)+R(NEXT)
          DO 19 L=2,JMAX-1
            DO 17 J=2,JMAX-1
              AA(J-1)=A(J,L)
              BB(J-1)=B(J,L)+R(N)
              CC(J-1)=C(J,L)
              RR(J-1)=PSI(J,L)-G(J,L)
17          CONTINUE
            CALL TRIDAG(AA,BB,CC,RR,UU,JMAX-2)
            DO 18 J=2,JMAX-1
              PSI(J,L)=-PSI(J,L)+TWO*R(N)*UU(J-1)
18          CONTINUE
19        CONTINUE
          DO 23 J=2,JMAX-1
            DO 21 L=2,JMAX-1
              AA(L-1)=D(J,L)
              BB(L-1)=E(J,L)+R(N)
              CC(L-1)=F(J,L)
              RR(L-1)=PSI(J,L)
21          CONTINUE
            CALL TRIDAG(AA,BB,CC,RR,UU,JMAX-2)
            DO 22 L=2,JMAX-1
              U(J,L)=UU(L-1)
              PSI(J,L)=-PSI(J,L)+RFACT*UU(L-1)
22          CONTINUE
23        CONTINUE
24      CONTINUE
        ANORM=ZERO
        DO 26 J=2,JMAX-1
          DO 25 L=2,JMAX-1
            RESID=A(J,L)*U(J-1,L)+(B(J,L)+E(J,L))*U(J,L)
     *          +C(J,L)*U(J+1,L)+D(J,L)*U(J,L-1)
     *          +F(J,L)*U(J,L+1)+G(J,L)
            ANORM=ANORM+ABS(RESID)
25        CONTINUE
26      CONTINUE
        IF(ANORM.LT.EPS*ANORMG)RETURN
27    CONTINUE
      PAUSE 'MAXITS exceeded'
      END
