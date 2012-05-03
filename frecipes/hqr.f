      SUBROUTINE HQR(A,N,NP,WR,WI)
      DIMENSION A(NP,NP),WR(NP),WI(NP)
      ANORM=ABS(A(1,1))
      DO 12 I=2,N
        DO 11 J=I-1,N
          ANORM=ANORM+ABS(A(I,J))
11      CONTINUE
12    CONTINUE
      NN=N
      T=0.
1     IF(NN.GE.1)THEN
        ITS=0
2       DO 13 L=NN,2,-1
          S=ABS(A(L-1,L-1))+ABS(A(L,L))
          IF(S.EQ.0.)S=ANORM
          IF(ABS(A(L,L-1))+S.EQ.S)GO TO 3
13      CONTINUE
        L=1
3       X=A(NN,NN)
        IF(L.EQ.NN)THEN
          WR(NN)=X+T
          WI(NN)=0.
          NN=NN-1
        ELSE
          Y=A(NN-1,NN-1)
          W=A(NN,NN-1)*A(NN-1,NN)
          IF(L.EQ.NN-1)THEN
            P=0.5*(Y-X)
            Q=P**2+W
            Z=SQRT(ABS(Q))
            X=X+T
            IF(Q.GE.0.)THEN
              Z=P+SIGN(Z,P)
              WR(NN)=X+Z
              WR(NN-1)=WR(NN)
              IF(Z.NE.0.)WR(NN)=X-W/Z
              WI(NN)=0.
              WI(NN-1)=0.
            ELSE
              WR(NN)=X+P
              WR(NN-1)=WR(NN)
              WI(NN)=Z
              WI(NN-1)=-Z
            ENDIF
            NN=NN-2
          ELSE
            IF(ITS.EQ.30)PAUSE 'too many iterations'
            IF(ITS.EQ.10.OR.ITS.EQ.20)THEN
              T=T+X
              DO 14 I=1,NN
                A(I,I)=A(I,I)-X
14            CONTINUE
              S=ABS(A(NN,NN-1))+ABS(A(NN-1,NN-2))
              X=0.75*S
              Y=X
              W=-0.4375*S**2
            ENDIF
            ITS=ITS+1
            DO 15 M=NN-2,L,-1
              Z=A(M,M)
              R=X-Z
              S=Y-Z
              P=(R*S-W)/A(M+1,M)+A(M,M+1)
              Q=A(M+1,M+1)-Z-R-S
              R=A(M+2,M+1)
              S=ABS(P)+ABS(Q)+ABS(R)
              P=P/S
              Q=Q/S
              R=R/S
              IF(M.EQ.L)GO TO 4
              U=ABS(A(M,M-1))*(ABS(Q)+ABS(R))
              V=ABS(P)*(ABS(A(M-1,M-1))+ABS(Z)+ABS(A(M+1,M+1)))
              IF(U+V.EQ.V)GO TO 4
15          CONTINUE
4           DO 16 I=M+2,NN
              A(I,I-2)=0.
              IF (I.NE.M+2) A(I,I-3)=0.
16          CONTINUE
            DO 19 K=M,NN-1
              IF(K.NE.M)THEN
                P=A(K,K-1)
                Q=A(K+1,K-1)
                R=0.
                IF(K.NE.NN-1)R=A(K+2,K-1)
                X=ABS(P)+ABS(Q)+ABS(R)
                IF(X.NE.0.)THEN
                  P=P/X
                  Q=Q/X
                  R=R/X
                ENDIF
              ENDIF
              S=SIGN(SQRT(P**2+Q**2+R**2),P)
              IF(S.NE.0.)THEN
                IF(K.EQ.M)THEN
                  IF(L.NE.M)A(K,K-1)=-A(K,K-1)
                ELSE
                  A(K,K-1)=-S*X
                ENDIF
                P=P+S
                X=P/S
                Y=Q/S
                Z=R/S
                Q=Q/P
                R=R/P
                DO 17 J=K,NN
                  P=A(K,J)+Q*A(K+1,J)
                  IF(K.NE.NN-1)THEN
                    P=P+R*A(K+2,J)
                    A(K+2,J)=A(K+2,J)-P*Z
                  ENDIF
                  A(K+1,J)=A(K+1,J)-P*Y
                  A(K,J)=A(K,J)-P*X
17              CONTINUE
                DO 18 I=L,MIN(NN,K+3)
                  P=X*A(I,K)+Y*A(I,K+1)
                  IF(K.NE.NN-1)THEN
                    P=P+Z*A(I,K+2)
                    A(I,K+2)=A(I,K+2)-P*R
                  ENDIF
                  A(I,K+1)=A(I,K+1)-P*Q
                  A(I,K)=A(I,K)-P
18              CONTINUE
              ENDIF
19          CONTINUE
            GO TO 2
          ENDIF
        ENDIF
      GO TO 1
      ENDIF
      RETURN
      END
