      SUBROUTINE REVCST(X,Y,IORDER,NCITY,N,DE)
C Rev. 12/13/85
      DIMENSION X(NCITY),Y(NCITY),IORDER(NCITY),N(6),XX(4),YY(4)
      ALEN(X1,X2,Y1,Y2)=SQRT((X2-X1)**2+(Y2-Y1)**2)
      N(3)=1+MOD((N(1)+NCITY-2),NCITY)
      N(4)=1+MOD(N(2),NCITY)
      DO 11 J=1,4
          II=IORDER(N(J))
          XX(J)=X(II)
          YY(J)=Y(II)
11    CONTINUE
      DE=-ALEN(XX(1),XX(3),YY(1),YY(3))-ALEN(XX(2),XX(4),YY(2),YY(4))
     *     +ALEN(XX(1),XX(4),YY(1),YY(4))+ALEN(XX(2),XX(3),YY(2),YY(3))
      RETURN
      END

      SUBROUTINE REVERS(IORDER,NCITY,N)
      DIMENSION IORDER(NCITY),N(6)
      NN=(1+MOD(N(2)-N(1)+NCITY,NCITY))/2
      DO 11 J=1,NN
          K=1+MOD((N(1)+J-2),NCITY)
          L=1+MOD((N(2)-J+NCITY),NCITY)    
          ITMP=IORDER(K)
          IORDER(K)=IORDER(L)
          IORDER(L)=ITMP
11    CONTINUE
      RETURN
      END

      SUBROUTINE TRNCST(X,Y,IORDER,NCITY,N,DE)
      DIMENSION X(NCITY),Y(NCITY),IORDER(NCITY),N(6),XX(6),YY(6)
      ALEN(X1,X2,Y1,Y2)=SQRT((X2-X1)**2+(Y2-Y1)**2)
      N(4)=1+MOD(N(3),NCITY)
      N(5)=1+MOD((N(1)+NCITY-2),NCITY)
      N(6)=1+MOD(N(2),NCITY)
      DO 11 J=1,6
          II=IORDER(N(J))
          XX(J)=X(II)
          YY(J)=Y(II)
11    CONTINUE
      DE=-ALEN(XX(2),XX(6),YY(2),YY(6))-ALEN(XX(1),XX(5),YY(1),YY(5))
     *     -ALEN(XX(3),XX(4),YY(3),YY(4))+ALEN(XX(1),XX(3),YY(1),YY(3))
     *     +ALEN(XX(2),XX(4),YY(2),YY(4))+ALEN(XX(5),XX(6),YY(5),YY(6))
      RETURN
      END    

      SUBROUTINE TRNSPT(IORDER,NCITY,N)
      PARAMETER(MAXCIT=1000)
      DIMENSION IORDER(NCITY),JORDER(MAXCIT),N(6)
      M1=1+MOD((N(2)-N(1)+NCITY),NCITY)
      M2=1+MOD((N(5)-N(4)+NCITY),NCITY)
      M3=1+MOD((N(3)-N(6)+NCITY),NCITY)
      NN=1  
      DO 11 J=1,M1
          JJ=1+MOD((J+N(1)-2),NCITY)
          JORDER(NN)=IORDER(JJ)
          NN=NN+1
11    CONTINUE
      IF (M2.GT.0) THEN
          DO 12 J=1,M2
              JJ=1+MOD((J+N(4)-2),NCITY)
              JORDER(NN)=IORDER(JJ)
              NN=NN+1
12        CONTINUE
      ENDIF
      IF (M3.GT.0) THEN
          DO 13 J=1,M3
              JJ=1+MOD((J+N(6)-2),NCITY)
              JORDER(NN)=IORDER(JJ)
              NN=NN+1
13        CONTINUE
      ENDIF
      DO 14 J=1,NCITY
          IORDER(J)=JORDER(J)
14    CONTINUE
      RETURN
      END

      SUBROUTINE METROP(DE,T,ANS)
      PARAMETER(JDUM=1)
      LOGICAL ANS
      ANS=(DE.LT.0.0).OR.(RAN3(JDUM).LT.EXP(-DE/T))
      RETURN
      END
