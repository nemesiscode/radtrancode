      SUBROUTINE MEMCOF(DATA,N,M,PM,COF,WK1,WK2,WKM)
      DIMENSION DATA(N),COF(M),WK1(N),WK2(N),WKM(M)
      P=0.
      DO 11 J=1,N
        P=P+DATA(J)**2
11    CONTINUE
      PM=P/N
      WK1(1)=DATA(1)
      WK2(N-1)=DATA(N)
      DO 12 J=2,N-1
        WK1(J)=DATA(J)
        WK2(J-1)=DATA(J)
12    CONTINUE
      DO 17 K=1,M
        PNEUM=0.
        DENOM=0.
        DO 13 J=1,N-K
          PNEUM=PNEUM+WK1(J)*WK2(J)
          DENOM=DENOM+WK1(J)**2+WK2(J)**2
13      CONTINUE
        COF(K)=2.*PNEUM/DENOM
        PM=PM*(1.-COF(K)**2)
        IF(K.NE.1)THEN
          DO 14 I=1,K-1
            COF(I)=WKM(I)-COF(K)*WKM(K-I)
14        CONTINUE
        ENDIF
        IF(K.EQ.M)RETURN
        DO 15 I=1,K
          WKM(I)=COF(I)
15      CONTINUE
        DO 16 J=1,N-K-1
          WK1(J)=WK1(J)-WKM(K)*WK2(J)
          WK2(J)=WK2(J+1)-WKM(K)*WK1(J+1)
16      CONTINUE
17    CONTINUE
      PAUSE 'never get here'
      END
