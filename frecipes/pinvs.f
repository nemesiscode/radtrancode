      SUBROUTINE PINVS(IE1,IE2,JE1,JSF,JC1,K,C,NCI,NCJ,NCK,S,NSI,NSJ)
      PARAMETER (ZERO=0.,ONE=1.,NMAX=10)
      DIMENSION C(NCI,NCJ,NCK),S(NSI,NSJ),PSCL(NMAX),INDXR(NMAX)
      JE2=JE1+IE2-IE1
      JS1=JE2+1
      DO 12 I=IE1,IE2
        BIG=ZERO
        DO 11 J=JE1,JE2
          IF(ABS(S(I,J)).GT.BIG) BIG=ABS(S(I,J))
11      CONTINUE
        IF(BIG.EQ.ZERO) PAUSE 'Singular matrix, row all 0'
        PSCL(I)=ONE/BIG
        INDXR(I)=0
12    CONTINUE
      DO 18 ID=IE1,IE2
        PIV=ZERO
        DO 14 I=IE1,IE2
          IF(INDXR(I).EQ.0) THEN
            BIG=ZERO
            DO 13 J=JE1,JE2
              IF(ABS(S(I,J)).GT.BIG) THEN
                JP=J
                BIG=ABS(S(I,J))
              ENDIF
13          CONTINUE
            IF(BIG*PSCL(I).GT.PIV) THEN
              IPIV=I
              JPIV=JP
              PIV=BIG*PSCL(I)
            ENDIF
          ENDIF
14      CONTINUE
        IF(S(IPIV,JPIV).EQ.ZERO) PAUSE 'Singular matrix'
        INDXR(IPIV)=JPIV
        PIVINV=ONE/S(IPIV,JPIV)
        DO 15 J=JE1,JSF
          S(IPIV,J)=S(IPIV,J)*PIVINV
15      CONTINUE
        S(IPIV,JPIV)=ONE
        DO 17 I=IE1,IE2
          IF(INDXR(I).NE.JPIV) THEN
            IF(S(I,JPIV).NE.ZERO) THEN
              DUM=S(I,JPIV)
              DO 16 J=JE1,JSF
                S(I,J)=S(I,J)-DUM*S(IPIV,J)
16            CONTINUE
              S(I,JPIV)=ZERO
            ENDIF
          ENDIF
17      CONTINUE
18    CONTINUE
      JCOFF=JC1-JS1
      ICOFF=IE1-JE1
      DO 21 I=IE1,IE2
        IROW=INDXR(I)+ICOFF
        DO 19 J=JS1,JSF
          C(IROW,J+JCOFF,K)=S(I,J)
19      CONTINUE
21    CONTINUE
      RETURN
      END
