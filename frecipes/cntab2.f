      SUBROUTINE CNTAB2(NN,NI,NJ,H,HX,HY,HYGX,HXGY,UYGX,UXGY,UXY)
      PARAMETER (MAXI=100,MAXJ=100,TINY=1.E-30)
      DIMENSION NN(NI,NJ),SUMI(MAXI),SUMJ(MAXJ)
      SUM=0
      DO 12 I=1,NI
        SUMI(I)=0.0
        DO 11 J=1,NJ
          SUMI(I)=SUMI(I)+NN(I,J)
          SUM=SUM+NN(I,J)
11      CONTINUE
12    CONTINUE
      DO 14 J=1,NJ
        SUMJ(J)=0.
        DO 13 I=1,NI
          SUMJ(J)=SUMJ(J)+NN(I,J)
13      CONTINUE
14    CONTINUE
      HX=0.
      DO 15 I=1,NI
        IF(SUMI(I).NE.0.)THEN
          P=SUMI(I)/SUM
          HX=HX-P*ALOG(P)
        ENDIF
15    CONTINUE
      HY=0.
      DO 16 J=1,NJ
        IF(SUMJ(J).NE.0.)THEN
          P=SUMJ(J)/SUM
          HY=HY-P*ALOG(P)
        ENDIF
16    CONTINUE
      H=0.
      DO 18 I=1,NI
        DO 17 J=1,NJ
          IF(NN(I,J).NE.0)THEN
            P=NN(I,J)/SUM
            H=H-P*ALOG(P)
          ENDIF
17      CONTINUE
18    CONTINUE
      HYGX=H-HX
      HXGY=H-HY
      UYGX=(HY-HYGX)/(HY+TINY)
      UXGY=(HX-HXGY)/(HX+TINY)
      UXY=2.*(HX+HY-H)/(HX+HY+TINY)
      RETURN
      END
