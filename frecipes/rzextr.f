      SUBROUTINE RZEXTR(IEST,XEST,YEST,YZ,DY,NV,NUSE)
      PARAMETER (IMAX=11,NMAX=10,NCOL=7)
      DIMENSION X(IMAX),YEST(NV),YZ(NV),DY(NV),D(NMAX,NCOL),FX(NCOL)
      X(IEST)=XEST
      IF(IEST.EQ.1) THEN
        DO 11 J=1,NV
          YZ(J)=YEST(J)
          D(J,1)=YEST(J)
          DY(J)=YEST(J)
11      CONTINUE
      ELSE
        M1=MIN(IEST,NUSE)
        DO 12 K=1,M1-1
          FX(K+1)=X(IEST-K)/XEST
12      CONTINUE
        DO 14 J=1,NV
          YY=YEST(J)
          V=D(J,1)
          C=YY
          D(J,1)=YY
          DO 13 K=2,M1
            B1=FX(K)*V
            B=B1-C
            IF(B.NE.0.) THEN
              B=(C-V)/B
              DDY=C*B
              C=B1*B
            ELSE
              DDY=V
            ENDIF
            V=D(J,K)
            D(J,K)=DDY
            YY=YY+DDY
13        CONTINUE
          DY(J)=DDY
          YZ(J)=YY
14      CONTINUE
      ENDIF
      RETURN
      END
