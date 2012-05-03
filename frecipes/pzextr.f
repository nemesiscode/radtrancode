      SUBROUTINE PZEXTR(IEST,XEST,YEST,YZ,DY,NV,NUSE)
      PARAMETER (IMAX=11,NCOL=7,NMAX=10)
      DIMENSION X(IMAX),YEST(NV),YZ(NV),DY(NV),QCOL(NMAX,NCOL),D(NMAX)
      X(IEST)=XEST
      DO 11 J=1,NV
        DY(J)=YEST(J)
        YZ(J)=YEST(J)
11    CONTINUE
      IF(IEST.EQ.1) THEN
        DO 12 J=1,NV
          QCOL(J,1)=YEST(J)
12      CONTINUE
      ELSE
        M1=MIN(IEST,NUSE)
        DO 13 J=1,NV
          D(J)=YEST(J)
13      CONTINUE
        DO 15 K1=1,M1-1
          DELTA=1./(X(IEST-K1)-XEST)
          F1=XEST*DELTA
          F2=X(IEST-K1)*DELTA
          DO 14 J=1,NV
            Q=QCOL(J,K1)
            QCOL(J,K1)=DY(J)
            DELTA=D(J)-Q
            DY(J)=F1*DELTA
            D(J)=F2*DELTA
            YZ(J)=YZ(J)+DY(J)
14        CONTINUE
15      CONTINUE
        DO 16 J=1,NV
          QCOL(J,M1)=DY(J)
16      CONTINUE
      ENDIF
      RETURN
      END
