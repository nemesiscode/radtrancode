      SUBROUTINE SCRSHO(FX)
      PARAMETER (ISCR=60,JSCR=21)
      CHARACTER*1 SCR(ISCR,JSCR),BLANK,ZERO,YY,XX,FF
      DIMENSION Y(ISCR)
      DATA BLANK,ZERO,YY,XX,FF/' ','-','l','-','x'/
1     CONTINUE
      WRITE (*,*) ' Enter X1,X2 (= to stop)'
      READ (*,*) X1,X2
      IF(X1.EQ.X2) RETURN
      DO 11 J=1,JSCR
        SCR(1,J)=YY
        SCR(ISCR,J)=YY
11    CONTINUE
      DO 13 I=2,ISCR-1
        SCR(I,1)=XX
        SCR(I,JSCR)=XX
        DO 12 J=2,JSCR-1
          SCR(I,J)=BLANK
12      CONTINUE
13    CONTINUE
      DX=(X2-X1)/(ISCR-1)
      X=X1
      YBIG=0.
      YSML=YBIG
      DO 14 I=1,ISCR
        Y(I)=FX(X)
        IF(Y(I).LT.YSML) YSML=Y(I)
        IF(Y(I).GT.YBIG) YBIG=Y(I)
        X=X+DX
14    CONTINUE
      IF(YBIG.EQ.YSML) YBIG=YSML+1.
      DYJ=(JSCR-1)/(YBIG-YSML)
      JZ=1-YSML*DYJ
      DO 15 I=1,ISCR
        SCR(I,JZ)=ZERO
        J=1+(Y(I)-YSML)*DYJ
        SCR(I,J)=FF
15    CONTINUE
      WRITE (*,'(1X,1PE10.3,1X,80A1)') YBIG,(SCR(I,JSCR),I=1,ISCR)
      DO 16 J=JSCR-1,2,-1
        WRITE (*,'(12X,80A1)') (SCR(I,J),I=1,ISCR)
16    CONTINUE
      WRITE (*,'(1X,1PE10.3,1X,80A1)') YSML,(SCR(I,1),I=1,ISCR)
      WRITE (*,'(12X,1PE10.3,40X,E10.3)') X1,X2
      GOTO 1
      END
