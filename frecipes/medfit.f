      SUBROUTINE MEDFIT(X,Y,NDATA,A,B,ABDEV)
      PARAMETER (NMAX=1000)
      EXTERNAL ROFUNC
      COMMON /ARRAYS/ NDATAT,XT(NMAX),YT(NMAX),ARR(NMAX),AA,ABDEVT
      DIMENSION X(NDATA),Y(NDATA)
      SX=0.
      SY=0.
      SXY=0.
      SXX=0.
      DO 11 J=1,NDATA
        XT(J)=X(J)
        YT(J)=Y(J)
        SX=SX+X(J)
        SY=SY+Y(J)
        SXY=SXY+X(J)*Y(J)
        SXX=SXX+X(J)**2
11    CONTINUE
      NDATAT=NDATA
      DEL=NDATA*SXX-SX**2
      AA=(SXX*SY-SX*SXY)/DEL
      BB=(NDATA*SXY-SX*SY)/DEL
      CHISQ=0.
      DO 12 J=1,NDATA
        CHISQ=CHISQ+(Y(J)-(AA+BB*X(J)))**2
12    CONTINUE
      SIGB=SQRT(CHISQ/DEL)
      B1=BB
      F1=ROFUNC(B1)
      B2=BB+SIGN(3.*SIGB,F1)
      F2=ROFUNC(B2)
1     IF(F1*F2.GT.0.)THEN
        BB=2.*B2-B1
        B1=B2
        F1=F2
        B2=BB
        F2=ROFUNC(B2)
        GOTO 1
      ENDIF
      SIGB=0.01*SIGB
2     IF(ABS(B2-B1).GT.SIGB)THEN
        BB=0.5*(B1+B2)
        IF(BB.EQ.B1.OR.BB.EQ.B2)GOTO 3
        F=ROFUNC(BB)
        IF(F*F1.GE.0.)THEN
          F1=F
          B1=BB
        ELSE
          F2=F
          B2=BB
        ENDIF
        GOTO 2
      ENDIF
3     A=AA
      B=BB
      ABDEV=ABDEVT/NDATA
      RETURN
      END
