      SUBROUTINE SPARSE(B,N,ASUB,ATSUB,X,RSQ)
      PARAMETER (NMAX=500,EPS=1.E-6)
      DIMENSION B(N),X(N),G(NMAX),H(NMAX),XI(NMAX),XJ(NMAX)
      EPS2=N*EPS**2
      IRST=0
1     IRST=IRST+1
      CALL ASUB(X,XI)
      RP=0.
      BSQ=0.
      DO 11 J=1,N
        BSQ=BSQ+B(J)**2
        XI(J)=XI(J)-B(J)
        RP=RP+XI(J)**2
11    CONTINUE
      CALL ATSUB(XI,G)
      DO 12 J=1,N
        G(J)=-G(J)
        H(J)=G(J)
12    CONTINUE
      DO 19 ITER=1,10*N
        CALL ASUB(H,XI)
        ANUM=0.
        ADEN=0.
        DO 13 J=1,N
          ANUM=ANUM+G(J)*H(J)
          ADEN=ADEN+XI(J)**2
13      CONTINUE
        IF(ADEN.EQ.0.)PAUSE 'very singular matrix'
        ANUM=ANUM/ADEN
        DO 14 J=1,N
          XI(J)=X(J)
          X(J)=X(J)+ANUM*H(J)
14      CONTINUE
        CALL ASUB(X,XJ)
        RSQ=0.
        DO 15 J=1,N
          XJ(J)=XJ(J)-B(J)
          RSQ=RSQ+XJ(J)**2
15      CONTINUE
        IF(RSQ.EQ.RP.OR.RSQ.LE.BSQ*EPS2)RETURN
        IF(RSQ.GT.RP)THEN
          DO 16 J=1,N
            X(J)=XI(J)
16        CONTINUE
          IF(IRST.GE.3)RETURN
          GO TO 1
        ENDIF
        RP=RSQ
        CALL ATSUB(XJ,XI)
        GG=0.
        DGG=0.
        DO 17 J=1,N
          GG=GG+G(J)**2
          DGG=DGG+(XI(J)+G(J))*XI(J)
17      CONTINUE
        IF(GG.EQ.0.)RETURN
        GAM=DGG/GG
        DO 18 J=1,N
          G(J)=-XI(J)
          H(J)=G(J)+GAM*H(J)
18      CONTINUE
19    CONTINUE
      PAUSE 'too many iterations'
      RETURN
      END
