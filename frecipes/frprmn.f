      SUBROUTINE FRPRMN(P,N,FTOL,ITER,FRET)
      PARAMETER (NMAX=50,ITMAX=200,EPS=1.E-10)
      DIMENSION P(N),G(NMAX),H(NMAX),XI(NMAX)
      FP=FUNC(P)
      CALL DFUNC(P,XI)
      DO 11 J=1,N
        G(J)=-XI(J)
        H(J)=G(J)
        XI(J)=H(J)
11    CONTINUE
      DO 14 ITS=1,ITMAX
        ITER=ITS
        CALL LINMIN(P,XI,N,FRET)
        IF(2.*ABS(FRET-FP).LE.FTOL*(ABS(FRET)+ABS(FP)+EPS))RETURN
        FP=FUNC(P)
        CALL DFUNC(P,XI)
        GG=0.
        DGG=0.
        DO 12 J=1,N
          GG=GG+G(J)**2
C         DGG=DGG+XI(J)**2
          DGG=DGG+(XI(J)+G(J))*XI(J)
12      CONTINUE
        IF(GG.EQ.0.)RETURN
        GAM=DGG/GG
        DO 13 J=1,N
          G(J)=-XI(J)
          H(J)=G(J)+GAM*H(J)
          XI(J)=H(J)
13      CONTINUE
14    CONTINUE
      PAUSE 'FRPR maximum iterations exceeded'
      RETURN
      END
