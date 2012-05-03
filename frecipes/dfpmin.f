      SUBROUTINE DFPMIN(P,N,FTOL,ITER,FRET)
      PARAMETER (NMAX=50,ITMAX=200,EPS=1.E-10)
      DIMENSION P(N),HESSIN(NMAX,NMAX),XI(NMAX),G(NMAX),DG(NMAX),
     *HDG(NMAX)
      FP=FUNC(P)
      CALL DFUNC(P,G)
      DO 12 I=1,N
        DO 11 J=1,N
          HESSIN(I,J)=0.
11      CONTINUE
        HESSIN(I,I)=1.
        XI(I)=-G(I)
12    CONTINUE
      DO 24 ITS=1,ITMAX
        ITER=ITS
        CALL LINMIN(P,XI,N,FRET)
        IF(2.*ABS(FRET-FP).LE.FTOL*(ABS(FRET)+ABS(FP)+EPS))RETURN
        FP=FRET
        DO 13 I=1,N
          DG(I)=G(I)
13      CONTINUE
        FRET=FUNC(P)
        CALL DFUNC(P,G)
        DO 14 I=1,N
          DG(I)=G(I)-DG(I)
14      CONTINUE
        DO 16 I=1,N
          HDG(I)=0.
          DO 15 J=1,N
            HDG(I)=HDG(I)+HESSIN(I,J)*DG(J)
15        CONTINUE
16      CONTINUE
        FAC=0.
        FAE=0.
        DO 17 I=1,N
          FAC=FAC+DG(I)*XI(I)
          FAE=FAE+DG(I)*HDG(I)
17      CONTINUE
        FAC=1./FAC
        FAD=1./FAE
        DO 18 I=1,N
          DG(I)=FAC*XI(I)-FAD*HDG(I)
18      CONTINUE
        DO 21 I=1,N
          DO 19 J=1,N
            HESSIN(I,J)=HESSIN(I,J)+FAC*XI(I)*XI(J)
     *        -FAD*HDG(I)*HDG(J)+FAE*DG(I)*DG(J)
19        CONTINUE
21      CONTINUE
        DO 23 I=1,N
          XI(I)=0.
          DO 22 J=1,N
            XI(I)=XI(I)-HESSIN(I,J)*G(J)
22        CONTINUE
23      CONTINUE
24    CONTINUE
      PAUSE 'too many iterations in DFPMIN'
      RETURN
      END
