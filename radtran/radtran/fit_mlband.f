      subroutine fit_mlband(X,Y,SIG,NDATA,SL,BL)
C     ***********************************************************************
C     Subroutine fit a Goody Voigt Band to an input transmission curve
C
C     Input variables
C	X(20)	REAL	X array of curve
C	Y(20)	REAL	Y array of curve
C	SIG(20)	REAL	s.d. array of curve
C 	NDATA 	INTEGER	Number of points in curve
C
C     Input/Output variables
C	knu	REAL	GV Band Parameters. At initiation values are first
C	delta	"	guess parameters. Contain fitted parameters on
C	y1	"	output
C
C     ***********************************************************************
C     Parameters for the least squares fit
      PARAMETER (MA=10,MFIT=10,MCA=10,MDATA=20,MMAX=20)
      DIMENSION X(MDATA),Y(MDATA),SIG(MDATA),A(MA),LISTA(MA),
     *  COVAR(MCA,MCA),ALPHA(MCA,MCA),ATRY(MMAX),BETA(MMAX),DA(MMAX)
      LOGICAL INFO
      DOUBLE PRECISION SL,BL
      INTEGER ERR

      NA=2
      NFIT=NA
      NCA=NA

      DO I=1,NA
       LISTA(I)=I
      END DO

      A(1)=SNGL(SL)
      A(2)=SNGL(BL)

C      print*,'fit_mlband. SL,BL (initial) = ',SL,BL
C      print*,'LISTA : ',(LISTA(I),I=1,NA)
C      print*,'A : ',(A(I),I=1,NA)
      ALAMDA=-1
      NITER=0

20    INFO=.FALSE.
      CALL MRQMINM(X,Y,SIG,NDATA,A,NA,LISTA,NFIT,
     *    COVAR,ALPHA,NCA,OCHISQ,CHISQ,ALAMDA,INFO,ERR)

      sigma=sqrt(chisq*sig(1)*sig(1)/na)
C      print*,'niter,alamda,chisq,sigma',niter,alamda,chisq,sigma
      NITER=NITER+1
      IF((sigma.gt.1.e-3.and.NITER.lt.100).OR.INFO)GOTO 20
      
      ALAMDA=0.

      CALL MRQMINM(X,Y,SIG,NDATA,A,NA,LISTA,NFIT,
     *    COVAR,ALPHA,NCA,OCHISQ,CHISQ,ALAMDA,INFO,ERR)

      IF(ERR.NE.-1)THEN
       SL=DBLE(A(1))
       BL=DBLE(A(2))
      ENDIF

      RETURN

      END

      SUBROUTINE MRQMINM(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,
     *    COVAR,ALPHA,NCA,OCHISQ,CHISQ,ALAMDA,INFO,ERR)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MA),
     *  COVAR(NCA,NCA),ALPHA(NCA,NCA),ATRY(MMAX),BETA(MMAX),DA(MMAX)
      LOGICAL INFO
      INTEGER ERR

C      print*,'MRQMINM A: ',(A(I),I=1,MA)
C      print*,'MRQMINM LISTA: ',(LISTA(I),I=1,MA)
C      print*,'ALAMBDA : ',alamda
      IF(ALAMDA.LT.0.)THEN
        KK=MFIT+1
        DO 12 J=1,MA
          IHIT=0
          DO 11 K=1,MFIT
            IF(LISTA(K).EQ.J)IHIT=IHIT+1
11        CONTINUE
          IF (IHIT.EQ.0) THEN
            LISTA(KK)=J
            KK=KK+1
          ELSE IF (IHIT.GT.1) THEN
            PRINT*,'Improper permutation in LISTA - IHIT.GT.1'
	    STOP
          ENDIF
12      CONTINUE
        IF (KK.NE.(MA+1)) THEN
         PRINT*, 'Improper permutation in LISTA - KK.NE.(MA+1)'
         STOP
        ENDIF
        ALAMDA=2000.0
        CALL MRQCOFM(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NCA,CHISQ)
        OCHISQ=CHISQ
C        print*,'OCHISQ = ',OCHISQ
        DO 13 J=1,MA
          ATRY(J)=A(J)
13      CONTINUE
      ENDIF
C      print*,'OK to here'
C      print*,'MFIT = ',MFIT
      DO 15 J=1,MFIT
        DO 14 K=1,MFIT
          COVAR(J,K)=ALPHA(J,K)
14      CONTINUE
        COVAR(J,J)=ALPHA(J,J)*(1.+ALAMDA)
C        print*,J,BETA(J)
        DA(J)=BETA(J)
15    CONTINUE
C      print*,'Calling GAUSSJ1'
      CALL GAUSSJ1(COVAR,MFIT,NCA,DA,1,1,ERR)
      IF(ALAMDA.EQ.0.)THEN
        CALL COVSRT(COVAR,NCA,MA,LISTA,MFIT)
        RETURN
      ENDIF
      DO 16 J=1,MFIT
        ATRY(LISTA(J))=A(LISTA(J))+DA(J)
        IF(ATRY(LISTA(J)).LE.0.0)THEN
C         PRINT*,'fit_band info: A(I) < 0. Set to 0.5*A(I)'
C         PRINT*,'J,LISTA(J),ATRY(LISTA(J))',J,LISTA(J),
C     &    A(LISTA(J)),ATRY(LISTA(J))
         ATRY(LISTA(J))=0.5*A(LISTA(J))
C         PRINT*,'ALAMDA, new ATRY = ',ALAMDA,ATRY(LISTA(J))
         INFO=.TRUE.
        ENDIF
16    CONTINUE
      CALL MRQCOFM(X,Y,SIG,NDATA,ATRY,MA,LISTA,MFIT,COVAR,DA,NCA,CHISQ)
C      print*,'OCHISQ,CHISQ,ALAMBDA',OCHISQ,CHISQ,ALAMDA
      IF(CHISQ.LE.OCHISQ)THEN
        ALAMDA=0.5*ALAMDA
        OCHISQ=CHISQ
        DO 18 J=1,MFIT
          DO 17 K=1,MFIT
            ALPHA(J,K)=COVAR(J,K)
17        CONTINUE
          BETA(J)=DA(J)
          A(LISTA(J))=ATRY(LISTA(J))
18      CONTINUE
      ELSE
        ALAMDA=5.*ALAMDA
        CHISQ=OCHISQ
      ENDIF
      RETURN
      END


      SUBROUTINE MRQCOFM(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NALP,
     *CHISQ)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),ALPHA(NALP,NALP),BETA(MA),
     *    DYDA(MMAX),LISTA(MFIT),A(MA)
      REAL U
      DO 12 J=1,MFIT
        DO 11 K=1,J
          ALPHA(J,K)=0.
11      CONTINUE
        BETA(J)=0.
12    CONTINUE
      CHISQ=0.
      DO 15 I=1,NDATA
        U=X(I)
        CALL TAU_ML(U,A,MA,YMOD,DYDA)
        SIG2I=1./(SIG(I)*SIG(I))
C        print*,'MRQCOF',I,Y(I),YMOD,SIG2I,DYDA(1),DYDA(2)
        DY=Y(I)-YMOD
        DO 14 J=1,MFIT
          WT=DYDA(LISTA(J))*SIG2I
C          print*,J,LISTA(J),WT
          DO 13 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+WT*DYDA(LISTA(K))
13        CONTINUE
          BETA(J)=BETA(J)+DY*WT
C          print*,J,BETA(J)
14      CONTINUE
        CHISQ=CHISQ+DY*DY*SIG2I
15    CONTINUE
      DO 17 J=2,MFIT
        DO 16 K=1,J-1
          ALPHA(K,J)=ALPHA(J,K)
16      CONTINUE
17    CONTINUE
      RETURN
      END


      SUBROUTINE TAU_ML(U,A,NA,YMOD,DYDA)
      INTEGER NA,MMAX
      PARAMETER(MMAX=20)
      REAL A(NA),DYDA(MMAX),A1(2),U,YMOD,Y1,Y2,DA

      DO 5 I=1,2
       A1(I)=A(I)
5     CONTINUE
      YMOD=TRANSMITM(A1,U)
C      print*,'YMOD = ',YMPD
C      PRINT*,'A1 : ',A1
      DO 10 I=1,NA
       DA=A(I)*0.01
C       print*,I,A(I),DA
       A1(I)=A(I)+DA
C       print*,U,A1
       Y2=TRANSMITM(A1,U)
C       print*,Y2
       A1(I)=A(I)-DA
C       print*,U,A1
       Y1=TRANSMITM(A1,U)
C       print*,Y1
       DY=Y2-Y1
C       IF(ABS(DY).LT.1e-6)THEN
C        DYDA(I)=0.0
C       ELSE
C        DYDA(I)=(Y2-Y1)/(2.*DA)
C       ENDIF
        DYDA(I)=(Y2-Y1)/(2.*DA)
C       print*,Y2,Y1,DA,DYDA(I)
       A1(I)=A(I)
10    CONTINUE

C      print*,'tau_ml: ymod, dyda',ymod,(dyda(i),i=1,2)

      RETURN

      END


      REAL FUNCTION TRANSMITM(A1,U)
      REAL A1(2),U,PI,dexp
      DOUBLE PRECISION SL,BL
      PARAMETER(PI=3.1415927)

      SL=DBLE(A1(1))
      BL=DBLE(A1(2))


      TRANSMITM=SNGL(dexp(-0.5*PI*BL*
     & (DSQRT(1.+4.*SL*U/(PI*BL))-1.)))

      return
      end


      SUBROUTINE GAUSSJ1(A,N,NP,B,M,MP,ERR)
      IMPLICIT REAL (A-H,O-Z)
      INTEGER ERR
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      ERR=0
C      print*,'gaussj1 called'
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PRINT*,'Gaussj1: Singular matrix'
                PRINT*,'IPIV(K) > 1'
                ERR=-1
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) THEN
         ERR=-1
         PRINT*,'GAUSSJ1: A(ICOL,ICOL)=0'
         PRINT*,'Singular matrix.'
         PIVINV=1.
        ELSE
         PIVINV=1./A(ICOL,ICOL)
        END IF
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
C      print*,'gaussj1 OK'
      RETURN
      END


