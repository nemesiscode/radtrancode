      subroutine calc_esum2(fknu,delad,y1,T,lcalch,g_ord,del_g,k_g,
     1ng,miter)
C     ************************************************************************
C     Subroutine to calculate an exponential-sum approximation to a 
C     transmission curve defined by the Goody-voigt parameters.
C
C     Input variables
C	fknu		real 	knu at temperature of layer
C       delad		real 	mean line spacing/(ad0*sqrt(T)
C	y1		real 	aL/aD at pressure and temperature of layer
C 	lcalch		integer	Switch to determine how transmission curve
C				is calculated.
C				1  = Goody-Voigt
C				2  = Malkmus-Lorentz
C				3  = K-distribution
C	g_ord(maxg)	real	G-ordinates of required k-distribution
C	del_g(maxg)	real	Weights of required k-distribution
C	ng		integer	Number of ordinates in k-distribution
C	miter		integer	Max number of iterations to find k_g
C
C     Output variable:
C	k_g(maxg)		real	Calculated k-distribution. Also contains 
C				input transmission function parameters if 
C				lcalch=3
C
C     Pat Irwin		1/3/95
C
C     ***********************************************************************
C
C
C     Variables -------------------------------------------------------------
      implicit DOUBLE PRECISION (a-h,o-z)
      include '../includes/arrdef.f'
      integer miter,ng
      real yv,g_ord(maxg),gin,kout
      real fknu,delad,y1,T,del_g(maxg),k_g(maxg)
      logical info
      integer lcalch,err
      double precision SL,BL
      parameter (ma=10,mca=10,mdata=20)
      dimension x(mdata),y(mdata),sig(mdata),a(ma),lista(ma),
     *  covar(mca,mca),alpha(mca,mca),w(ma)
      real x1(mdata),yy1(mdata),sig1(mdata)
C     -----------------------------------------------------------------------

C     Check input parameters and return if zero.
      if(fknu.eq.0.or.delad.eq.0.or.y1.eq.0)then
       do j=1,ng
        k_g(j)=0.
       end do
       return
      end if

      if(lcalch.lt.1.or.lcalch.gt.3)then
       print*,'Calc_esum2. lcalch not valid. Setting to 1'
       lcalch=1
       return
      end if

      ndata=20

C     Calculate transmission curve to fit to

      call tran_curve(fknu,delad,y1,T,lcalch,k_g,del_g,ng,
     1ndata,x1,yy1,sig1)

      do loop=1,ndata
       x(loop)=dble(x1(loop))
       y(loop)=dble(yy1(loop))
       sig(loop)=dble(sig1(loop))
      end do

      na=ng
      nfit=na
      nca=na

C     Calculate first guess k-distribution using ML band
      yv=0.25*y1/delad
      call ml_lac1(fknu,yv,SL,BL)
      do i=1,ng
        gin=g_ord(i)
        call kml(SL,BL,gin,kout)
        a(i)=dexp(dble(-kout))
      end do

      do i=1,na
       w(i)=dble(del_g(i))
       lista(i)=i
      end do     

      alamda=-1.
      niter=0

20    info=.false.
      call mrqminp1(x,y,sig,ndata,a,na,lista,nfit,
     *    covar,alpha,nca,chisq,alamda,w,info,err)
      sigma=sqrt(chisq*sig(1)*sig(1)/ndata)
      niter=niter+1

C     Print out parameters if fitting has gone wrong
      if(err.lt.0)then
       print*,'Error in fitting'
       print*,' '
       print*,'fknu,delad,y1,T'
       print*,fknu,delad,y1,T
       print*,' '
       print*,'x  y  sig'
       do i=1,ndata
        print*,x(i),y(i),sig(i)
       end do
       print*,' '
       print*,'k-guess  attempted fit'
       do i=1,ng
        gin=g_ord(i)
        call kml(SL,BL,gin,kout)
        print*,kout,sngl(-dlog(a(i)))
       end do
      end if

      if(sigma.gt.1.e-3.and.niter.lt.miter)goto 20
      print*,'Number of iterations = ',niter
      print*,'Successful convergence = ',.not.(info)
      print*,'Esum : sd of fit to T(m) (%) = ',sigma*100.
      alamda=0.

      call mrqminp1(x,y,sig,ndata,a,na,lista,nfit,
     *    covar,alpha,nca,chisq,alamda,w,info,err)

C     Print out parameters if fitting has gone wrong
      if(err.lt.0)then
       print*,'Error in fitting'
       print*,' '
       print*,'fknu,delad,y1,T'
       print*,fknu,delad,y1,T
       print*,' '
       print*,'x  y  sig'
       do i=1,ndata
        print*,x(i),y(i),sig(i)
       end do
       print*,' '
       print*,'k-guess  attempted fit'
       do i=1,ng
        gin=g_ord(i)
        call kml(SL,BL,gin,kout)
        print*,kout,sngl(-dlog(a(i)))
       end do
      end if


      do i=1,na
       k_g(i)=sngl(-dlog(a(i)))
      end do

      return
      end


      subroutine ml_lac1(knu,yv,SL,BL)
C
***********************************************************************
C     
C     Convert EKS Malkmus-Lorentz Parameters to those used by Lacis and
C     Oinas (1991)
C
C     Pat Irwin 3/6/94
C
C
***********************************************************************
      implicit none
      real knu,yv,pi
      double precision A,B
      double precision SL,BL
      parameter (pi=3.1415927)

      A=dble(knu)
      B=dble(pi*A*yv)

      SL = A
      BL = 4.*B/(pi*A)

      return
      end

      SUBROUTINE MRQMINP1(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,
     *    COVAR,ALPHA,NCA,CHISQ,ALAMDA,W,INFO,ERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MA),
     *  COVAR(NCA,NCA),ALPHA(NCA,NCA),ATRY(MMAX),BETA(MMAX),DA(MMAX),
     *   W(MA)
      LOGICAL INFO
      INTEGER ERR
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
            PRINT*, 'Improper permutation in LISTA - IHIT>1'
            STOP
          ENDIF
12      CONTINUE
        IF (KK.NE.(MA+1)) THEN
         PRINT*, 'Improper permutation in LISTA - KK.NE.(MA+1)'
         STOP
        ENDIF
        ALAMDA=0.001
        CALL MRQCOFP1(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NCA,
     1 CHISQ,W)
        OCHISQ=CHISQ
        DO 13 J=1,MA
          ATRY(J)=A(J)
13      CONTINUE
      ENDIF
      DO 15 J=1,MFIT
        DO 14 K=1,MFIT
          COVAR(J,K)=ALPHA(J,K)
14      CONTINUE
        COVAR(J,J)=ALPHA(J,J)*(1.+ALAMDA)
        DA(J)=BETA(J)
15    CONTINUE
      CALL GAUSSJD(COVAR,MFIT,NCA,DA,1,1,ERR)
      IF(ERR.EQ.-1)THEN
       print*,err
       CHISQ=0.
       DO J=1,MA
        A(J)=1.
       END DO
       RETURN
      END IF
      IF(ALAMDA.EQ.0.)THEN
        CALL COVSRTD(COVAR,NCA,MA,LISTA,MFIT)
        RETURN
      ENDIF
      DO 16 J=1,MFIT
        ATRY(LISTA(J))=A(LISTA(J))+DA(J)
        IF(ATRY(LISTA(J)).LE.0.0)THEN
C         print*,'Calc_Esum2 Info: ATRY < 0. Set to A/2'
         INFO=.TRUE.
         ATRY(LISTA(J))=0.5*A(LISTA(J))
        END IF
        IF(ATRY(LISTA(J)).GT.1.0)THEN
         INFO=.TRUE.
C         print*,'Calc_Esum2 Info: ATRY > 1. Set to (1+A)/2'
         ATRY(LISTA(J))=0.5*(1+A(LISTA(J)))
        END IF
16    CONTINUE
23    TEMP1=0.
      DO 77 J=1,MFIT-1
       IF(ATRY(LISTA(J)).LT.ATRY(LISTA(J+1)))THEN
C        print*,'Calc_Esum2 Info: ATRY non-monotonic. Reset'
        INFO=.TRUE.
        TEMP1=ATRY(LISTA(J))
        ATRY(LISTA(J))=ATRY(LISTA(J+1))
        ATRY(LISTA(J+1))=TEMP1
       END IF
77    CONTINUE
      IF(TEMP1.GT.0.)GOTO 23

      CALL MRQCOFP1(X,Y,SIG,NDATA,ATRY,MA,LISTA,MFIT,COVAR,DA,NCA,
     1 CHISQ,W)
       IF(CHISQ.LT.OCHISQ)THEN
        ALAMDA=0.2*ALAMDA
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
         IF(ALAMDA.GT.1)ALAMDA=1.
         CHISQ=OCHISQ
       ENDIF
      RETURN
      END


      SUBROUTINE MRQCOFP1(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NALP,
     *CHISQ,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),ALPHA(NALP,NALP),BETA(MA),
     *    DYDA(MMAX),LISTA(MFIT),A(MA),W(MA)

      DO 12 J=1,MFIT
        DO 11 K=1,J
          ALPHA(J,K)=0.
11      CONTINUE
        BETA(J)=0.
12    CONTINUE
      CHISQ=0.
      DO 15 I=1,NDATA
        CALL EXPSUM1(X(I),A,YMOD,DYDA,W,MA)
        SIG2I=1./(SIG(I)*SIG(I))
        DY=Y(I)-YMOD
        DO 14 J=1,MFIT
          WT=DYDA(LISTA(J))*SIG2I
          DO 13 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+WT*DYDA(LISTA(K))
13        CONTINUE
          BETA(J)=BETA(J)+DY*WT
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

      SUBROUTINE EXPSUM1(X,A,Y,DYDA,W,NA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(NA),A(NA),DYDA(NA)
      Y=0.
      DO 10 I=1,NA
       Y=Y+W(I)*(A(I)**X)
       DYDA(I)=X*W(I)*(A(I)**(X-1))
10    CONTINUE
      RETURN
      END

      SUBROUTINE COVSRTD(COVAR,NCVM,MA,LISTA,MFIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
      END

      SUBROUTINE GAUSSJD(A,N,NP,B,M,MP,ERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ERR
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      ERR=0
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
                PRINT*,'GAUSSJD: Singular matrix'
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
         PRINT*,'GAUSSJD: A(ICOL,ICOL)=0'
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
      RETURN
      END

