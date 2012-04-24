      subroutine calc_esum5(ndata,x1,yy1,sig1,ng,del_g,k_g,errk,alamda1)
C     ************************************************************************
C     Fits and exponential sum series to an input transmission curve using
C     the Levenburg-Marquedt method non-lnear least squares method 
C     described in Numerical recipes.
C
C     Routine solves monotonicity requirement by defining:
C	k(i) = k(i-1) + a(i)	where a(i) >= 0.
C
C     Input Variables:
C	ndata	     	integer	Number of points in transmission curve
C	x1(mdata)	real	Amounts
C	yy1(mdata)	real	Transmissions
C	sig1(mdata)	real	Transmission errors
C	ng		integer	Number of k-coefficients
C	del_g(maxg)	real	Weights of coefficients
C	k_g(maxg)		real	First guess k-coefficients
C	alamda1		real	First value of lamda for MRQ iterations
C
C     Output variables
C	k_g(maxg)		real	Fitted coefficients
C	errk(maxg)	real	Estimated fitting errors
C     ***********************************************************************
      implicit DOUBLE PRECISION (a-h,o-z)
      include '../includes/arrdef.f'
      integer miter,ng,ndata,nerr
      real pi,yv,g_ord(maxg),gin,kout,errk(maxg),k1,kold,alamda1
      parameter (pi=3.1415927)
      real fknu,delad,y1,T,del_g(maxg),k_g(maxg),tau_goody_voigt2
      logical info
      integer lcalch,err
      double precision SL,BL
      parameter (ma=10,mfit=10,mca=10,mdata=20,mmax=20,miter=100)
      dimension x(mdata),y(mdata),sig(mdata),a(ma),lista(ma),
     *  covar(mca,mca),alpha(mca,mca),atry(mmax),beta(mmax),da(mmax),
     *  w(ma),dyda(mdata)
      real x1(mdata),yy1(mdata),sig1(mdata),tmp
      character*1,ans
C     -----------------------------------------------------------------------

      print*,'Now in esum5'
      alamdaset=dble(alamda1)
 
      do loop=1,ndata
       x(loop)=dble(x1(loop))
       y(loop)=dble(yy1(loop))
       sig(loop)=dble(sig1(loop))
      end do

      print*,'flag 1'
      na=ng
      nfit=na
      nca=na


C     Calculate first guess
      kold=0.
      do i=1,ng
        kout=k_g(i)
        a(i)=kout-kold
        kold=kout
      end do
    
      do i=1,na
       w(i)=dble(del_g(i))
       lista(i)=i
      end do     

      alamda=-1.
      niter=0
      nerr = 0
      ochisq=1e10
      print*,'flag 2'

20    info=.false.
      call mrqminp5(x,y,sig,ndata,a,na,lista,nfit,
     *    covar,alpha,nca,ochisq,chisq,alamda,w,info,err,alamdaset)
C      print*,'mrqminp5 OK'

      sigma=sqrt(chisq*sig(1)*sig(1)/ndata)
      print*,niter,alamda,ochisq,chisq,sigma

      niter=niter+1

      if((sigma.gt.1e-3.and.niter.lt.miter))goto 20

      print*,'Successful convergence. Alamda : ',alamda

c      print*,'Calculating COVAR. niter = ',niter
      alamda=0.

      call mrqminp5(x,y,sig,ndata,a,na,lista,nfit,
     *    covar,alpha,nca,ochisq,chisq,alamda,w,info,err,alamdaset)

      sum=0.
      sumerr=0.
      do i=1,na
	if (covar(I,I).ge.1.e30) covar(I,I) = 1.e30
	tmp = sum + a(I)
       	sum=min(1.e30,tmp)
       k_g(i)=sum
       if(covar(i,i).lt.0.and.nerr.lt.200)then
        print*,'Esum5. Convergence error, covar < 0. nerr = ',nerr
        print*,'i,covar(i,i) ',i,covar(i,i)
        alamda=alamdaset
        nerr = nerr+1
        goto 20
       else
	tmp = abs(covar(I,I))
	tmp = sumerr+sngl(abs(covar(i,i)))
        sumerr=min(1.e30,tmp)
        e1=sqrt(sumerr)
        errk(i) = 100.*e1/k_g(i)
       end if
      end do

      return
      end

      SUBROUTINE MRQMINP5(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,
     *    COVAR,ALPHA,NCA,OCHISQ,CHISQ,ALAMDA,W,INFO,ERR,ALAMDASET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MA),
     *  COVAR(NCA,NCA),ALPHA(NCA,NCA),ATRY(MMAX),BETA(MMAX),DA(MMAX),
     *   W(MA)
      LOGICAL INFO
      INTEGER ERR
C      print*,'In mrqminp5. Alamda = ',ALAMDA
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
            PAUSE 'Improper permutation in LISTA'
          ENDIF
12      CONTINUE
        IF (KK.NE.(MA+1)) PAUSE 'Improper permutation in LISTA'
        ALAMDA=ALAMDASET
        CALL MRQCOFP5(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NCA,
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
C      print*,'calling gaussjd'
      CALL GAUSSJD(COVAR,MFIT,NCA,DA,1,1,ERR)
C      print*,'gaussjd OK'
      IF(ERR.EQ.-1)THEN
c       print*,'Calc_esum5. Error = ',err
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
c          print*,'Calc_Esum5 Info: ATRY < 0. Set to A'
          INFO=.TRUE.
C          ALAMDA=ALAMDASET
          ATRY(LISTA(J)) = A(LISTA(J))
        END IF
16    CONTINUE

      CALL MRQCOFP5(X,Y,SIG,NDATA,ATRY,MA,LISTA,MFIT,COVAR,DA,NCA,
     1 CHISQ,W)
       IF(CHISQ.LT.OCHISQ)THEN
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
         ALAMDA=2.0*ALAMDA
c         IF(ALAMDA.GT.1.0D290)ALAMDA = 1.0D290
         IF(ALAMDA.GT.1.0D150)ALAMDA = 1.0D150
         CHISQ=OCHISQ
       ENDIF
      RETURN
      END


      SUBROUTINE MRQCOFP5(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NALP,
     *CHISQ,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),ALPHA(NALP,NALP),BETA(MA),
     *    DYDA(MMAX),LISTA(MFIT),A(MA),W(MA)
      
C      print*,'Now in MRQCOFP5'
      DO 12 J=1,MFIT
        DO 11 K=1,J
          ALPHA(J,K)=0.
11      CONTINUE
        BETA(J)=0.
12    CONTINUE
      CHISQ=0.
      DO 15 I=1,NDATA
        CALL EXPSUM3(X(I),A,YMOD,DYDA,W,MA)
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

      SUBROUTINE EXPSUM3(X,A,Y,DYDA,W,NA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../includes/arrdef.f'
      DIMENSION W(NA),A(NA),DYDA(NA)
      DIMENSION TEMP(MAXG)

      ARG=0.
      DO 10 I=1,NA
       ARG=ARG+A(I)
       TEMP(I)=DEXP(-ARG*X)
10    CONTINUE
      
      Y=0.
      DO 20 I=1,NA
       Y=Y+W(I)*TEMP(I)
       DYDA(I)=0.
c       DO 30 J=1,I				! Think this is right
       DO 30 J=I,NA				! But this works better???
        DYDA(J)=DYDA(J) - X*W(J)*TEMP(J)
30     CONTINUE
20    CONTINUE

      RETURN
      END


