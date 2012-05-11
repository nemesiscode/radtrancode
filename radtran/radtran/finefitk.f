      subroutine finefitk(imethod,ndata,x1,yy1,sig1,ng,
     &  del_g,k_g,alamda1,icalc)
C     ************************************************************************
C     Fits and exponential sum series to an input transmission curve using
C     the Levenburg-Marquedt method non-lnear least squares method 
C     described in Numerical recipes.
C
C     Routine solves for log(k) to ensure positive numbers.
C     Routine shuffles coefficients if they become non-monotonic
C
C     Input Variables:
C	imethod		integer	Desired method of fitting k-coefficients
C				0 : a(i)=alog(k_g(i))
C				1 : a(i)=k_g(i)-k_g(i-1)
C	ndata	     	integer	Number of points in transmission curve
C	x1(mdata)	real	Amounts
C	yy1(mdata)	real	Transmissions
C	sig1(mdata)	real	Transmission errors
C	ng		integer	Number of k-coefficients
C	del_g(10)	real	Weights of coefficients
C	k_g(10)		real	First guess k-coefficients
C	alamda1		real	First value of lamda for MRQ iterations
C
C     Output variables
C	k_g(10)		real	Fitted coefficients
C	icalc		integer	Quality indicator:
C				0 : successful fit
C				-1 : Fitting errors encountered during fit
C				-2 : overflow in alambda
C				-3 : Fitting errors encountered at end
C
C     Pat Irwin
C       Original July 2003
C	Revised and documented	9/1/04
C	Revised again		24/1/04
C
C     ***********************************************************************
      implicit double precision (a-h,o-z)
      include '../includes/arrdef.f'
      integer miter,ng,ndata,imethod
      real alamda1,kout,kold
      real del_g(maxg),k_g(maxg),chistore(4),t(4)
      logical aover
      integer err,icalc,istore
      parameter (ma=10,mca=10,mdata=20,miter=100)
      dimension x(mdata),y(mdata),sig(mdata),a(ma),
     *  covar(mca,mca),alpha(mca,mca),
     *  w(ma),yfit(mdata)
      real x1(mdata),yy1(mdata),sig1(mdata),tmp
C     -----------------------------------------------------------------------

      alamdaset=dble(alamda1)
      icalc=0
      do loop=1,ndata
       x(loop)=dble(x1(loop))
       y(loop)=dble(yy1(loop))
       sig(loop)=dble(sig1(loop))
      end do

      na=ng
      nca=na

C     Calculate first guess
      kold=0.  
      do i=1,ng
        if(imethod.eq.0)then
         a(i)=alog(k_g(i))
        else
         kout=k_g(i)   
         a(i)=kout-kold
         kold=kout
        endif
C        print*,a(i)        
      end do
    
      do i=1,na
       w(i)=dble(del_g(i))
      end do     

      alamda=-1.
      niter=0
      ochisq=1e10
      istore=1
      ifull = 0

20    continue
C      print*,'niter = ',niter
      call mrqminp7(x,y,yfit,sig,ndata,a,na,
     *    covar,alpha,nca,ochisq,chisq,alamda,w,err,aover,
     &    alamdaset,imethod)
    

      chistore(istore)=sngl(chisq)
      istore=istore+1
      if(istore.eq.5)then
       istore = 1
       ifull = 1
      endif

      if(err.eq.-1)then
       print*,'Finefitk: Singular matrix'
       icalc = -1
      endif 

      if(aover)then
       print*,'ALAMBDA overflowing. Fit aborted'
       icalc = -2
      endif

      sigma = 0.0
      do i=1,ndata
       dy = abs(y(i)-yfit(i))
       if(dy.gt.sigma)sigma=dy
C       print*,y(i),yfit(i)   
      enddo

      niter=niter+1

      if(ifull.eq.1)then
       do i=1,4
        j=istore-5+i
        if(j.lt.1)then
         j=j+4
        endif
        t(i)=chistore(j)
C        print*,i,t(i)
       enddo
       if((t(2).le.t(1)).and.(t(3).lt.t(2)).and.(t(4).lt.t(3)))then
        dx = 100.0*(t(3)-t(4))/t(3)
        if(dx.lt.0.0001)then
         iconverge=1
        else
         iconverge=0
        endif
       endif
      endif

C      print*,iconverge
      if(iconverge.eq.1)goto 21

C      print*,sigma,icalc
      if(sigma.gt.0.001.and.niter.lt.miter.and.icalc.ne.-2)goto 20

21    alamda=0.
C      print*,'niter,dx,iconverge = ',niter,dx,iconverge

      call mrqminp7(x,y,yfit,sig,ndata,a,na,
     &    covar,alpha,nca,ochisq,chisq,alamda,w,err,aover,
     &    alamdaset,imethod)

      if(err.eq.-1)then
       print*,'Finefitk: Singular matrix'
       print*,'Aborting'
       icalc = -3
      endif

      if(imethod.eq.0)then
       do i=1,na
        k_g(i)=sngl(exp(a(I)))
       end do
      else
       sum=0.0
       do i=1,na
        tmp = sngl(sum + a(I))
        sum=min(1.e30,tmp)
        k_g(i)=sngl(sum)
       enddo
      endif

      return

      end

      SUBROUTINE MRQMINP7(X,Y,YFIT,SIG,NDATA,A,MA,
     *    COVAR,ALPHA,NCA,OCHISQ,CHISQ,ALAMDA,W,ERR,AOVER,
     *    ALAMDASET,IMETHOD)
C     ***************************************************************
C     Version of Marqardt-Levenburg fitting routine tailored for fitting 
C     exponential sums to transmission curve data
C
C     Input variables
C	X(NDATA)	double	x-coeffients of data to be fitted (amounts)
C	Y(NDATA)	double	transmissions
C	SIG(NDATA)	double	transmission errors
C	NDATA		integer Number of points in transmission curves
C	A(MA)		double	Variables to vary to fit transmission curve
C	MA		integer Array length of A()
C	NCA		integer	defines size of arrays COVAR and ALPHA
C	ALAMDA		double	Alambda variable M-L parameter
C	W(MA)		double	K-G weights
C	ALAMDASET	double	Initial set ALAMDA
C	IMETHOD		integer	Required representation of K_G() by A()
C
C     Output variables
C	YFIT(NDATA)	double	Fitted transmissions
C	COVAR(NCA,NCA)	double	Can't remember
C	ALPHA(NCA,NCA)	double	Can't remember
C	CHISQ		double	Chi^2 between Y() and calculated YFIT()
C	ERR		integer	Error flag from matrix inverter GAUSSJD
C	AOVER		logical	Flag indicates if ALAMBDA has overflown
C	
C     Pat Irwin		9/1/04
C	
C     ***************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),YFIT(NDATA),SIG(NDATA),A(MA),
     &  COVAR(NCA,NCA),ALPHA(NCA,NCA),ATRY(MMAX),LISTA(MMAX),
     &  BETA(MMAX),DA(MMAX,1), W(MA)
      INTEGER ERR,IMETHOD,JMOD,ICHECK
      LOGICAL AOVER

      AOVER=.FALSE.
      IF(ALAMDA.LT.0.)THEN
        ALAMDA=ALAMDASET
C        PRINT*,'SET: ALAMBDA,A',ALAMDA,(A(I),I=1,MA)
        CALL MRQCOFP7(X,Y,YFIT,SIG,NDATA,A,MA,ALPHA,
     1   BETA,NCA,CHISQ,W,IMETHOD)
C         PRINT*,'ALAMBDA,CHISQ',ALAMBDA,CHISQ
        OCHISQ=CHISQ
        DO 13 J=1,MA
          ATRY(J)=A(J)
13      CONTINUE
      ENDIF

      DO 15 J=1,MA
        DO 14 K=1,MA
          COVAR(J,K)=ALPHA(J,K)
14      CONTINUE
        COVAR(J,J)=ALPHA(J,J)*(1.+ALAMDA)
        LISTA(J)=J
        DA(J,1)=BETA(J)
C        print*,'beta',j,beta(j)
15    CONTINUE

C      print*,'COVAR'
C      open(34,file='cov.dat',status='unknown')
C       do i=1,ma
C        print*,(covar(i,j),j=1,ma)
C        write(34,*)(covar(i,j),j=1,ma)
C       enddo
C      close(34)

      JMOD=0
      ICHECK=0
      TEST = 1.0
C      print*,'Calling DINVMARQ'
      CALL DINVMARQ(JMOD,ICHECK,COVAR,MA,NCA,DA,ERR)
C      print*,'ERR = ',ERR
      IF(ERR.EQ.-1)THEN
        TEST = 0.0
        ERR = 0
        print*,'Error in inversion'
        print*,'ALAMDA = ',ALAMDA
        print*,'TEST = ',TEST
      END IF
      IF(ALAMDA.EQ.0.)THEN
        RETURN
      ENDIF
C      print*,'TEST1',TEST
      IF(TEST.GT.0.0)THEN
        DO 16 J=1,MA
C        print*,A(J),DA(J,1)
        ATRY(J)=A(J)+DA(J,1)        
        IF(ATRY(J).LT.-100.0.OR.ATRY(J).GT.75)THEN
         PRINT*,'Gone pear-shaped'
         TEST=0.0
         ATRY(J)=A(J)
        ENDIF
C        print*,ATRY(J)
16     CONTINUE
      ELSE
       DO J=1,MA
        ATRY(J)=A(J)
       ENDDO
      ENDIF
C      print*,'TEST2',TEST

      IF(IMETHOD.EQ.0)THEN
C        print*,'calling sorta'
        CALL SORTA(MA,ATRY)
C        print*,'sorta OK'
        IF(ATRY(MA).GT.75.0)THEN
          TEST=0.0
          print*,'Finefitk. overflow',atry(ma)
        ENDIF
      ELSE
       DO J=1,MA
        IF(ATRY(J).LE.0.0)THEN
         print*,'Finefitk Info: ATRY < 0. Reset to A'
         TEST=0.0
        END IF
       ENDDO
      ENDIF

C      print*,'TEST3',TEST
      IF(TEST.LT.1.0)THEN
       PRINT*,'Error in Finefitk. K has gone too big.'
       CHISQ = 2*OCHISQ
      ELSE
C       PRINT*,'ALAMBDA,ATRY',ALAMDA,(ATRY(I),I=1,MA)
       CALL MRQCOFP7(X,Y,YFIT,SIG,NDATA,ATRY,MA,COVAR,
     1  DA,NCA,CHISQ,W,IMETHOD)
C       PRINT*,'AA: ALAMBDA,CHISQ',ALAMDA,CHISQ
      ENDIF

C      PRINT*,ALAMDA,OCHISQ,CHISQ
      IF(CHISQ.LE.OCHISQ)THEN
        ALAMDA=ALAMDA/3.0
        OCHISQ=CHISQ
        DO 18 J=1,MA
          DO 17 K=1,MA
            ALPHA(J,K)=COVAR(J,K)
17        CONTINUE
          BETA(J)=DA(J,1)
          A(J)=ATRY(J)
18      CONTINUE
      ELSE
         ALAMDA=10.0*ALAMDA
         IF(ALAMDA.GT.1.0D20)THEN
           AOVER=.TRUE.
           RETURN
         ENDIF 
         CHISQ=OCHISQ
      ENDIF

      RETURN
      END


      SUBROUTINE MRQCOFP7(X,Y,YFIT,SIG,NDATA,A,MA,
     & ALPHA,BETA,NALP,CHISQ,W,IMETHOD)
C     ************************************************************
C     Calculates YFIT,ALPHA and BETA for next iteration of the fitting
C     model
C
C     Input variables
C	X(NDATA)	DOUBLE	Amounts
C	Y(NDATA)	DOUBLE  Transmissions
C	SIG(NDATA)	DOUBLE	Errors
C	NDATA		integer Number of points in transmission curves
C	A(MA)		DOUBLE	Variables to vary to fit transmission curve
C	MA		integer Array length of A()
C	NALP		integer	defines size of arrays ALPHA and BETA
C	W(MA)		DOUBLE	Weights of k_g
C	IMETHOD		integer	Required representation of K_G() by A()
C
C     Output variables
C	YFIT(NDATA)	DOUBLE	Fitted transmissions
C	ALPHA(NALP,NALP)	DOUBLE	See Press et al., p522
C	BETA(MA)	DOUBLE	See Press et al., p522
C	CHISQ		DOUBLE	Chi^2 between Y() and calculated YFIT()
C
C     Pat Irwin		9/1/04
C
C     ************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMAX=20)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),ALPHA(NALP,NALP),BETA(MA),
     *    DYDA(MMAX),A(MA),W(MA), YFIT(NDATA)
      INTEGER IMETHOD

      DO 12 J=1,MA
        DO 11 K=1,J
          ALPHA(J,K)=0.
11      CONTINUE
        BETA(J)=0.
C        print*,'cof',j,a(j)
12    CONTINUE
      CHISQ=0.
      DO 15 I=1,NDATA
        IF(IMETHOD.EQ.0)THEN
         CALL EXPSUM(X(I),A,YMOD,DYDA,W,MA)
        ELSE
         CALL EXPSUM3(X(I),A,YMOD,DYDA,W,MA)
        ENDIF

        SIG2I=1./(SIG(I)*SIG(I))
        YFIT(I)=YMOD
        DY=Y(I)-YMOD
C        print*,'cof',I,Y(I),YFIT(I)
        DO 14 J=1,MA
          WT=DYDA(J)*SIG2I
          DO 13 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+WT*DYDA(K)
13        CONTINUE
          BETA(J)=BETA(J)+DY*WT
14      CONTINUE
        CHISQ=CHISQ+DY*DY*SIG2I
15    CONTINUE
      DO 17 J=2,MA
        DO 16 K=1,J-1
          ALPHA(K,J)=ALPHA(J,K)
16      CONTINUE
17    CONTINUE
      RETURN
      END

      SUBROUTINE EXPSUM(X,A,Y,DYDA,W,NA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(NA),A(NA),DYDA(NA)

      Y=0.
      DO 20 I=1,NA
       U = X*DEXP(A(I))
       Y=Y+W(I)*DEXP(-U)
       DYDA(I)=-W(I)*U*DEXP(-U)
C       IF(DYDA(I).EQ.0.0)THEN
C        PRINT*,'Error in EXPSUM. DYDA too small'
C        PRINT*,I,X,U,W(I)*U*DEXP(-U)
C       ENDIF
20    CONTINUE

      RETURN
      END


      SUBROUTINE EXPSUM3(X,A,Y,DYDA,W,NA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(NA),A(NA),DYDA(NA)
      DIMENSION TEMP(10)
        
      ARG=0.
      DO 10 I=1,NA
       ARG=ARG+A(I)
       TEMP(I)=DEXP(-ARG*X)
10    CONTINUE

      Y=0.
      DO 20 I=1,NA
       Y=Y+W(I)*TEMP(I)
       DYDA(I)=0.
c       DO 30 J=1,I                             ! Think this is right
       DO 30 J=I,NA                             ! But this works better???
        DYDA(J)=DYDA(J) - X*W(J)*TEMP(J)
30     CONTINUE 
20    CONTINUE

      RETURN
      END


      SUBROUTINE SORTA(MA,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      DIMENSION A(MA)
 
10    IREPEAT=0
      DO I=2,MA  
       IF(A(I).LT.A(I-1))THEN
        TMP = A(I)
        A(I)=A(I-1)
        A(I-1)=TMP
        IREPEAT=1
       ENDIF
      ENDDO

      IF(IREPEAT.EQ.1)GOTO 10

      RETURN

      END            
