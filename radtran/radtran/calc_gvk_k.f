      subroutine calc_gvk_k(fknu,delad,y,T,umean,lcalch,g_ord,k_g,ng)
C ***************************************************************************
C
C     Calculates the mean value of k(g) in a g-interval using the 
C     Goody-Voigt model
C
C     Input variables
C	knu		real 	knu at temperature of layer
C	y		real 	aL/aD at pressure and temperature of layer
C       delad		real 	mean line spacing/ad0
C  	T	        real 	layer temperature
C       g2              real    Minimum value of g in interval
C       g1              real	Maximum value of g in interval
C
C     Output variable.
C	calc_gvk	real	mean k-coefficient from gmin to gmax
C
C	Pat Irwin	?/?/??	Original
C	Pat Irwin	26/4/12	Commented.
C
C ***************************************************************************
C
C -Fits set of analytical functions to transmittance curve
C -Performs inverse Laplacian transform to obtain f(k)
C -Finds Gibbs error and makes new limits of f(k)
C -Removes Gibbs error again and calculates g(k) and h(k,u)
C -Fits function to g(k) and h(k) in order to calculate mean transmission
C  between g1 and g2
C
C Set number of points in function 
C Set number of fitting functions
C Twenty seems to be about ok (from trial and error!)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include '../includes/arrdef.f'
      parameter(ndata=20,ma=20,mfit=20,ncvm=20)
      double precision x(ndata),ytau(ndata),offset(ma)
      double precision rk(ndata),fk(ndata),gk(ndata),hk(ndata)
      double precision lrk(ndata)
C
C Input parameters
C
      integer ng
      real fknu,delad,y,T,g_ord(maxg),k_g(maxg),tau_goody_voigt2
      real umean
      logical lcalch
C
C
C Parameters for the least squares fit
C
      double precision a(ma),covar(ncvm,ncvm),sig(ndata),yfit(ndata)
      integer lista(ma)
C
C Parameters used in fitting h(g) and k(g)
C
      double precision dk2(ndata),dh2(ndata)
C
C Set limits of U, the abcissa in absorber space 
C This has to be done carefully to minimise error in transform 
C 
C Limits are set so that ytau(U1) = 0.99 and ytau(U2)=0.01
C Firstly find where this is the case by linear interpolation
C in log10(U) space
C
C Set first approximate limits for U
C

      if(fknu.eq.0.or.delad.eq.0.or.y.eq.0.or.umean.eq.0.)then
       do j=1,ng-1
        k_g(j)=0.
       end do
       return
      end if

      U1=-4.0-log10(fknu)
      U2=U1+10.0+log10(delad)

C
C Now find where ytau(U)=0.99 and 0.01
C
      imin=2
      imax=ndata
      lmin=0
      lmax=0
C
      do 15 i=1,ndata
        UL=U1 + (U2-U1)*real(i-1)/real(ndata-1)
        U=10.0**UL
        x(i)=U
        ytau(i)=exp(-tau_goody_voigt2(fknu,delad,y,T,SNGL(U)))
C
        if(ytau(i).lt.0.99.and.lmin.eq.0)then
          imin=i
          lmin=1
        endif
        if(ytau(i).lt.0.01.and.lmax.eq.0)then
          imax=i
          lmax=1
        endif    
C
15    continue
C
      dxmin=log10(x(imin))-log10(x(imin-1))
      dymin=ytau(imin)-ytau(imin-1)
      dxmax=log10(x(imax))-log10(x(imax-1))
      dymax=ytau(imax)-ytau(imax-1)
C
C Now set limits in absorber space using simple linear interpolation 
C
      U1=log10(x(imin-1))+dxmin*(0.99-ytau(imin-1))/dymin
      U2=log10(x(imax-1))+dxmax*(0.01-ytau(imax-1))/dymax
C
C Set arbitrary limits for inverse space abcissa
C
      F1=-U2-1.0
      F2=-U1-1.0
C
C Set offsets for fitted function so that limits are U1 to U2
C Offsets are linearly spaced in Log_10(U)
C
      do 20 i=1,ma
        UL=U1 + (U2-U1)*real(i-1)/real(ma-1)
        offset(i)=10.0**UL
20    continue
C
      do 25 i=1,ndata
        UL=U1 + (U2-U1)*real(i-1)/real(ndata-1)
        U=10.0**UL
        x(i)=U
        ytau(i)=exp(-tau_goody_voigt2(fknu,delad,y,T,SNGL(U)))
C
        FL=F1 + real(i-1)*(F2-F1)/real(ndata-1)
        F=10.0**FL
        rk(i)=F
        lrk(i)=FL
25    continue

C      print*,'data selected'

C               __
C Fit a series  \  i!(a_i)/(p+c)^i+1 where i=1,2,3...n
C               /_
C Initialise errors in y and indices of parameters to be adjusted
C
      do 30 i=1,ndata
        sig(i)=1.0
30    continue
      do 31 i=1,ma
        a(i)=1.0
        lista(i)=i
31    continue
C
C The fitting routine from Numerical recipes (Press etal, 1986)
C
      call lfit_m(x,ytau,sig,ndata,a,ma,lista,mfit,covar,ncvm,
     +          chisq,offset)
C
C Now transform function: fk=SUM(exp(-at)) ...the PDF
C

      do 50 i=1,ndata
        fk(i)=0.0
        do 52 j=1,ma          
          fk(i)=fk(i)+a(j)*exp(-offset(j)*rk(i))
52       continue
50    continue
C
C Find wings of function fk(rk)
C
      call wingf(rk,fk,F1,F2)
C
C Obtain new array of rk so that it is stretched over domain {F1,F2}
C
      do 54 i=1,ndata
        FL=F1 + real(i-1)*(F2-F1)/real(ndata-1)
        F=10.0**FL
        rk(i)=F
54    continue
C
C Using new rk array, calculate fk again as sum of exponentials
C
      do 56 i=1,ndata
        fk(i)=0.0
        do 57 j=1,ma
          fk(i)=fk(i)+a(j)*exp(-offset(j)*rk(i))
57      continue
56    continue
C
C Find wings once more, and now remove them
C
      call wingr(rk,fk)
C
C Explicitly integrate fk to find gk ...the cumulative PDF
C
      do 60 i=2,ndata
        fint=0.0
        do 65 j=2,i
          efk1=fk(j-1)
          efk2=fk(j)
          trapez=(rk(j)-rk(j-1)) * (efk1+efk2)/2.0
          fint=fint+trapez
65      continue
        gk(i)=fint
60    continue
      gk(1)=0.0
C
C Normalise the fitted functions fk and gk using assumption gk(ndata)=1
C
      fac=gk(ndata)
      do 70 i=1,ndata
        gk(i)=gk(i)/fac
        fk(i)=fk(i)/fac
70    continue
C
C Calculate h(k)
C
      do 90 i=2,ndata
        fint=0.0
        do 100 j=2,i
          efk1=fk(j-1)*exp(-umean*rk(j-1))
          efk2=fk(j)*exp(-umean*rk(j))
          trapez=(rk(j)-rk(j-1)) * (efk1+efk2)/2.0
          fint=fint+trapez
100     continue
        hk(i)=fint
90    continue
      hk(1)=0.0
C
C Given g1 and g2, calculate mean k in interval
C Do this by cubic spline interpolation to find h(g) 
C Use linear interpolation (in log10(k) space) to get k_mean(g)
C
C Get limits of g(k) to use with interpolation
C
      if(lcalch)then
        if(gk(1).eq.gk(2))then
          print*,'gk(1) = gk(2) = ',gk(1)
          stop
        end if
        dhdg1=(hk(2)-hk(1))/(gk(2)-gk(1))
        mpoint=ndata
120     continue
        if(abs(gk(mpoint-1)-1D0).lt.1.D-13)then
           mpoint=mpoint-1
           goto 120
        end if 
        dhdgn=(hk(mpoint)-hk(mpoint-1))/(gk(mpoint)-gk(mpoint-1))
        call spline_m(gk,hk,mpoint,dhdg1,dhdgn,dh2)
C
        do 150 i=1,ng-1
          gmind=DBLE(g_ord(i))
          gmaxd=DBLE(g_ord(i+1))
          call splint_m(gk,hk,dh2,mpoint,gmind,rhmin)
          call splint_m(gk,hk,dh2,mpoint,gmaxd,rhmax)
          deltah=rhmax-rhmin
          deltag=gmaxd-gmind
C
C Check for deltah=0, corresponding to no transmission
C
          if(deltah.le.1.0E-15)then
            print*,'DELTA-H ZERO: NO TRANSMISSION'
            k_g(i)=0.0
          else
            k_g(i)=SNGL(-log(deltah/deltag)/umean)
          endif
150     continue
C
      else
        if(gk(1).eq.gk(2))then
           print*,'gk(1) = gk(2) = ',gk(1)
           stop
        end if
        mpoint=ndata
125     continue
        if(abs(gk(mpoint-1)-1D0).lt.1D-13)then
           mpoint=mpoint-1
           goto 125
        end if
        do 160 i=1,ng-1
          gval=DBLE(g_ord(i)+g_ord(i+1))/2.0
          lmin=0
          imin=1
          do 180 j=1,ndata
          if(gk(j).gt.gval.and.lmin.eq.0)then
            imin=j-1
            lmin=1
          endif
180       continue
          dlkdg=(lrk(imin+1)-lrk(imin))/(gk(imin+1)-gk(imin))
          rlk=lrk(imin)+dlkdg*(gval-gk(imin))
          k_g(i)=SNGL(10.0**rlk)
160     continue
      endif
C
      return
      end
C
      SUBROUTINE LFIT_M(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,NCVM,
     *CHISQ,offset)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMAX=50)
      DOUBLE PRECISION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),
     *    COVAR(NCVM,NCVM),BETA(MMAX),AFUNC(MMAX),offset(ma)
      INTEGER LISTA(MA)
      KK=MFIT+1
      DO 12 J=1,MA
        IHIT=0
        DO 11 K=1,MFIT
          IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
        IF (IHIT.EQ.0) THEN
          LISTA(KK)=J
          KK=KK+1
        ELSE IF (IHIT.GT.1) THEN
          PRINT*, 'Improper set in LISTA - IHIT.GT.1'
          STOP
        ENDIF
12    CONTINUE
      IF (KK.NE.(MA+1)) THEN
        PRINT*, 'Improper set in LISTA'
        STOP
      ENDIF
      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE
      DO 18 I=1,NDATA
        CALL FUFIT(X(I),AFUNC,MA,offset)
        YM=Y(I)
        IF(MFIT.LT.MA) THEN
          DO 15 J=MFIT+1,MA
            YM=YM-A(LISTA(J))*AFUNC(LISTA(J))
15        CONTINUE
        ENDIF
        SIG2I=1./SIG(I)**2
        DO 17 J=1,MFIT
          WT=AFUNC(LISTA(J))*SIG2I
          DO 16 K=1,J
            COVAR(J,K)=COVAR(J,K)+WT*AFUNC(LISTA(K))
16        CONTINUE
          BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18    CONTINUE
      IF (MFIT.GT.1) THEN
        DO 21 J=2,MFIT
          DO 19 K=1,J-1
            COVAR(K,J)=COVAR(J,K)
19        CONTINUE
21      CONTINUE
      ENDIF
      CALL GAUSSJ_M(COVAR,MFIT,NCVM,BETA,1,1)
      DO 22 J=1,MFIT
        A(LISTA(J))=BETA(J)
22    CONTINUE
      CHISQ=0.
      DO 24 I=1,NDATA
        CALL FUFIT(X(I),AFUNC,MA,offset)
        SUM=0.
        DO 23 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
23      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
24    CONTINUE
      RETURN
      END

      SUBROUTINE GAUSSJ_M(A,N,NP,B,M,MP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=50)
      DOUBLE PRECISION A(NP,NP),B(NP,MP),IPIV(NMAX)
      INTEGER INDXR(NMAX),INDXC(NMAX)
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
                PRINT*, 'Singular matrix'
		STOP
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
          PRINT*,'Singular matrix.'
	  STOP
	ENDIF
        PIVINV=1./A(ICOL,ICOL)
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

      subroutine wingr(rk,fk)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      parameter (ndata=20)
      double precision rk(ndata),fk(ndata)
C
C Find wings of function fk(rk), and remove them
C Firstly find maximum of fk by simple search
C
      idum=1
      do 10 i=2,ndata
        if(fk(i).gt.fk(idum))idum=i
10    continue
C
C Now eliminate wings from function
C
      iwing=0
      do 20 i=idum,1,-1
        if(fk(i).le.0.0.or.iwing.eq.1)then
          iwing=1
          fk(i)=1.0e-10
        endif
20    continue
C
      iwing=0
      do 30 i=idum+1,ndata
        if(fk(i).le.0.0.or.iwing.eq.1)then
          iwing=1
          fk(i)=1.0e-10
        endif
30    continue
C
      return
      end


      subroutine wingf(rk,fk,F1,F2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      parameter (ndata=20)
      double precision rk(ndata),fk(ndata)
C
C Find wings of function fk(rk)
C Firstly find maximum of fk by simple search
C
      idum=1
      do 10 i=2,ndata
        if(fk(i).gt.fk(idum))idum=i 
10    continue
C
C Now eliminate wings from function and find new limits of fk
C
      iwing=0
      do 20 i=idum,1,-1
        if(fk(i).le.0.0.and.iwing.eq.0)then
          iwing=1
          F1=log10(rk(i))
        endif
20    continue

      iwing=0
      do 30 i=idum+1,ndata
        if(fk(i).le.0.0.and.iwing.eq.0)then
          iwing=1
          F2=log10(rk(i))
        endif
30    continue
C
      return
      end

      subroutine fufit(x,afunc,ma,offset)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                         
C evaluate a/(p^2+a^2) where a=1,2,3...ma
C
      parameter(mmax=50)
      double precision afunc(mmax),offset(ma)
C
      do 200 j=1,ma
        afunc(j)=1.0/(x+offset(j))
200   continue
C
      return
      end

      real function tau_goody_voigt2(knu,delad,y,T,U)
C     $Id: calc_gvk_k.f,v 1.3 2011-06-17 15:40:25 irwin Exp $
C     *********************************************************************
C
C     Version 1.1
C     Calculate the transmission of a path using the the formulation of the
C     Goody-Voigt random band model given in Strong (1992) D.Phil.Thesis
C
C     Integration of the function over x is performed in two stages using 
C     Simpsons rule with 101 points over the ranges:
C	1)	x = 0 to 30*y
C 	2)	x = 30*y to 3030*y
C     This is the same method used in tran_eks.f by K.Strong.
C
C     ----------------------------------------------------------------------
C     
C     Version 1.2
C     Version 1.1 gives spurious results if the coefficient of V(x,y) in the
C     denominator of the integral is large. We really want to set the first
C     limit of x to be where the integrand is 1% of the integrand when x=0.
C
C     The integrand may be written:
C	B(x,y) =         V(x,y)
C 		    ---------------
C	  	      1 + A*V(x,y)
C
C     Put C = 0.01*Bmax = 0.01*B(0,y). When the function B is 1% of max
C     then:
C		V1 = V(x,y) = C/(1-AC)
C
C     In the wings, the Voigt line shape looks like a Lorentz line. Thus:
C  	 	(nu - nu0) is proportional to 1/sqrt(k)
C     Now when x=30y, V(x,y) is always less than 1% of Vmax = V(0,y). Thus
C     Putting D = V(0,y), when the factor A is large, the integration limit
C     must be scaled to:
C		x_int = 30*y*sqrt(0.01*D/V1) 
C
C     and the integration is then more appropriate. Further checks to ensure
C     numerical accuracy are in progress.
C 
C     ------------------------------------------------------------------------
C     Input Variables
C
C	delad	real	line spacing / alpha_D0 (K**-0.5)
C	P	real	Mean pressure of path (atm)
C	T	real	Mean temperature of path (K)
C	U	real	Absorber amount in path (*1e20 (molec/cm2))
C	q	real	Fractional abundance of active gas 
C
C     Output variable
C
C 	tau_goody_voigt2	real	Transmission of path
C     ------------------------------------------------------------------------
C
C	Pat Irwin	26/9/94
C
C     *********************************************************************
      implicit none
      real knu,delad,U,y,T,C1,T0,P,P0
      integer i,idim,ndim
      parameter (P0 = 1.,T0=296.0, C1=1.439,idim=500)
      real SQT0
      parameter (SQT0 = 17.20465053)		! Sqrt(296.0)
      real x,dx,s(idim),A1,A2,humlic,simp_int_m,Aconst,lx,Cconst,Vconst
      real Bconst,Dconst,Ratio,Tconst

      if(knu.eq.0.or.U.eq.0.)then
       tau_goody_voigt2=0.
      else
C
      Dconst = humlic(0.,y)
      Aconst = knu*U*delad/sqrt(T)
      Bconst = Dconst/(1. + Aconst*Dconst)
      Cconst = 0.01*Bconst      

      Vconst = Cconst / (1. - Cconst*Aconst)

      Ratio = sqrt(0.01*Dconst/Vconst)
     
      lx = 30.*y*Ratio

777   Tconst = humlic(lx,y)/(1. + Aconst*humlic(lx,y))
      if(Tconst.gt.Cconst)then
       lx = lx*1.5
       goto 777
      endif

      do 10 i=1,101
       x=lx*(i-1)/100.
       s(i)=humlic(x,y)/(1. + Aconst*humlic(x,y))
10    continue

      ndim=101
      dx = lx/100.

      A1=simp_int_m(s,idim,ndim,dx)
 
      do 20 i=1,101
       x=lx + (i-1)*lx
       s(i)=humlic(x,y)/(1. + Aconst*humlic(x,y))
20    continue      

      ndim=101
      dx = lx

      A2=simp_int_m(s,idim,ndim,dx)

      tau_goody_voigt2=2.*knu*U*(A1+A2)
  
      end if

      return
      end



      SUBROUTINE SPLINE_M(X,Y,N,YP1,YPN,Y2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=100)
      DOUBLE PRECISION X(N),Y(N),Y2(N),U(NMAX)
C      print*,yp1,ypn
C      do i=1,n
C       print*,x(i),y(i)
C      end do
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

      SUBROUTINE SPLINT_M(XA,YA,Y2A,N,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) THEN
       PRINT*, 'Bad XA input.'
       STOP
      ENDIF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END

       FUNCTION SIMP_INT_M(X,NDIM,N,FH)
       IMPLICIT REAL (A-H,O-Z)
C
       DIMENSION X(NDIM)
       IF(0.5*N.EQ.INT(0.5*N))GOTO 999
       SUM=X(1)+X(N)
       DO 10 I=2,N-1,2
        SUM=SUM+4*X(I)
10     CONTINUE
       DO 20 I=3,N-2,2
        SUM=SUM+2*X(I)
20     CONTINUE
       SIMP_INT_M=FH*SUM/3.0
       RETURN
999    STOP 'Simpson: Integer N must be ODD'
       END

