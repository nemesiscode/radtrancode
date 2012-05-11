      subroutine esumset(knu0,delad,y0,EL,SFB,c1,c2,qrot,P,T,q,x,y,
     1 sig,ndata,gabsc,k_g,ng,eml)
C     **************************************************************
C     Calculates first approximation to the k-distribution, which is then 
C     improved upon by calc_esum5.f
C
C     Pat Irwin	?/?/??	Original
C     Pat Irwin	26/4/12	Commented.
C
C     **************************************************************

      implicit none
      include '../includes/arrdef.f'
      real knu0,delad,y0,EL,SFB,c1,c2,qrot,P,T,q,x(20),y(20)
      real sig(20), xcorr, tautmp,dpexp
      real gabsc(maxg),k_g(maxg),eml

      real U,tau_eks,T1,U1,U2,DU,E_TRAN,yv,tau_goody_voigt1
      real gin,kout,sds,tml,pi,p0,t0
      parameter (pi=3.1415927,p0=1.,t0=296.)
      integer NDATA,I,ng
      double precision SL,BL



      	U=1.
122	if (c1.ne.0.) then
		xcorr = c1 * u * (p/p0)*(t/t0)**c2
	else	
		xcorr = 0
	endif
	tautmp = xcorr + TAU_GOODY_VOIGT1(KNU0,DELAD,Y0,EL,SFB,
     1		QROT,P,T,U,Q)

	T1=DPEXP(-tautmp)
      	IF(T1.LT.0.5)THEN
       		U=U/2.
       		GOTO 122
      	END IF

123	if (c1.ne.0.) then
		xcorr = c1 * u * (p/p0)*(t/t0)**c2
	else	
		xcorr = 0
	endif
	tautmp = xcorr + TAU_GOODY_VOIGT1(KNU0,DELAD,Y0,EL,SFB,
     1		QROT,P,T,U,Q)
	T1=DPEXP(-tautmp)
      	IF(T1.GT.0.5)THEN
       		U=U*2.
       		GOTO 123
      	END IF

C     find upper limit of curve
      	U2=U
222	if (c1.ne.0.) then
		xcorr = c1 * u2 * (p/p0)*(t/t0)**c2
	else	
		xcorr = 0
	endif
	tautmp = xcorr + TAU_GOODY_VOIGT1(KNU0,DELAD,Y0,EL,SFB,
     1		QROT,P,T,U2,Q)
	T1=DPEXP(-tautmp)
      	IF(T1.GT.0.001)THEN
       		U2=U2*2.
       		GOTO 222
      	END IF

C     find lower part
      	U1=U
322	if (c1.ne.0.) then
		xcorr = c1 * u1 * (p/p0)*(t/t0)**c2
	else	
		xcorr = 0
	endif
	tautmp = xcorr + TAU_GOODY_VOIGT1(KNU0,DELAD,Y0,EL,SFB,
     1		QROT,P,T,U1,Q)
	T1=DPEXP(-tautmp)
      	IF(T1.LT.0.999)THEN
       		U1=U1/2.
       		GOTO 322
      	END IF

      	U1=LOG(U1)
      	U2=LOG(U2)


      	DU=(U2-U1)/(NDATA-1)

      	E_TRAN=0.0001

      	DO 101 I=1,NDATA
       		U=U1+DU*(I-1)
       		U=DPEXP(U)
		if (c1.ne.0.) then
			xcorr = c1 * u * (p/p0)*(t/t0)**c2
		else	
			xcorr = 0
		endif
		tautmp = xcorr + TAU_GOODY_VOIGT1(KNU0,DELAD,Y0,EL,SFB,
     1			QROT,P,T,U,Q)
		Y(I) = DPEXP(-tautmp)
       		X(I)=U
       		SIG(I)=E_TRAN
101   	CONTINUE

      yv=0.25*y0/delad

      print*,'knu0,yv,EL,SFB,c1,c2,qrot,P,T,q,SL,BL',knu0,yv,EL,SFB,
     1		c1,c2,qrot,P,T,q,SL,BL

      CALL ml_lacis(knu0,yv,EL,SFB,qrot,P,T,q,SL,BL)


      CALL fit_mlband(X,Y,SIG,NDATA,SL,BL)


      sds = 0.
      do i=1,ndata
       U=x(i)
       tml = sngl(dexp(-0.5*PI*BL*(DSQRT(1.+4.*SL*U/(PI*BL))-1.)))
       sds = sds + (y(i) - tml)**2
      end do

      sds = sds/ndata
      sds = 100*sqrt(sds)

      eml = sds


      do i=1,ng
       gin=gabsc(i)
       call kml(SL,BL,gin,kout)
       k_g(i)=kout
      end do

      return

      end
