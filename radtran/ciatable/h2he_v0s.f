      subroutine h2he_v0s(norm,temp1,fnumin,fnumax,dnu,nf1,
     1freq1,alfatot)
C     *************************************************************
C     Subroutine to calculate an H2-He collision induced absorption spectrum
C     in the v0 rotational-translational band.
C
C     Input variables:
C       norm            integer 0 = equilibrium hydrogen
C                               1 = 3:1 ortho/para hydrogen
C       temp            double  Temperature (K)
C       fnumin          double  Lowest wavenumber in spectrum (cm-1)
C       fnumax          double  Highest wavenumber in spectrum (cm-1)
C       dnu             double  Wavenumber step (cm-1)
C
C     Output variables
C       nf1             integer Number of points in spectrum
C       freq1(601)      double  Frequency grid (cm-1)
C       alfatot(601)    double  Absorption coefficient (cm-1 amagat-2)
C
C     Pat Irwin         31/1/96
C     Pat Irwin		2/3/12	Updated for Radtrans2.0
C
C     *************************************************************


C ============================================================================
C Copyright (C) 1993  Aleksandra Borysow  (April 1993)
C ============================================================================
C 
C Copyright Notice:
C You may use this program  for your scienftific applications,
C but please do not  distribute it yourself.
C Direct all your inquires to the  author: e-mail aborysow@phy.mtu.edu
C Final request: If you publish  your work which benefited from this
C program, please acknowledge using this program and QUOTE 
C the original paper describing the procedure: 
c	
c     Jacek Borysow, L.Frommhold, G.Birnbaum, 
c	Astrophys.J. vol. 326, p.509-515, (1988)

C ============================================================================

C	implicit double precision (a-h,o-z)
C
c 	===============    R/T H2-He CIA  40 - 3000K =========================

C     THESE DATA SERVE AS THE INPUT TO THE MAIN PROGRAM ADDEM.
C     TEMP = TEMPERATURE IN KELVIN, SHOULD BE BETWEEN 40. AND 3000.
C     FNUMIN = LOWEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU)
C     FNUMAX = HIGHEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU)
C     DNU = FREQUENCY INCREMENT IN CM-1. DNU SHOULD BE CHOSEN SO
C           THAT NOT MORE THAN 600 STEPS ARE NEEDED TO GO FROM
C           FNUMIN TO FNUMAX (ELSE ARRAY DIMENSIONS OF FREQ,ABSCOEF
C           MUST BE ADJUSTED IN ADDEM).
C     NORMAL=1 ADJUSTS THE H2 STAT.WEIGHTS SO THAT THE ORTHO/PARA RATIO
C     EQUALS 3:1. SET NORMAL=0 FOR EQUILIBRIUM HYDROGEN.
C
C     THIS PROGRAM GENERATES THE H2-He ROTATIONAL/TRANSLATIONAL 
C     CIA SPECTRA. IT IS BASED ON QUANTUM LINE SHAPE COMPUTATIONS AND
C     THE AB INITIO DIPOLE DATA BY W. MEYER. 
C
C     File out.h2he IS OUTPUT: HEADER PLUS ABSORPTION COEFF. ALPHA(NU)
C
	implicit double precision (a-h,o-z)
	dimension freq1(601),alfatot(601)
      COMMON /H2PARTB/ Q,WH2(2),B0,D0,JRANGE1,NORMAL
      COMMON /RESULT1/ NF
      common /result2/ FREQ(601),ABSCOEF(601)
      Y(X,A,B,C,D)=A*DEXP(B*X + C*X*X + D*X*X*X)
C

      TEMP=TEMP1
      NORMAL=NORM

c initialise the weights for H2
      WH2(1)=1.0
      WH2(2)=3.0
c this is needed otherwise they are not resetted after a normal calculation

      NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1
      NF1=NF
      IF (NF.GT.601) NF=601
      if(temp.le.40.d0 .or. temp.gt.3000.d0) then
       print*,'h2he_v0s: Warning'
       print*,'Temperature should be 40 < T < 3000'
      end if

      CALL PARTSUM3 (TEMP)
C
C     THE H2-He SPECTRA: 40K -- 3000K
C     =================
C
      X=DLOG(TEMP)
      DO 10 I=1,NF
      FREQ(I)=FNUMIN+DFLOAT(I-1)*DNU
	alfatot(i)=0.d0
	freq1(i)=freq(i)
   10 ABSCOEF(I)=0.
C
C     THE LAMBDA,L = 23 COMPONENT:
C     (QUADRUPOLE PLUS OVERLAP INDUCED COMPONENT)
      XM0=Y(X,.262d-62, -0.71411d0, 0.14587d0, -0.00352d0 )
      XM1=Y(X, 0.269d-49, -0.98315d0, 0.21988d0, -0.00729d0)
      XM2=Y(X, 0.406d-35, -2.25664d0, 0.50098d0, -0.01925d0)
	call BCparam(temp, XM0,XM1,XM2, t1, t2)
	ifun = 1
      CALL ADDSPEC3(XM0,T1,T2,TEMP,0,0,-1,-1,2,3,0,0,ifun)
	do 111 i=1, nf
C	write(3,*) freq(i), abscoef(i)
111	alfatot(i) = alfatot(i) + abscoef(i)

C     PARAMETERS FOR 21 (PURE OVERLAP) COMPONENT
      XM0=Y(X, 0.424d-64, -0.37075d0, 0.17473d0, -0.00489d0)
      XM1=Y(X, 0.174d-49, -1.89232d0, 0.44399d0, -0.02029d0)
      XM2=Y(X, 0.874d-35, -3.27717d0, 0.69166d0, -0.02865d0)
	call K0param(temp, XM0, XM1, XM2, t1, t2)
	ifun = 0
      CALL  ADDSPEC3(XM0,T1,T2,TEMP,0,0,-1,-1,2,1,0,0,ifun)
	do 311 i=1, nf
C	write(3,*) freq(i), abscoef(i)
311	alfatot(i) = alfatot(i) + abscoef(i)
C
C     PARAMETERS FOR 01 (PURE OVERLAP) COMPONENT
      XM0=Y(X, 0.223d-61, -1.89198d0, 0.45505d0, -0.02238d0)
      XM1=Y(X, 0.625d-48, -1.96486d0, 0.47059d0, -0.02402d0)
      XM2=Y(X, 0.316d-33, -3.40400d0, 0.72793d0, -0.03277d0)
 	S=XM0
	call K0param(temp, XM0, XM1, XM2, t1, t2)
	ifun = 0
      CALL  ADDSPEC3(XM0, T1,T2,TEMP,0,0,-1,-1,0,1,0,0, ifun)
	do 319 i=1, nf
C	write(3,*) freq(i), abscoef(i)
319	alfatot(i) = alfatot(i) + abscoef(i)

C      WRITE (3,140) FREQ(1),FREQ(NF),FREQ(2)-FREQ(1)
C      WRITE (3,150) (Alfatot(I),I=1,NF)
c	write(3, 160) (freq(i), alfatot(i)*1.d7, i=1,nf)
160	format( f7.1, e12.4)


  140 FORMAT (/,' Total ABSORPTION COEFFICIENT ALPHA(fnu)',/,
     3  ' FROM',F8.1,' CM-1 TO',F8.1,' CM-1, AT',F6.2,
     3  ' CM-1 INCREMENTS',' IN UNITS OF CM-1 AMAGAT-2',/)
  150 FORMAT (  6e12.4)

C
   30 FORMAT ( 34H1ABSORPTION SPECTRA OF HYDROGEN AT,F8.1,  2H K,/1X,43(
     11H=),/, 11H MIN.FREQ.=,F8.1,  5H CM-1,10X, 10HMAX.FREQ.=,F8.1,  5H
     2 CM-1,10X, 15HFREQ.INCREMENT=,F8.2,  5H CM-1,5X,  2HIN,I5,  6H STE
     3PS,//)
   40 FORMAT ( 51H1ABSORPTION SPECTRUM OF HYDROGEN-HELIUM MIXTURES AT,F8
     1.1,  2H K,/1X,59(1H=)/3F15.1,I10//)

      RETURN

      END


      SUBROUTINE ADDSPEC3 (G0,TAU1,TAU2,TEMP,MP,LIKE,LAMBDA1
     1,LAMBDA2,LAMBDA,LVALUE,NVIB1,NVIB2,ifun)
C
C     THIS PROGRAM GENERATES A LISTING OF THE CIA TR ALFA(OMEGA)
C     IF BOTH LAMBDA1 AND LAMBDA2 ARE NEGATIVE: SINGLE TRANSITIONS;
C     DOUBLE TRANSITIONS ARE ASSUMED OTHERWISE.
C     MP=1 GIVES LISTINGS OF INTERMEDIATE RESULTS.
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.
C
	implicit double precision (a-h,o-z)
      COMMON /H2PARTB/ Q,WH2(2),B0,D0,JRANGE1,NORMAL
      COMMON /RESULT1/ NF
      common /result2/ FREQ(601),ABSCOEF(601)
      DATA CLOSCHM,BOLTZWN/2.68675484E19,.6950304/
      DATA HBAR,PI,CLIGHT/1.054588757d-27,3.1415926535898,2.9979250E10/
      EH2(N,I)=4395.34*(DFLOAT(N)+0.5)-117.905*(DFLOAT(N)+0.5)**2
     2  +(60.809-
     1  2.993*(DFLOAT(N)+0.5)+.025*(DFLOAT(N)+.5)**2)*DFLOAT(I)-
     3 (.04648-.00134*(DFLOAT(N)+.5))*DFLOAT(I*I)
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*
     1  DEXP(-1.4387859/T*EH2(0,J*(J+1)))

	do 777 i=1, nf
	abscoef(i) = 0.d0
777	continue

      TWOPIC=2.*PI*CLIGHT
c      IF (MP.NE.1) MP=0
c      IF (LIKE.NE.1) LIKE=0
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2
      CALIB=CALIB/DFLOAT(1+LIKE)
      BETA=1./(BOLTZWN*TEMP)
      LIST=NF
C
C     ROTATIONAL SPECTRUM FOR THE DETAILED LISTING   *******************
C
C      WRITE (3,110) LAMBDA,LVALUE,G0,TAU1,TAU2
  110 FORMAT (/, ' LAMBDA,LVALUE=',2I3,'  COMPONENT INCLUDED.',/,
     1  ' LINE SHAPE PARAMETERS:',3E12.3 ,/ )
C
      IF ((LAMBDA1.LT.0).AND.(LAMBDA2.LT.0)) GO TO 60
      JPLUSL=JRANGE1+MAX0(LAMBDA1,LAMBDA2)
      DO 50 I1=1,JRANGE1
         J1=I1-1
      DO 50 IP1=1,JPLUSL
         JP1=IP1-1
         CG1S=CLEBSQR3(J1,LAMBDA1,JP1)
         IF (CG1S) 50,50,10
   10    P1=PH2(J1,TEMP)/Q
         IF (P1.LT.0.001) GO TO 50
         OMEGA1=EH2(NVIB1,JP1*IP1)-EH2(0,J1*I1)
         DO 40 I2=1,JRANGE1
            J2=I2-1
         DO 40 IP2=1,JPLUSL
            JP2=IP2-1
            CG2S=CLEBSQR3(J2,LAMBDA2,JP2)
            IF (CG2S) 40,40,20
   20       P2=PH2(J2,TEMP)/Q
            IF (P2.LT.0.001) GO TO 40
            OMEGA2=EH2(NVIB2,JP2*IP2)-EH2(0,J2*I2)
            FAC=CALIB*P1*P2*CG1S*CG2S
            DO 30 I=1,LIST
               FRQ=FREQ(I)-OMEGA1-OMEGA2
               WKI=FREQ(I)*(1.-DEXP(-BETA*FREQ(I)))
               WKF=WKI*FAC
               XBG=G0*BGAMA3(FRQ,TAU1,TAU2,TEMP,ifun)
               ABSCOEF(I)=ABSCOEF(I)+XBG*WKF
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      GO TO 100
c	Single transitions here:
   60 JPLUSL=JRANGE1+LAMBDA
      DO 90 I=1,JRANGE1
         J=I-1
      DO 90 IP=1,JPLUSL
         JP=IP-1
         CGS=CLEBSQR3(J,LAMBDA,JP)
         IF (CGS) 90,90,70
   70    P=PH2(J,TEMP)/Q
         IF (P.LT.0.001) GO TO 90
         OMEGA1=EH2(NVIB1,JP*IP)-EH2(0,J*I)
         FAC=CALIB*P*CGS
         DO 80 IQ=1,LIST
            FRQ=FREQ(IQ)-OMEGA1
            WKI=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))
            WKF=WKI*FAC
            XBG=G0*BGAMA3(FRQ,TAU1,TAU2,TEMP,ifun)
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF
   80    CONTINUE
   90 CONTINUE
  100 CONTINUE
C      WRITE (3,140) FREQ(1),FREQ(NF),FREQ(2)-FREQ(1)
C      WRITE (3,150) (ABSCOEF(I),I=1,LIST)
      RETURN
C
  140 FORMAT (/,' ABSORPTION COEFFICIENT ALPHA(fnu)',/,
     3  ' FROM',F8.1,' CM-1 TO',F8.1,' CM-1, AT',F6.2,
     3  ' CM-1 INCREMENTS',' IN UNITS OF CM-1 AMAGAT-2',/)
  150 FORMAT (  6E12.4)
C
      END

      FUNCTION CLEBSQR3 (L,LAMBDA,LP)
C
C     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)
C     FOR INTEGER ARGUMENTS ONLY
C     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.
C
	implicit double precision (a-h,o-z)
      FC=DFLOAT(2*LP+1)
      GO TO 10
C
      ENTRY THREEJ23
C
C     THIS ENTRY RETURNS THE SQUARED 3-J SYMBOL   L LAMBDA LP
C                                                 0    0    0
C     INSTEAD OF THE CLEBSCH-GORDAN COEFFICIENT
C     (LIMITATION TO INTEGER ARGUMENTS ONLY)
C
C     NOTE THAT THE THREE-J SYMBOLS ARE COMPLETELY SYMMETRIC IN THE
C     ARGUMENTS. IT WOULD BE ADVANTAGEOUS TO REORDER THE INPUT ARGUMENT
C     LIST SO THAT LAMBDA BECOMES THE SMALLEST OF THE 3 ARGUMENTS.
C
      FC=1.
   10 CLEBSQR3=0.
      IF (((L+LAMBDA).LT.LP).OR.((LAMBDA+LP).LT.L).OR.((L+LP).LT.LAMBDA)
     1) RETURN
      IF (MOD(L+LP+LAMBDA,2).NE.0) RETURN
      IF ((L.LT.0).OR.(LP.LT.0).OR.(LAMBDA.LT.0)) RETURN
      F=1./DFLOAT(L+LP+1-LAMBDA)
      IF (LAMBDA.EQ.0) GO TO 30
      I1=(L+LP+LAMBDA)/2
      I0=(L+LP-LAMBDA)/2+1
      DO 20 I=I0,I1
   20 F=F*DFLOAT(I)/DFLOAT(2*(2*I+1))
   30 P=FC*F*FCTL3(LAMBDA+L-LP)*FCTL3(LAMBDA+LP-L)
      CLEBSQR3=P/(FCTL3((LAMBDA+L-LP)/2)*FCTL3((LAMBDA+LP-L)/2))**2
      RETURN
C
      END
      FUNCTION FCTL3 (N)
	implicit double precision (a-h,o-z)
      P(Z)=((((-2.294720936d-4)/Z-(2.681327160d-3))/Z+(3.472222222d-3))/
     1Z+(8.333333333d-2))/Z+1.
      FCTL3=1.
      IF (N.LE.1) RETURN
      IF (N.GT.15) GO TO 20
      J=1
      DO 10 I=2,N
   10 J=J*I
      FCTL3=DFLOAT(J)
      RETURN
   20 Z=DFLOAT(N+1)
      FCTL3=(DEXP(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)
      RETURN
      END


      SUBROUTINE PARTSUM3 (TEMP)
C
C     H2 ROTATIONAL PARTITION SUM Q = Q(T).
C
	implicit double precision (a-h,o-z)
      COMMON /H2PARTB/ Q,WH2(2),B0,D0,JRANGE1,NORMAL
C
C     NORMAL=0 INITIATES EQUILIBRIUM HYDROGEN;
C     NORMAL=1 ADJUSTS FOR ORTHO/PARA HYDROGEN RATIO OF 3:1
C
      DATA B0,D0,WH2(1),WH2(2)/59.3392,0.04599,1.,3./
      EH2(N,I)=4395.34*(DFLOAT(N)+.5)-117.905*(DFLOAT(N)+.5)**2
     2  +(60.809-2.993*(DFLOAT(N)+.5)+.025*(DFLOAT(N)+.5)**2)*
     2  DFLOAT(I)- (.04648-.00134*(DFLOAT(N)+.5))*DFLOAT(I*I)
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*DEXP(-1.4387859*
     1  EH2(0,J*(J+1)
     1)/T)
C
C     Q,B0,D0,WH2 - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR H2
C
      Q=0.
      J=0
   10 DQ=PH2(J,TEMP)
      Q=Q+DQ
      J=J+1
      IF (DQ.GT.Q/900.) GO TO 10
      JRANGE1=J
      IF (NORMAL) 20,20,30
   20 CONTINUE
C     WRITE (3,60) Q,JRANGE1,TEMP
      RETURN
   30 J=-1
      S=0.
      SEV=0.
   40 J=J+1
      DS=PH2(J,TEMP)
      S=S+DS
      SEV=SEV+DS
      J=J+1
      S=S+PH2(J,TEMP)
      IF (J.GT.JRANGE1) GO TO 50
      GO TO 40
   50 SODD=S-SEV
      WH2(2)=WH2(2)*3.*SEV/SODD
      Q=4.*SEV
C      WRITE (3,70) Q,JRANGE1,TEMP,WH2(1),WH2(2)
      RETURN
C
   60 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF H2: Q=,E12.4,10X,
     1 7HJ MAX =,I3,10X, 26HTHERMAL EQUILIBRIUM H2; T=,F10.1,  2H K,/)
   70 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF H2: Q=,E12.4,4X,
     16HJ MAX=,I3,4X, 25HN O R M A L  HYDROGEN, T=,F6.1,  2H K,4X,  6HW1
     2,W2=,2F6.2/)
C
      END

	subroutine BCparam(temp, G0, G1, G2, tau1, tau2)
	implicit double precision (a-h,o-z)
c	G0,G1,G2 are tree lowest translational moments	
      data HBAR/1.05458875D-27  /
      data BOLTZK/1.38054D-16/
	T = temp
      TAU0=HBAR/(2.*BOLTZK*T)
      TTA=G0*TAU0/G1
      TAU1=DSQRT((G2*TTA-G0*(1.+TAU0**2/TTA))/(G0*(TAU0/TTA)**2))
      TAU2=TTA/TAU1
C	write(3,10)  G0, G1, G2
10	format(' M0, M1, M2:', 3e12.4)
	return
	end

	subroutine K0param(temp, G0,G1,G2, tau1, tau2)
	implicit double precision (a-h,o-z)
c	G0,G1,G2 are tree lowest translational moments	
      data HBAR/1.05458875D-27  /
      data BOLTZK/1.38054D-16/
	T = temp
      TAU0=HBAR/(2.*BOLTZK*T)
      DELT=(TAU0*G1)**2 -4.*(G1*G1/G0+G1/TAU0-G2)*TAU0*TAU0*G0
      if(delt.le.0.d0) go to 88
      TAU1=(-DSQRT(DELT)-TAU0*G1)/(2.*(G1*G1/G0+G1/TAU0-G2))
      if(tau1.le.0.d0) go to 88
      TAU1=DSQRT(TAU1)
      TAU2=TAU0*TAU1*G0/(G1*TAU1*TAU1-TAU0*G0)
C	write(3,10)  G0, G1, G2
10	format(' M0, M1, M2:', 3e12.4)
	return
88    continue
C     write (6,177) delt, tau10
177   format(' Problem: one of the following is negative for K0, ',
     1  /,' delt, tau1:', 2e13.4)
      stop 9955   
	end

      FUNCTION BGAMA03(FNU,TAU5,TAU6,TEMP)
C     K0 LINE SHAPE MODEL
C     NORMALIZATION SUCH THAT 0TH MOMENT EQUALS UNITY
C     FNU IS THE FREQUENCY IN CM-1; TEMP IN KELVIN.
      implicit double precision (a-h,o-z)
      DATA TWOPIC,BKW,HBOK,PI/1.883651568d11,.6950304256d0,
     1 7.638280918d-12, 3.141592654d0/

      TAU4=dSQRT(TAU5*TAU5+(HBOK/(TEMP*2.d0))**2)
      OMEGA=TWOPIC*FNU
      XNOM=1.d0/(TAU6*TAU6)+OMEGA*OMEGA
      X=TAU4*DSQRT(XNOM)
      TAU56=TAU5/TAU6
      TAU56=DMIN1(TAU56,430.d0)
      BGAMA03=(TAU5/PI)*DEXP(TAU56+FNU/(2.d0*BKW*TEMP))*
     1 XK03(X)
      RETURN
      END

      FUNCTION BGAMA13(FNU,TAU1,TAU2,TEMP)
C     BIRNBAUM S CIA LINE SHAPE MODEL (K1)
C     NORMALIZATION SUCH THAT 0TH MOMENT EQUALS UNITY
C     FNU IS THE FREQUENCY IN CM-1; TEMP IN KELVIN.
       implicit double precision (a-h,o-z)
      DATA TWOPIC,BKW,HBOK,PI/1.883651568d11,.6950304256d0,
     1 7.638280918d-12, 3.141592654d0/
      TAU3=dSQRT(TAU2*TAU2+(HBOK/(TEMP*2.))**2)
      OMEGA=TWOPIC*FNU
      DENOM=1.d0+(OMEGA*TAU1)**2
      X=(TAU3/TAU1)*dSQRT(DENOM)
      AAA=TAU2/TAU1
      AAA=dMIN1(AAA,430.d0)
      BGAMA13=(TAU1/PI)*dEXP(AAA+FNU/(2.d0*BKW*TEMP))*
     1 XK13(X)/DENOM
      RETURN
      END

      FUNCTION BGAMA3 (FNU,T1,T2,TEMP,ifun)
      implicit double precision (a-h,o-z)
	if(ifun.eq.0) bgama3=bgama03(fnu,t1,t2,temp)
	if(ifun.eq.1) bgama3=bgama13(fnu,t1,t2,temp)
      RETURN
      END


      FUNCTION XK03(X)
C     MODIFIED BESSEL FUNCTION K0(X)
C     ABRAMOWITZ AND STEGUN P.379
       implicit double precision (a-h,o-z)
      IF(X-2.d0) 10,10,20
   10 T=(X/3.75d0)**2
      FI0=(((((.0045813*T+.0360768)*T+.2659732)*T
     1 +1.2067492)*T+3.0899424)*T+3.5156229)*T+1.
      T=(X/2.)**2
      P=(((((.00000740*T+.00010750)*T+.00262698)*T
     1 +.03488590)*T+.23069756)*T+.42278420)*T+
     2 (-.57721566)
      X=DABS(X)
      XK03=-DLOG(X/2.)*FI0+P
      RETURN
   20 T=(2./X)
      P=(((((.00053208*T-.00251540)*T+.00587872)*T
     1 -.01062446)*T+.02189568)*T-.07832358)*T+
     2 1.25331414
      X=DMIN1(X,330.d0)
      XK03=DEXP(-X)*P/DSQRT(X)
      RETURN
      END

 
      FUNCTION XK13(X)
C     MODIFIED BESSEL FUNCTION K1(X) TIMES X
C     PRECISION IS BETTER THAN 2.2e-7 EVERYWHERE.
C     ABRAMOWITZ AND S,TEGUN, P.379; TABLES P.417.
       implicit double precision (a-h,o-z)
      IF(X-2.) 10,10,20
   10 T=(X/3.75)**2
      FI1=X*((((((.00032411*T+.00301532)*T+.02658733)*T+.15084934)
     1 *T+.51498869)*T+.87890594)*T+.5)
      T=(X/2.)**2
      P=(((((-.00004686*T-.00110404)*T-.01919402)*T-.18156897)*T-
     1 .67278579)*T+.15443144)*T+1.
      XK13=X*dLOG(X/2)*FI1+P
      RETURN
   20  T=2./X
      P=(((((-.00068245*T+.00325614)*T-.00780353)*T+.01504268)*T-
     1 .03655620)*T+.23498619)*T+1.25331414
      X=dMIN1(X,330.d0)
      XK13=dSQRT(X)*dEXP(-X)*P
      RETURN
      END

