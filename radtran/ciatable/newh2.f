      subroutine h2h2_v0s(norm,temp,fnumin,fnumax,dnu,nf1,
     1freq1,alfatot)
cC     *************************************************************
C     Subroutine to calculate an H2-H2 collision induced absorption spectrum
C     in the v0 rotational-translational band.
C
C     Input variables:
C	norm		integer	0 = equilibrium hydrogen
C				1 = 3:1 ortho/para hydrogen
C	temp		double	Temperature (K)
C	fnumin		double	Lowest wavenumber in spectrum (cm-1)
C	fnumax		double	Highest wavenumber in spectrum (cm-1)
C	dnu		double	Wavenumber step (cm-1)
C
C     Output variables
C	nf1		integer	Number of points in spectrum
C	freq1(601)	double	Frequency grid (cm-1)
C	alfatot(601)	double	Absorption coefficient (cm-1 amagat-2)
C
C     Pat Irwin		31/1/96
C     Pat Irwin		2/3/12	Updated for Nemesis2.0
C
C ============================================================================
C Copyright (C) 1991  Aleksandra Borysow and Lothar Frommhold
C ============================================================================
C 
C Copyright Notice:
C You may use this program  for your scienftific applications,
C but please do not  distribute it yourself.
C Direct all your inquires to the  author: e-mail aborysow@phy.mtu.edu
C Final request: If you publish  your work which benefited from this
C program, please acknowledge using this program and QUOTE 
C the original paper describing the procedure: 
c	 The Astrophysical Journal
c      vol. 296, p. 644--654, (1985) by J. Borysow, L. Trafton,
C      L. Frommhold and G. Birnbaum.
C ============================================================================

	
C
c	redone from h2h2RT.cyber on 27 Feb. 1991
c	checked with exp/model data - works excellent!
C
C     THIS PROGRAM GENERATES THE H2-H2 TRANSLATIONAL/ROTATIONL
C     CIA SPECTRA. IT IS BASED ON QUANTUM LINE SHAPE COMPUTATIONS AND
C     THE AB INITIO DIPOLE DATA BE W. MEYER. DIMER FINE STRUCTURES ARE
C     SUPPRESSED. THE PROGRAM WAS WRITTEN BY ALEKSANDRA BORYSOW AND
C     LOTHAR FROMMHOLD. THIS IS THE NOVEMBER 1985 VERSION
C
C     H2-H2 COMPUTATIONS REFERENCE: MEYER, FROMMHOLD AND BIRNBAUM,
C     TO BE PUBLISHED IN PHYS.REV.A IN 1985;
C     THE H2-H2 MODELING USED HERE IS BEING PUBLISHED: J.BORYSOW,
C     L.TRAFTON, L.FROMMHOLD, G.BIRNBAUM, AP.J. (1985)
C
C     TAPE3 IS OUTPUT: HEADER PLUS ABSORPTION COEFF. ALPHA(NU)
C     h2part use to be auxil(5)
	implicit double precision (a-h,o-z)
	dimension alfatot(601),FREQ1(601)
      COMMON /H2PART1/  Q,WH2(2),B0,D0
      common /h2part2/ jrange1,NORMAL      
      COMMON /RESULT1/ NF
      common /result2/ FREQ(601),ABSCOEF(601)
      Y(X,A,B,C)=A*DEXP((C*X+B)*X)
C
      NORMAL=NORM

      IF(TEMP.LT.40.OR.TEMP.GT.300)THEN
       PRINT*,'h2h2_s: Warning this is a new routine'
       PRINT*,'Temperature should be 40 < T < 300'
      END IF

      NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1
      IF (NF.GT.601) NF=601
      NF1=NF
      print*, temp,fnumax,fnumin,nf,dnu
        CALL PARTSUM (TEMP)
C
C     THE H2-H2 SPECTRA
C     =================
C
      X=DLOG(TEMP)
      DO 10 I=1,NF
       FREQ(I)=FNUMIN+DFLOAT(I-1)*DNU
        FREQ1(i)=FREQ(I)      
	alfatot(i)=0.d0
   10 ABSCOEF(I)=0.
C
C     THE LAMBDA1,LAMBDA2,LAMBDA,L = 2023 AND 0223 COMPONENTS:
C     (QUADRUPOLE PLUS OVERLAP INDUCED COMPONENT)
C	sum of 2023+0223:
      S=Y(X,2.881d-61,-1.1005d0,0.1310d0)
      E=Y(X,7.3485d0,-1.3874d0,0.1660d0)
      T1=Y(X,7.883d-13,-.3652d0,-.0271d0)
      T2=Y(X,3.803d-13,-.4048d0,-.0091d0)
      T3=Y(X,1.0922d-12,-.4810d0,-.0127d0)
      T4=Y(X,5.174d-12,-.9841d0,0.0483d0)
	write(*,*)x,S,E,T1,t2,T3,T4
      CALL  ADDSPEC (S,E,T1,T2,T3,T4,TEMP,0,1,0,2,2,3,0,0,1.)
      do 1000 i=1,nf
 1000    print*,'quadrupole',i,freq(i),abscoef(i), alfatot(i)
	do 111 i=1, nf
111	alfatot(i) = alfatot(i) + abscoef(i)

C     PARAMETERS FOR 4045 AND 0445 (PURE HEXADECAPOLE) COMPONENTS
      S=Y(X,2.404d-64,-1.4161d0,0.1847d0)
      E=Y(X,-.8033d0,-.4474d0,-.0235d0)
      T1=Y(X,3.873d-13,-.4226d0,-.0183d0)
      T2=Y(X,2.743d-13,-.3566d0,-.0140d0)
      T3=Y(X,4.171d-13,-.5223d0,0.0097d0)
      T4=Y(X,2.2725d-12,-1.1056d0,0.0139d0)
	write(*,*)x,S,E,T1,t2,T3,T4
      CALL ADDSPEC(S,E,T1,T2,T3,T4,TEMP,0,1,4,0,4,5,0,0,1.)
      do 2000 i=1,nf
 2000    print*,'hexadecapole',i,freq(i),abscoef(i),alfatot(i)
	do 211 i=1, nf
211	alfatot(i) = alfatot(i) + abscoef(i)

C     PARAMETERS FOR 0221 AND 2021 (PURE OVERLAP) COMPONENTS
      S=Y(X,6.393d-63,-1.5964d0,0.2359d0)
      E=Y(X,21.414d0,-1.2511d0,0.1178d0)
      T1=Y(X,1.876d-13,-.4615d0,-.0012d0)
      T2=Y(X,4.839d-13,-.5158d0,0.0075d0)
      T3=Y(X,4.550d-13,-.5507d0,0.0095d0)
      T4=Y(X,2.045d-12,-.5266d0,-.0240d0)
	write(*,*)x,S,E,T1,t2,T3,T4
      CALL  ADDSPEC(S,E,T1,T2,T3,T4,TEMP,0,1,0,2,2,1,0,0,1.)
      do 3000 i=1,nf
 3000    print*,'overlap',i,freq(i),abscoef(i),alfatot(i)
	do 311 i=1, nf
311	alfatot(i) = alfatot(i) + abscoef(i)
C
C     PARAMETERS FOR 2233 QUADRUPOLE INDUCED DOUBLE TRANSITIONS
      S=Y(X,5.965d-63,-1.0394d0,0.1184d0)
      E=Y(X,6.674d0    ,-.9459d0,0.1124d0)
      T1=Y(X,4.764d-13,-.1725d0,-.0450d0)
      T2=Y(X,4.016d-13,-.3802d0,-.0134d0)
      T3=Y(X,1.0752d-12,-.4617d0,-.0085d0)
      T4=Y(X,1.1405d-11,-1.2991d0,0.0729d0)
	write(*,*)x,S,E,T1,t2,T3,T4
      CALL ADDSPEC(S,E,T1,T2,T3,T4,TEMP,0,1,2,2,3,3,0,0,1.)
      do 4000 i=1,nf
 4000    print*,'double t',i,freq(i),abscoef(i),alfatot(i)
	do 411 i=1, nf
411	alfatot(i) = alfatot(i) + abscoef(i)
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
	stop
      END


      SUBROUTINE ADDSPEC (G0,EP,TAU1,TAU2,TAU5,TAU6,TEMP,MP,LIKE,LAMBDA1
     1,LAMBDA2,LAMBDA,LVALUE,NVIB1,NVIB2,FACTOR)
C
C     THIS PROGRAM GENERATES A LISTING OF THE CIA TR ALFA(OMEGA)
C     IF BOTH LAMBDA1 AND LAMBDA2 ARE NEGATIVE: SINGLE TRANSITIONS;
C     DOUBLE TRANSITIONS ARE ASSUMED OTHERWISE.
C     MP=1 GIVES LISTINGS OF INTERMEDIATE RESULTS.
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.
C
	implicit double precision (a-h,o-z)
      COMMON /H2PART1/ Q,WH2(2),B0,D0
      common /h2part2/ JRANGE1,NORMAL
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
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2
      CALIB=CALIB/DFLOAT(1+LIKE)
      BETA=1./(BOLTZWN*TEMP)
      LIST=NF
C
C     ROTATIONAL SPECTRUM FOR THE DETAILED LISTING   *******************
      print*,'tau1',tau1
      print*,'tau2',tau2
      print*,'tau5',tau5
      print*,'tau6',tau6
      print*,'lambdas',lambda1,lambda2,lambda
      print*,'lval go ep',lvalue,go,ep
      print*,'temp',temp

C
      IF ((LAMBDA1.LT.0).AND.(LAMBDA2.LT.0)) GO TO 60
      JPLUSL=JRANGE1+MAX0(LAMBDA1,LAMBDA2)
      DO 50 I1=1,JRANGE1
         J1=I1-1
      DO 50 IP1=1,JPLUSL
         JP1=IP1-1
         CG1S=CLEBSQR(J1,LAMBDA1,JP1)
         IF (CG1S) 50,50,10
   10    P1=PH2(J1,TEMP)/Q
         IF (P1.LT.0.001) GO TO 50
         OMEGA1=EH2(NVIB1,JP1*IP1)-EH2(0,J1*I1)
         DO 40 I2=1,JRANGE1
            J2=I2-1
         DO 40 IP2=1,JPLUSL
            JP2=IP2-1
            CG2S=CLEBSQR(J2,LAMBDA2,JP2)
            IF (CG2S) 40,40,20
   20       P2=PH2(J2,TEMP)/Q
            IF (P2.LT.0.001) GO TO 40
            OMEGA2=EH2(NVIB2,JP2*IP2)-EH2(0,J2*I2)
            FAC=CALIB*P1*P2*CG1S*CG2S
            DO 30 I=1,LIST
               FRQ=FREQ(I)-OMEGA1-OMEGA2
               WKI=FREQ(I)*(1.-DEXP(-BETA*FREQ(I)))
               WKF=WKI*FAC
               XBG=G0*BGAMA(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)
               ABSCOEF(I)=ABSCOEF(I)+XBG*WKF
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      GO TO 100
   60 JPLUSL=JRANGE1+LAMBDA
      DO 90 I=1,JRANGE1
         J=I-1
      DO 90 IP=1,JPLUSL
         JP=IP-1
         CGS=CLEBSQR(J,LAMBDA,JP)
         IF (CGS) 90,90,70
   70    P=PH2(J,TEMP)/Q
         IF (P.LT.0.001) GO TO 90
         OMEGA1=EH2(NVIB1,JP*IP)-EH2(0,J*I)
         FAC=CALIB*P*CGS
         DO 80 IQ=1,LIST
            FRQ=FREQ(IQ)-OMEGA1
            WKI=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))
            WKF=WKI*FAC
            XBG=G0*BGAMA(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF
   80    CONTINUE
   90 CONTINUE
  100 CONTINUE
      RETURN
C
  110 FORMAT (/, 32H LAMBDA1,LAMBDA2, LAMBDA,LVALUE=,2I3,2X,2I3, 20H COM
     1PONENT INCLUDED.,/15X, 22HLINE SHAPE PARAMETERS:,6E12.3,5X/ 
     25HG(0)=,E12.3/)
  140 FORMAT (/,' ABSORPTION COEFFICIENT ALPHA(fnu)',/,
     3  ' FROM',F8.1,' CM-1 TO',F8.1,' CM-1, AT',F6.2,
     3  ' CM-1 INCREMENTS',' IN UNITS OF CM-1 AMAGAT-2',/)
  150 FORMAT (  6E12.4)
C
      END
      FUNCTION CLEBSQR (L,LAMBDA,LP)
C
C     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)
C     FOR INTEGER ARGUMENTS ONLY
C     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.
C
	implicit double precision (a-h,o-z)
      FC=DFLOAT(2*LP+1)
      GO TO 10
C
      ENTRY THREEJ2
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
   10 CLEBSQR=0.
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
   30 P=FC*F*FCTL(LAMBDA+L-LP)*FCTL(LAMBDA+LP-L)
      CLEBSQR=P/(FCTL((LAMBDA+L-LP)/2)*FCTL((LAMBDA+LP-L)/2))**2
      RETURN
C
      END
      FUNCTION FCTL (N)
	implicit double precision (a-h,o-z)
      P(Z)=((((-2.294720936d-4)/Z-(2.681327160d-3))/Z+(3.472222222d-3))/
     1Z+(8.333333333d-2))/Z+1.
      FCTL=1.
      IF (N.LE.1) RETURN
      IF (N.GT.15) GO TO 20
      J=1
      DO 10 I=2,N
   10 J=J*I
      FCTL=DFLOAT(J)
      RETURN
   20 Z=DFLOAT(N+1)
      FCTL=(DEXP(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)
      RETURN
      END

      FUNCTION BGAMA (FNU,T1,T2,EPS,T3,T4,TEMP)
C
C     EQUATION 13, SO-CALLED EBC MODEL, OF BORYSOW,TRAFTON,FROMMHOLD,
C     AND BIRNBAUM, ASTROPHYS.J., TO BE PUBLISHED (1985)
C     NOTE THAT BGAMA REDUCES TO THE BC PROFILE FOR EPS=0.
C
	implicit double precision (a-h,o-z)
      REAL K0
      DATA PI,CLIGHT/3.1415926535898,2.99792458E10/
      DATA HBAR,BOLTZ/1.0545887d-27,1.380662d-16/
      P1(X)=((((((.0045813*X+.0360768)*X+.2659732)*X+1.2067492)*X+3.0899
     1424)*X+3.5156229)*X+1.)
      P2(X)=((((((.00000740*X+.00010750)*X+.00262698)*X+.03488590)*X+.23
     1069756)*X+.42278420)*X-.57721566)
      P3(X)=((((((.00032411*X+.00301532)*X+.02658733)*X+.15084934)*X+.51
     1498869)*X+.87890594)*X+.5)
      P4(X)=((((((-.00004686*X-.00110404)*X-.01919402)*X-.18156897)*X-.6
     17278579)*X+.15443144)*X+1.)
      P5(X)=((((((.00053208*X-.00251540)*X+.00587872)*X-.01062446)*X+.02
     1189568)*X-.07832358)*X+1.25331414)
      P6(X)=((((((-.00068245*X+.00325614)*X-.00780353)*X+.01504268)*X-.0
     13655620)*X+.23498619)*X+1.25331414)
C
      OMEGA=2.*PI*CLIGHT*FNU
      T0=HBAR/(2.*BOLTZ*TEMP)
      Z=SQRT((1.+(OMEGA*T1)**2)*(T2*T2+T0*T0))/T1
      IF (Z-2.) 10,10,20
   10 XK1=Z*Z*DLOG(Z/2.)*P3((Z/3.75)**2)+P4((Z/2.)**2)
      GO TO 30
   20 XK1=SQRT(Z)*DEXP(-Z)*P6(2./Z)
   30 IF (EPS.EQ.0.) GO TO 70
      ZP=SQRT((1.+(OMEGA*T4)**2)*(T3*T3+T0*T0))/T4
      IF (ZP-2.) 40,40,50
   40 K0=-DLOG(ZP/2.)*P1((ZP/3.75)**2)+P2((ZP/2.)**2)
      GO TO 60
   50 K0=DEXP(-ZP)*P5(2./ZP)/SQRT(ZP)
   60 BGAMA=((T1/PI)*DEXP(T2/T1+T0*OMEGA)*
     2 XK1/(1.+(T1*OMEGA)**2)+EPS*(T3/
     1PI)*DEXP(T3/T4+T0*OMEGA)*K0)/(1.+EPS)
      RETURN
   70 BGAMA=(T1/PI)*DEXP(T2/T1+T0*OMEGA)*XK1/(1.+(OMEGA*T1)**2)
      RETURN
C
      END
      SUBROUTINE PARTSUM (TEMP)
C
C     H2 ROTATIONAL PARTITION SUM Q = Q(T).
C
	implicit double precision (a-h,o-z)
      COMMON /H2PART1/ Q,WH2(2),B0,D0
      common /h2part2/ JRANGE1,NORMAL
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
   20 WRITE (3,60) Q,JRANGE1,TEMP
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
      WRITE (3,70) Q,JRANGE1,TEMP,WH2(1),WH2(2)
      RETURN
C
   60 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF H2: Q=,E12.4,10X,
     1 7HJ MAX =,I3,10X, 26HTHERMAL EQUILIBRIUM H2; T=,F10.1,  2H K,/)
   70 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF H2: Q=,E12.4,4X,
     16HJ MAX=,I3,4X, 25HN O R M A L  HYDROGEN, T=,F6.1,  2H K,4X,  6HW1
     2,W2=,2F6.2/)
C
      END
