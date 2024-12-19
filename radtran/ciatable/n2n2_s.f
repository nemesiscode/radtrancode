      subroutine n2n2_s(temp,fnumin,fnumax,dnu,nf,freq,alfatot,slit1)

C       ****************************************************************
C     subroutine adapted from Borysow code to calculate CIA absorption
C     of N2-N2 collisions (for Titan atmosphere).
C
C     Input variables:
C       temp            double  Temperature (K)
C       fnumin          double  Lowest wavenumber in spectrum (cm-1)
C       fnumax          double  Highest wavenumber in spectrum (cm-1)
C       dnu             double  Wavenumber step (cm-1)
C       slit1           double  slit width in cm-1
C
C     Output variables
C       nf              integer Number of points in spectrum
C       freq(601)       double  Frequency grid (cm-1)
C       alfatot(601)    double  Absorption coefficient (cm-1 amagat-2)
C
C     C. Nixon 30-1-97
C	Nick Teanby		17/6/04	Bug Fix
c                                   -so that this subroutine and the 
c                                   orginal program by Borysow give the 
c                                   same output (alphatot)
c                                   
c                       *1.e80 change to 1.d80 in:
c                       RSILO(I)=dlog(RSI(I)*1.d80) for comatability with
c                       double precision variables
c                       
c                       *extra ,1 removed from the ldelvis(54) array 
c                       (according to conors/ a. borysow's paper)
c                       
c                       *nsri and wnrmax3 set back to values in the original
c                       program (this puts the dimer spikes in the right 
c                       positions and removes some spurious features in the
c                       absorbtion, like a spike at 135cm-1)
c
c		NB. I would not recommend using this with slit<4.3cm-1, as
c		some assumptions used in making the original code may be violated
c		(specificly the Dimer spikes are poorly constrained by lab data)
c                       
c	Pat Irwin	2/3/12	Updated for Radtrans2.0                       
c                       
CC       =========================================
C       PROGRAM PREPARED BY ALEKSANDRA BORYSOW (APRIL'1987)
C       (UNIVERSITY OF TEXAS AT AUSTIN, PHYSICS DEPARTMENT)
c	  Original version: written on Cyber    
c
C       PROGRAM GENERATES N2-N2 COLLISION-INDUCED SPECTRA AT
C       TEMPERATURES BETWEEN 50 TO 300 K.
C       CIA SPECTRA MODELED AFTER PAPER (*)
C       ALEKSANDRA BORYSOW AND LOTHAR FROMMHOLD,
C       ASTROPHYSICAL JOURNAL, VOL. 311, PAGES 1043-1057, (1986)
c
C	  Revised by Glenn Orton (1989) - to work on Sun workstations
c	  and on the VAX machines (Fortran-77)
c	  Passes standard test on Sun, at 200K (July 1992)
c	  Ready for the distribution.
C       ====================================================
c
c          Also in revision: double precision for all f.p. variables
c
      implicit double precision (a-h,o-z)
      character*5 lgas
c
c      COMMON/APP3/SLIT,DX,NSRI,WNRMAX3,NS,NSRIUP
      COMMON/APP3a/SLIT,DX,WNRMAX3
      COMMON/APP3b/NSRI,NS,NSRIUP
      COMMON/RSILO/RSILO(201)
      COMMON/BOU43/INITB
c     COMMON/BB/OMEG,RSI,RSIGG,NSOL,ALFA
      COMMON/BBa/OMEG,RSI,RSIGG,ALFA
      COMMON/BBc/NSOL
      COMMON/BF/G0BF,DELBF,OM0
      COMMON/STATT/ QW3(51,3),QW4(51,4)
      COMMON/LIKEB/LIKE,LGAS
      COMMON/K1K0/IK1K0
      COMMON/BBB/IBOUND
      COMMON/ENERG/EB(41,6), NIV(6)
      common/dimer/nlines
      DIMENSION FREQ(601),ABSCOEF(601),ALFATOT(601)
      DIMENSION RSI(201),RSIGG(201),TT(2),SS(1),OMEG(201)
      Y(X,A,B,C)=A*dexp((C*X+B)*X)
c
C     TEMP = TEMPERATURE IN KELVIN, SHOULD BE BETWEEN 50. AND 300.
C     FNUMIN = LOWEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU)
C     FNUMAX = HIGHEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU)
C       LINE SHAPE MODELLING WILL BE MOST ACCURATE WITHIN RANGE OF
C       R-T SPECTRAL INTENSITIES AS 1:100.
C     DNU = FREQUENCY INCREMENT IN CM-1. DNU SHOULD BE CHOSEN SO
C           THAT NOT MORE THAN 600 STEPS ARE NEEDED TO GO FROM
C           FNUMIN TO FNUMAX (ELSE ARRAY DIMENSIONS OF FREQ,ABSCOEF
C           MUST BE ADJUSTED IN ADDEM).
c
      DATA NIV/ 41, 32, 24, 17, 9, 2/
      DATA LGAS /'N2-N2'/

      if (temp .lt. 50 .or. temp .gt. 300) then
         print*, 'n2n2_s: Warning'
         print*, 'Temperature should be in range 50 < T < 300 K'
      end if

      LIKE=1
      slit = slit1 !(to avoid problems with variable declaration in
C                    common block conflicting with parameter passing!)

C I'm being a bit lazy here (partly because I'm not entirely sure
C how to fix it).
      if (slit .lt. 0.50) then
         print*, 'n2n2_s: Warning'
         print*, 'Slit size approaches minimum grid scale'
         print*, 'resolution of integration: 0.22 cm-1'
      end if

C
      NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1
C save the true value of nf for passing back 
      nf1 =nf
      IF (NF.GT.601) NF=601

      FNUMAX=FNUMIN+FLOAT(NF-1)*DNU
c
      CALL PARTSUM7(TEMP)
C
C     THE N2-N2 SPECTRA   FOR 50-300K
C     ================================
C
      X=dlog(TEMP)
      DO 10 I=1,NF
      FREQ(I)=FNUMIN+FLOAT(I-1)*DNU
      ALFATOT(I)=0.0
10    ABSCOEF(I)=0.
C    ===========================================
      data (eb(1,j),j=1,6) /-54.99996d0,-31.89437d0,-16.05019d0,
     $ -6.49343d0,-1.76583d0,-0.17133d0/
      data (eb(2,j),j=1,6) /-54.86228d0,-31.77215d0,-15.94640d0,
     $ -6.41131d0,-1.70887d0,-0.14341d0/
      data (eb(3,j),j=1,6) /-54.58697d0,-31.52779d0,-15.73896d0,
     $ -6.24732d0,-1.59552d0,0.d0/
      data (eb(4,j),j=1,6) /-54.17413d0,-31.16143d0,-15.42815d0,
     $ -6.00202d0,-1.42700d0,0.d0/
      data (eb(5,j),j=1,6) /-53.62391d0,-30.67334d0,-15.01440d0,
     $ -5.67623d0,-1.20523d0,0.d0/
      data (eb(6,j),j=1,6) /-52.93648d0,-30.06382d0,-14.49830d0,
     $ -5.27111d0,-0.93302d0,0.d0/
      data (eb(7,j),j=1,6) /-52.11211d0,-29.33328d0,-13.88057d0,
     $ -4.78813d0,-0.61434d0,0.d0/
      data (eb(8,j),j=1,6) /-51.15108d0,-28.48222d0,-13.16213d0,
     $ -4.22919d0,-0.25504d0,0.d0/
      data (eb(9,j),j=1,6) /-50.05374d0,-27.51123d0,-12.34407d0,
     $ -3.59665d0,0.13641d0,0.d0/
      data (eb(10,j),j=1,6) /-48.82049d0,-26.42099d0,-11.42771d0,
     $ -2.89345d0,2*0.d0/
      data (eb(11,j),j=1,6) /-47.45179d0,-25.21229d0,-10.41455d0,
     $ -2.12325d0,2*0.d0/
      data (eb(12,j),j=1,6) /-45.94815d0,-23.88603d0,-9.30639d0,
     $ -1.29074d0,2*0.d0/
      data (eb(13,j),j=1,6) /-44.31014d0,-22.44322d0,-8.10531d0,
     $ -0.40202d0,2*0.d0/
      data (eb(14,j),j=1,6) /-42.53841d0,-20.88502d0,-6.81376d0,
     $ 0.53450d0,2*0.d0/
      data (eb(15,j),j=1,6) /-40.63365d0,-19.21272d0,-5.43459d0,
     $ 1.50455d0,2*0.d0/
      data (eb(16,j),j=1,6) /-38.59665d0,-17.42777d0,-3.97121d0,
     $ 2.48212d0,2*0.d0/
      data (eb(17,j),j=1,6) /-36.42824d0,-15.53182d0,-2.42768d0,
     $ 3.46665d0,2*0.d0/
      data (eb(18,j),j=1,6) /-34.12937d0,-13.52669d0,-0.80899d0,
     $ 3*0.d0/
      data (eb(19,j),j=1,6) /-31.70105d0,-11.41446d0,0.87859d0,3*0.d0/
      data (eb(20,j),j=1,6) /-29.14439d0,-9.19750d0,2.62689d0,3*0.d0/
      data (eb(21,j),j=1,6) /-26.46061d0,-6.87848d0,4.42334d0,3*0.d0/
      data (eb(22,j),j=1,6) /-23.65103d0,-4.46049d0,6.24733d0,3*0.d0/
      data (eb(23,j),j=1,6) /-20.71709d0,-1.94714d0,8.06983d0,3*0.d0/
      data (eb(24,j),j=1,6) /-17.66041d0, 0.65736d0,9.90464d0,3*0.d0/
      data (eb(25,j),j=1,6) /-14.48271d0, 3.34788d0,4*0.d0/
      data (eb(26,j),j=1,6) /-11.18593d0, 6.11816d0,4*0.d0/
      data (eb(27,j),j=1,6) /-7.77221d0, 8.95978d0,4*0.d0/
      data (eb(28,j),j=1,6) /-4.24393d0,11.86130d0,4*0.d0/
      data (eb(29,j),j=1,6) /-0.60374d0,14.80383d0,4*0.d0/
      data (eb(30,j),j=1,6) /3.14531d0,17.75924d0,4*0.d0/
      data (eb(31,j),j=1,6) /6.99978d0,20.71774d0,4*0.d0/
      data (eb(32,j),j=1,6) /10.95566d0,23.71589d0,4*0.d0/
      data (eb(33,j),j=1,6) /15.00818d0,5*0.d0/
      data (eb(34,j),j=1,6) /19.15136d0,5*0.d0/
      data (eb(35,j),j=1,6) /23.37787d0,5*0.d0/
      data (eb(36,j),j=1,6) /27.67681d0,5*0.d0/
      data (eb(37,j),j=1,6) /32.03237d0,5*0.d0/
      data (eb(38,j),j=1,6) /36.42278d0,5*0.d0/
      data (eb(39,j),j=1,6) /40.83668d0,5*0.d0/
      data (eb(40,j),j=1,6) /45.29436d0,5*0.d0/
      data (eb(41,j),j=1,6) /49.79246d0,5*0.d0/
        JJ=1
 442    continue
 1023   FORMAT(5x,10f12.5)
        JJ=JJ+1
        IF(JJ.EQ.42) GO TO 444
        GO TO 442
C       EB(JJ,IV) JJ-ROTATIONAL LEVEL "L", IV- VIBRATIONAL LEVEL "V";
444     CONTINUE
c
C     =================================
C     QUADRUPOLAR INDUCTION: (50-300K)
C     =================================
      EPS=1.d-5
      TT(1)=10.d0
      CALL BBOUND32(TEMP,RSI,NSOL)
	ij = 0
      DO 88 I=1,NSOL
	ij = ij + 1
c	mod can be only 0 or 1 or 2
	if ( mod(ij,3).eq.0) rslow1 = 1.5d-60
	if ( mod(ij,3).eq.1) rslow1 = 1.7d-60
	if ( mod(ij,3).eq.2) rslow1 = 1.6d-60
      if (rsi(i).lt.1.d-60) rsi(i)=rslow1
      RSILO(I)=dlog(RSI(I)*1.d80)
      OMEG(I)=dfloat(i-1)*DX
88    continue
c
9991    FORMAT(1X, 10(10E12.5,/))
c
      CALL SPLINE7(NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)
c
        IK1K0=1
        IBOUND=1
C       B-C LINESHAPE HERE
C       THESE VALUES (S,T1,T2) REPLACE VALUES GIVEN IN PAPER (*):
C       PUBLISHED IN AN ERRATUM, ASTROPHYSICAL JOURNAL, VOL.320, P.437
C       (1987)
        S=Y(X,0.29723d+1, -0.99569d0, 0.09464d0)
        T1=Y(X,0.12962d-11, -0.13048d0, -0.03128d0)
        T2=Y(X,0.37969d-13, 1.03681d0, -0.14336d0)
        E=0.d0
        T3=0.d0
        T4=0.d0
c
      CALL ADDSPEC7(S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,2,0,2,3)
      DO 20 I=1,NF
c	quad(i) = abscoef(i)
20    ALFATOT(I)=ABSCOEF(I)+ALFATOT(I)
c
C     ===========================================================
C     PARAMETERS FOR 4045  (HEXADECAPOLE) COMPONENTS:
C     ================================================================
c      CALL BBOUND54(TEMP,RSI,NSOL)
      CALL BBOUND54(TEMP,RSI,NSOL)
	ij = 0
      DO 111 I=1,NSOL
	ij = ij + 1
c	mod can be only 0 or 1 or 2
	if ( mod(ij,3).eq.0) rslow1 = 1.5d-60
	if ( mod(ij,3).eq.1) rslow1 = 1.7d-60
	if ( mod(ij,3).eq.2) rslow1 = 1.6d-60

      if (rsi(i).lt.1.d-60) rsi(i)=rslow1
      RSILO(I)=dlog(RSI(I)*1.d80)
111   OMEG(I)=dfloat(I-1)*DX
      CALL SPLINE7(NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)
C     ==============================
C       TEMPERATURES 50-140K
C     ===========================
        IF(TEMP.GE.140.) GO TO 333
c
      S=Y(X,1.80926d-1,-1.69153d0,0.18605d0)
      E=Y(X,0.3d0, 0.d0,0.d0)
      T1=Y(X,0.66017d-15,2.59982d0,-0.31831d0)
      T2=Y(X,0.12481d-11,-0.57028d0, 0.05983d0)
      T3=Y(X,0.52681d-12,-0.24719d0,0.00519d0)
      T4=Y(X,0.27518d16,-25.38969d0,2.46542d0)
c
      IK1K0=0
      IBOUND=1
      CALL ADDSPEC7(S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,4,0,4,5)
      DO 50 I=1,NF
c	hexa(i) = abscoef(i)
50    ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)
      GO TO 334
C       ================================================
C       TEMPERATURES 140 K-300 K
C       ============================================
333   IK1K0=0
      IBOUND=1
cx    S=Y(X,8.25299d-62,-1.25562,0.12981)
      S=Y(X,8.25299d-2,-1.25562d0,0.12981d0)
      E=Y(X,0.3d0, 0.d0,0.d0)
      T1=Y(X,0.36611d-14,1.47688d0,-0.16537d0)
      T2=Y(X,0.61264d-9, -2.25011d0,0.15289d0)
      T3=Y(X,0.79820d-9,-2.76152d0, 0.21847d0)
      T4=Y(X,0.52868d-21,7.66253d0,-0.77527d0)
c
      CALL ADDSPEC7(S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,4,0,4,5)
      DO 550 I=1,NF
c	hexa(i) = abscoef(i)
550    ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)
c
C       ============================================================
C       DOUBLE TRANSITIONS
C       ==============================================================
C       TEMPERATURES 50-300K
C       ==========================================
334     CONTINUE
        IK1K0=1
        IBOUND=0
cx      S=Y(X,1.19261d-58, -3.78587,0.34024)
        S=Y(X,1.19261d+02, -3.78587d0,0.34024d0)
        T1=Y(X,0.93777d-11, -0.66548d0, 0.0033d0)
        T2=Y(X,0.30395d-12,0.24728d0, -0.06607d0)
        T3=0.
        T4=0.
      CALL ADDSPEC7(S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,2,2,3,3)
      DO 650 I=1,NF
650    ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)
c
C       ==========================================
C       1202 OVERLAP NEGLECTED (EXTREMELY SMALL)
C       =============================================
c
      return
      end


C ************************** SUBROUTINES ****************************
C *******************************************************************
      SUBROUTINE ADDSPEC7(G0,EP,TAU1,TAU2,TAU5,TAU6,TEMP,NF,FREQ,
     $ ABSCOEF,MP,LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE)
C
C     THIS PROGRAM GENERATES LISTING OF R-T CIA  ALFA(OMEGA)
C     IF EITHER LAMBDA1 OR LAMBDA2 EQUAL TO ZERO - SINGLE TRANSITIONS;
C     DOUBLE TRANSITIONS ARE ASSUMED OTHERWISE.
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.
C
      implicit double precision (a-h,o-z)
      COMMON/BBa/OMEG(201),RSI(201),RSIGG(201),BETA
      COMMON/BBc/NSOL
      COMMON/RSILO/RSILO(201)
      COMMON/APP3a/SLIT,DX,WNRMAX3
      COMMON/APP3b/NSRI,NS,NSRIUP
      COMMON/BBB/IBOUND
      COMMON/N2PARTC/Q1,WN2(2),B01,D01,JRANGE2
      COMMON/ENERG/EB(41,6), NIV(6)
      DIMENSION ABSCOEF(601),FREQ(601)
      DATA CLOSCHM,BOLTZWN/2.68675484d19,.6950304d0/
      DATA HBAR,PI,CLIGHT/1.054588757d-27,3.1415926535898d0,
     $ 2.9979250d10/
c
      EN2(I)=(B01-dfloat(I)*D01)*dfloat(I)
      PN2(J,T)=dfloat(2*J+1)*WN2(1+MOD(J,2))*dexp(-1.4387859/T*EN2(J*(J
     1 +1)))
      TWOPIC=2.*PI*CLIGHT
      IF(LIKE.NE.1)LIKE=0
c        Take care of factor of 1.e-60 here.
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*(1.d-30*CLOSCHM)**2
      CALIB=CALIB/FLOAT(1+LIKE)
      BETA=1.d0/(BOLTZWN*TEMP)
      LIST=NF
      DO 88 I=1,LIST
88    ABSCOEF(I)=0.0d0
C
      IF((LAMBDA1.EQ.0).OR.(LAMBDA2.EQ.0))GOTO 152
      JPLUSL=JRANGE2+MAX0(LAMBDA1,LAMBDA2)
c
      jsum=0
      DO 150 I1=1,JRANGE2
      J1=I1-1
      DO 150 IP1=1,JPLUSL
      JP1=IP1-1
      CG1S=CLEBSQR7(J1,LAMBDA1,JP1)
      IF(CG1S)150,150,130
130   P1=PN2(J1,TEMP)/Q1
      jsum=jsum+1
      OMEGA1=EN2(JP1*IP1)-EN2(J1*I1)
      DO 148 I2=1,JRANGE2
      J2=I2-1
      DO 148 IP2=1,JPLUSL
      JP2=IP2-1
      CG2S=CLEBSQR7(J2,LAMBDA2,JP2)
      IF(CG2S)148,148,132
132   P2=PN2(J2,TEMP)/Q1
      OMEGA2=EN2(JP2*IP2)-EN2(J2*I2)
      FAC=CALIB*P1*P2*CG1S*CG2S
      DO 146 I=1,LIST
      FRQ=FREQ(I)-OMEGA1-OMEGA2
      WKI=FREQ(I)*(1.-dexp(-BETA*FREQ(I)))
      WKF=WKI*FAC
      XBG=G0*BGAMA7(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)
      IF(IBOUND.EQ.0) GO TO 555
      if (dabs(frq).le.wnrmax3)
     $ XBG = XBG + SPECFCT7(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)
555   ABSCOEF(I)=ABSCOEF(I)+XBG*WKF
146   CONTINUE
148   CONTINUE
150   CONTINUE
c
        GO TO 2222
c
C     SINGLE TRANSITIONS AT NITROGEN'S ROTATIONAL FREQUENCIES
C     =======================================================
152   JPLUSL=JRANGE2+LAMBDA
      DO 200 I=1,JRANGE2
      J=I-1
      DO 200 IP=1,JPLUSL
      JP=IP-1
      CGS=CLEBSQR7(J,LAMBDA,JP)
      IF(CGS)200,200,210
210   P=PN2(J,TEMP)/Q1
      jsum=jsum+1
      OMEGA1=EN2(JP*IP)-EN2(J*I)
      FAC=CALIB*P*CGS
      DO 199 IQ=1,LIST
      FRQ=FREQ(IQ)-OMEGA1
cx    WKI=FREQ(IQ)*(1.-EXP(-BETA*FREQ(IQ)))
      WKI=FREQ(IQ)*(1.-dexp(-BETA*FREQ(IQ)))
      WKF=WKI*FAC
      XBG=G0*BGAMA7(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)
      IF(IBOUND.EQ.0) GO TO 444
      if (dabs(frq).le.wnrmax3)
     $ XBG = XBG + SPECFCT7(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)
444   ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF
199   CONTINUE
200   CONTINUE
2222  continue !PRINT 44,(ABSCOEF(I),I=1,LIST)
C44    FORMAT((1X,10E12.4,/))
      RETURN
      END
      SUBROUTINE PARTSUM7(TEMP)
C     N2 ROTATIONAL PARTITION SUM Q = Q(T).
c
      implicit double precision (a-h,o-z)
c
      COMMON/N2PARTC/Q,WN2(2),B0,D0,JRANGE1
      DATA B0,D0,WN2(1),WN2(2)/1.98957, 0.58d-5, 2., 1./
      EN2(I)=(B0-dfloat(I)*D0)*dfloat(I)
      PN2(J,T)=dfloat(2*J+1)*WN2(1+MOD(J,2))*dexp(-1.4387859*EN2
     1 (J*(J+1))/T)
c
C     Q,B0,D0,WN2 - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR N2
      Q=0.
      J=0
50    DQ=PN2(J,TEMP)
      Q=Q+DQ
      J=J+1
      IF(DQ.GT.Q/900.) GOTO 50
      JRANGE1=J
30    FORMAT(/' ROTATIONAL PARTITION FUNCTION OF N2: Q=',F8.2,
     1   10X,'J MAX =',I3/)
      RETURN
      END
c
      SUBROUTINE PROFILE7(X,Y)
C     A TRIANGULAR SLIT FUNCTION IS USED.
c
      implicit double precision (a-h,o-z)
c
      COMMON/BL3/RSI(401)
c     COMMON/APP3/SLIT,DX,NSRI,WNRMAX3,NS,NSRIUP
      COMMON/APP3a/SLIT,DX,WNRMAX3
      COMMON/APP3b/NSRI,NS,NSRIUP
      COMMON/BBBB/ IDELV,IV,IVP,IDELL,IL,ILP
c
      IF(Y)105,106,1
1     CONTINUE
      X0=(NSRI+1.)+X/DX
      NC=X0
      N1=NC+1
      SLOPE=Y/SLIT
      NU=X0-NS
      IF(NU.LT.1)NU=1
      IF(NU.GT.NSRIUP)RETURN
      NO=X0+NS
      IF(NO.GT.NSRIUP)NO=NSRIUP
      IF(NO.LT.1)RETURN
      IF(NC.GT.NSRIUP)NC=NSRIUP
      IF(NC.LE.1)GOTO 101
      DO 100 I=NU,NC
      XI=(I-1.)*DX-WNRMAX3
      DR=SLOPE*(XI-(X-SLIT))
      IF(DR.LE.0.)GOTO 100
      RSI(I)=RSI(I)+DR
  100 CONTINUE
  101 CONTINUE
c
        IF(NC.GE.NSRIUP)RETURN
      IF(N1.LT.1)N1=1
      DO 102 I=N1,NO
      XI=(I-1.)*DX-WNRMAX3
      DR=Y-SLOPE*(XI-X)
      IF(DR.LE.0.)GOTO 102
      RSI(I)=RSI(I)+DR
  102 CONTINUE
      RETURN
  105 continue !PRINT 10,SLIT
   10 FORMAT(/' A TRIANGULAR SLIT FUNCTION OF',F6.3,' CM-1 HALFWIDTH IS
     ' USED'/)
106   continue
      RETURN
      END
c
cx    FUNCTION SPECFCT7(FREQ,OMEGA,PHI,PHI2,N,RTEMP)
      double precision function specfct7(freq,omega,phi,phi2,n,rtemp)
C     THIS INTERPOLATES THE SPECTRAL FUNCTION PHI(FREQ) DEFINED AT
C     OMEGA(N) AS PHI(N). PHI2 IS THE SECOND DERIVATIVE AT OMEGA
C     WHICH MUST BE OBTAINED FIRST (USE SPLINE7 FOR THAT PURPOSE).
C     RTEMP IS THE RECIPROCAL TEMPERATURE IN CM-1 UNITS.
C     NOTE THAT WE INTERPOLATE 1.E80 TIMES THE LOGARITHM OF PHI(OMEGA)
c         Note that in GSO's revision, this factor is removed.
c           (revision modified)
c
      implicit double precision (a-h,o-z)
c
      DIMENSION PHI(N),PHI2(N),OMEGA(N)
c  ** need to dimension F and GP or compiler complains about rank mismatch
      DIMENSION F(2),GP(1)
      TFAC=0.
      F(1)=FREQ
      IF(F(1))10,20,20
10    F(1)=ABS(F(1))
      TFAC=(-RTEMP*F(1))
20    IF(F(1).LE.OMEGA(N))GOTO 30
      SPECFCT7=EXP(-(PHI(N-1)-PHI(N))*(F(1)-OMEGA(N))/
     $(OMEGA(N)-OMEGA(N-1))+PHI(N)+TFAC)*(1.d-80)
cx    SPECFCT7=dexp(-(PHI(N-1)-PHI(N))*(F-OMEGA(N))/
cx   $(OMEGA(N)-OMEGA(N-1))+PHI(N)+TFAC)
      RETURN
30    CALL IXPOLAT7(N,1,0,1.d-6,OMEGA,PHI,F,GP,SI,NR,PHI2)
      SPECFCT7=EXP(TFAC+GP(1))*(1.d-80)
cx    SPECFCT7=dexp(TFAC+GP)
      RETURN
      END

C *******************************************************************
      SUBROUTINE BBOUND32 (TEMP,RSI,NSOL)
C
      implicit double precision (a-h,o-z)
c
      COMMON/ENERG/ EB(41,6), NIV(6)
c     COMMON/APP3/SLIT,DX,NSRI,WNRMAX3,NS,NSRIUP
      COMMON/APP3a/SLIT,DX,WNRMAX3
      COMMON/APP3b/NSRI,NS,NSRIUP
      COMMON/BL3/RSIBB(401)
      common/dimer/nlines
c
        COMMON/BBBB/ LDELVI,IVI,IVIP,LDELEL,LL,LLP
        DIMENSION RSI(201)
c
c          Stored values
      dimension ldelvis(63),ivis(63),ivips(63),ldelels(63),as(63),bs(63)
      data ldelvis/0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $ 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,
     $ 4,4,4,4,4,4/
      data ivis/0,0,1,1,2,2,3,3,4,4,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,
     $ 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,1,2,0,1,2,0,1,0,1,0,1,2,0,1,
     $ 0,1,0,1,0,1/
      data ivips/0,0,1,1,2,2,3,3,4,4,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,
     $ 4,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,3,4,5,3,4,5,3,4,3,4,3,4,5,4,
     $ 5,4,5,4,5,4,5/
      data ldelels/1,3,1,3,1,3,1,3,1,3,-3,-1,1,3,-3,-1,1,3,-3,-1,1,3,3,
     $ 1,-1,-3,-3,-1,1,3,-3,-1,1,3,-3,-1,1,3,-1,-3,1,3,1,1,1,-1,-1,-1,
     $ -3,-3,3,3,-3,-3,-3,-3,-3,3,3,1,1,-1,-1/
      data as /0.44844d-39,0.44356d-39,0.29345d-39,0.28850d-39,
     $ 0.16441d-39,0.15899d-39,0.72882d-40,0.67748d-40,0.10378d-40,
     $ 0.13041d-41,0.15006d-40,0.15370d-40,0.16139d-40,0.17143d-40,
     $ 0.19985d-40,0.20169d-40,0.20994d-40,0.22094d-40,0.16360d-40,
     $ 0.16281d-40,0.16714d-40,0.17326d-40,0.80425d-41,0.80862d-41,
     $ 0.80093d-41,0.81366d-41,0.24471d-41,0.25406d-41,0.26629d-41,
     $ 0.28064d-41,0.46227d-41,0.47150d-41,0.48513d-41,0.50133d-41,
     $ 0.39968d-41,0.39840d-41,0.39810d-41,0.39687d-41,0.11806d-41,
     $ 0.13458d-41,0.38746d-41,0.39219d-41,0.73334d-42,0.13390d-41,
     $ 0.13041d-41,0.71401d-42,0.13461d-41,0.65776d-42,0.69847d-42,
     $ 0.13517d-41,0.75545d-42,0.13268d-41,0.69847d-42,0.13517d-41,
     $ 0.74640d-42,0.21322d-42,0.26037d-42,0.20823d-42,0.20632d-42,
     $ 0.21067d-42,0.20531d-42,0.21218d-42,0.23006d-42/
      data bs /0.000430d0,0.000460d0,0.000830d0,0.000890d0,0.001700d0,
     $ 0.001860d0,0.004100d0,0.004570d0,0.000000d0,0.000000d0,
     $ 0.000999d0,0.000523d0,0.000149d0,-0.000168d0,0.001837d0,
     $ 0.001153d0,0.000660d0,0.000254d0,0.003603d0,0.002677d0,
     $ 0.002101d0,0.001738d0,0.005950d0,0.006843d0,0.000000d0,
     $ 0.007035d0,0.001025d0,0.000642d0,0.000254d0,-0.000164d0,
     $ 0.002342d0,0.001975d0,0.001640d0,0.001328d0,0.004943d0,
     $ 0.004999d0,0.005461d0,0.006839d0,0.000000d0,0.010993d0,
     $ 0.000000d0,0.000000d0,0.001367d0,0.005262d0,0.000000d0,
     $ 0.001601d0,0.004510d0,0.000000d0,0.001828d0,0.004175d0,
     $ 0.048160d0,0.007033d0,0.001828d0,0.004175d0,0.009338d0,
     $ 0.003733d0,0.008508d0,0.006979d0,0.000000d0,0.005035d0,
     $ 0.000000d0,0.004169d0,0.000000d0/
c
        DATA TWOPIC/1.88365183d11/
        DATA  PI/ 3.141592654d0/
c
C     EB(I,K) - BOUND ENERGIES
C     AM =  (MTX.EL. (L,BETA,L') )**2
C     M(L,L',V,V') OF PAPER (*) ARE TO BE CORRECTED:
C     M(L,L',V,V') = AM * (2L+1) * C(L,3,L')**2
C     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN
C     WNRMAX3 - THE FREQUENCY RANGE OF B-B CONTRIBUTION
C     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH
C     B-B SPECTRUM, IN [CM-1].
C       A,B, COEFFICIENTS , EQ. 7, A.BORYSOW, L.FROMMHOLD,
C       AP. J. VOL.311, 1043-1057, (1986)
c

C this is the number of integration points (max)
c      NSRI=200
c changed back to original (NT)
	NSRI=190

C this is the limits of the wavenumber integration range
C calculated by adding the (variable) slit size to the maximum 
C Delta E of the transition.
c      WNRMAX3= 43.0 + SLIT
c changed back to original (NT)
      WNRMAX3= 45.

      NSRIUP=2*NSRI+1
      DX=WNRMAX3/dfloat(NSRI)
C note that there is a minimum value of DX = 43.0/200 = 2.15
      NS=INT(SLIT/DX)
c
      DO 300 I=1,401
300   RSIBB(I)=0.0
c
      ALFA=1./(0.69519*TEMP)
      RM= 13.90 * 1.672649d-24
C     RM - REDUCED MASS FOR N2-N2
c
      PF=((RM*(1.380662d-16)*TEMP*2.*PI/(6.626176d-27)**2)**1.5)
c       Scale down by 1.d60
      pf=1.d-60*pf
c
        NR=0
 555    continue
        NR=NR+1
        ldelvi=ldelvis(nr)
        ivi=ivis(nr)
        ivip=ivips(nr)
        ldelel=ldelels(nr)
        a=as(nr)
        b=bs(nr)
        IV=IVI+1
        IVP=IVIP+1
C       THESE ARE ENERGY COLUMNS (V, V')
c
        NNII=NIV(IV)
c
        DO 20 L=1,NNII
C       LOOP OVER INITIAL L-VALUES...
C
        AM = A*EXP( -B*dfloat( (L-1)*L ) )
C       L - NUMBER OF ROW, (L-1) - ROTATIONAL LEVEL
        LP=L+LDELEL
        LL=L-1
        LLP=LP-1
C       LL,LLP ARE INITIAL L AND FINAL L' ANGULAR MOMENTUM QUANTUM NUMBERS
c
        IF(LP.GT.NIV(IVP) .OR. LP.LT. 1) GO TO 20
        IF( EB(LP,IVP).EQ.0.0) GO TO 20
        IF( EB(L,IV).EQ.0.0) GO TO 20
        STOKE = EB(LP,IVP) - EB(L,IV)
c
        STOKI=AM * dexp(-ALFA*EB(L,IV))/PF*dfloat(2*LL+1)*
     1  CLEBSQR7(LL,3,LLP)
c
        CALL PROFILE7(STOKE,STOKI)
          if (stoki.gt.0.d0) nlines=nlines+1
c
        STOKIP=AM * dexp(-ALFA*EB(LP,IVP))/PF*dfloat(2*LLP+1)*
     1  CLEBSQR7(LLP,3,LL)
c
        CALL PROFILE7(-STOKE,STOKIP)
          if (stokip.gt.0.d0) nlines=nlines+1
c
20    CONTINUE
      IF(NR.EQ.63) GO TO 56
      GO TO 555
56    CONTINUE
C     32 ENTRIES FOR (3220)+(3202) IN TABLE 6, 63 IN ALL:
C     DATA EXPANDED AND INCLUDE NOW ALL POSSIBLE B-B TRANSITIONS
c
      DO 90 N=1,NSRIUP
90    RSIBB(N)=RSIBB(N)/TWOPIC/SLIT
c
      NSOL=NSRI+1
      NSOL2=NSOL+1
C     RSI - CONTRIBUTION FOR POSITIVE FREQUENCY SHIFTS
c
       DO 22 I=1,NSOL
22     RSI(I) = RSIBB(NSOL-1+I)
c
999     FORMAT(50(10E12.5,/))
c
      RETURN
      END
c
C ******************************************************************
      SUBROUTINE BBOUND54(TEMP,RSI,NSOL)
c
      implicit double precision (a-h,o-z)
c
      COMMON/ENERG/ EB(41,6), NIV(6)
      COMMON/APP3a/SLIT,DX,WNRMAX3
      COMMON/APP3b/NSRI,NS,NSRIUP
      COMMON/BL3/RSIBB(401)
      common/dimer/nlines
      DIMENSION RSI(201)
c          Stored values
      dimension ldelvis(54),ivis(54),ivips(54),ldelels(54),as(54),bs(54)
      data ldelvis/0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $ 1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
      data ivis/0,0,0,1,1,1,2,2,2,3,3,3,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,
     $ 2,2,2,3,3,3,3,3,3,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2/
      data ivips/0,0,0,1,1,1,2,2,2,3,3,3,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,
     $ 3,3,3,4,4,4,4,4,4,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4/
      data ldelels/1,3,5,1,3,5,1,3,5,1,3,5,-5,-3,-1,1,3,5,-5,-3,-1,1,
     $ 3,5,-5,-3,-1,1,3,5,-5,-3,-1,1,3,5,-5,-3,-1,1,3,5,-5,-3,-1,1,
     $ 3,5,-5,-3,-1,1,3,5/
      data as/0.79332d-41,0.78207d-41,0.77235d-41,0.45815d-41,
     $ 0.44834d-41,0.44059d-41,0.21730d-41,0.20824d-41,0.20250d-41,
     $ 0.77222d-42,0.70351d-42,0.66815d-42,0.49611d-42,0.52232d-42,
     $ 0.52979d-42,0.54652d-42,0.56827d-42,0.59277d-42,0.57330d-42,
     $ 0.60620d-42,0.60862d-42,0.62104d-42,0.63809d-42,0.65698d-42,
     $ 0.39501d-42,0.41599d-42,0.41033d-42,0.41097d-42,0.41339d-42,
     $ 0.41530d-42,0.15858d-42,0.15976d-42,0.15478d-42,0.15066d-42,
     $ 0.14554d-42,0.13848d-42,0.99241d-43,0.10109d-42,0.10396d-42,
     $ 0.10758d-42,0.11176d-42,0.11636d-42,0.16460d-42,0.16470d-42,
     $ 0.16617d-42,0.16837d-42,0.17085d-42,0.17327d-42,0.11797d-42,
     $ 0.11593d-42,0.11405d-42,0.11174d-42,0.10853d-42,0.10401d-42/
      data bs/0.000612d0,0.000635d0,0.000677d0,0.001137d0,0.001201d0,
     $ 0.001341d0,0.002290d0,0.002449d0,0.002870d0,0.005426d0,
     $ 0.005876d0,0.007450d0,0.001000d0,0.000883d0,0.000609d0,
     $ 0.000392d0,0.000207d0,0.000037d0,0.001625d0,0.001624d0,
     $ 0.001305d0,0.001084d0,0.000927d0,0.000821d0,0.002978d0,
     $ 0.003273d0,0.002994d0,0.002954d0,0.003153d0,0.003668d0,
     $ 0.005799d0,0.006423d0,0.006733d0,0.007960d0,0.010937d0,
     $ 0.019179d0,0.001229d0,0.000993d0,0.000767d0,0.000543d0,
     $ 0.000309d0,0.000051d0,0.002456d0,0.002300d0,0.002210d0,
     $ 0.002193d0,0.002273d0,0.002506d0,0.004556d0,0.004825d0,
     $ 0.005454d0,0.006725d0,0.009431d0,0.016672d0/
      DATA TWOPIC/1.88365183d11/
      DATA  PI/ 3.141592654/
c
C the number of integration points (200-max)
c      NSRI=200
c changed back to original (NT)
      NSRI=190

C the integration range
c      WNRMAX3=44.0 + SLIT
c changed back to original (NT)
	WNRMAX3=47.
      
      NSRIUP=2*NSRI+1
      DX=WNRMAX3/dfloat(NSRI)
C minimum DX is 44/200 = 0.22

      NS=INT(SLIT/DX)
c
      DO 300 I=1,401
300   RSIBB(I)=0.0
      ALFA=1./(0.69519*TEMP)
      RM= 13.90 * 1.672649d-24
C     RM - REDUCED MASS FOR N2-N2
c
      PF=((RM*(1.380662d-16)*TEMP*2.*PI/(6.626176d-27)**2)**1.5)
c       Scale downward by 1.d60
      pf=1.d-60*pf
c
        NR=0
 555    continue
        NR=NR+1
        ldelvi=ldelvis(nr)
        ivi=ivis(nr)
        ivip=ivips(nr)
        ldelel=ldelels(nr)
        a=as(nr)
        b=bs(nr)

        IV=IVI+1
        IVP=IVIP+1
        NNII=NIV(IV)
        DO 20 L=1,NNII
        AM = A*EXP( -B*dfloat( L*(L+1) ) )
        LP=L+LDELEL
        LL=L-1
        LLP=LP-1
        IF(LP.GT.NIV(IVP) .OR. LP.LT. 1) GO TO 20
        IF( EB(LP,IVP).EQ.0.0) GO TO 20
        IF( EB(L,IV).EQ.0.0) GO TO 20
c
        STOKE = EB(LP,IVP) - EB(L,IV)
        STOKI= AM * dexp(-ALFA*EB(L,IV))/PF*dfloat(2*LL+1)*
     1  CLEBSQR7(LL,5,LLP)
        CALL PROFILE7( STOKE,STOKI)
          if (stoki.gt.0.d0) nlines=nlines+1
        STOKIP= AM * dexp(-ALFA*EB(LP,IVP))/PF*dfloat(2*LLP+1)*
     1  CLEBSQR7(LLP,5,LL)
        CALL PROFILE7(-STOKE,STOKIP)
          if (stokip.gt.0.d0) nlines=nlines+1
20      CONTINUE
        IF(NR.EQ.54) GO TO 56
c
        GO TO 555
56      CONTINUE
C       54 ENTRIES FOR (5440)=(5404) IN TABLE 6
c
      DO 90 N=1,NSRIUP
90    RSIBB(N)=RSIBB(N)/TWOPIC/SLIT
c
      NSOL=NSRI+1
      NSOL2=NSOL+1
c
      DO 22 I=1,NSOL
22    RSI(I)=RSIBB(NSOL-1+I)
c
999     FORMAT(50(10E12.5,/))
c
      RETURN
      END
      DOUBLE PRECISION FUNCTION CLEBSQR7(L,LAMBDA,LP)
C     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)
C     FOR INTEGER ARGUMENTS ONLY
C     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.
c
      implicit double precision (a-h,o-z)
c
      FC=dfloat(2*LP+1)
      GOTO 2
C
      ENTRY THREEJ27
C     THIS ENTRY RETURNS THE SQUARED 3-J SYMBOL   L LAMBDA LP
C                                                 0    0    0
C     INSTEAD OF THE CLEBSCH-GORDAN COEFFICIENT
C     (LIMITATION TO INTEGER ARGUMENTS ONLY)
C
C     NOTE THAT THE THREE-J SYMBOLS ARE COMPLETELY SYMMETRIC IN THE
C     ARGUMENTS. IT WOULD BE ADVANTAGEOUS TO REORDER THE INPUT ARGUMENT
C     LIST SO THAT LAMBDA BECOMES THE SMALLEST OF THE 3 ARGUMENTS.
      FC=1.
2     CLEBSQR7=0.
      IF(((L+LAMBDA).LT.LP).OR.((LAMBDA+LP).LT.L).OR.((L+LP).LT.
     $LAMBDA))RETURN
      IF(MOD(L+LP+LAMBDA,2).NE.0)RETURN
      IF((L.LT.0).OR.(LP.LT.0).OR.(LAMBDA.LT.0))RETURN
      F=1./dfloat(L+LP+1-LAMBDA)
      IF(LAMBDA.EQ.0)GOTO 22
      I1=(L+LP+LAMBDA)/2
      I0=(L+LP-LAMBDA)/2+1
      DO 20 I=I0,I1
20    F=F*dfloat(I)/dfloat(2*(2*I+1))
22    P=FC*F*FCTL7(LAMBDA+L-LP)*FCTL7(LAMBDA+LP-L)
      CLEBSQR7=P/(FCTL7((LAMBDA+L-LP)/2)*FCTL7((LAMBDA+LP-L)/2))**2
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCTL7(N)
c
      implicit double precision (a-h,o-z)
c
      P(Z)=((((-2.294720936d-4)/Z-(2.681327160d-3))/Z+(3.472222222d-3))/
     $Z+(8.333333333d-2))/Z+1.
      FCTL7=1.
      IF(N.LE.1)RETURN
      IF(N.GT.15)GOTO 20
      J=1
      DO 10 I=2,N
10    J=J*I
      FCTL7=dfloat(J)
      RETURN
20    Z=dfloat(N+1)
      FCTL7=(dexp(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)
      RETURN
      END
      double precision FUNCTION BGAMA7(FNU,T1,T2,EPS,T3,T4,TEMP)
C     SPECTRAL FUNCTION "EBC", FOR REFERENCE:
C     SEE "PHENOMENA INDUCED BY INTERMOLECULAR INTERACTIONS",
C     ED. G. BIRNBAUM; J. BORYSOW AND L. FROMMHOLD, P.67, (1985)
C       ============================================
c
      implicit double precision (a-h,o-z)
      double precision K0
c
        COMMON/K1K0/IK1K0
C       IF IK1K0=1 ONLY B-C; EBC OTHERWISE
      DATA PI,CLIGHT/3.1415926535898,2.99792458E10/
      DATA HBAR,BOLTZ/1.0545887d-27,1.380662d-16/
      P1(X)=((((((.0045813d0*X+.0360768d0)*X+.2659732d0)*X
     1 +1.2067492d0)*X+3.0899424d0)*X+3.5156229d0)*X+1.d0)
      P2(X)=((((((.00000740d0*X+.00010750d0)*X+.00262698d0)*X
     2 +.03488590d0)*X+.23069756d0)*X+.42278420d0)*X-.57721566d0)
      P3(X)=((((((.00032411d0*X+.00301532d0)*X+.02658733d0)*X
     3 +.15084934d0)*X+.51498869d0)*X+.87890594d0)*X+.5)
      P4(X)=((((((-.00004686d0*X-.00110404d0)*X-.01919402d0)*X
     4 -.18156897d0)*X-.67278579d0)*X+.15443144d0)*X+1.d0)
      P5(X)=((((((.00053208d0*X-.00251540d0)*X+.00587872d0)*X
     5 -.01062446d0)*X+.02189568d0)*X-.07832358d0)*X+1.25331414d0)
      P6(X)=((((((-.00068245d0*X+.00325614d0)*X-.00780353d0)*X
     6 +.01504268d0)*X-.03655620d0)*X+.23498619d0)*X+1.25331414d0)
C
      OMEGA=2.d0*PI*CLIGHT*FNU
      T0=HBAR/(2.d0*BOLTZ*TEMP)
      Z=dsqrt((1.d0+(OMEGA*T1)**2)*(T2*T2+T0*T0))/T1
      IF(Z-2.d0)10,10,12
10    XK1=Z*Z*dlog(Z/2.d0)*P3((Z/3.75d0)**2)+P4((Z/2.d0)**2)
      GOTO 20
12    XK1=dsqrt(Z)*dexp(-Z)*P6(2.d0/Z)
20    CONTINUE
      BGAMBC=(T1/PI)*dexp(T2/T1+T0*OMEGA)*XK1/(1.d0+(T1*OMEGA)**2)
      IF (IK1K0.EQ.1) GO TO 55
      ZP=SQRT((1.d0+(OMEGA*T4)**2)*(T3*T3+T0*T0))/T4
      IF(ZP-2.d0)22,22,24
22    K0=-dlog(ZP/2.d0)*P1((ZP/3.75d0)**2)+P2((ZP/2.d0)**2)
      GO TO 30
24    K0=dexp(-ZP)*P5(2.d0/ZP)/dsqrt(ZP)
30    BGAMA7 =(BGAMBC+EPS*(T3/PI)*dexp(T3/T4+T0*OMEGA)*K0)/(1.d0+EPS)
      GO TO 66
55    BGAMA7=BGAMBC
66    continue
      RETURN
      END
      SUBROUTINE SPLINE7(L,M,K,EPS,X,Y,T,SS,SI,NR,S2)
c
      implicit double precision (a-h,o-z)
c
C     SPLINE INTERPOLATION AND QUADRATURE, THIRD ORDER AFTER GREVILLE.
C     INPUT ARGUMENTS L...Y, OUTPUT SS...NR.
C     L DATA POINTS X(1), Y(1) ... X(L),Y(L)
C     EPS=ERROR CRITERION, TYPICALLY EPS=1.d-5 FOR 5 DECI. PLACES ACCURA
C     M ARGUMENTS T(1)..T(M) FOR WHICH FUNCTION VALUES SS(1)..SS(M), FOR
C     K=0; OR FIRST OR SECOND DERIVATIVE FOR K=1 OR -1, RESPECTIVELY.
C     NOTE THAT M HAS TO BE AT LEAST EQUAL TO 1.
C     SI=INTEGRAL (OVER WHOLE INTERVAL) FOR K=2 ONLY.
C     FOR 'NATURAL' SPLINE FUNCTIONS, S2(1)=S2(L)=0. MUST BE INPUT*NOTE*
C     N0 INDICATES THE NUMBER OF OUT-OF-RANGE CALLS. X(1)<T(I)<X(L)
C     EXTRAPOLATE WITH CAUTION. (ASSUMPTION D2Y/DX2 = 0.)
C     S2(I) IS THE 2ND DERIVATIVE AT X=X(I) AND IS COMPUTED INTERNALLY.
c     DIMENSION X(L),Y(L),T(M),SS(M),S2(L)
      DIMENSION X(201),Y(201),T(2),SS(1),S2(201)
      DELY(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
      B(I)=(X(I)-X(I-1))*0.5d0/(X(I+1)-X(I-1))
      C(I)=3.d0*(DELY(I)-DELY(I-1))/(X(I+1)-X(I-1))
c
      N=L
      N1=N-1
      NR=0
4     DO 52 I=2,N1
52    S2(I)=C(I)/1.5d0
      OMEGA=1.0717968d0
      IC=0
C     'NATURAL' SPLINE FUNCTIONS OF THIRD ORDER.
      S2(1)=0.d0
      S2(N)=0.d0
5     ETA=0.d0
      IC=IC+1
      SM=dabs(S2(1))
      DO 25 I=2,N
      IF (dabs(S2(I)).GT.SM) SM=dabs(S2(I))
25    CONTINUE
      EPSI=EPS*SM
6     DO 10 I=2,N1
7     W=(C(I)-B(I)*S2(I-1)-(0.5d0-B(I))*S2(I+1)-S2(I))*OMEGA
8     IF(dabs(W)-ETA)10,10,9
9     ETA=dabs(W)
10    S2(I)=S2(I)+W
13    IF(ETA-EPSI)14,5,5
c
c     ENTRY IXPOLAT7
      entry ixpolat7(L,M,K,EPS,X,Y,T,SS,SI,NR,S2)
C     THIS ENTRY USEFUL WHEN ITERATION PREVIOUSLY COMPLETED
      N=L
      N1=N-1
      NR=0
      IC=-1
14    IF(K-2) 15,20,15
15    DO 61 J=1,M
16    I=1
54    IF(T(J)-X(1))58,17,55
55    IF(T(J)-X(N))57,59,158
56    IF(T(J)-X(I))60,17,57
57    I=I+1
      GOTO 56
58    NR=NR+1
      HT1=T(J)-X(1)
      HT2=T(J)-X(2)
      YP1=DELY(1)+(X(1)-X(2))*(2.*S2(1)+S2(2))/6.d0
      IF(K)72,70,71
71    SS(J)=YP1+HT1*S2(1)
      GOTO 61
70    SS(J)=Y(1)+YP1*HT1+S2(1)*HT1*HT1/2.d0
      GOTO 61
72    SS(J)=S2(I)
      GOTO 61
158   HT2=T(J)-X(N)
      HT1=T(J)-X(N1)
      NR=NR+1
      YPN=DELY(N1)+(X(N)-X(N1))*(S2(N1)+2.*S2(N))/6.d0
      IF(K)82,80,81
81    SS(J)=YPN+HT2*S2(N)
      GOTO 61
80    SS(J)=Y(N)+YPN*HT2+S2(N)*HT2*HT2/2.d0
      GOTO 61
82    SS(J)=S2(N)
      GOTO 61
59    I=N
60    I=I-1
c
17    HT1=T(J)-X(I)
      HT2=T(J)-X(I+1)
      PROD=HT1*HT2
      S3=(S2(I+1)-S2(I))/(X(I+1)-X(I))
      SS2=S2(I)+HT1*S3
      DELSQS=(S2(I)+S2(I+1)+SS2)/6.d0
      IF(K)43,41,42
41    SS(J)=Y(I)+HT1*DELY(I)+PROD*DELSQS
      GOTO 61
42    SS(J)=DELY(I)+(HT1+HT2)*DELSQS+PROD*S3/6.d0
      GOTO 61
43    SS(J)=SS2
61    CONTINUE
20    SI=0.
      DO 62 I=1,N1
      H=X(I+1)-X(I)
62    SI=SI+0.5d0*H*(Y(I)+Y(I+1))-H**3*(S2(I)+S2(I+1))/24.d0
      IF(K.EQ.2)NR=IC
      RETURN
      END
