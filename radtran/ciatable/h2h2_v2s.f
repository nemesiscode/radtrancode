      subroutine h2h2_v2s(ntype1,temp1,freqlo1,freqmax,deltastep,list1,
     1freq1,alfatot1)
C     *************************************************************
C     Subroutine to calculate an H2-H2 collision induced absorption spectrum
C     in the v2 first overtone band.
C
C     Input variables:
C       ntype1          integer 0 = equilibrium hydrogen
C                               1 = 3:1 ortho/para hydrogen
C       temp1           double  Temperature (K)
C       freqlo1         double  Lowest wavenumber in spectrum (cm-1)
C       freqmax         double  Highest wavenumber in spectrum (cm-1)
C       deltastep       double  Wavenumber step (cm-1)
C
C     Output variables
C       list1           integer Number of points in spectrum
C       freq1(601)      double  Frequency grid (cm-1)
C       alfatot1(601)   double  Absorption coefficient (cm-1 amagat-2)
C
C     Pat Irwin         31/1/96
C     Pat Irwin		2/3/12	Updated for Radtrans2.0
C
C     *************************************************************

c      This program computes RV CIA for H2_H2 pairs
c        in the first overtone band
c      at temperatures between 20 and 500K.
C ============================================================================
C Copyright (C) 1994  Chunguang Zheng and Aleksandra Borysow
C ============================================================================
       IMPLICIT double precision (A-H,O-Z)
       CHARACTER*10  MARKER,lab(18),labeta
       dimension ibgama(18),LA1(18),LA2(18),LA(18),LL(18)
       dimension alfatot(0:18,601),fit(18,6,0:3)
       dimension alfatot1(601),freq1(601)
       common/ integ / rlo, rup, frac, Ninteg
       common/temp/temp
       COMMON /PRTSUM/ JRANGE1,MARKER
       common/prtsum1/ q1, w1(2)
       common /wagib/ w(2)
       common/ola/  xxm0
       common/sumrul/ g0, g1
       COMMON /SPECIT/ LIST
       common/ specit1/ FREQ(601), ABSCOEF(601)
       COMMON/SPEKTRM/LIKE,nvib_F(2,2),nvib_I,FREQLO,DFREQ,NTYPE
       common/ read/ go
       common/om/omega_0
       common/uncor/ xm0,tau1,tau2
       common/corr/ omega1,omega2,tauh
       COMMON/YMAX/YMAXXR,YMAXXL
       common /BB/ J1, Jp1, J2, Jp2, labeta

       DATA LIKE, nvib_F, nvib_I / 1, 0,1, 2,1, 0 /
c      nvib_F - upper state, nvib_I - lower state
       DATA W(1), W(2) /1., 3./   ! rotational weights;

       data nf/ 601 /
       data ibgama/1,1,0,0,0,1,1,1,1,1,1,0,0,0,1,1,1,1/
       data LA1/2,0,0,2,0,2,2,4,0,2,0,0,2,0,2,2,4,0/
       data LA2/0,2,0,0,2,2,2,0,4,0,2,0,0,2,2,2,0,4/
       data LA /2,2,0,2,2,3,1,4,4,2,2,0,2,2,3,1,4,4/
       data LL /3,3,1,1,1,3,1,5,5,3,3,1,1,1,3,1,5,5/
        DATA lab/'00|2023|02','00|0223|02','00|0001|02'
     &          ,'00|2021|02','00|0221|02','00|2233|02'
     &          ,'00|2211|02','00|4045|02','00|0445|02'
     &          ,'00|2023|11','00|0223|11','00|0001|11'
     &          ,'00|2021|11','00|0221|11','00|2233|11'
     &          ,'00|2211|11','00|4045|11','00|0445|11'/
        data ((fit( 1,i,j),j=0,3),i=1,6) /
     &         -63.3046,-2.90099,1.11488,-0.123286,
     &         -50.2794,-3.1198,1.20418,-0.125843,
     &         -36.1987,-4.3048,1.80452,-0.179853,
     &         1.23371,-1.36526,0.680679,-0.0627319,
     &         0.911984,0.389325,-0.430318,0.133539,
     &         -12.1644,-0.82163,0.208866,-0.0410492/
        data ((fit( 2,i,j),j=0,3),i=1,6) /
     &         -63.0276,-2.96191,1.18582,-0.150176,
     &         -50.2065,-3.10239,1.25046,-0.157528,
     &         -36.2885,-4.52031,2.05348,-0.248609,
     &         1.09598,-1.76608,0.988936,-0.131329,
     &         1.18926,-0.124895,-0.158453,0.0790217,
     &         -12.1563,-0.716906,0.144517,-0.0226413/
        data ((fit( 3,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit( 4,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit( 5,i,j),j=0,3),i=1,6) /
     &         -63.7765,-3.92644,1.62045,-0.178311,
     &         -50.4974,-3.907,1.61054,-0.177028,
     &         -36.4676,-4.38734,1.74076,-0.149647,
     &         1.4131,-0.625082,0.213482,0.0144708,
     &         0.955977,0.400642,-0.399686,0.123877,
     &         -12.2838,-0.691115,0.109554,-0.0194189/
        data ((fit( 6,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit( 7,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit( 8,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit( 9,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit(10,i,j),j=0,3),i=1,6) /
     &         -62.6761,-2.97135,1.1584,-0.131847,
     &         -49.7147,-3.07092,1.18876,-0.126675,
     &         -35.6231,-4.29454,1.80539,-0.182219,
     &         1.05338,-1.08379,0.517464,-0.0340486,
     &         0.771492,0.524109,-0.517489,0.154438,
     &         -12.1218,-0.868801,0.235473,-0.0464067/
        data ((fit(11,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit(12,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit(13,i,j),j=0,3),i=1,6) /
     &         -63.4485,-4.23085,1.79698,-0.209452,
     &         -50.2283,-4.04232,1.71151,-0.201378,
     &         -36.2593,-4.40465,1.76462,-0.159138,
     &         1.26484,-0.369052,0.0645006,0.0401403,
     &         0.853682,0.516928,-0.475046,0.140707,
     &         -12.2553,-0.729886,0.125514,-0.0208309/
        data ((fit(14,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit(15,i,j),j=0,3),i=1,6) /
     &         -63.8781,-2.83143,1.09649,-0.128031,
     &         -51.0455,-2.88536,1.094,-0.115844,
     &         -37.014,-4.24996,1.83075,-0.194018,
     &         0.893777,-1.13806,0.572189,-0.0444672,
     &         0.770947,0.551549,-0.552963,0.16236,
     &         -12.0842,-0.916309,0.274018,-0.0547391/
        data ((fit(16,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit(17,i,j),j=0,3),i=1,6) /24*0.0/

        data ((fit(18,i,j),j=0,3),i=1,6) /24*0.0/

        temp=temp1
        ntype=ntype1
        freqlo = freqlo1

c initialise the weights for H2
        W(1)=1.0
        W(2)=3.0
c this is needed otherwise they are not reset after a normal calculation


       if(temp.lt.20.d0 .or. temp.gt.500.d0) then
        print*,'h2h2_v2s: Warning'
        print*,'Temperature should be 20 < T < 500'
       endif

      


       list = 1+ INT(FREQMAX - FREQLO)/DELTASTEP
       list1=list
       if(list.gt.601) then
       write(*,*) 'too many steps.... Currently the max. number is: 601'
       stop 555
       endif

      DO 10 I=1,list
      FREQ(I)=FREQLO+DBLE(I-1)*DELTASTEP
      FREQ1(I)=FREQ(I)
      do j = 0,18
       alfatot(j,i)=0.
      end do
10    continue           !ALFATOT(I)=0.

C+++++++++++++       Compute the partition function Q1: +++++++++++++++++++

      CALL PARTFCT1 (TEMP, Q1, NTYPE, JRANGE1, MARKER, W )
       w1(1) = w(1)
       w1(2) = w(2)
        ALFA =1./(0.69519*TEMP)
        hbar = 1.054588757D-27
       boltzk = 1.38054d-16
       tau0 = hbar/(2.*boltzk*temp)
c       *****************************************************************
       do 11 n=1,18
        if (fit(n,1,0).eq.0.0) goto 11
       if (n.le.9) OMEGA_0 = H2ELEV1(2,0) - H2ELEV1(0,0)
       else OMEGA_0 = (H2ELEV1(1,0) - H2ELEV1(0,0))*2.0
        labeta=lab(n)

        xt=dlog10(temp)
        xm0=10.0**
     &  (fit(n,1,0)+fit(n,1,1)*xt+fit(n,1,2)*xt*xt+fit(n,1,3)*xt**3)
        xm1=10.0**
     &  (fit(n,2,0)+fit(n,2,1)*xt+fit(n,2,2)*xt*xt+fit(n,2,3)*xt**3)
        xm2=10.0**
     &  (fit(n,3,0)+fit(n,3,1)*xt+fit(n,3,2)*xt*xt+fit(n,3,3)*xt**3)
        omega1=10.0**
     &  (fit(n,4,0)+fit(n,4,1)*xt+fit(n,4,2)*xt*xt+fit(n,4,3)*xt**3)
        omega2=10.0**
     &  (fit(n,5,0)+fit(n,5,1)*xt+fit(n,5,2)*xt*xt+fit(n,5,3)*xt**3)
        TAUH=10.0**
     &  (fit(n,6,0)+fit(n,6,1)*xt+fit(n,6,2)*xt*xt+fit(n,6,3)*xt**3)
c       if(xm0.eq.0.0) goto 11

        if (ibgama(n).eq.0) then
        DELT=(TAU0*XM0)**2-4.*(XM1*XM1/XM0+XM1/TAU0-XM2)*TAU0*TAU0*XM0
        TAU1=(-DSQRT(DELT)-TAU0*xm1)/(2.*(XM1*XM1/XM0+XM1/TAU0-XM2))
       TAU1=DSQRT(TAU1)
       TAU2=TAU0*TAU1*xm0/(xm1*TAU1*TAU1-TAU0*xm0)
       else
       TTA=XM0*TAU0/XM1
       XXX = (XM2*TTA-XM0*(1.d0 + TAU0**2/TTA))/(XM0*(TAU0/TTA)**2)
C       print*,XXX
       TAU1=DSQRT(ABS(XXX))
C       print*,TTA,TAU1
       TAU2=TTA/TAU1
       end if

       call range1(ibgama(n))

        CALL ADDSPEC2(TEMP,LIKE, LA1(n), LA2(n), LA(n), LL(n), 
     1   freqcut, ibgama(n),XM0,tau1,tau2,omega1,omega2,tauH,
     2   1+int(n/9))
c       * * * * ** * * * * * * * * * * ** * * * * * * * * * * ** * * * * *

       do 198 l=1,list
        if (n.eq.15) ABSCOEF(l)=ABSCOEF(l)/2
        alfatot(0,l)=alfatot(0,l)+ABSCOEF(l)*2
        alfatot(n,l)=ABSCOEF(l)*2
198    continue
11     continue


       do iloop=1,list
        alfatot1(iloop)=alfatot(0,iloop)
       end do

       RETURN

       END
       
      SUBROUTINE ADDSPEC2( TEMP, LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE, 
     2  CUT, ibgama,  G0, tau1,tau2,  omega11, omega22, tauH ,m)
C

C     This subroutine generates a detailed listing of CIA ALFA(OMEGA)
c       IBGAMA=0 means uses K0, IBGAMA=1 (uses BC)
C     LIKE=1 for like systems (as H2-H2) 
C
      IMPLICIT double precision (A-H,O-Z)
      CHARACTER*10 MESSAGE, MARKER
      CHARACTER*10 labeta
      character*6 note
      COMMON /PRTSUM/ JRANGE1,MARKER
      COMMON/SPEKTRM/LKE,nvib_F(2,2),nvib_I,FREQLO,DFREQ,NTYPE
      common/sumrul/ gg0, gg1
      common/prtsum1/ q1, w1(2)
      double precision scalex,x
      COMMON /SPECIT/ LIST
      common/specit1/ FREQ(601),ABSCOEF(601)
      common/BB/ J1, Jp1, J2, Jp2, labeta
      common/ola/  XM00
      common/gam/ gamma(3)
      common/om/omeg_0
      common/ integ / rlo, rup, frac, ninteg
      COMMON/YMAX/YMAXXR, YMAXXL
      data nf /601 /

      DATA  SCALEF/1.D80/,  BOLTZWN/ 0.6950304D0/
      DATA HBAR,PI,CLIGHT/1.054588757D-27, 3.141592653589797D0,
     1   2.9979250D10/
      DATA Boltzk / 1.380658d-16 /
      DO 10 I=1, nf
   10 ABSCOEF(I)=0.D0
      TWOPIC=2.*PI*CLIGHT
      MESSAGE='NON-COMUL.'

      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*(2.68676D19)**2
      CALIB=CALIB/DBLE(1+LIKE)
      FREQCUT = CUT              ! do not alter the outside parameter CUT

      BETA1 = 1.D0/(BOLTZWN*TEMP)

        JPLUSL=JRANGE1+MAX0(LAMBDA1,LAMBDA2)
      DO 160 I1=1,JRANGE1
         J1=I1-1
      DO 160 IP1=1,JPLUSL
         JP1=IP1-1
         CG1S=CLEBSQR2(J1,LAMBDA1,JP1)
         IF (CG1S) 160,160,110
  110    P1=H2POPL1(nvib_I,J1,TEMP,W1)/Q1
         OMEGA1=H2ELEV1(nvib_F(m,1),JP1)-H2ELEV1(nvib_I,J1)
         DO 150 I2=1,JRANGE1
            J2=I2-1
         DO 150 IP2=1,JPLUSL
            JP2=IP2-1
            CG2S=CLEBSQR2(J2,LAMBDA2,JP2)
            IF (CG2S) 150,150,120
  120       P2=H2POPL1(nvib_I,J2,TEMP,W1)/Q1
            OMEGA2=H2ELEV1(nvib_F(m,2),JP2)-H2ELEV1(NVIB_i,J2)

       if(j1.eq.jp1) note = 'Single'
       if(j1.ne.jp1) note = 'Double'

        FAC=CALIB*P1*P2*CG1S*CG2S   
        X=0.0
        FAC=FAC*scalex(X)

            DO 140 I=1,LIST
       FRQ=FREQ(I)-OMEGA1-OMEGA2
       if (frq.gt.ymaxxR) go to 140
       if (frq.lt.ymaxxL) go to 140
  130   WKF=FREQ(I) * (1.-DEXP(- BETA1 * FREQ(I))) * FAC

       if (ibgama.eq.0) xx = g0 * bgama01(frq,tau1,tau2,temp)
       if (ibgama.eq.1) xx = g0 * bgama1a(frq,tau1,tau2,temp)
       yy = g0/2. * (bgama_h1(frq-(omega11+omega22),tauH,temp) -
     1     bgama_h1(frq-(omega22-omega11), tauH, temp))
       gg = xx + yy

        ABSCOEF(I) = ABSCOEF(I) + gg * WKF
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE

      RETURN
      END

      SUBROUTINE PARTFCT1 (TEMP, Q, NTYPE, JRANGE, MARKER, W)
C
C     NTYPE = 0 for equilibrium hydrogen 
C     NTYPE = 1 redefines the weights W for "normal" hydrogen with a
C               ratio of ortho:para H2 of 3:1. 
C
      IMPLICIT double precision (A-H,O-Z)
      CHARACTER*10  MARKER
      DIMENSION W(2)
      DATA  BOLTZWN /.6950304D0/
C
      Q=0.
      J=0
   30 DQ=H2POPL1(0,J,TEMP,W) + h2popl1(1,j,temp,w) + 
     1   h2popl1(2,j,temp,w) + h2popl1(3,j,temp,w) 

      Q = Q + DQ
      J = J + 1
      if ( (w(2).eq.0.).and.(mod(j,2).eq.1)) go to 30 ! odd J: do not compare
      IF (DQ.GT.Q/990.) GO TO 30

      J=-1
      S=0.

      IF (NTYPE.EQ.1) GO TO 60
   40 J=J+1
       DDD = H2POPL1(0, J, TEMP, W) + h2popl1(1,j,temp,w) + 
     1   h2popl1(2,j,temp,w) + h2popl1(3,j,temp,w) 

      S= S + DDD/Q
      IF (S.LE.0.999D0) GO TO 40
      JRANGE=J+2
      MARKER='EQUILIBRM '

        GO TO 90
C
C     "NORMAL" HYDROGEN REQUIRES REDEFINITION OF W(2):
C
   60 SEV=0.
   70 J=J+1
      DS=H2POPL1(0, J, TEMP, W)+ h2popl1(1,j,temp,w) + 
     1   h2popl1(2,j,temp,w) + h2popl1(3,j,temp,w) 

      S=S+DS
      SEV=SEV+DS
      J=J+1
      S=S+H2POPL1(0, J, TEMP, W)+ h2popl1(1,j,temp,w) + 
     1   h2popl1(2,j,temp,w) + h2popl1(3,j,temp,w) 

      IF (S.LE.0.999*Q) GO TO 70
      JRANGE=J+2
      SODD=S-SEV
      W(2)=W(2)*(3.*SEV)/SODD
C
C     DEFINITION OF "NORMAL" HYDROGEN: 3*S(EVEN) = S(ODD)
C
      Q=4.*SEV
      MARKER='NORMAL (!)'
   90  continue
      RETURN
      END


      FUNCTION H2POPL1 (NVIB, J, T, W)
C     UNNORMALIZED (!) Boltzmann factor of the roto- vibrational levels
C          of the H2 molecule.  (To be normalized by Q = partition sum)
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION W(2)
      DATA  BOLTZWN /.6950304D0/
C
20       DDD = H2ELEV1(nvib,J)
        H2POPL1 = DBLE(2*J+1)*W(1+MOD(J,2))*DEXP(-DDD/(BOLTZWN*T))
      RETURN
      END

      FUNCTION H2ELEV1 (NVIB, JROT)
C     Roto- vibrational energy levels of the H2 molecule
C     Needed for the computation of the partition sum Q of H2 ******
C
      IMPLICIT double precision (A-H,O-Z)
C
c       The following results are from LEVELS (levels.lev file) of
c       Kolos et all. (HSINGLX.FOR function), fitted as function of J(J+1)
c       Fitting program FITLEVELS.FOR (.for), below only for v=0, 1, 2
C       NOTE: ROTATIONAL CONSTANTS B AND D HAVE BEEN FIXED WHILE FITTING THE 
C       ROTATIONAL LEVELS, AFTER S.L.BRAGG, J.W.BRAULT,W.H.SMITH, AP.J.,VOL.263,
C       P.999-1004, (1982).

         EH2(I, A, B, D, H, XL, G, P, R)= A + B*DBLE(I) -
     1   D*1.D-2*DBLE(I)**2 + H * 1.D-5 * DBLE(I)**3 -
     1   XL*1.D-8*DBLE(I)**4 + G * 1.D-11*DBLE(I)**5 -
     1   P * 1.D-14 * DBLE(I)**6 + R*1.D-17*DBLE(I)**7
C       A - THE ROTATIONAL ENERGY (NVIB, J=0)

c       |E(v=0,J=0) - E(v=1,J=0)| = 4162.14 cm-1

       if(nvib.eq.0) h2elev1=eh2(jrot*(jrot+1), -0.361132D5, 
     1  59.33451D0, 4.56510D0, 4.568777D0, 4.366951D0, 2.863721D0, 
     1  1.051977D0, 0.156903D0)  + 36113.d0
c       add 36113.d0 to each energy to prevent h2elev1/(kT) being too
c       big at small temperatures
c       ======================================================       
       if(nvib.eq.1) h2elev1=eh2(jrot*(jrot+1), -0.319511D5, 
     1  56.3742D0, 4.4050D0, 4.515341D0, 4.614924D0, 3.301549D0, 
     1  1.32896D0, 0.212922D0)   + 36113.d0

       if(nvib.eq.2) h2elev1=eh2(jrot*(jrot+1),-0.280238D5, 53.482D0, 
     1  4.28D0, 4.552623D0, 4.957454D0, 3.695921D0, 1.469646D0, 
     1  0.199078D0)   + 36113.d0

       if(nvib.eq.3) h2elev1=eh2(jrot*(jrot+1), -.243268d+5, 
     1 0.5050d+02, 0.3818d+01, 0.2393d+01, -.8313d+00, -.4652d+01, 
     1 -.4660d+01, -.1628d+01) + 36113.d0

       if(nvib.eq.4) h2elev1=eh2(jrot*(jrot+1), -.208566d+05, 
     1  0.4762d+02, 0.3605d+01, 0.1511d+01, -.3982d+01, -.1106d+02, 
     1 -.1108d+02, -.4167d+01) + 36113.d0

       if(nvib.eq.5) h2elev1=eh2(jrot*(jrot+1), -.176133d+05, 
     1  0.4481d+02, 0.3511d+01, 0.1099d+01, -.6643d+01, -.1913d+02, 
     1 -.2174d+02, -.9433d+01) + 36113.d0

       if(nvib.gt.5) stop 888
C
      RETURN
      END


      FUNCTION XK0A(X)
C     MODIFIED BESSEL FUNCTION K0(X)
C     ABRAMOWITZ AND STEGUN P.379
       implicit double precision (a-h,o-z)
      IF(X-2.) 10,10,20
   10 T=(X/3.75d0)**2
      FI0=(((((.0045813*T+.0360768)*T+.2659732)*T
     1 +1.2067492)*T+3.0899424)*T+3.5156229)*T+1.
      T=(X/2.)**2
      P=(((((.00000740*T+.00010750)*T+.00262698)*T
     1 +.03488590)*T+.23069756)*T+.42278420)*T+
     2 (-.57721566)
      X=ABS(X)
      XK0A=-dLOG(X/2.)*FI0+P
      RETURN
   20 T=(2./X)
      P=(((((.00053208*T-.00251540)*T+.00587872)*T
     1 -.01062446)*T+.02189568)*T-.07832358)*T+
     2 1.25331414
      X=dMIN1(X,330.d0)
      XK0A=dEXP(-X)*P/dSQRT(X)
      RETURN
      END

      FUNCTION XK1A(X)
C     MODIFIED BESSEL FUNCTION K1(X) TIMES X
C     PRECISION IS BETTER THAN 2.2E-7 EVERYWHERE.
C     ABRAMOWITZ AND S,TEGUN, P.379; TABLES P.417.
       implicit double precision (a-h,o-z)
      IF(X-2.) 10,10,20
   10 T=(X/3.75)**2
      FI1=X*((((((.00032411*T+.00301532)*T+.02658733)*T+.15084934)
     1 *T+.51498869)*T+.87890594)*T+.5)
      T=(X/2.)**2
      P=(((((-.00004686*T-.00110404)*T-.01919402)*T-.18156897)*T-
     1 .67278579)*T+.15443144)*T+1.
      XK1A=X*dLOG(X/2)*FI1+P
      RETURN
   20  T=2./X
      P=(((((-.00068245*T+.00325614)*T-.00780353)*T+.01504268)*T-
     1 .03655620)*T+.23498619)*T+1.25331414
      X=dMIN1(X,330.d0)
      XK1A=dSQRT(X)*dEXP(-X)*P
      RETURN
      END


      FUNCTION BGAMA01(FNU,TAU1,TAU2,TEMP)
C     K0 LINE SHAPE MODEL
C     NORMALIZATION SUCH THAT 0TH MOMENT EQUALS UNITY
C     FNU IS THE FREQUENCY IN CM-1; TEMP IN KELVIN.
       implicit double precision (a-h,o-z)
      DATA TWOPIC,BKW,HBOK,PI/1.883651568d11,.6950304256d0,
     1 7.638280918d-12, 3.141592654d0/
       TAU0 = HBOK /(2.*TEMP)
      TAU4=DSQRT(TAU1*TAU1+(TAU0)**2)
      OMEGA=TWOPIC*FNU
      XNOM=1.d0/(TAU2*TAU2)+OMEGA*OMEGA
      X=TAU4*dSQRT(XNOM)
      TAU12=TAU1/TAU2
      TAU12=DMIN1(TAU12,430.d0)
      BGAMA01=(TAU1/PI)*DEXP(TAU12+FNU/(2.d0*BKW*TEMP))*
     1 XK0A(X)
       end


      FUNCTION CLEBSQR2 (L,LAMBDA,LP)
C     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)
C     FOR INTEGER ARGUMENTS ONLY
C     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.
C
      IMPLICIT double precision (A-H,O-Z)
      FC=DBLE(2*LP+1)
      GO TO 10
C
      ENTRY THREEJ2B
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
   10 CLEBSQR2=0.
      IF (((L+LAMBDA).LT.LP).OR.((LAMBDA+LP).LT.L).OR.((L+LP).LT.LAMBDA)
     1) RETURN
      IF (MOD(L+LP+LAMBDA,2).NE.0) RETURN
      IF ((L.LT.0).OR.(LP.LT.0).OR.(LAMBDA.LT.0)) RETURN
      F=1./DBLE(L+LP+1-LAMBDA)
      IF (LAMBDA.EQ.0) GO TO 30 
      I1=(L+LP+LAMBDA)/2        
      I0=(L+LP-LAMBDA)/2+1
      DO 20 I=I0,I1       
   20 F=F*DBLE(I)/DBLE(2*(2*I+1))
   30 P=FC*F*FCTL2(LAMBDA+L-LP)*FCTL2(LAMBDA+LP-L)
      CLEBSQR2=P/(FCTL2((LAMBDA+L-LP)/2)*FCTL2((LAMBDA+LP-L)/2))**2
      RETURN
      END

      FUNCTION FCTL2 (N)
C     Simple FACTORIALS program
      IMPLICIT double precision (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      P(Z)=((((-2.294720936D-4)/Z-(2.681327160D-3))/Z+(3.472222222D-3))/
     1Z+(8.3333333333D-2))/Z+1.
      FCTL2=1            
      IF (N.LE.1) RETURN
      IF (N.GT.15) GO TO 20
      J=1           
      DO 10 I=2,N   
   10 J=J*I         
      FCTL2=DBLE(J)
      RETURN        
   20 Z=DBLE(N+1) 
      FCTL2=(DEXP(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)
      RETURN
      END

      FUNCTION BGAMA1A(FNU,TAU1,TAU2,TEMP)
C     BIRNBAUM CIA LINE SHAPE MODEL (K1)
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
      BGAMA1A=(TAU1/PI)*dEXP(AAA+FNU/(2.d0*BKW*TEMP))*
     1 XK1A(X)/DENOM
      RETURN
      END

      FUNCTION BGAMA_H1(FNU,tau ,TEMP)
c       "K1" one parameter model lineshape (14-Jul-1989 )
c       with Egelstaff time; unique for "H" function
C     NORMALIZATION SUCH THAT 0TH MOMENT EQUALS UNITY
C     FNU IS THE FREQUENCY IN CM-1; TEMP IN KELVIN.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA TWOPIC,BKW,HBOK,PI/1.883651568E11,.6950304256,
     1 7.638280918D-12, 3.141592654/

      TAU3=DSQRT( 3./2.*TAU**2 + (HBOK/(TEMP*2.))**2 )
      OMEGA=TWOPIC*FNU
      ARG = DABS(OMEGA)*TAU3
       IF(ARG.EQ.0.D0) GO TO 5
      BGAMA_H1 =1./PI * DEXP(FNU/(2.*BKW*TEMP))* XK1A(ARG)
     1      * (TAU/TAU3)**2 * TAU * (3./2.)**(3./2.)
       RETURN
5       BGAMA_H1 =1./PI * DEXP(FNU/(2.*BKW*TEMP))
     1      * (TAU/TAU3)**2 * TAU * (3./2.)**(3./2.)
      RETURN
      END

       SUBROUTINE RANGE1(IBGAMA)
       implicit double precision (a-h,o-z)
      COMMON/YMAX/YMAXXR, YMAXXL
       COMMON/UNCOR/ G0, TAU1, TAU2
       COMMON/CORR/ omega1, omega2, tauH 
       common/temp/temp

c       choosing YMAXX here:

       if(ibgama.eq.1) gmax =g0 * bgama1a(0.d0,tau1,tau2,temp)
       if(ibgama.eq.0) gmax =g0 * bgama01(0.d0,tau1,tau2,temp)

       STEP = 100.D0
       XS = 100.D0
       omegapl = omega1 + omega2
       omegamn = omega2 - omega1

4441       if(ibgama.eq.1) gx = g0 * bgama1a(xs,tau1,tau2,temp)
       if(ibgama.eq.0) gx =g0 * bgama01(xs,tau1,tau2,temp)

       gy =  g0/2. * ( bgama_h1(xs- omegapl,tauH,temp) -
     1       bgama_h1(xs - omegamn, tauH, temp) )
       GLOWR = GX + GY
       XS = XS + STEP
       if(xs.gt.5000. ) go to 6565
       IF(GLOWR. GT. GMAX/5000.D0) GO TO 4441
6565       YMAXXR = XS-STEP
       if(ymaxxr.lt.1000.) ymaxxr = 2000.

       XS = -100.D0
       if(ibgama.eq.1) gx =g0 * bgama1a(0.d0,tau1,tau2,temp)
       if(ibgama.eq.0) gx =g0 * bgama01(0.d0,tau1,tau2,temp)

6441       XS = XS - STEP
       if(xs.le.-5000.) go to 7766
       if(ibgama.eq.1) gx =g0 * bgama1a(xs,tau1,tau2,temp)
       if(ibgama.eq.0) gx =g0 * bgama01(xs,tau1,tau2,temp)

       gy =  g0/2. * ( bgama_h1(xs- omegapl,tauH,temp) -
     1       bgama_h1(xs - omegamn, tauH, temp) )
       if((gx+gy).le.0.d0) go to 7766
       GLOWL = GX + GY
       IF(GLOWL. GT. GMAX/5000.D0) GO TO 6441
7766       YMAXXL = XS + STEP
       return
       end
        
        function scalex(X)
        IMPLICIT double precision (A-H,O-Z)
        external va1, beta_1
        CHARACTER*10 lab
        COMMON /BB/ J1, Jp1, J2, Jp2, LAB
        common /temp/temp
        DATA ANGAU, FOURPI, HBAR/ 0.52917706, 12.566370614359,
     &                            1.05458876D-34/
        T= TEMP * 1.380662D-23
        FMP=2.0*1.67265D-27
        F2=HBAR*HBAR/(12.*FMP*T*T*ANGAU**2*1.D-20)
        h=0.05
        bb=betacom1(0.0,0)
        g0=0.0d0
        g1=0.0d0
        do 100 i=60,220
        r=h*dble(i)
        dr=r/100.0d0
        CALL DERIVAS1(VA1,R,DR,V,DV,DDV,Z3,Z4,JFLAG)
        G = DEXP(-V/T)*R*R
        Gp=G*(1.+F2*((DV/T-4./R)*DV-2.*DDV))
        g0=g0+gp*( BETA_1(R,0) )**2
        g1=g1+gp*( BETA_1(R,1) )**2
 100    continue
        g0=FOURPI*G0*h*(ANGAU*1.D-8)**3
        g1=FOURPI*G1*h*(ANGAU*1.D-8)**3
        scalex=g1/g0
        return
        end

      FUNCTION VA1(R)
C       tested, OK  24-JUL-1989 18:02:16 
c       double potential; initial state: H2(v=0)-H2(v=0);
c       final state: H2(v=0) - H2(v=2);                17-JUL-1989 
C     ==============================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*10 NOTE(2), LGAS, LABEL
      COMMON /EPOT/ EPS(2),RV(2),FM,RMIN,RMAX,HQ
      COMMON /COMUNIC/ NOTE,LGAS,LABEL,FMP,EPSP(2),RM(2),LVALUE
      DIMENSION A(2), GAMA(2), ALPHA(2), D(2), RMIN1(2), E(2),
     1  C6(2), C8(2), C10(2), C6S(2), C8S(2), C10S(2)

        FACTOR = (4.803250D-10)**2/(0.52916607D-8) *(1.D-7)
        FMP=2.D0
        LGAS =  'H2-H2 GAS '

C       POTENTIAL H2(v=0) - H2(v=0)
        NOTE(1)='00 - 00 H2'
        RM(1)  = 6.5*0.52917706D0
       RMIN1(1)= 6.5D0
       E(1)=    0.11163D-3
        EPSP(1)= 0.11163D-3 * FACTOR
       A(1)=    0.37996D7
       GAMA(1)= 2.39588D0
       ALPHA(1)=14.84169D0
       D(1) =   1.02354D0
       C6(1)=  -12.134D0
       C8(1)=  -221.26D0
       C10(1)= -0.50358D4

C       POTENTIAL H2(v=2) - H2(v=0)
        NOTE(2)='02 - 02 H2'
        RM(2)  = 6.5*0.52917706D0
        EPSP(2)= 0.11658D-3 * FACTOR
       RMIN1(2)=6.5D0
       E(2)=   0.11658D-3
       A(2) = 0.38000D7
       GAMA(2)= 2.96022
       ALPHA(2)= 14.93242
       D(2)=      1.33860
       C6(2)=   -13.52281
       C8(2)=   -228.55879
       C10(2)=  -5525.0086

       DO 10 I=1,2
       C6S(I) = C6(I) /( E(I) * RMIN1(I)**6 )
       C8S(I) = C8(I) /( E(I) * RMIN1(I)**8 )
10       C10S(I) = C10(I) /( E(I) * RMIN1(I)**10 )
        
        K1=1
        X = R/6.5      !/RV(K1)

        V_REP = A(K1) * X**GAMA(K1)  * DEXP(- ALPHA(K1)*X)
        F = 1.D0
        IF(X .LE. D(K1) ) F = DEXP(- ( D(K1)/X-1.)**2)
        V_DISP = F * (C6S(K1)/X**6 + C8S(K1)/X**8 + C10S(K1)/X**10)
        FF = V_REP + V_DISP
        VA1 = FF * EPSp(K1)
        RETURN
        END

       FUNCTION BETA_1(R,kk)

CZ  ** Changed for scaling  Feb. 3, 1994 ****************

      IMPLICIT double precision (A-H,O-Z)
      CHARACTER*10 NOTE,LGAS,LABEL,LAB
      dimension N(18), BN(18), B0(18), A(18), B(18), LL(18),
     &                 BN1(18),B01(18),A1(18),B1(18),
     &                 BN2(18),B02(18),A2(18),B2(18),
     &                 BN3(18),B03(18),A3(18),B3(18),
     &                 BN4(18),B04(18),A4(18),B4(18)
      COMMON /BB/ J1, J1p, J2, J2p, LAB          ! , XM0
      COMMON /COMUNIC/ NOTE(2),LGAS,LABEL,FMP,EPSP(2),RMP(2),LVALUE
      DATA DEBYE / 2.541769713d-18 /
      DATA LL /3,3,1,1,1,3,1,5,5,3,3,1,1,1,3,1,5,5/
      DATA N  /4,4,7,0,0,4,0,6,6,4,4,0,7,7,4,0,6,6/

      DATA BN /-0.597d-1,0.949d-1,-0.540d1,0.0,0.0,
     &         0.102d-1,0.0,-0.0562,-0.102,
     &         0.114, -0.114, 0.0, -1.00, 1.00, -0.375d-1,
     &         0.0,0.36,-0.36/
      DATA B0 /-0.629d-6,-0.15d-4,0.662d-5,0.42d-5,0.156d-4,
     &         -0.327d-7,0.203d-7,-0.543d-6,-0.523d-5,
     &         0.635d-5, -0.635d-5, 0.0, -0.152d-4, 0.152d-4,
     &         0.479d-5,-0.143d-5,-0.68d-8,0.68d-8/
      DATA A  /-2.179,-1.415,-2.818,-1.800,-1.298,
     &         -3.122,-3.277,-1.755,-1.451,
     &         -1.081,-1.081,0.0,-1.498,-1.498, -1.203,
     &         -1.468,-6.462, -6.462/
      DATA B  /-0.373,0.004,-0.309,-0.08,-0.011,
     &         -0.254,-0.321,-0.175,-0.015,
     &         -0.399, -0.399, 0.0, -0.147, -0.147, -0.01,
     &         -0.081,-0.959, -0.959/

      DATA BN1 /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &          0.210d-2,-0.282d-2,0.0,0.0,0.0,
     &          5.54d-4,0.0,2.69d-3,-1.12d-2/
      DATA B01 /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &          -0.324d-7,-0.373d-6,0.0,-0.143d-6,0.695d-6,
     &          -0.388d-7,0.197d-7,-0.632d-9,0.0/
      DATA A1  /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &          -2.01,-1.486,0.0,-1.385,-1.568,
     &          -1.454,-1.417,-3.966,0.0/
      DATA B1  /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &          0.0,-0.029,0.0,-0.107,-0.026,
     &          -0.015,-0.031,-0.318,0.0/

      DATA BN2 /6.05d-5,-0.107d-2,0.0,0.0,0.0,
     &          -3.12d-4,0.0,-2.43d-3,-1.4d-2,
     &          0.282d-2,-0.210d-2,0.0,0.0,0.0,
     &          -5.54d-4,0.0,1.12d-2,-2.69d-3/
      DATA B02 /-0.113d-6,-0.417d-6,-2.25d-6,0.954d-7,0.7d-6,
     &          0.287d-7,-0.131d-7,0.0,0.0,
     &          0.373d-6,0.324d-7,0.0,-0.695d-6,0.143d-6,
     &          0.388d-7,-0.197d-7,0.0,0.632d-9/
      DATA A2  /-1.933,-1.515,-1.615,-2.106,-1.552,
     &          -1.448,-1.521,0.0,0.0,
     &          -1.486,-2.010,0.0,-1.568,-1.385,
     &          -1.454,-1.417,0.0,-3.966/
      DATA B2  /-0.095,-0.022,-0.079,-0.088,-0.031,
     &          -0.038,-0.096,0.0,0.0,
     &          -0.029,0.0,0.0,-0.026,-0.107,
     &          -0.015,-0.031,0.0,-0.318/

      DATA BN3 /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &          -2.06d-3,2.79d-3,0.0,0.0,0.0,
     &          -5.43d-4,0.0,-2.62d-3,1.11d-2/
      DATA B03 /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &          3.35d-8,3.71d-7,0.0,0.137d-6,-0.689d-6,
     &          0.377d-7,-0.191d-7,0.786d-9,0.0/
      DATA A3  /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &          -1.995,-1.488,0.0,-1.384,-1.570,
     &          -1.457,-1.426,-3.784,0.0/
      DATA B3  /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &          0.0,-0.029,0.0,-0.108,-0.026,
     &          -0.014,-0.033,-0.285,0.0/

      DATA BN4 /1.5d-4,8.62d-4,0.0,0.0,0.0,
     &          2.84d-4,0.0,2.79d-3,1.4d-2,
     &          -2.79d-3,2.06d-3,0.0,0.0,0.0,
     &          5.43d-4,0.0,-1.11d-2,2.62d-3/
      DATA B04 /1.25d-7,4.24d-7,2.3d-6,-0.113d-6,-0.704d-6,
     &          -0.277d-7,0.128d-7,0.0,0.0,
     &          -3.71d-7,-3.35d-8,0.0,0.689d-6,-0.137d-6,
     &          -0.377d-7,0.191d-7,0.0,-0.786d-9/
      DATA A4  /-1.89,-1.518,-1.61,-2.115,-1.548,
     &          -1.45,-1.509,0.0,0.0,
     &          -1.488,-1.995,0.0,-1.57,-1.384,
     &          -1.457,-1.426,0.0,-3.784/
      DATA B4  /-0.086,-0.022,-0.078,-0.105,-0.031,
     &          -0.043,-0.101,0.0,0.0,
     &          -0.029,0.0,0.0,-0.026,-0.108,
     &          -0.014,-0.033,0.0,-0.285/

      x=R-6.0
      betap=( (B01(i)*DEXP(A1(i)*X+B1(i)*X**2)+BN1(i)/R**N(i))*
     &         j1*(j1+1)
     &      + (B02(i)*DEXP(A2(i)*X+B2(i)*X**2)+BN2(i)/R**N(i))*
     &         j2*(j2+1)
     &      + (B03(i)*DEXP(A3(i)*X+B3(i)*X**2)+BN3(i)/R**N(i))*
     &         j1p*(j1p+1)
     &      + (B04(i)*DEXP(A4(i)*X+B4(i)*X**2)+BN4(i)/R**N(i))*
     &         j2p*(j2p+1)  )*debye
      beta_1=(B0(i)*DEXP(A(i)*X+B(i)*X**2)+BN(i)/R**N(i))*debye
      if (kk.eq.1) beta_1=beta_1+betap
      RETURN

      ENTRY BETACOM1               ! the same for all terms
      i = 0
      if (LAB .eq. '00|2023|02') I=1
      if (LAB .eq. '00|0223|02') I=2
      if (LAB .eq. '00|0001|02') I=3
      if (LAB .eq. '00|2021|02') I=4
      if (LAB .eq. '00|0221|02') I=5
      if (LAB .eq. '00|2233|02') I=6
      if (LAB .eq. '00|2211|02') I=7
      if (LAB .eq. '00|4045|02') I=8
      if (LAB .eq. '00|0445|02') I=9
      if (LAB .eq. '00|2023|11') I=10
      if (LAB .eq. '00|0223|11') I=11
      if (LAB .eq. '00|0001|11') I=12
      if (LAB .eq. '00|2021|11') I=13
      if (LAB .eq. '00|0221|11') I=14
      if (LAB .eq. '00|2233|11') I=15
      if (LAB .eq. '00|2211|11') I=16
      if (LAB .eq. '00|4045|11') I=17
      if (LAB .eq. '00|0445|11') I=18

      if ( i.eq.0 ) stop 99
      LVALUE=LL(i)
      BETA_1 = 0.D0
      RETURN
      END
      SUBROUTINE DERIVAS1(FCT,X0,DX1, Y,Y1,Y2,Y3,Y4,IFLAG)
C     COMPUTES OTH..4TH DERIVATIVES OF Y=FCT(X) AT X=X0 USING FIVE
C     ABSCISSAE AND CENTRAL DIFFERENCES. FCT MUST BE DECLARED EXTERNAL
C     IN CALLING ROUTINE. NOTE THAT FCT(X) MUST BE DEFINED FOR
C     X0-2DX<=X<=X0+2DX
C     A TEST OF THE DEFINITION OF THE FIRST DERIVATIVE: IFLAG=1 OUTPUT
C     SEE ABRAMOWITZ AND STEGUN, P. 914
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION  F(7)
9     X=X0 
      H=DX1 
      HSQ=H*H 
      IFLAG=0
      DO 10 I=1,7
10    F(I)=FCT(X-DBLE(4-I)*H)
      Y=F(4)
      A=F(1)+F(7)   
      A1=F(1)-F(7)
      B=F(2)+F(6)   
      B1=F(2)-F(6)
      C=F(3)+F(5)   
      C1=F(3)-F(5)

C     THESE ARE HIGHEST ORDER DERIVATIVES
      Y1 = (-A1+9.*B1-45.*C1)/(60.*H)
      Y2 = (A-13.5*B+135.*C-245.*Y)/(90.*HSQ)
      Y3 = (A1-8.*B1+13.*C1)/(8.*HSQ*H)
      Y4 = (-A+12.*B-39.*C+56.*Y)/(6.*HSQ*HSQ)
      RETURN
      END

