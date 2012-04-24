      subroutine h2he_v1s(ntype1,temp1,freqlo1,freqmax,dfreq1,nf1,
     1freq1,alfatot1)
C     *************************************************************
C     Subroutine to calculate an H2-He collision induced absorption spectrum
C     in the v1 fundamental band.
C
C     Input variables:
C       ntype1          integer 0 = equilibrium hydrogen
C                               1 = 3:1 ortho/para hydrogen
C       temp1           double  Temperature (K)
C       freqlo1         double  Lowest wavenumber in spectrum (cm-1)
C       freqmax         double  Highest wavenumber in spectrum (cm-1)
C       dfreq1          double  Wavenumber step (cm-1)
C
C     Output variables
C       nf1             integer Number of points in spectrum
C       freq1(601)      double  Frequency grid (cm-1)
C       alfatot1(601)   double  Absorption coefficient (cm-1 amagat-2)
C
C     Pat Irwin         31/1/96
C     Pat Irwin         2/3/12	Updated for Radtrans2.0
C
C     *************************************************************

 
C     To add several contributions of spectral functions
C     and to generate  detailed listings of the roto-vibrational spectra
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*10 MESSAGE, MARKER
 
      COMMON /PRTSUM14/ Q1,W1(2),Q2,W2(2)
      common /prtsum24/ JRANGE1,JRANGE2,MARKER
      COMMON /SPECIT14/ LIST
      dimension alfatot1(601),freq1(601)
      common /specit24/ FREQ(601),ABSCOEF(601),ALFATOT(601)
      COMMON/SPEKTRM14/ NF,LIKE, NVIB1,NVIB2 
      COMMON /SPEKTRM24/ FREQLO, DFREQ,fmax
      COMMON/TEMP4/TEMP

      DATA LIKE,NVIB1,NVIB2 / 0, 1, 0/
      DATA FMAX / 3500.d0/
 
c       Data for "01" term from fitting the Q.M. lineshapes
 
c        t1 "X"    r(%)  1.563
                data A01A, A01B, A01C, A01D, A01E /0.33040D+03,
     1  -.21649D+01, 0.52749D+00, -.67195D-01, 0.30651D-02/
c        t2 "X"    r(%)  3.364
                data B01A, B01B, B01C, B01D, B01E /0.54319D+01,
     1   0.20547D+01,-.84063D+00, 0.10838D+00, -.48703D-02/
c        G0 total  r(%)  0.000
                data C01A, C01B, C01C, C01D, C01E /0.83228D+02,
     1  -.15051D+01, 0.37391D+00, -.17139D-01, -.93260D-04/
c        epsilon   r(%)  0.000
                data D01A, D01B, D01C, D01D, D01E /0.40478D-02,
     1   0.18367D+01, -.42343D+00, 0.41328D-01, -.16901D-02/
c        t1 "Y"    r(%)  0.641
                data E01A, E01B, E01C, E01D, E01E / 0.49231D-03,
     1   0.56766D+01, -.11506D+01, 0.83468D-01, -.18655D-02/
c        t2 "Y"    r(%)  0.000
                data F01A, F01B, F01C, F01D, F01E /0.18547D+04,
     1   -.12944D+01, -.20196D+00, 0.60123D-01, -.35471D-02/
c       ================================================
c       These are results of fitted "23" lineshape parameters;
c
c        t1 "X"    r(%)  0.501
        data A23A, A23B, A23C, A23D, A23E /0.12235D+03,
     1  -0.61810, 0.18320D-01, -0.61509D-02, 0.43208D-03/
C        t2 "X"    r(%)  1.148
        data B23A, B23B, B23C, B23D, B23E /0.26987D+02,
     1  -0.70068D-01, -0.14427D+00, 0.20516D-01, -0.10215D-02/
C        G0 total  r(%)  0.171
        data C23A, C23B, C23C, C23D, C23E / 0.95270D+01,
     1  -0.61553D+00, 0.11515D+00, 0.41458D-02, -0.50238D-03/
C       epsilon   r(%)  2.671
        data D23A, D23B, D23C, D23D, D23E /0.43532D+01,
     1   -0.20351D+01, 0.24114, -0.48096D-03 , -0.99351D-03/
C        t1 "Y"    r(%)  6.618
        data E23A, E23B, E23C, E23D, E23E / 0.53493D+04,
     1   -0.34232D+01, 0.47824D+00 , -0.24050D-01, 0.25984D-04/
C        t2 "Y"    r(%)  3.638
        data F23A, F23B, F23C, F23D, F23E / 0.12477D+02,
     1   0.26576D+00, -0.35979D-01, -0.11942D-01, 0.92788D-03/
 
c       "Y" FOR  "21" TERM
        data a21_k0,b21_k0,c21_k0,d21_k0,E21_k0,F21_k0/0.65849D-66 ,
     1  0.14311D+00,-0.20553D+00, 0.79330D-01,-0.79669D-02, 0.24910D-03/
        data a211,b211,c211,d211,e211,f211/0.65884d-12, -1.10380d0,
     1  0.12066d0, 0.00997d0, -0.00491d0, 0.00035d0/
c       data for tau1 ("Y", "21")
        data a212,b212,c212,d212,e212,f212/0.54418d-12, -0.83587d0,
     1   0.20806d0, -0.05608d0, 0.00713d0, -0.00034d0/
c       data for tau2 ("Y", 21 )
 
c       "X"  FOR  "21" TERM
        data a21_F0,b21_F0,c21_F0,d21_F0,E21_F0,F21_F0/0.12524D-64,
     1  -0.20424D+00,-0.12240D+00, 0.76295D-01,-0.85239D-02,0.30576D-03/
        data a21_F1,b21_F1,c21_F1,d21_F1,E21_F1,F21_F1/0.29219D-51,
     1  -0.19336D+00,-0.11755D+00, 0.73404D-01,-0.81282D-02,0.28816D-03/
        data a21_F2,b21_F2,c21_F2,d21_F2,E21_F2,F21_F2/0.90029D-39 ,
     1  0.28128D+01, -0.13524D+01,0.30663D+00,-0.27084D-01,0.86170D-03/
C       ===============================================================
        S(T,A,B,C,D,E,F) = A*dexp(b*dlog(T) + C*(dlog(T))**2 +
     1   D*(dlog(T))**3 + E*(dlog(T))**4 + F*(dlog(T))**5 )
C       ===============================================================

      temp=temp1
      ntype=ntype1
      freqlo = freqlo1
      dfreq=dfreq1

      if(temp.lt.18.or.temp.gt.7000)then
       print*,'h2he_v1s: Warning'
       print*,'Temperature should be  18 < T < 7000'
      end if

      NF = 1 + INT((FREQMAX-FREQLO)/DFREQ)
      nf1=nf

      LIST = NF
      DO 10 I=1,NF
      FREQ(I)=FREQLO+DFLOAT(I-1)*DFREQ
      freq1(i)=freq(i)
      alfatot1(i)=0.
10    ALFATOT(I)=0.
 
 
        CALL PARTFCT4 (TEMP, Q1, NTYPE, JRANGE1, MARKER, W1)
C       PARTITION FUNCTION at temperature TEMP
        ALFA=1./(0.69519*TEMP)
c       *****************************************************************
c       For "21" use K0 + Y
        ibgama=0
c       moments for K0 (==X) first, for "21" term
        GF0=S(temp,A21_F0,B21_F0,C21_F0,D21_F0,E21_F0,F21_F0)
        GF1=S(temp,A21_F1,B21_F1,C21_F1,D21_F1,E21_F1,F21_F1)
        GF2=S(temp,A21_F2,B21_F2,C21_F2,D21_F2,E21_F2,F21_F2)
 
        call paramk0(gF0,gF1,gF2,temp,tauF1,tauF2)
 
        GK0=S(temp,A21_K0,B21_K0,C21_K0,D21_K0,E21_K0,F21_K0)
        tauk1=s(temp,a211,b211,c211,d211,e211,f211)
        tauk2=s(temp,a212,b212,c212,d212,e212,f212)
 
C       IBGAMA tellS in ADDSPEC whever to call BGAMA1 (BC) or BGAMA0 (K0)
        CALL ADDSPEC4(ibgama,TEMP,LIKE,GF0,tauF1,tauF2,GK0,tauK1,tauK2,
     1   -1,-1,2,1, FMAX,NVIB1,NVIB2)
        DO 20 I=1,NF
20      ALFATOT(I)=ABSCOEF(I)
C
c       For "01" use K0 + Y
        ibgama=0
 
c       Use fitted parameters (A,B,C,...)
        TAUF1=S(TEMP, A01A, A01B, A01C, A01D, A01E,0.D0)*1.D-14
        TAUF2=S(TEMP, B01A, B01B, B01C, B01D, B01E,0.D0)*1.D-14
        TAUK1=S(TEMP, E01A, E01B, E01C, E01D, E01E,0.D0)*1.D-14
        TAUK2=S(TEMP, F01A, F01B, F01C, F01D, F01E,0.D0)*1.D-14
        EPSK =S(TEMP, D01A, D01B, D01C, D01D, D01E,0.D0)
        GTOTAL=S(TEMP, C01A, C01B, C01C, C01D, C01E,0.D0)*1.D-65
        GF0=GTOTAL/(1.D0+EPSK)
        GK0=GTOTAL* EPSK/(1.D0+EPSK)
 
        ymaxx=6000.d0
        if(temp.lt.50.d0) ymaxx=2800.d0
C       IBGAMA tellS in ADDSPEC whever to call BGAMA1 (BC) or BGAMA0 (K0)
 
        CALL ADDSPEC4(ibgama,TEMP,LIKE,GF0,tauF1,tauF2,GK0,tauK1,tauK2,
     1   -1,-1,0,1, ymaxx ,NVIB1,NVIB2)
 
        DO 22 I=1,NF
22      ALFATOT(I)=ALFATOT(I) + ABSCOEF(I)
 
C       for "23" use BC + Y
 
C       FITTED PARAMETERS USED ("BC+Y")
        TAUF1=S(TEMP, A23A, A23B, A23C, A23D, A23E,0.D0)*1.D-14
        TAUF2=S(TEMP, B23A, B23B, B23C, B23D, B23E,0.D0)*1.D-14
        TAUK1=S(TEMP, E23A, E23B, E23C, E23D, E23E,0.D0)*1.D-14
        TAUK2=S(TEMP, F23A, F23B, F23C, F23D, F23E,0.D0)*1.D-14
        EPSK=S(TEMP, D23A, D23B, D23C, D23D, D23E,0.D0)
        GTOTAL=S(TEMP, C23A, C23B, C23C, C23D, C23E,0.D0)*1.D-65
        GF0=GTOTAL/(1.D0+EPSK)
        GK0=GTOTAL* EPSK/(1.D0+EPSK)
 
C       IBGAMA tellS in ADDSPEC whever to call BGAMA1 (BC) or BGAMA0 (K0)
        IBGAMA=1
        CALL ADDSPEC4(ibgama,TEMP,LIKE,GF0,tauF1,tauF2,GK0,tauK1,tauK2,
     1   -1,-1,2,3, FMAX,NVIB1,NVIB2)
 
      DO 24 I=1,NF
24    ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)
C
c       "45" TERM, NEGLIGIBLE
c
        MESSAGE = 'GRANDTOTAL'
        CALL ALFAMOM4 (ALFATOT,TEMP,MESSAGE,3,LIKE,-1,-1,-1,-1)

       do i=1,nf
        alfatot1(i)=alfatot(i)
       end do

11     continue
       
       RETURN

       END
 
      SUBROUTINE ADDSPEC4(ibgama,temp,like,gF0,tauF1,tauF2,
     1  GK0,tauK1,tauK2,LAMBDA1,LAMBDA2,LAMBDA,LVALUE,
     1  FREQCUT,NVIB1,NVIB2)
C
C     This program generates detailed listing of CIA ALFA(OMEGA).
C     FREQCUT serves two purposes: cut-off for moment integrals, and
C     to skip interpolation of G(OMEGA) for OMEGA > FREQCUT.
C     LIKE=1 for like systems (as H2-H2) (zero else)
c     IBGAMA=0 (X:==K0); IBGAMA=1 (X:==BC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*10 MESSAGE, MARKER
      CHARACTER*20 NHOL
      COMMON /PRTSUM14/ Q1,W1(2),Q2,W2(2)
      common /prtsum24/ JRANGE1,JRANGE2,MARKER
      COMMON /SPECIT14/ LIST
      common /specit24/ FREQ(601),ABSCOEF(601),ALFATOT(601)
      DIMENSION  BETA(3)
      DATA  SCALEF/1.D80/,  BOLTZWN/.6950304D0/
      DATA HBAR,PI,CLIGHT/1.054588757D-27,3.141592653589797D0,
     1   2.9979250D10/
C
        if ( (IBGAMA.eq.0).or.(ibgama.eq.1)) go to 1113
1112    continue
C       write(6,1124)
1124    format(/, ' IBGAMA undefined, STOP 123',/)
        stop 123
1113    continue
 
      DO 10 I=1,LIST
   10 ABSCOEF(I)=0.D0
      TWOPIC=2.D0*PI*CLIGHT
      MESSAGE='NON-COMUL.'
      IF (LIKE.NE.1) LIKE=0
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*(2.68676D19)**2
      CALIB=CALIB/FLOAT(1+LIKE)
      BETA(1)=1.D0/(BOLTZWN*TEMP)
 
C        WRITE (6,270) LAMBDA1,LAMBDA2,LAMBDA,LVALUE,JRANGE1,MARKER,
C     1   W1(1),W1(2)
270      FORMAT ('1 LAMBDA1,LAMBDA2, LAMBDA,LVALUE=',2I3,2X,2I3,8X,
     1 'RANGE1=',I2,8X, 'OUTPUT ADDSPEC',10X,A10,2F10.2,/,
     2 1X,45(1H=) )
 
C     ROTATIONAL SPECTRUM FOR THE DETAILED LISTING   ************
 
170    JPLUSL=JRANGE1+LAMBDA
C
C     Molecule interacting with an atom (but not another molecule):
C
      DO 210 I=1,JRANGE1
         J=I-1
      DO 210 IP=1,JPLUSL
         JP=IP-1
         CGS=CLEBSQR4(J,LAMBDA,JP)
         IF (CGS) 210,210,180
  180    P=H2POPL4(0,J,TEMP,W1)/Q1
         OMEGA1=H2ELEV4(NVIB1,JP)-H2ELEV4(0,J)
 
         FAC=CALIB*P*CGS
         DO 200 IQ=1,LIST
            FRQ=FREQ(IQ)-OMEGA1
            IF (DABS(FRQ)-FREQCUT) 190,200,200
  190       WKF=FREQ(IQ)*(1.-DEXP(-BETA(1)*FREQ(IQ)))*FAC
                if(Ibgama.eq.0) XX=GF0*BGAMA04(frq,tauF1,tauF2,temp)
                if(Ibgama.eq.1) XX=GF0*BGAMA14(frq,tauF1,tauF2,temp)
       IF(GK0.NE.0.0) XX = XX + GK0 * BGAMAY(FRQ,tauK1,tauK2,temp)
                ABSCOEF(IQ)=ABSCOEF(IQ) + XX*WKF
  200    CONTINUE
  210 CONTINUE
C
220   CALL ALFAMOM4 (ABSCOEF,TEMP,MESSAGE,6,LIKE,
     1              LAMBDA1,LAMBDA2,LAMBDA,LVALUE)
      RETURN
C
  290 FORMAT (1X,10E13.4)
  300 FORMAT (1X,10F13.1)
      END
 
      SUBROUTINE PARTFCT4 (TEMP, Q, NTYPE, JRANGE, MARKER, W)
C
C     NTYPE = 0 for equilibrium hydrogen (or other gas in therm. eq.)
C     NTYPE = 1 redefines the weights W for "normal" hydrogen with a
C               ratio of ortho:para H2 of 3:1.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*10  MARKER
      DIMENSION W(2)
      DATA  BOLTZWN /.6950304D0/
C
      Q=0.
      J=0
   30 DQ=H2POPL4(0,J,TEMP,W) + h2popl4(1,j,temp,w) +
     1   h2popl4(2,j,temp,w) + h2popl4(3,j,temp,w) +
     1   h2popl4(4,j,temp,w) + h2popl4(5,j,temp,w)
 
      Q=Q+DQ
      J=J+1
      IF (DQ.GT.Q/1000.D0) GO TO 30
      J=-1
      S=0.
      IF (NTYPE.EQ.1) GO TO 60
   40 J=J+1
        ddd = H2POPL4(0, J, TEMP, W) + h2popl4(1,j,temp,w) +
     1   h2popl4(2,j,temp,w) + h2popl4(3,j,temp,w) +
     1   h2popl4(4,j,temp,w) + h2popl4(5,j,temp,w)
 
      S=S+ddd/Q
      IF (S.LE.0.996D0) GO TO 40
      JRANGE=J+3
      MARKER='EQUILIBRM.'
      GO TO 90
C
C     "NORMAL" HYDROGEN REQUIRES REDEFINITION OF W(2):
C
   60 SEV=0.
   70 J=J+1
      DS=H2POPL4(0, J, TEMP, W)+ h2popl4(1,j,temp,w) +
     1   h2popl4(2,j,temp,w) + h2popl4(3,j,temp,w) +
     1   h2popl4(4,j,temp,w) + h2popl4(5,j,temp,w)
 
      S=S+DS
      SEV=SEV+DS
      J=J+1
      S=S+H2POPL4(0, J, TEMP, W)+ h2popl4(1,j,temp,w) +
     1   h2popl4(2,j,temp,w) + h2popl4(3,j,temp,w) +
     1   h2popl4(4,j,temp,w) + h2popl4(5,j,temp,w)
 
      IF (S.LE.0.996*Q) GO TO 70
      JRANGE=J+3
      SODD=S-SEV
      W(2)=W(2)*(3.*SEV)/SODD
C
C     DEFINITION OF "NORMAL" HYDROGEN: 3*S(EVEN) = S(ODD)
C
      Q=4.*SEV
      MARKER='NORMAL (!)'
C
   90 continue
C     WRITE (6,310) TEMP,Q,W(1),W(2)
  310 FORMAT (' PARTITION SUM OF H2 AT T=',F8.2,
     1   8H K IS Q=,E12.4,20X,2F10.3)
      RETURN
      END
 
      FUNCTION H2POPL4 (NVIB, J, T, W)
C     UNNORMALIZED (!) Boltzmann factor of the roto- vibrational levels
C          of the H2 molecule.  (To be normalized by Q = partition sum)
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION W(2)
      DATA  BOLTZWN /.6950304D0/,   INIT /1/
C
      IF (INIT) 10,20,10
   10 INIT=0
      W(1) = 1.
      W(2) = 3.
C     set W(2) = 0. as A TEMPORARY MODIFICATION TO SIMULATE para-H2
C
20      ddd = H2ELEV4(nvib,J)
        H2POPL4 = DFLOAT(2*J+1)*W(1+MOD(J,2))*DEXP(-ddd/(BOLTZWN*T))
c       print 777, ddd, h2popl,j
777     format(' In H2POPL: H2ELEV=',e10.3,' H2POPL=',e10.3,' J=',i5)
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION H2ELEV4 (NVIB, JROT)
C     Roto- vibrational energy levels of the H2 molecule
C     Needed for the computation of the partition sum Q of H2 ******
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
c       The following VALUES ARE THE RESULT OF FITTING THE BOUND STATE
C       LEVELS OF H-H BASED ON THE KOLOS et al. POTENTIAL FUNCTION;
C       fitted as function of J(J+1) FOR EACH VIBRATIONAL LEVELS 'V'
C       RESULTS BELOW below ARE for v=0,1,2,3,4,5.
 
         EH2(I, A, B, D, H, XL, G, P, R)= A + B*DFLOAT(I) -
     1   D*1.D-2*DFLOAT(I)**2 + H * 1.D-5 * DFLOAT(I)**3 -
     1   XL*1.D-8*DFLOAT(I)**4 + G * 1.D-11*DFLOAT(I)**5 -
     1   P * 1.D-14 * DFLOAT(I)**6 + R*1.D-17*DFLOAT(I)**7
C       A - THE ROTATIONAL ENERGY FOR (NVIB, J=0)
 
        if(nvib.eq.0) h2elev4=eh2(jrot*(jrot+1), -0.361132D5,
     1  59.33451D0, 4.56510D0, 4.568777D0, 4.366951D0, 2.863721D0,
     1  1.051977D0, 0.156903D0)  + 36113.d0
c       add 36113.d0 to each energy to prevent h2elev/(kT) being too
c       big at small temperatures; DOES NOT AFFECT THE RESULTANT POPULATION
C       FACTORS P(J)
c       ======================================================
C       ACCURACY OF FIT: 0.04% !!!
        if(nvib.eq.1) h2elev4=eh2(jrot*(jrot+1), -0.319511D5,
     1  56.3742D0, 4.4050D0, 4.515341D0, 4.614924D0, 3.301549D0,
     1  1.32896D0, 0.212922D0)   + 36113.d0
 
C       ACCURACY OF FIT: 0.03% !!!
        if(nvib.eq.2) h2elev4=eh2(jrot*(jrot+1),-0.280238D5, 53.482D0,
     1  4.28D0, 4.552623D0, 4.957454D0, 3.695921D0, 1.469646D0,
     1  0.199078D0)   + 36113.d0
 
C       ACCURACY OF FIT: 0.034% !!!
        if(nvib.eq.3) h2elev4=eh2(jrot*(jrot+1), -.243268d+05,
     1 0.5050d+02, 0.3818d+01, 0.2393d+01, -.8313d+00, -.4652d+01,
     1  -.4660d+01, -.1628d+01) + 36113.d0
 
        if(nvib.eq.4) h2elev4=eh2(jrot*(jrot+1), -.208566d+05,
     1  0.4762d+02, 0.3605d+01, 0.1511d+01, -.3982d+01, -.1106d+02,
     1 -.1108d+02, -.4167d+01) + 36113.d0
 
        if(nvib.eq.5) h2elev4=eh2(jrot*(jrot+1), -.176133d+05,
     1  0.4481d+02, 0.3511d+01, 0.1099d+01, -.6643d+01, -.1913d+02,
     1  -.2174d+02, -.9433d+01) + 36113.d0
 
        if(nvib.gt.5) stop 888
C
      RETURN
      END
 
      SUBROUTINE ALFAMOM4(ALFATOT,TEMP,MESSAGE,K,LIKE,
     1                   LAMBDA1,LAMBDA2,LAMBDA,LVALUE)
C
C     Writes out the ALFA (absorption) array and computes moments
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*10 MESSAGE
      COMMON /SPECIT14/ LIST
      common /specit24/ FREQ(601),A1(601),A2(601)
      DIMENSION ABS2(601), XI(2), FX(2), ABSP(601), ALFATOT(601)
      dimension freq1(601)
      integer LIST1,MS,KS
      DATA EPS, XI(1) /1.D-8, 1.D0/, PI /3.141592653589797D0/
      DATA CLIGHT,CLOSCHM /2.9979250D10, 2.68675484D19/
      DATA HBAR,BK / 1.054588757D-27, 1.38066244D-16/
C
      TWOPIC= 2.D0*PI*CLIGHT
      CALIB = TWOPIC*((4.D0*PI**2)/(3.D0*HBAR*CLIGHT))*CLOSCHM**2
     1        /DFLOAT(1+LIKE)
C      WRITE (K,330) FREQ(1),FREQ(LIST),FREQ(2)-FREQ(1),TEMP,
C     1   LAMBDA1,LAMBDA2,LAMBDA,LVALUE,MESSAGE
330   FORMAT ('1 ABSORPTION FROM',F10.1,' CM-1 TO',F7.1,' CM-1.',
     1 ' every',F8.2, ' CM-1',5x, ' Temperature:', f8.2,' K',/,
     1  2X,4I2,5X,' IN UNITS OF CM-1 AMAGAT-2.',4X,A10,/,
     1 1x,45(1h=),/)
C      WRITE (K,340) (ALFATOT(I),I=1,LIST)
340   FORMAT ( 90(10E13.5,/))

      LIST1=LIST
      DO I=1,LIST
       FREQ1(I)=FREQ(I)
      END DO
      MS=1
      KS=2

C
C     MOMENTS OF THE TRANSL/ROTATIONAL CIA SPECTRUM
C
      CALL SPLINE(LIST1,MS,KS,EPS,FREQ1,ALFATOT,XI,FX,ALFA1,NR,ABS2)
      ALFA1=ALFA1*TWOPIC/CLOSCHM**2
      DO 250 I=2,LIST
         E2=DEXP(HBAR*TWOPIC*FREQ(I)/(BK*TEMP))
         DIV=FREQ(I)*(E2-1.)/(E2+1.)
         ABSP(I)=ALFATOT(I)/DIV
  250 CONTINUE
      FR=(FREQ(2)/FREQ(3))**2
      ABSP(1)=(ABSP(2)-FR*ABSP(3))/(1.-FR)
      CALL SPLINE(LIST1,MS,KS,EPS,FREQ1,ABSP,XI,FX,GAMA1,NR,ABS2)
      GAMA0=GAMA1*TWOPIC/CALIB
      GAMA1=GAMA1*HBAR/(2.*BK*TEMP*CLOSCHM**2)
      DO 260 I=1,LIST
  260 ABSP(I)=ABSP(I)*FREQ(I)**2
      CALL SPLINE(LIST1,MS,KS,EPS,FREQ1,ABSP,XI,FX,DELTA1,NR,ABS2)
      DELTA1=DELTA1*TWOPIC**3/CALIB
      CORR=-1.
      RETURN
      END
 
      FUNCTION CLEBSQR4 (L,LAMBDA,LP)
C     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)
C     FOR INTEGER ARGUMENTS ONLY
C     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FC=DFLOAT(2*LP+1)
      GO TO 10
C
      ENTRY THREEJ24
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
   10 CLEBSQR4=0.
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
   30 P=FC*F*FCTL4(LAMBDA+L-LP)*FCTL4(LAMBDA+LP-L)
      CLEBSQR4=P/(FCTL4((LAMBDA+L-LP)/2)*FCTL4((LAMBDA+LP-L)/2))**2
      RETURN
      END
 
      FUNCTION FCTL4 (N)
C     Simple FACTORIALS program
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      P(Z)=((((-2.294720936D-4)/Z-(2.681327160D-3))/Z+(3.472222222D-3))/
     1Z+(8.333333333D-2))/Z+1.
      FCTL4=1
      IF (N.LE.1) RETURN
      IF (N.GT.15) GO TO 20
      J=1
      DO 10 I=2,N
   10 J=J*I
      FCTL4=DFLOAT(J)
      RETURN
   20 Z=DFLOAT(N+1)
      FCTL4=(DEXP(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)
      RETURN
      END
 
        SUBROUTINE PARAMK0(G0,G1,G2,TEMP,TAU1,TAU2)
c       *****************************************************************
c                       K0
c       ******************************************************************
c       Computes the parameters for K0 function from the three first moments
        implicit double precision (a-h,o-z)
        DATA HBAR, BOLTZK / 1.05458875D-27, 1.380662D-16 /
        T=TEMP
        TAU0 = HBAR/(2.D0 * BOLTZK * T)
        DELT=(TAU0*G1)**2 -4.*(G1*G1/G0+G1/TAU0-G2)*TAU0*TAU0*G0
        if(delt.le.0.d0) go to 999
        tau1=(-DSQRT(DELT)-TAU0*G1)/(2.*(G1*G1/G0+G1/TAU0-G2))
        if(tau1.le.0.d0) go to 999
        tau1=DSQRT(tau1)
        tau2=TAU0*tau1*G0/(G1*tau1*tau1-TAU0*G0)
        go to 889
999     continue
C       write (6,1177) delt, tau1
1177    format(' One of the following is negative for K0, ',
     1  'K0 skipped ',/,' delt, tau1:', 2e13.4,/,' STOP 111')
        stop 111
889     return
        end
 
        SUBROUTINE PARAMK(G0,G1,G2,TEMP,TAU1,TAU2)
c       *****************************************************************
c                        Y function (BGAMAY)
c       ******************************************************************
        implicit double precision (a-h,o-z)
        DATA HBAR, BOLTZK / 1.05458875D-27, 1.380662D-16 /
c       G0, G1  G2 moments for Y  function
        AA=g1/g0
        BB=g2/g0
        T = TEMP
        TAU0 = HBAR/(2.D0 * BOLTZK * T)
 
        T0 = TAU0
        X = AA - 1.d0/t0
        cc =  4.d0* BB - 3.d0*AA**2 - 6.d0*AA/t0 + 9.d0/t0**2
        IF(CC.LT.0.0) GO TO 999
        GO TO 998
999     G0=0.0
C       SET PARAMETERS FOR ZERO-VALUED Y_FUNCTION (CAN NOT BE DETERMINED)
        TAU1=1.E-14
        TAU2=1.E-14
c       WRITE(6,222)
222     FORMAT(/, ' FUNCTION "Y" CAN NOT BE DETERMINED, ',
     1  ' TAKEN AS ZERO ',/)
        RETURN
998     Y = DSQRT(CC)
        iset=0
c       Set #1:
        if((-X-Y).lt.0.d0) go to 333
        iset=1
        tau2=dsqrt(2.d0*t0/(-X-Y))
        tau1=2.d0*t0/(3.d0*X+Y)/dabs(tau2)
c       Set #2:
333     if( (-X+Y) .lt.0.d0 ) go to 335
        iset=2
        tau2=dsqrt(2.d0*t0/(-X+Y))
        tau1= 2.d0*t0/(3.d0*X-Y)/dabs(tau2)
 
335     continue
C       write(6,1999) iset, tau1,tau2
1999    format(' Parameters for Y function (set #',i2,'):',2e12.5)
        if (iset.eq.1 .or. iset.eq.2) Return
        if(iset .eq.0)  write(6,994)
994     format(/,' Parameters for Y function not found, STOP 155',/)
        stop 155
        end
 
        SUBROUTINE PARAMBC(g0,g1,g2,temp,tau1,tau2)
c       ************************************************************************
c                       BC
c       ************************************************************************
        implicit double precision (a-h,o-z)
        DATA HBAR, BOLTZK / 1.05458875D-27, 1.380662D-16 /
        T=temp
        TAU0 = HBAR/(2.D0 * BOLTZK * T)
        TTA = G0 * TAU0/ G1
        XXX = (G2*TTA-G0*(1.+TAU0**2/TTA))/(G0*(TAU0/TTA)**2)
        if(xxx.le.0.d0) go to 99
        go to 98
99      continue
C       write (6,122)
        stop 777
122     format (/, ' xxx<0 in PARAMBC, STOP 777',/)
98      TAU1=DSQRT(xxx)
        TAU2=TTA/TAU1
        return
        end
 
      FUNCTION XK04(X)
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
      XK04 = - DLOG(X/2.)*FI0+P
      RETURN
   20 T=(2./X)
      P=(((((.00053208*T-.00251540)*T+.00587872)*T
     1 -.01062446)*T+.02189568)*T-.07832358)*T+
     2 1.25331414
      X=dMIN1(X,330.d0)
      XK04=DEXP(-X)*P/DSQRT(X)
      RETURN
      END
 
      FUNCTION XK14(X)
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
      XK14=X*dLOG(X/2)*FI1+P
      RETURN
   20  T=2./X
      P=(((((-.00068245*T+.00325614)*T-.00780353)*T+.01504268)*T-
     1 .03655620)*T+.23498619)*T+1.25331414
      X=dMIN1(X,330.d0)
      XK14=DSQRT(X)*DEXP(-X)*P
      RETURN
      END
 
      FUNCTION BGAMAY(FNU,TAU1,TAU2,TEMP)
C     NORMALIZATION SUCH THAT 0TH MOMENT EQUALS UNITY, OK
C     FNU IS THE FREQUENCY IN CM-1; TEMP IN KELVIN.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      DATA TWOPIC,BKW,HBOK,PI/1.883651568D11,.6950304256D0,
     1 7.638280918D-12, 3.141592654D0/
 
        B1 =  BGAMA14(FNU,TAU1,TAU2,TEMP)
c       B1 is a value of a B-C lineshape at FNU
        omega = twopic * fnu
        tau0 = hbok /(2.*Temp)
        BGAMAY = DEXP( -(tau2-dabs(tau2))/tau1) * tau1 * abs(tau2)/tau0
     1           * omega * B1
c       has units of BGAMA1 (BC)
        RETURN
        END
 
      FUNCTION BGAMA14(FNU,TAU1,TAU2,TEMP)
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
      AAA=DMIN1(AAA,430.d0)
      BGAMA14=(TAU1/PI)*dEXP(AAA+FNU/(2.d0*BKW*TEMP))*
     1 XK14(X)/DENOM
      RETURN
      END
 
      FUNCTION BGAMA04(FNU,TAU5,TAU6,TEMP)
C     K0 LINE SHAPE MODEL
C     NORMALIZATION SUCH THAT 0TH MOMENT EQUALS UNITY
C     FNU IS THE FREQUENCY IN CM-1; TEMP IN KELVIN.
        implicit double precision (a-h,o-z)
      DATA TWOPIC,BKW,HBOK,PI/1.883651568d11,.6950304256d0,
     1 7.638280918d-12, 3.141592654d0/
 
      TAU4=DSQRT(TAU5*TAU5+(HBOK/(TEMP*2.d0))**2)
      OMEGA=TWOPIC*FNU
      XNOM=1.d0/(TAU6*TAU6)+OMEGA*OMEGA
      X=TAU4*dSQRT(XNOM)
      TAU56=TAU5/TAU6
      TAU56=dMIN1(TAU56,430.d0)
      BGAMA04=(TAU5/PI)*dEXP(TAU56+FNU/(2.d0*BKW*TEMP))*
     1 XK04(X)
      RETURN
      END
 
      SUBROUTINE SPLINE(L,M,K,EPS,X,Y,T,SS,SI,NR,S2)
C
C     Spline interpolation and quadrature, third order after Greville.
C     Input arguments L...Y, output SS...NR. (Except NR=-1 allows to
C           read in S2(0) and S2(L) if desired)
C     L data points X(1), Y(1) ... X(L),Y(L)
C     EPS=error criterion, typically EPS=1.E-5 for 5 deci. places accuracy
C     M arguments T(1)..T(M) for which function values SS(1)..SS(M), for
C           K=0; or first or second derivative for K=1 or -1, respectively.
C     Note that M has to be at least equal to 1.
C     SI=integral (over whole interval) for K=2 only.
C     'NATURAL' spline functions (if NR.NE.-1) with S2(0)=S2(L)=0 UNLESS
C           NR was set initially to -1; in that case, read in S2(0),S2(L).
C     NR indicates the number of out-of-range calls. X(1)@T(I)@X(L)
C     Extrapolate with caution. (ASSUMPTION D2Y/DX2 = 0.)
C     S2(I) is the 2nd derivative at X=X(I) and is computed internally.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(601), Y(601), T(2), SS(2), S2(601)
      DELY(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
      B(I)=(X(I)-X(I-1))*0.5/(X(I+1)-X(I-1))
      C(I)=3.*(DELY(I)-DELY(I-1))/(X(I+1)-X(I-1))
      N=L
      N1=N-1
      DO 10 I=2,N1
   10 S2(I)=C(I)/1.5
      OMEGA=1.0717968
      IC=0
C
C     'NATURAL' spline functions of third order unless NR=-1:
C
      IF (NR.EQ.-1) GOTO 20
      S2(N)=0.
      S2(1)=0.
   20 ETA=0.
      IC=IC+1
      SM=DABS(S2(1))
      DO 30 I=2,N
         IF (DABS(S2(I)).GT.SM) SM=DABS(S2(I))
   30 CONTINUE
      EPSI=EPS*SM
      DO 50 I=2,N1
         W=(C(I)-B(I)*S2(I-1)-(0.5-B(I))*S2(I+1)-S2(I))*OMEGA
         IF (DABS(W)-ETA) 50,50,40
   40    ETA=DABS(W)
   50 S2(I)=S2(I)+W
      IF (ETA-EPSI) 60,20,20
C      ENTRY IXPOLAT
C
C      THIS ENTRY USEFUL WHEN ITERATION PREVIOUSLY COMPLETED
C
C      N=L
C      N1=N-1
C      IC=-1
   60 IF (K.EQ.2) GO TO 260
      NR=0
      DO 250 J=1,M
         I=1
         IF (T(J)-X(1)) 110,210,80
   80    IF (T(J)-X(N)) 100,190,150
   90    IF (T(J)-X(I)) 200,210,100
  100    I=I+1
         GO TO 90
  110    NR=NR+1
         HT1=T(J)-X(1)
         HT2=T(J)-X(2)
         YP1=DELY(1)+(X(1)-X(2))*(2.*S2(1)+S2(2))/6.
         IF (K) 140,130,120
  120    SS(J)=YP1+HT1*S2(1)
         GO TO 250
  130    SS(J)=Y(1)+YP1*HT1+S2(1)*HT1*HT1/2.
         GO TO 250
  140    SS(J)=S2(I)
         GO TO 250
  150    HT2=T(J)-X(N)
         HT1=T(J)-X(N1)
         NR=NR+1
         YPN=DELY(N1)+(X(N)-X(N1))*(S2(N1)+2.*S2(N))/6.
         IF (K) 180,170,160
  160    SS(J)=YPN+HT2*S2(N)
         GO TO 250
  170    SS(J)=Y(N)+YPN*HT2+S2(N)*HT2*HT2/2.
         GO TO 250
  180    SS(J)=S2(N)
         GO TO 250
  190    I=N
  200    I=I-1
  210    HT1=T(J)-X(I)
         HT2=T(J)-X(I+1)
         PROD=HT1*HT2
         S3=(S2(I+1)-S2(I))/(X(I+1)-X(I))
         SS2=S2(I)+HT1*S3
         DELSQS=(S2(I)+S2(I+1)+SS2)/6.
         IF (K) 240,220,230
  220    SS(J)=Y(I)+HT1*DELY(I)+PROD*DELSQS
         GO TO 250
  230    SS(J)=DELY(I)+(HT1+HT2)*DELSQS+PROD*S3/6.
         GO TO 250
  240    SS(J)=SS2
  250 CONTINUE
  260 SI=0.
      DO 270 I=1,N1
         H=X(I+1)-X(I)
  270 SI=SI+0.5*H*(Y(I)+Y(I+1))-H**3*(S2(I)+S2(I+1))/24.
      IF (K.EQ.2) NR=IC
      RETURN
      END
