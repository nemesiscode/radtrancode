      subroutine h2h2_v1s(ntype1,temp1,freqlo1,freqmax,deltastep,list1,
     1freq1,alfatot1)
C     *************************************************************
C     Subroutine to calculate an H2-H2 collision induced absorption spectrum
C     in the v1 fundamental band
C
C     Input variables:
C       ntype		integer 0 = equilibrium hydrogen
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
C     *******************************************************************
c     An error has been found in the  program modeling
c       the collision induced spectra
c       of hydrogen pairs in the fundamental band.
c       The slight mistake affects mostly the spectra at higher
c       temperatures (T>200K), mainly near the centers of the roto-translational lines.
c       This version of the program (ver. 2) accounts for bound states as
c       well as for free states, and fully reconstructs the quantum mechanical lineshape
c       computations reported before (Meyer,Borysow,Frommhold, PRA,40,6931,1990).
c       The temperature range has been extended to include T=350K (max), to facilitate
c       modeling of Jupiter's atmosphere.
c       With my deepest apologies;	Aleksandra Borysow (Oct 28, 1993)

c       COPYRIGHT ALEKSANDRA BORYSOW 1991 and 1993
c       Refer and quote a paper:
c       A. Borysow, Icarus, vol.92, 273 (1991).
c       This is an updated version (erratum will be published, 1993)

c       ***********************************************************
c       Program written by: Aleksandra Borysow
c       Physics Department
c       Michigan Technological University
c       1400 Townsend Drive
c       Houghton, MI 49931 - 1295
c       e-mail: ABORYSOW@phy.mtu.edu
c       phone: (906) 487 - 2198
c       ***********************************************************
c       Program to model H2-H2 CIA spectra at low temperatures
c       for planetary interest (20-300K)
c       NOTE: low resolution model, bound states included (1993)
c       ================================================================
c       Program asks interactively for:
c       TEMPERATURE [K], choose from 20 - 300K
c       min. frequency (FREQLO) (in cm^[-1])
c       max. frequency (FREQMAX);
c       frequency interval (in cm^{-1})
c       Max. number of points: 601 (look for "NF" in the program)
c
c       Program computes absorption spectrum at temperatures
c       between  20 and 300K H2-H2 fundamental band
c       J-dependence accounted for
c       For details, see paper by W.Meyer,A.Borysow,L.Frommhold,
c       on H2-H2 fundamental band, Phys.Rev.A, vol.40, 6931, 1990.
c       Program uses model lineshapes for J1=J1'=J2=J2'==0
c       and rescales their magnitude by their M0 moments.
c

c       COPYRIGHT A.BORYSOW 1991
C     Program generates  a detailed listings of the roto-vibrational spectra
c       of H2-H2 (free-free transitions only) in the fundamental band
c     (0-->1 vibrational transition) at temperatures less than 300K
c       The model reproduces exactly the quantum mechanical computations as
c       published in: W. Meyer; A. Borysow and L. Frommhold;
c       Physical Review A40, 6931 (1989)
c
c       Program uses new RV lineshape, as described in the paper
c       by G. Birnbaum and A. Borysow, Molecular Physics,vol.73,no.1,57-68 (1991)
c
c       For this program:	Refer and quote a paper:
c       A. Borysow, Icarus, vol.92, 273 (1991).
c
        IMPLICIT double precision (A-H,O-Z)
        CHARACTER*10 MESSAGE, MARKER
        CHARACTER*11 labeta
c        PARAMETER NF = 601 ! determines max. # of points;
c       Adjust NF in all subroutines the same; take dimension of the
c       arrays equal to number of needed frequency points;
c
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet


C     NTYPE=0 for equilibrium hydrogen; NTYPE=1 for normal H2



        common/ integ / rlo, rup, frac, Ninteg
        common/temp/temp
      COMMON /PRTSUM/ JRANGE1,MARKER
        common/prtsum1/ q1, w1(2)
        common /wagia/ w(2)
        common/ bb1/ labeta
        common /bb2/ j1, j1p, j2, j2p
        common/ola/  xxm0
        common/sumrul/ g0, g1
      COMMON /SPECIT/ LIST
        dimension freq1(601),alfatot1(601)
        common/ specit1/ FREQ(601), ABSCOEF(601), ALFATOT(601)
      COMMON/SPEKTRM1/ LIKE, nvib_F, nvib_I, NTYPE
       common /spektrm2/ FREQLO,DFREQ,FMAX
        common/ read/ go
        common/om/omega_0
        common/ uncor/ xm0, tau1, tau2
        common/corr/ omega1, omega2, tauh
      COMMON/YMAX/YMAXXR, YMAXXL
      DATA LIKE, nvib_F, nvib_I / 1, 1, 0 /
c       		nvib_F - upper state, nvib_I - lower state

      DATA W /1., 3./   ! rotational weights;
      DATA FMAX / 3000.d0/       ! to be kept this way
c
        data nf/ 601 /

        Z(T,a,b,c,d) = A * DEXP( B * dlog(T) + C * (dlog(T))**2 +
     1  D * (dlog(T))**3 )

        temp=temp1
        freqlo=freqlo1
        ntype=ntype1

c initialise the weights for H2
        W(1)=1.0
        W(2)=3.0
c this is needed otherwise they are not reset after a normal calculation
      
	if ( (temp. lt. 20.) .or. (temp.gt.300.)) then
          if(idiag.gt.0)then
           print*,'h2h2_v1s: Warning'
	   print*,'Temperatures should be in range 20 < T < 300'
          endif
        end if

C	CALCULATE OMEGA_0
C	DETERMINING THE PURELY VIBRATIONAL TRANSITION FREQUENCY:(CM^-1)
	JJ=0
        OMEGA_0 = H2ELEV(NVIB_F,JJ) - H2ELEV(NVIB_I,JJ)


	list = 1 + INT (FREQMAX - FREQLO)/DELTASTEP

        list1=list

	IF(list.GT.nf) STOP 963
C	==>	CHANGE ' character nf = ' OR ADJUST DELTASTEP FOR A LARGER ONE!

C	write(*, 1233) FREQLO, FREQMAX, DELTASTEP, list
1233	FORMAT (' THE SPECTRUM WILL BE COMPUTED AT FREQUENCIES',/,
     1  ' FROM', F12.3,' TO ', F12.3, ' EVERY', F10.1, ' IN ',I4,
     1  ' STEPS ')

      DO 10 I=1,list
      FREQ(I)=FREQLO+DFLOAT(I-1)*DELTASTEP
      FREQ1(I)=FREQ(I)
10    ALFATOT(I)=0.

C      print*,'A',alfatot(100)
C	write(10, 1212) list, freqlo, deltastep
1212	format(' Rotovibrational CIA spectrum of H2-H2 ',/,
     2  ' Fundamental band',/,
     1  ' Program version 2, revised October 1993 (A. Borysow)',/,
     2 ' Absorption given at  NF=',i4,' points',/,
     1  ' From ', f10.2, ' cm-1, every ', f10.2, ' cm-1',/)

c	===============================================================
c	Compute the partition function Q1:
        CALL PARTFCT (TEMP, Q1, NTYPE, JRANGE1, MARKER, W )
	w1(1) = w(1)
	w1(2) = w(2)
      ALFA =1./(0.69519*TEMP)
      hbar = 1.054588757D-27
	boltzk = 1.38054d-16
	tau0 = hbar/(2.*boltzk*temp)
c	*****************************************************************
c	compute radial distribution function g(R) at temperature TEMP
	ldist=1
	labeta ='  call g(R)'

	rlo = 1.d0
	rup = 12.d0
	Ninteg = 601
	frac = 100.d0
	CALL MOMTS (ldist, labeta, RLO,RUP,FRAC,Ninteg,TEMP,XXM0)
c	*****************************************************************
	labeta = '<0|B0001|1>'
c	ver. 1993
	XM0 = Z(temp, 0.66747d-61, -0.34043d+1, 0.62353d00,-0.27926d-01)
	XM1 = Z(temp, 0.17096d-47, -0.33461d+1, 0.61460d00,-0.28201d-01)
	g0 = xm0
	g1 = xm1
	XM2 = Z(temp, 0.94242d-34,-0.31984d+1, 0.52278d00,-0.12498d-01)

c	calculate tau1 and tau2 of K0 from the three above moments
	DELT = (TAU0*XM0)**2-4.*(XM1*XM1/XM0+XM1/TAU0-XM2)*TAU0*TAU0*XM0
      TAU1=(-DSQRT(DELT)-TAU0*xm1)/(2.*(XM1*XM1/XM0+XM1/TAU0-XM2))
      TAU1=DSQRT(TAU1)
      TAU2=TAU0*TAU1*xm0/(xm1*TAU1*TAU1-TAU0*xm0)

c	ver. 1993
	omega1 = Z(temp,0.27279d+2,-0.69452d00, 0.11661d00,-0.61575d-3)
	omega2 = Z(temp,0.11077d+2,0.46570d-2,-0.49445d-1, 0.10357d-1)
	tauH =   Z(temp,0.39182d-12,-0.52643d00, 0.39936d-5, 0.95064d-3)


	ibgama = 0		! means, X==K0
	freqcut = 4000.d0
	call range(ibgama)

      CALL ADDSPEC1(TEMP,LIKE, 0, 0, 0, 1, freqcut,nvib_F,nvib_I,
     2  ibgama,  XM0, tau1,tau2, omega1, omega2,  tauH )
        DO 920 I=1,list
920     ALFATOT(I) = 2.* ABSCOEF(I)
	Labeta='<0|B0221|1>'


c	Moments fitted     00|0221|01
c	1993 ver:
	XM0 = Z(temp,  0.74632d-64, -0.49090d00, 0.59048d-1,0.89534d-2)
	XM1 = Z(temp, 0.21555d-50,  -0.51456d00, 0.70540d-1,0.72538d-2)
	g0 = xm0
	g1 = xm1
	XM2 = z(temp, 0.93317d-37,-0.138170d00,-0.75397d-1,0.26852d-1)

c	calculate tau1 and tau2 of K0 from the three above moments
	DELT = (TAU0*XM0)**2-4.*(XM1*XM1/XM0+XM1/TAU0-XM2)*TAU0*TAU0*XM0
      TAU1=(-DSQRT(DELT)-TAU0*xm1)/(2.*(XM1*XM1/XM0+XM1/TAU0-XM2))
      TAU1=DSQRT(TAU1)
c	1993 ver:
c	parameters for correction term:
	omega1 = z(temp, 0.19834d+2, -0.44580d00,0.60392d-1, 0.35041d-2)
	omega2 = z(temp,0.87918d+1,0.19156d00,-0.92862d-1,0.13604d-1)
	tauH = z(temp, 0.37724d-12,-0.51442d00,-0.21899d-2, 0.10429d-2)

	ibgama = 0		! means, ==K0
	call range(0)
      CALL ADDSPEC1(TEMP,LIKE, 0, 2, 2, 1, FMAX, nvib_F,nvib_I,
     2  ibgama,  XM0, tau1,tau2, omega1, omega2,  tauH )

      DO 24 I=1,list
24    ALFATOT(I)= ALFATOT(I)+ 2.* ABSCOEF(I)
c	* * * * * * * * * * * * * ** * * * * * * * * * * *
	labeta='<0|B0223|1>'


C	Moments fitted     00|0223|01
c	ver. 1993:
 	XM0 = Z(temp,0.34819d-61,-0.22398d+01,0.36696d+00,-0.15267d-01)
	XM1 = Z(TEMP,0.44881d-48,-0.24665d+01,0.40935d+00,-0.15843d-01)
	XM2 = Z(temp,0.37817d-34,-0.32010d+01,0.55644d+00,-0.17697d-01)
	g0 = xm0
	g1 = xm1

c	find tau1 and tau2 of BC...

      TTA=XM0*TAU0/XM1
	XXX = (XM2*TTA-XM0*(1.d0 + TAU0**2/TTA))/(XM0*(TAU0/TTA)**2)
      TAU1=DSQRT(XXX)
      TAU2=TTA/TAU1

c	ver. 1993:
	omega1 = Z(temp,0.27405d+2,-0.17864d+1,0.39631d00,-0.19679d-1)
	omega2 = Z(temp,0.10116d+2,-0.74263d-2,-0.63547d-1,0.12726d-1)
	tauH = Z(temp,0.55383d-12,-0.70431d00,0.54877d-01,-0.42722d-2)


	ibgama = 1		! means, ==BC
	call range(1)
      CALL ADDSPEC1(TEMP,LIKE, 0, 2, 2, 3, FMAX, nvib_F,nvib_I,
     2  ibgama,  XM0, tau1,tau2, omega1, omega2,  tauH )

        DO 324 I=1,list
 324    ALFATOT(I)= ALFATOT(I)+ 2.* ABSCOEF(I)
c	* * * * * * * * * * * ** * * * * * * * * * * *
	labeta='<0|B2023|1>'


C	Moments fitted     00|2023|01
c	ver. 1993
	XM0 = Z(temp,0.14553d-61,-0.22892d1,0.39167d00,-0.19424d-1)
	XM1 = Z(temp,0.12817d-48,-0.24553d1,0.42439d00,-0.20329d-1)
	XM2 = Z(temp,0.93224d-35,-0.34191d1,0.64311d00,-0.27519d-1)
	g0 = xm0
	g1 = xm1

c	find tau1 and tau2 of BC...
      TTA=XM0*TAU0/XM1
	XXX = (XM2*TTA-XM0*(1.d0 + TAU0**2/TTA))/(XM0*(TAU0/TTA)**2)
      TAU1=DSQRT(XXX)
      TAU2=TTA/TAU1

c	ver. 1993
	omega1 = Z(temp,  0.21178d2,  -0.20090d1,0.46603d00,-0.25764d-1)
	omega2 = Z(temp, 0.96463d1, 0.74521d-2,-0.74189d-1,0.13590d-1)
	tauH = Z(temp,  0.65879d-12, -0.76779d00,0.71105d-1,-0.51727d-2)


	ibgama = 1		! means, ==BC
	call range(1)
      CALL ADDSPEC1( TEMP,LIKE, 2, 0, 2, 3, FMAX, nvib_F,nvib_I,
     2  ibgama,  XM0, tau1,tau2, omega1, omega2,  tauH )

        DO 724 I=1,list
724     ALFATOT(I)= ALFATOT(I)+ 2.* ABSCOEF(I)
c	* * * * * * * * * * * ** * * * * * * * * * * *

	labeta='<0|B2233|1>'


C	Moments fitted     00|2223|01
c	ver. 1993:
 	XM0 = Z(temp,0.86252d-63, -0.22504d+01, 0.38397d00,-0.19256d-1)
	XM1 = Z(TEMP, 0.62723d-50, -0.23286d+1, 0.39342d00,-0.18206d-1)
	XM2 = Z(temp,  0.43604d-36,-0.33521d+1, 0.63337d00,-0.27082d-1)
	g0 = xm0
	g1 = xm1

c	find tau1 and tau2 of BC...

      TTA=XM0*TAU0/XM1
        XXX = (XM2*TTA-XM0*(1.d0 + TAU0**2/TTA))/(XM0*(TAU0/TTA)**2)
      TAU1=DSQRT(XXX)
      TAU2=TTA/TAU1

c      ver. 1993:
       omega1 = Z(temp,0.17717d+2,  -0.20101d+1,0.47054d00, -0.26063d-1)
       omega2 = Z(temp,0.88423d+1,   0.10092d00,-0.10153d00, 0.15948d-1)
       tauH = Z(temp,  0.72376d-12, -0.84037d00,0.92159d-1, -0.70604d-2)

	ibgama = 1		! means, ==BC
	call range(1)
      CALL ADDSPEC1(TEMP,LIKE, 2, 2, 3, 3, FMAX, nvib_F,nvib_I,
     2  ibgama,  XM0, tau1,tau2, omega1, omega2,  tauH )

        DO 824 I=1,list
 824    ALFATOT(I)= ALFATOT(I)+ 2.* ABSCOEF(I)


c	* * * * * * * * * * * ** * * * * * * * * * * *
c	Very weak terms like 2211, 0445, 4045, 2021, 2243, 2245 are neglected
c	* * * * * * * * * * * ** * * * * * * * * * * ** * * * * * * * * * * *
        MESSAGE = 'GRANDTOTAL'
        CALL ALFAMOM (TEMP,MESSAGE,10,LIKE,-1,-1,-1,-1)

       do i=1,list
        alfatot1(i)=alfatot(i)
       end do


       RETURN

       END

      SUBROUTINE ADDSPEC1( TEMP,
     $     LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE, CUT, nvib_F, nvib_I,
     2  ibgama,  G0, tau1,tau2,  omega11, omega22, tauH )
C

C     This subroutine generates a detailed listing of CIA ALFA(OMEGA)
c	IBGAMA=0 means uses K0, IBGAMA=1 (uses BC)
C     LIKE=1 for like systems (as H2-H2)
C
      IMPLICIT double precision (A-H,O-Z)
      CHARACTER*10 MESSAGE, MARKER
      CHARACTER*11 labeta
      character*6 note
      COMMON /PRTSUM/ JRANGE1,MARKER
	common/sumrul/ gg0, gg1
	common/prtsum1/ q1, w1(2)
      COMMON /SPECIT/ LIST
	common/specit1/ FREQ(601),ABSCOEF(601), ALFATOT(601)
	common /BB1/ labeta
        common /bb2/ J1, Jp1, J2, Jp2
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
      CALIB=CALIB/DFLOAT(1+LIKE)
      FREQCUT = CUT		! do not alter the outside parameter CUT

      BETA1 = 1.D0/(BOLTZWN*TEMP)

270      FORMAT (//,' LAMBDA1,LAMBDA2, LAMBDA,LVALUE=',2I3,2X,2I3,8X,
     1   'RANGE1=',I2,/, ' OUTPUT ADDSPEC',10X,A10,'HYDROGEN',2F10.2,/,
     2   1X,45(1H=) )

C      Molecule interacting with molecule (as opposed to atom):

	XM00 =  G0	! pure quantum result

        JPLUSL=JRANGE1+MAX0(LAMBDA1,LAMBDA2)
      DO 160 I1=1,JRANGE1
         J1=I1-1
      DO 160 IP1=1,JPLUSL
         JP1=IP1-1
         CG1S=CLEBSQR1(J1,LAMBDA1,JP1)
         IF (CG1S) 160,160,110
  110    P1=H2POPL(nvib_I,J1,TEMP,W1)/Q1
         OMEGA1=H2ELEV(nvib_I,JP1)-H2ELEV(nvib_I,J1)
C	"FIRST" H2 does not undergo vibrational transition

         DO 150 I2=1,JRANGE1
            J2=I2-1
         DO 150 IP2=1,JPLUSL
            JP2=IP2-1
            CG2S=CLEBSQR1(J2,LAMBDA2,JP2)
            IF (CG2S) 150,150,120
  120       P2=H2POPL(nvib_I,J2,TEMP,W1)/Q1
            OMEGA2=H2ELEV(nvib_F,JP2)-H2ELEV(NVIB_i,J2)
c	"second" H2 does undergo vibrational transition

	ldist = 0
	CALL MOMTS (ldist, labeta, RLO,RUP,FRAC,Ninteg,TEMP, XM0j )

	if(j1.eq.jp1) note = 'Single'
	if(j1.ne.jp1) note = 'Double'

c	Below a correction made for different M0(J1,J1',J2,J2'):
        FAC=CALIB*P1*P2*CG1S*CG2S  *  ( Xm0j/Xm00)

            DO 140 I=1,LIST
       FRQ=FREQ(I)-OMEGA1-OMEGA2
	if (frq.gt.ymaxxR) go to 140
	if (frq.lt.ymaxxL) go to 140
  130   WKF=FREQ(I) * (1.-DEXP(- BETA1 * FREQ(I))) * FAC

	if (ibgama.eq.0) xx = g0 * bgama0(frq,tau1,tau2,temp)
	if (ibgama.eq.1) xx = g0 * bgama1(frq,tau1,tau2,temp)
	yy = g0/2. * (bgama_h(frq-(omega11+omega22),tauH,temp) -
     1     bgama_h(frq-(omega22-omega11), tauH, temp))
	gg = xx + yy

        ABSCOEF(I) = ABSCOEF(I) + gg * WKF
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE

      CALL ALFAMOM(TEMP,MESSAGE,10,LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE)
      RETURN
      END

      SUBROUTINE PARTFCT (TEMP, Q, NTYPE, JRANGE, MARKER, W)
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
   30 DQ=H2POPL(0,J,TEMP,W) + h2popl(1,j,temp,w) +
     1   h2popl(2,j,temp,w) + h2popl(3,j,temp,w)

      Q = Q + DQ
      J = J + 1
      if ( (w(2).eq.0.).and.(mod(j,2).eq.1)) go to 30 ! odd J: do not compare
      IF (DQ.GT.Q/990.) GO TO 30

      J=-1
      S=0.

      IF (NTYPE.EQ.1) GO TO 60
   40 J=J+1
       DDD = H2POPL(0, J, TEMP, W) + h2popl(1,j,temp,w) +
     1   h2popl(2,j,temp,w) + h2popl(3,j,temp,w)

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
      DS=H2POPL(0, J, TEMP, W)+ h2popl(1,j,temp,w) +
     1   h2popl(2,j,temp,w) + h2popl(3,j,temp,w)

      S=S+DS
      SEV=SEV+DS
      J=J+1
      S=S+H2POPL(0, J, TEMP, W)+ h2popl(1,j,temp,w) +
     1   h2popl(2,j,temp,w) + h2popl(3,j,temp,w)

      IF (S.LE.0.999*Q) GO TO 70
      JRANGE=J+2
      SODD=S-SEV
      W(2)=W(2)*(3.*SEV)/SODD
C
C     DEFINITION OF "NORMAL" HYDROGEN: 3*S(EVEN) = S(ODD)
C
      Q=4.*SEV
      MARKER='NORMAL (!)'
C
   90  continue

      RETURN
      END


      FUNCTION H2POPL (NVIB, J, T, W)
C     UNNORMALIZED (!) Boltzmann factor of the roto- vibrational levels
C          of the H2 molecule.  (To be normalized by Q = partition sum)
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION W(2)
      DATA  BOLTZWN /.6950304D0/
C
20	DDD = H2ELEV(nvib,J)
        H2POPL = DFLOAT(2*J+1)*W(1+MOD(J,2))*DEXP(-DDD/(BOLTZWN*T))
      RETURN
      END

      FUNCTION H2ELEV (NVIB, JROT)
C     Roto- vibrational energy levels of the H2 molecule
C     Needed for the computation of the partition sum Q of H2 ******
C
      IMPLICIT double precision (A-H,O-Z)
C
c	The following results are from LEVELS (levels.lev file) of
c	Kolos et all. (HSINGLX.FOR function), fitted as function of J(J+1)
c	Fitting program FITLEVELS.FOR (.for), below only for v=0, 1, 2
C	NOTE: ROTATIONAL CONSTANTS B AND D HAVE BEEN FIXED WHILE FITTING THE
C	ROTATIONAL LEVELS, AFTER S.L.BRAGG, J.W.BRAULT,W.H.SMITH, AP.J.,VOL.263,
C	P.999-1004, (1982).

         EH2(I, A, B, D, H, XL, G, P, R)= A + B*DFLOAT(I) -
     1   D*1.D-2*DFLOAT(I)**2 + H * 1.D-5 * DFLOAT(I)**3 -
     1   XL*1.D-8*DFLOAT(I)**4 + G * 1.D-11*DFLOAT(I)**5 -
     1   P * 1.D-14 * DFLOAT(I)**6 + R*1.D-17*DFLOAT(I)**7
C	A - THE ROTATIONAL ENERGY (NVIB, J=0)

c	|E(v=0,J=0) - E(v=1,J=0)| = 4162.14 cm-1

	if(nvib.eq.0) h2elev=eh2(jrot*(jrot+1), -0.361132D5,
     1  59.33451D0, 4.56510D0, 4.568777D0, 4.366951D0, 2.863721D0,
     1  1.051977D0, 0.156903D0)  + 36113.d0
c	add 36113.d0 to each energy to prevent h2elev/(kT) being too
c	big at small temperatures
c	======================================================
	if(nvib.eq.1) h2elev=eh2(jrot*(jrot+1), -0.319511D5,
     1  56.3742D0, 4.4050D0, 4.515341D0, 4.614924D0, 3.301549D0,
     1  1.32896D0, 0.212922D0)   + 36113.d0

	if(nvib.eq.2) h2elev=eh2(jrot*(jrot+1),-0.280238D5, 53.482D0,
     1  4.28D0, 4.552623D0, 4.957454D0, 3.695921D0, 1.469646D0,
     1  0.199078D0)   + 36113.d0

	if(nvib.eq.3) h2elev=eh2(jrot*(jrot+1), -.243268d+5,
     1 0.5050d+02, 0.3818d+01, 0.2393d+01, -.8313d+00, -.4652d+01,
     1 -.4660d+01, -.1628d+01) + 36113.d0

	if(nvib.eq.4) h2elev=eh2(jrot*(jrot+1), -.208566d+05,
     1  0.4762d+02, 0.3605d+01, 0.1511d+01, -.3982d+01, -.1106d+02,
     1 -.1108d+02, -.4167d+01) + 36113.d0

	if(nvib.eq.5) h2elev=eh2(jrot*(jrot+1), -.176133d+05,
     1  0.4481d+02, 0.3511d+01, 0.1099d+01, -.6643d+01, -.1913d+02,
     1 -.2174d+02, -.9433d+01) + 36113.d0

	if(nvib.gt.5) stop 888
C
      RETURN
      END

      SUBROUTINE ALFAMOM(TEMP,MESSAGE,K,LIKE,
     1                   LAMBDA1,LAMBDA2,LAMBDA,LVALUE)
C
C     Writes out the ALFA (absorption) array
C
      IMPLICIT double precision (A-H,O-Z)
      CHARACTER*10 MESSAGE
c      PARAMETER NF = 601
      COMMON /SPECIT/ LIST
c      common/specit1/ FREQ(nf),Abscoef(nf),Alfatot(nf)
      common/specit1/ FREQ(601),Abscoef(601),Alfatot(601)
      DIMENSION  XI(1),  aa(601), abs2(601), absp(601), fx(1)
      DATA EPS, XI(1) /1.D-6, 3000.D0/, PI /3.141592653589797D0/
      DATA CLIGHT,CLOSCHM /2.9979250D10, 2.68675484D19/
      DATA HBAR,BK / 1.054588757D-27, 1.38066244D-16/
	data nf /601 /

	do 1555 i=1, list
	absp(i) = 0.d0
	abs2(i) = 0.d0
1555	continue

	IF(MESSAGE.EQ. 'GRANDTOTAL' ) go to 11
	if(message.eq. 'NON-COMUL.') go to 12
	stop 987
11	do 22 i=1, list
22	aa(i)=alfatot(i)
	go to 13
12	do 23 i=1, list
23	aa(i)=abscoef(i)
13	continue

      TWOPIC= 2.*PI*CLIGHT
      CALIB = TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2
     1        /DFLOAT(1+LIKE)
C      WRITE (K,330) FREQ(1),FREQ(LIST),FREQ(2)-FREQ(1),TEMP,
C     1   LAMBDA1,LAMBDA2,LAMBDA,LVALUE,MESSAGE
  330 FORMAT (/,' ABSORPTION FROM',F10.1,' CM-1 TO',F7.1,' CM-1.',
     1 ' every',F8.2, 10x, /, ' Temperature:', f8.2, /,
     1  2X,4I2,5X,' IN UNITS OF CM-1 AMAGAT-2.',4X,A10/)
C      WRITE (K,340) ( AA(I),I=1,LIST)
  340 FORMAT ( 110(6d12.4,/))
c	IF(MESSAGE.EQ. 'GRANDTOTAL' ) print 122, ( freq(i), aa(i), i=1, nf)
122	format(f10.2, 2x, e12.4 )

      RETURN
      END

      FUNCTION XK0(X)
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
      XK0=-dLOG(X/2.)*FI0+P
      RETURN
   20 T=(2./X)
      P=(((((.00053208*T-.00251540)*T+.00587872)*T
     1 -.01062446)*T+.02189568)*T-.07832358)*T+
     2 1.25331414
      X=dMIN1(X,330.d0)
      XK0=dEXP(-X)*P/dSQRT(X)
      RETURN
      END

      FUNCTION XK1(X)
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
      XK1=X*dLOG(X/2)*FI1+P
      RETURN
   20  T=2./X
      P=(((((-.00068245*T+.00325614)*T-.00780353)*T+.01504268)*T-
     1 .03655620)*T+.23498619)*T+1.25331414
      X=dMIN1(X,330.d0)
      XK1=dSQRT(X)*dEXP(-X)*P
      RETURN
      END


      FUNCTION BGAMA0(FNU,TAU1,TAU2,TEMP)
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
      BGAMA0=(TAU1/PI)*DEXP(TAU12+FNU/(2.d0*BKW*TEMP))*
     1 XK0(X)
	end

	SUBROUTINE MOMTS (IFIRST, labeta, RLO,RUP,FRAC,N,TEMP, G0)
c
c	If IFIRST = 1 then do only preliminary set-up; set ifirst=0
c	when calling moments,
c	G0 is the zeroth semiclassical moment

c	Recommended for very high temperatures: RLO=0.5, RUP=7.5 (Angstroms)
c	FRAC = 100., N = 400

C	COMPUTES 0-TH, SEMICLASSICAL TRANSLATIONAL SPECTRAL MOMENT
C	WORKS WITH "DOUBLE" POTENTIAL FUNCTIONS

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	CHARACTER*10 LGAS, LABEL, NOTE(2)
	CHARACTER*11 labeta

        COMMON /COMUNIC/ NOTE,LGAS,LABEL,FFMP,EPSP(2),RMN(2),
     1  LVALUE
        COMMON/EPOT/EPS(2),RM(2),FM,RMIN,RMAX,HHH
        DATA ANGAU / 0.529917706 /

C     RLO,RUP - LOWER AND UPPER LIMIT OF INTEGRATION IN ANG
C     N - NUMBER OF STEPS, FRAC REQUIRED BY DERIVAS, USUALLY SET=100.
C     FMP - REDUCED MASS IN HYDROGEN MASS, LVALUE = L
C     LGAS - 10 CHARACTER DESCRIPTION OF SYSTEM

	if (ifirst.eq.1) go to 222
	go to 333
cc222      BBB=BETACOM(0.D0)		! call BETA_JJ.for (Argument not needed - crashes compiler)
222      BBB=BETACOM	! call BETA_JJ.for
cc      VVV=VACOMUN(0.D0)			! call VA(R,k) (Argument not needed - crashes compiler)
      VVV=VACOMUN			! call VA(R,k)
      RM(1)=RMN(1)/ANGAU
      EPS(1)=EPSP(1)
      RM(2)=RMN(2)/ANGAU
      EPS(2)=EPSP(2)

10    FORMAT(3F10.5,I5)
11    FORMAT(F10.5)

      RLO1=RLO/ANGAU
      RUP1=RUP/ANGAU
      H=(RUP1-RLO1)/DFLOAT(N)
      T= TEMP * 1.380662D-23
	FMP=FFMP/2.d0
      	FMPI=2.*FMP*1.67265D-27
	return
333      CALL MOMENTS(RLO1,H,N,T,LVALUE,FRAC,FMPI,Gamma, gammap)
        G0=GAMMA
        RETURN
        END

      SUBROUTINE MOMENTS(RLO,H,N,T,LVALUE,FRAC,FMP,G0, g0p)
C     COMPUTES ZEROTH MOMENT (G0) FOR ANY L VALUE
c	G0p is the confirmation of the accuracy; lower order diff.
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ORDER/Y1P,Y2P,Y3P,Y4P
      EXTERNAL VTILDE,BETA
      DATA ANGAU, FOURPI, HBAR/ 0.52917706, 12.566370614359,
     1  1.05458876D-34/

      SP=0.
      SQ = 0.
      R=RLO-H
      F2=HBAR*HBAR/(12.*FMP*T*T*ANGAU**2*1.D-20)

        DO 10 I=1,N
        R=R+H
        DR=R/FRAC
	B = beta(R)
        CALL DERIVAS(VTILDE,R,DR,V,DV,DDV,Z3,Z4,JFLAG)
        DVP=Y1P
        DDVP=Y2P
        G = DEXP(-V/T)*R*R
        G2=G*(1.+F2*((DV/T-4./R)*DV-2.*DDV))
        G2P=G*(1.+F2*((DVP/T-4./R)*DVP-2.*DDVP))
        SP=SP+B*B*G2
        SQ=SQ+B*B*G2P
10      CONTINUE
      G0  = FOURPI*SP*H*(ANGAU*1.D-8)**3
      G0p = FOURPI*SQ*H*(ANGAU*1.D-8)**3
      RETURN
      END

      SUBROUTINE DERIVAS(FCT,X0,DX1, Y,Y1,Y2,Y3,Y4,IFLAG)
C     COMPUTES OTH..4TH DERIVATIVES OF Y=FCT(X) AT X=X0 USING FIVE
C     ABSCISSAE AND CENTRAL DIFFERENCES. FCT MUST BE DECLARED EXTERNAL
C     IN CALLING ROUTINE. NOTE THAT FCT(X) MUST BE DEFINED FOR
C     X0-2DX<=X<=X0+2DX
C     A TEST OF THE DEFINITION OF THE FIRST DERIVATIVE: IFLAG=1 OUTPUT
C     SEE ABRAMOWITZ AND STEGUN, P. 914
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ORDER/Y1P,Y2P,Y3P,Y4P
      DIMENSION  F(7)
9     X=X0
      H=DX1
      HSQ=H*H
      IFLAG=0
      DO 10 I=1,7
10    F(I)=FCT(X-DFLOAT(4-I)*H)
      Y=F(4)
      A=F(1)+F(7)
      A1=F(1)-F(7)
      B=F(2)+F(6)
      B1=F(2)-F(6)
      C=F(3)+F(5)
      C1=F(3)-F(5)

C     FOLLOWING ARE LOWER ORDER DERIVATIVES FOR TESTING
      Y1P=(B1-8.*C1)/(12.*H)
      Y2P=(16.*C-B-30.*Y)/(12.*HSQ)
      Y3P=(2.*C1-B1)/(2.*H*HSQ)
      Y4P=(B-4.*C+6.*Y)/(HSQ*HSQ)

C     THESE ARE HIGHEST ORDER DERIVATIVES
      Y1 = (-A1+9.*B1-45.*C1)/(60.*H)
      Y2 = (A-13.5*B+135.*C-245.*Y)/(90.*HSQ)
      Y3 = (A1-8.*B1+13.*C1)/(8.*HSQ*H)
      Y4 = (-A+12.*B-39.*C+56.*Y)/(6.*HSQ*HSQ)
      RETURN
      END

	FUNCTION VTILDE(R)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	K=0
	VTILDE = VA(R,K)
	RETURN
	END

      FUNCTION BETA(R)
      IMPLICIT double precision (A-H,O-Z)
      character*11 labeta
      CHARACTER*10 NOTE,LGAS,LABEL
      COMMON /BB1/ LABETA
      common /bb2/ J1, J1p, J2, J2p
	common/ ola/  XM0
      COMMON /COMUNIC/ NOTE(2),LGAS,LABEL,FMP,EPSP(2),RMP(2),LVALUE
      DATA DEBYE / 2.541769713d-18 /

      E(X) =  B0 * DEXP(A*X + B*X**2)
      E1(X) = B01 * DEXP(A1*X + B1*X**2)
      E2(X) = B02 * DEXP(A2*X + B2*X**2)

	if (labeta .eq. '<0|B0001|1>') go to 10
	if (labeta .eq. '<0|B2023|1>') go to 20
	if (labeta .eq. '<0|B0223|1>') go to 30
	if (labeta .eq. '<0|B2233|1>') go to 40
	if (labeta .eq. '<0|B2021|1>') go to 50
	if (labeta .eq. '<0|B0221|1>') go to 60
	if (labeta .eq. '<0|B2211|1>') go to 70
	if (labeta .eq. '<0|B4045|1>') go to 80
	if (labeta .eq. '<0|B0445|1>') go to 90
C	write(*, 555) labeta
555	format(' Unpredicted LABEL in BETA:', a11)
C	stop 999	! then, what term is that???

10	continue	! 0001 here
	LAMDA = 0
	LVALUE = 1
        LABEL = '00|0001|01'
	N = 7
	BN = -38.20D0
	B0 = 0.652d-3
	A = -1.521D0
	B = -0.033D0
	N1=1 ! really 0 but does not matter
	BN1=0.0d0
	B01 = 0.883d-6
	A1=-1.517D0
	B1 = -0.063D0
	N2=7
	BN2  =0.161d-1
	B02  =-0.318d-7
	A2  =-2.65D0
	B2 = -0.232D0
	go to 100
20	continue	! 2023 here
        LAMDA = 2
	LVALUE =3
	LABEL = '00|2023|01'
	N =4
	BN =-0.622D0
	B0 =0.811d-5
	A =-2.268D0
	B =-0.074D0
	N1= 4
	BN1= -0.152d-1
	B01= -0.645d-6
	A1= -1.307D0
	B1= -0.044D0
	N2= 4
	BN2= 0.148d-1
	B02= 0.667d-6
	A2= -1.339D0
	B2= -0.036D0
	go to 100
30	continue	! 0223
	LAMDA = 2
	LVALUE =3
	LABEL = '00|0223|01'
	N = 4
	BN = 0.853D0
	B0 = 0.118d-3
	A =  -1.479D0
	B  = -0.026D0
	N1 = 4
	BN1 = 0.154d-1
	B01 = 0.883d-6
	A1 = -1.402D0
	B1= -0.035D0
	N2 =  4
	BN2 = -0.149d-1
	B02 = -0.616d-6
	A2 = -1.311D0
	B2 = -0.036D0
	go to 100
40	continue	! 2233 here
        LAMDA = 2
	LVALUE  =3
        LABEL   = '00|2233|01'
	N  =4
	BN  =0.163D0
	B0  =-0.107d-4
	A  =-1.507D0
	B  =-0.028D0
        N1  =4
	BN1  =0.244d-2
	B01  =-0.122d-6
	A1  =-1.557D0
	B1 =-0.019D0
	N2  =4
	BN2  =-0.228d-2
	B02  =0.105d-6
	A2  =-1.576D0
	B2  =-0.016D0
	go to 100
50	continue ! 2021 here
        LAMDA =2
        LVALUE   = 1
        LABEL   = '00|2021|01'
         N  = 7
	BN  =15.5D0
	B0  =-0.148d-4
	A  =-1.219D0
	N1  = 1 ! really 0
	BN1  = 0.0D0
	B01  =0.152d-5
	A1  =-1.597D0
	B1  = -0.014D0
	N2  = 1
	BN2  =0.d0
	B02  =-0.150d-5
	A2  =-1.617D0
	B2  = -0.017D0
	go to 100
60	continue	! 0221 here
	LAMDA = 2
	LVALUE = 1
	LABEL =  '00|0221|01'
	N = 7
	BN = -0.161d2
	B0 = -0.153d-3
	A = -1.806D0
	B =-0.130D0
	N1 = 1 	! really 0
	BN1 = 0.d0
	B01 = -.19d-5
	A1 = -1.598D0
	B1 =-0.018D0
	N2 = 1
	BN2 = 0.d0
	B02 = 0.141d-5
	A2 = -1.629D0
	B2 = -0.017D0
	go to 100
70	continue	! 2211
	LAMDA = 2
	LVALUE =1
	LABEL = '00|2211|01'
	N =  7
	BN =  -0.105d1
	B0 =  0.862d-5
	A =  -1.431D0
	B = -0.002D0
	N1 =  1
	BN1 =  0.0d0
	B01 =  0.572d-7
	A1 =  -1.665D0
	B1 =-0.070D0
	N2 =  7
	BN2 =  0.534d-2
	B02 =  -0.679d-7
	A2 =  -1.569D0
	B2 = -0.022D0
	go to 100
80	continue	! 4045 here
	LAMDA =4
	LVALUE = 5
	LABEL ='00|4045|01'
	N = 6
	BN = -0.779D0
	B0 = 0.316d-6
	A = -3.755D0
	B = -0.341D0
	N1 =6
	BN1 = -0.234d-1
	B01 = 0.
	A1 = 0.
	B1= 0.
	N2 = 6
	BN2 = 0.176d-1
	B02 = 0.99d-7
	A2 = -1.91D0
	B2 = -0.332D0
	go to 100
90	continue	! 0445 here
	LAMDA = 4
	LVALUE = 5
	LABEL = '00|0445|01'
	N =  6
	BN =  0.220d1
	B0 =  0.295d-4
	A =  -1.430D0
	B =-0.131D0
	N1 =  6
	BN1 =  0.311d-1
	B01 =  0.
	A1 =  0.
	B1 = 0.
	N2 =  6
	BN2 =  -0.17d-1
	B02 =  -0.602d-7
	A2 =  -2.227D0
	B2 =-.497D0

100      BETAjj = (E1(R-6.D0)+ BN1/R**n1) * dfloat(J2*(J2+1))
     1      + (E2(R-6.D0) + BN2/R**N2) * dfloat(J2p*(J2p+1))
      BETA00 = E(R-6.D0) + BN/R**N
      beta = (betajj + beta00) * debye
      RETURN
      ENTRY BETACOM 		! the same for all terms
      BETA = 0.D0
      RETURN
      END

      FUNCTION VA(R,K)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*10 NOTE(2), LGAS, LABEL
      COMMON /EPOT/ EPS(2),RV(2),FM,RMIN,RMAX,HQ
      COMMON /COMUNIC/ NOTE,LGAS,LABEL,FMP,EPSP(2),RM(2),LVALUE
      DIMENSION A(2), GAMA(2), ALPHA(2), D(2), RMIN1(2), E(2),
     1  C6(2), C8(2), C10(2), C6S(2), C8S(2), C10S(2)

      K1=K+1
	X = R/RV(K1)

	V_REP = A(K1) * X**GAMA(K1)  * DEXP(- ALPHA(K1)*X)
	F = 1.D0
	IF(X .LE. D(K1) ) F = DEXP(- ( D(K1)/X-1.)**2)
	V_DISP = F * (C6S(K1)/X**6 + C8S(K1)/X**8 + C10S(K1)/X**10)
	FF = V_REP + V_DISP
      VA = FF * EPS(K1)
      RETURN

      ENTRY VACOMUN

        FACTOR = (4.803250D-10)**2/(0.52916607D-8) *(1.D-7)
        FMP=2.D0
        LGAS =  'H2-H2 GAS '

C	POTENTIAL H2(v=0) - H2(v=0)
      NOTE(1)='00 - 11 H2'
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

C	POTENTIAL H2(v=1) - H2(v=0)
      NOTE(2)='00 - 11 H2'
      RM(2)  = 6.5*0.52917706D0
      EPSP(2)= 0.11658D-3 * FACTOR
	RMIN1(2)=6.5D0
	E(2)=   0.11658D-3
	A(2) =  0.38000D7
	GAMA(2)=2.52498
	ALPHA(2)=14.76805
	D(2)=     1.03796D0
	C6(2)=   -12.396D0
	C8(2)=   -290.878D0
	C10(2)=  -4370.587D0

	DO 10 I=1,2
	C6S(I) = C6(I) /( E(I) * RMIN1(I)**6 )
	C8S(I) = C8(I) /( E(I) * RMIN1(I)**8 )
10	C10S(I) = C10(I) /( E(I) * RMIN1(I)**10 )
      VA=0.D0
      RETURN
      END

      FUNCTION CLEBSQR1 (L,LAMBDA,LP)
C     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)
C     FOR INTEGER ARGUMENTS ONLY
C     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.
C
      IMPLICIT double precision (A-H,O-Z)
      FC=DFLOAT(2*LP+1)
      GO TO 10
C
      ENTRY THREEJ2A
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
   10 CLEBSQR1=0.
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
   30 P=FC*F*FCTL1(LAMBDA+L-LP)*FCTL1(LAMBDA+LP-L)
      CLEBSQR1=P/(FCTL1((LAMBDA+L-LP)/2)*FCTL1((LAMBDA+LP-L)/2))**2
      RETURN
      END

      FUNCTION FCTL1 (N)
C     Simple FACTORIALS program
      IMPLICIT double precision (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      P(Z)=((((-2.294720936D-4)/Z-(2.681327160D-3))/Z+(3.472222222D-3))/
     1Z+(8.3333333333D-2))/Z+1.
      FCTL1=1
      IF (N.LE.1) RETURN
      IF (N.GT.15) GO TO 20
      J=1
      DO 10 I=2,N
   10 J=J*I
      FCTL1=DFLOAT(J)
      RETURN
   20 Z=DFLOAT(N+1)
      FCTL1=(DEXP(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)
      RETURN
      END

      FUNCTION BGAMA1(FNU,TAU1,TAU2,TEMP)
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
      BGAMA1=(TAU1/PI)*dEXP(AAA+FNU/(2.d0*BKW*TEMP))*
     1 XK1(X)/DENOM
      RETURN
      END

      FUNCTION BGAMA_H(FNU,tau ,TEMP)
c	"K1" one parameter model lineshape (14-Jul-1989 )
c	with Egelstaff time; unique for "H" function
C     NORMALIZATION SUCH THAT 0TH MOMENT EQUALS UNITY
C     FNU IS THE FREQUENCY IN CM-1; TEMP IN KELVIN.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA TWOPIC,BKW,HBOK,PI/1.883651568E11,.6950304256,
     1 7.638280918D-12, 3.141592654/

      TAU3=DSQRT( 3./2.*TAU**2 + (HBOK/(TEMP*2.))**2 )
      OMEGA=TWOPIC*FNU
      ARG = DABS(OMEGA)*TAU3
	IF(ARG.EQ.0.D0) GO TO 5
      BGAMA_H =1./PI * DEXP(FNU/(2.*BKW*TEMP))* XK1(ARG)
     1      * (TAU/TAU3)**2 * TAU * (3./2.)**(3./2.)
	RETURN
5	BGAMA_H =1./PI * DEXP(FNU/(2.*BKW*TEMP))
     1      * (TAU/TAU3)**2 * TAU * (3./2.)**(3./2.)
C	LIMIT K1(X) (X-->0) = 1/X; SO
C	LIMIT X*K1(X) (X-->0) = X/X = 1. EXACT
c	G1 = (beta *hbar)/(tau**2)
c	G2 = (5./3.)*beta*beta*hbar*hbar/(tau**4) + 2./(tau**2)
c	where 	beta = 1./(boltzk * temp) in cgs units
      RETURN
      END

	SUBROUTINE RANGE(IBGAMA)
	implicit double precision (a-h,o-z)
      COMMON/YMAX/YMAXXR, YMAXXL
	COMMON/UNCOR/ G0, TAU1, TAU2
	COMMON/CORR/ omega1, omega2, tauH
	common/temp/temp

c	choosing YMAXX here:

	if(ibgama.eq.1) gmax =g0 * bgama1(0.d0,tau1,tau2,temp)
	if(ibgama.eq.0) gmax =g0 * bgama0(0.d0,tau1,tau2,temp)

	STEP = 100.D0
	XS = 100.D0
	omegapl = omega1 + omega2
	omegamn = omega2 - omega1
4441	if(ibgama.eq.1) gx = g0 * bgama1(xs,tau1,tau2,temp)
	if(ibgama.eq.0) gx =g0 * bgama0(xs,tau1,tau2,temp)

	gy =  g0/2. * ( bgama_h(xs- omegapl,tauH,temp) -
     1	bgama_h(xs - omegamn, tauH, temp) )
	GLOWR = GX + GY
	XS = XS + STEP
	if(xs.gt.5000. ) go to 6565
	IF(GLOWR. GT. GMAX/5000.D0) GO TO 4441
6565	YMAXXR = XS-STEP
	if(ymaxxr.lt.1000.) ymaxxr = 2000.

	XS = -100.D0
	if(ibgama.eq.1) gx =g0 * bgama1(0.d0,tau1,tau2,temp)
	if(ibgama.eq.0) gx =g0 * bgama0(0.d0,tau1,tau2,temp)
6441	XS = XS - STEP
	if(xs.le.-5000.) go to 7766
	if(ibgama.eq.1) gx =g0 * bgama1(xs,tau1,tau2,temp)
	if(ibgama.eq.0) gx =g0 * bgama0(xs,tau1,tau2,temp)

	gy =  g0/2. * ( bgama_h(xs- omegapl,tauH,temp) -
     1	bgama_h(xs - omegamn, tauH, temp) )
	if((gx+gy).le.0.d0) go to 7766
	GLOWL = GX + GY
	IF(GLOWL. GT. GMAX/5000.D0) GO TO 6441
7766	YMAXXL = XS + STEP
	return
	end
