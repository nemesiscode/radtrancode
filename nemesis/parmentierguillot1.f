      SUBROUTINE PARMENTIERGUILLOT1(IPLANET,LATITUDE,NPRO,
     1 PIN,HIN,ALPHA,BETA,KIR,GAMMAV1,GAMMAV2,TSTAR,RSTAR,SDIST,TINT,
     2 TOUT,GRADTOUT)
C     **************************************************************   
C     Subroutine to compute double grey analytic TP profile of Parmentier and
C     Guillot (2014) and Line et al. (2013)

C     Input variables
C	IPLANET	INTEGER		Planet number (needed to compute gravity)
C	LATITUDE  REAL		Latitude (needed to compute gravity)
C	AMFORM	INTEGER		Profile type
C	NPROIN	INTEGER		Number of vertical levels
C	PIN(MAXPRO)	REAL	Pressure(atm)
C	HIN(MAXPRO)	REAL	Heights (km)
C	ALPHA		REAL	Parameter alpha. Weighting between two streams
C	BETA		REAL	Parameter beta. Albedo/emissivity weight.
C	KIR		REAL	Parameter K_IR. Thermal IR opacity parameter. (units: cm2/g)
C	GAMMAV1		REAL	Ratio of visible stream 1 opacity to thermal opacity
C	GAMMAV2		REAL	Ratio of visible stream 2 opacity to thermal opacity
C	TSTAR		REAL	Star temperature (K)
C	RSTAR		REAL	Star radius (km)
C	SDIST		REAL	Distance of planet from star (km).
C	TINT		REAL	T_internal (K)
C	
C
C     Output variables
C	TOUT(MAXPRO)	REAL	Output TP profile
C       GRADT(MAXPRO,5) REAL	Gradient of computed TP profile with
C 				respect to ALPHA,BETA,KIR,GAMMAV1 and 
C				GAMMAV2 respectively.
C
C     	Pat Irwin	25/1/21	Original
C			
C     **************************************************************   
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE '../radtran/includes/constdef.f'

      REAL LATITUDE,TOUT(MAXPRO)
      REAL PIN(MAXPRO),HIN(MAXPRO),ZETA,GRADTOUT(MAXPRO,5)
      REAL ALPHA,BETA,KIR,GAMMAV1,GAMMAV2,TAU,TINT
      INTEGER IPLANET,I,NPRO
      REAL RADIUS,G,TSTAR,RSTAR,SDIST,TIRR,C1,CX,X,G1
      REAL DCXDBETA,ZETA1,ZETA2,DZ1DTAU,DZ1DGA
      REAL DZ2DTAU,DZ2DGA,DXDALPHA,DXDBETA,DXDTAU
      REAL DXDKIR,DXDG1,DXDG2
      CHARACTER*8 PNAME


      TIRR=BETA*SQRT(0.5*RSTAR/SDIST)*TSTAR

      C1 = 3.0*(TINT**4)/4.0
      CX = 3.0*(TIRR**4)/4.0
      DCXDBETA = 3.0*(TIRR**3)*TIRR/BETA

      
C     Calc gravity at bottom of atmosphere (decide later if need at 1 bar, etc)
      CALL NEWGRAV(IPLANET,LATITUDE,HIN(1),RADIUS,G,PNAME)


C      print*,'G = (m s-2)',G
C     Test gravity: log(g)(cm s-2)=3.341
C      G1 = (10**3.341)/100.0
C      print*,'G1 = ',G1
C      G=G1

      DO I=1,NPRO

C      P/G first converted to SI units by multiplication by 1.013E5
C      P/G then converted from SI (kg m-2) to CGI (g cm-2) by dividing by 10.   
       TAU=KIR*(PIN(I)*1.013E4)/G

       X = C1*(2.0/3.0 + TAU)
       CALL CALCZETA(GAMMAV1,TAU,ZETA1,DZ1DTAU,DZ1DGA)
       CALL CALCZETA(GAMMAV2,TAU,ZETA2,DZ2DTAU,DZ2DGA)

       X = X + CX*((1.0-ALPHA)*ZETA1 + ALPHA*ZETA2)

       TOUT(I)=X**0.25  
       G1 = 0.25*X**(-0.75)

       DXDALPHA = -CX*ZETA1 +CX*ZETA2
       DXDBETA = ((1.0-ALPHA)*ZETA1+ALPHA*ZETA2)*DCXDBETA
       DXDTAU = CX*((1.0-ALPHA)*DZ1DTAU+ALPHA*DZ2DTAU)
       DXDTAU = DXDTAU+C1
       DXDKIR = DXDTAU*TAU/KIR
       DXDG1 = CX*(1.0-ALPHA)*DZ1DGA
       DXDG2 = CX*ALPHA*DZ2DGA

       GRADTOUT(I,1)=DXDALPHA*G1
       GRADTOUT(I,2)=DXDBETA*G1
       GRADTOUT(I,3)=DXDKIR*G1
       GRADTOUT(I,4)=DXDG1*G1
       GRADTOUT(I,5)=DXDG2*G1

      ENDDO


      RETURN

      END


      SUBROUTINE CALCZETA(GAMMA,TAU, ZETA, DZETADTAU, DZETADGA)
C     *********************************************************
C     Subroutine to calculate the zeta function of the double grey 
C     analytic TP profile of Parmentier and Guillot (2014) and 
C     Line et al. (2013)
C
C     Input variables
C	GAMMA	REAL	GAMMA parameter
C	TAU	REAL	TAU parameter
C    
C     Output variables
C       ZETA	REAL	Computed Zeta function
C	DZETADTAU REAL	Analytical partial derivative of ZETA 
C			   w.r.t. TAU
C	DZETADGA REAL	Analytical partial derivative of ZETA 
C			   w.r.t. GA
C
C     Pat Irwin		28/1/21
C
C     *********************************************************
      REAL X,EXPINT,GAMMA,TAU,Y,GRAD,SARG,C1,C2,C3,C4
      DOUBLE PRECISION ARG


      C0 = 2.0/3.0
      C1 = C0/GAMMA
      C2 = C1/GAMMA
      ARG=DBLE(GAMMA)*DBLE(TAU)
      SARG=GAMMA*TAU
      X=C0 + C1*(1.0 + (0.5*SARG-1)*EXP(-SARG))
      CALL E2(ARG,Y,GRAD)
      X=X+C0*GAMMA*(1-0.5*TAU**2)*Y

      ZETA=X
      
      DZETADTAU=0.0
      DZETADTAU=DZETADTAU+C1*(-(0.5*SARG-1)*GAMMA*EXP(-SARG))
      DZETADTAU=DZETADTAU+C1*0.5*GAMMA*EXP(-SARG)
      DZETADTAU=DZETADTAU+C0*GAMMA*(1-0.5*TAU**2)*GRAD*GAMMA
      DZETADTAU=DZETADTAU-C0*GAMMA*Y*TAU

      DZETADGA = C1*(-(0.5*SARG-1)*EXP(-SARG)*TAU + 
     &   0.5*TAU*EXP(-SARG))
      DZETADGA=DZETADGA-C2*(1.0 + (0.5*SARG-1)*EXP(-SARG))
      C3 = C0*(1-0.5*TAU**2)
      DZETADGA=DZETADGA+C3*(GAMMA*GRAD*TAU + Y)

      RETURN
 
      END


      SUBROUTINE E2(XIN,Y,GRAD)
C     ********************************************************
C     Routine for calculating the E2 exponential function and gradient
C     E2 function previously computed with IDL expint.pro procedure and
C     here tabulated to provide easy and quick interpolation
C
C     Input variables
C	XIN	REAL	Input X value
C
C     Output variables
C	Y	REAL	Interpolated E2 value
C	GRAD	REAL	Interpolated gradient d_E2/d_X
C
C     Pat Irwin		26/1/21 Original
C
C     ********************************************************
      INTEGER I,N
      PARAMETER(N=100)
      REAL YL(N),Y,GRAD,Z(N)
      DOUBLE PRECISION X(N),XIN,Z1,X1,X2,FX,G1
      DATA (YL(I),I=1,100) /0.00000,0.00000,0.00000,0.00000,0.00000,
     & 0.00000,0.00000,0.00000,0.00000,0.00000,-2.58860e-08,
     & -2.58860e-08, -2.58860e-08, -2.58860e-08, -5.17719e-08,
     & -5.17719e-08, -7.76579e-08, -7.76579e-08, -1.03544e-07,
     & -1.55316e-07, -2.07088e-07, -2.58860e-07, -3.36518e-07,
     & -4.40062e-07, -5.43606e-07, -7.24808e-07, -9.31896e-07,
     & -1.21664e-06, -1.55316e-06, -2.01911e-06, -2.61449e-06,
     & -3.39108e-06, -4.40064e-06, -5.66906e-06, -7.32579e-06,
     & -9.47437e-06, -1.22442e-05, -1.58166e-05, -2.04245e-05,
     & -2.63527e-05, -3.39637e-05, -4.37754e-05, -5.64092e-05,
     & -7.26162e-05, -9.34584e-05, -0.000120179, -0.000154463,
     & -0.000198384, -0.000254637, -0.000326597, -0.000418571,
     & -0.000536015, -0.000685872, -0.000876814, -0.00111995,
     & -0.00142903, -0.00182161, -0.00231951, -0.00295015,
     & -0.00374779, -0.00475518, -0.00602551, -0.00762473,
     & -0.00963471, -0.0121566, -0.0153151, -0.0192639,
     & -0.0241918, -0.0303303, -0.0379631, -0.0474377, -0.0591792,
     & -0.0737073, -0.0916590, -0.113815, -0.141135, -0.174802,
     & -0.216279, -0.267385, -0.330392, -0.408156, -0.504287,
     & -0.623379, -0.771314, -0.955672, -1.18628, -1.47592, -1.84134,
     & -2.30451, -2.89434, -3.64897, -4.61874, -5.87009, -7.49083,
     & -9.59684, -12.3411, -15.9257, -20.6167, -26.7657, -34.8358/

      DO I=1,100 
       Z(I)=-10.0+12.0*(I-1)/100.0
       X(I)=10.0**Z(I)
      ENDDO
 
      Z1=LOG10(XIN)
      IF(Z1.LT.-10.0)THEN
       Y=1.0
       GRAD=-20.0
       RETURN
      ENDIF
      IF(Z1.GT.1.89)THEN
       Y=0.0
       GRAD=0.0
       RETURN
      ENDIF

      I=1+INT((Z1+10.0)/0.12)

      IF(I.EQ.N)I=N-1

      X1=X(I)
      X2=X(I+1)
      FX=(XIN-X1)/(X2-X1)
  
      G1 = (YL(I+1)-YL(I))/(X2-X1)

      YLINT = (1.0-FX)*YL(I)+FX*YL(I+1)
      Y=10.0**YLINT

      GRAD=SNGL(Y*LOG(10.0)*G1)
C     Numerical overflow fix.
      IF(XIN.LT.1e-8)GRAD=-20.0

      RETURN

      END
