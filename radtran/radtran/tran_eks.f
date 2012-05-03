C-------------------------------------------------------------------------------
	FUNCTION TRAN_EKS(IMOD,CONDIT,LPAR,PAR)
C-------------------------------------------------------------------------------
C
C_TITLE: TRAN_EKS computes the transmission given the LPAR parameters PAR 
C	 for the model IMOD.
C
C_ARGS:	 IMOD		: INTEGER   : model number 
C	 CONDIT(MAXCOND) 	    : array of lab conditions
C          CONDIT(1)	: REAL      : number of molec/cm2 (*1.0E-24) 
C          CONDIT(2)	: REAL      : pressure (bar)
C          CONDIT(3)	: REAL      : temperature (Kelvin)
C          CONDIT(4)	: REAL      : mixing ratio of the absorbing gas
C        LPAR		: INTEGER   : number of parameters to initialise 
C        PAR(NPAR)	: REAL      : array of parameter values
C	 MAXCOND	: PARAMETER : maximum number of lab conditions
C
C_DESCR: Returns the transmission calculated from model parameters for 
C        various models.  The models implemented are:
C
C 	 IMOD =	1: Goody model with Lorentz lines.
C	
C		2: Goody model with Voigt lines, fitting only parameters 1,2.
C	
C  		3: Goody model with Voigt lines, fitting all 3 parameters.
C
C		4: Zachor-King model with 3 parameters.
C
C		5: Gibson-Pierluissi (fudged Zachor-King) with 5 parameters.
C
C		6: Smith model - simple empirical exponential.
C
C		7: Strong Limit.
C
C		8: Weak Limit of the Goody-Lorentz and Malkmus-Lorentz models.
C
C  		9: Goody model with Lorentz lines (as for IMOD=1) but with 
C     		   temperature dependence of K added as parameter 3.
C
C  	       10: Goody model with Voigt lines (as for IMOD=3) but with 
C      		   temperature dependence of K added as parameter 4.
C
C  	       11: Weak Limit (as for IMOD=8) but with temperature 
C		   dependence of K added as parameter 2.
C
C  	       12: Goody model with Voigt lines and variable PAR(3) 
C      		   (as for IMOD=3) but a foreign-broadening parameter added 
C		   as parameter 4.
C
C  	       13: The Goody model with Voigt lines (as for IMOD=3) but with 
C		   the temperature dependence of K added as parameter 4 (as 
C		   in IMOD=10), and a foreign-broadening parameter added as 
C		   parameter 5 (as in IMOD=12).
C
C	       14: Goody-Voigt model with K(T) and parameters 5 and 6 added 
C		   to account for far wing continuum contributions. 
C
C 	       15: Malkmus model with Lorentz lines.
C
C 	       16: Malkmus model with Voigt lines.
C
C 	       17: Malkmus model with Lorentz lines (as for IMOD=15) but with 
C     		   temperature dependence of K added as parameter 3.
C
C  	       18: Malkmus model with Voigt lines (as for IMOD=16) but with 
C      		   temperature dependence of K added as parameter 4.
C
C	       19: Malkmus-Voigt model with K(T) and parameters 5 and 6 added 
C		   to account for far wing continuum contributions. 
C
C	       20: Smith model (as for IMOD=6) with T dependence of K added 
C		   as parameter 5.
C
C	       21: Zachor-King model (as for IMOD=4) with Lorentz lineshape 
C	           explicitly used and T dependence of K added as parameter 4.
C
C	       22: Gibson-Pierluissi (as for IMOD=5) with Lorentz lineshape 
C	           explicitly used and T dependence of K added as parameter 6.
C
C  	       23: Goody model with Lorentz lines and k(T) (as for IMOD=9) 
C		   but with but a foreign-broadening parameter added 
C		   as parameter 4.
C
C  	       24: Malkmus model with Lorentz lines and k(T) (as for IMOD=17) 
C		   but with but a foreign-broadening parameter added 
C		   as parameter 4.
C
C  	       25: Malkmus model with Voigt lines and k(T) (as for IMOD=18) 
C		   but with but a foreign-broadening parameter added 
C		   as parameter 4.
C
C_FILES:  
C
C_CALLS: FUNCTION VOIGTF_EKS (a derivative of RJW's DRAYSN.FOR) returns the 
C	  		     Voigt line shape given standard parameters x and y
C        ROUTINE STRONG_EKS  computes the strong limit for Zachor-King model,
C         		     calls Numerical Recipes routines GAMMP and GAMMLN
C        ROUTINE STRONGLOR_EKS  computes the strong limit for Zachor-King 
C			     model explicitly using the Lorentz lineshape, 
C         		     calls Numerical Recipes routines GAMMP and GAMMLN
C
C_BUGS:  Must be LINKed with USER:[wells]recipes/lib on OXNIMS/ISAMS
C
C_HIST:  30sep85 SBC ORIGINAL VERSION
C        12feb86 SBC added code for weak and strong limits and the "1.0-f" 
C	 	 case, and added temperature dependence to the computation 
C		 of X in the Goody model
C	 16jan91 EKS began modifications from TRAN.FOR to TRAN_EKS.FOR
C	 25sep91 EKS began final revisions for thesis
C	 05oct94 EKS Modified from TRAN_EKS. Limits of integration for the
C		 Goody-Voigt Model were incorrect. Limits changed as described
C		 in code
C
C_END:
C
C---------------------------------------------------------------------------

	PARAMETER (MAXCOND=6)
	PARAMETER (PI=3.1415926)
	PARAMETER (P0=1.01325)
	PARAMETER (T0=296.0)

	INTEGER IMOD,LPAR

	REAL PAR(LPAR),WT(0:100),CONDIT(MAXCOND)
	REAL U,P,T,C
	REAL K,KU,KUDEL,PS,Y
	REAL DELX,DELX2,TAU,TAU2,X,X2,VT,VT2
	REAL ST,SQST,WEAKT,SQWT
	REAL TERM,P3

        DOUBLE PRECISION VOIGTF_EKS

C  	WT(0:100) contains the weighting factors for 101 point integration 
C  	using Simpson's Rule, used in evaluating the Goody-Voigt transmission.

	DATA WT/1.,4.,2.,4.,2.,4.,2.,4.,2.,4.,
     #	2.,4.,2.,4.,2.,4.,2.,4.,2.,4.,
     #	2.,4.,2.,4.,2.,4.,2.,4.,2.,4.,
     #	2.,4.,2.,4.,2.,4.,2.,4.,2.,4.,
     #	2.,4.,2.,4.,2.,4.,2.,4.,2.,4.,
     #	2.,4.,2.,4.,2.,4.,2.,4.,2.,4.,
     #	2.,4.,2.,4.,2.,4.,2.,4.,2.,4.,
     #	2.,4.,2.,4.,2.,4.,2.,4.,2.,4.,
     #	2.,4.,2.,4.,2.,4.,2.,4.,2.,4.,
     #	2.,4.,2.,4.,2.,4.,2.,4.,2.,4.,1./

C---------------------------------------------------------------------------

C  Rename the lab conditions to simplify the equations for TRAN_EKS.

	U = CONDIT(1)
	P = CONDIT(2)
	T = CONDIT(3)
	C = CONDIT(4)

C---------------------------------------------------------------------------
C
C  1: Goody model with Lorentz lines
C
C  The transmission is defined as:
C
C             { 		-1	            	    }
C  TRAN = exp {---------------------------------------------}
C	      {     [    1	             1	         ]  }
C	      { sqrt[ --------  +  --------------------- ]  }
C	      {     [ (K*U)**2     pi*K*Y*U*P*sqrt(T0/T) ]  }
C
C  where
C  PAR(1) = K(T0) = S0/d = absorption coefficient at T=T0
C  PAR(2) = y = a0(self)/P0*d = pressure coefficient
C
C  so
C             {		             -1	                                   }
C  TRAN = exp {------------------------------------------------------------}
C  	      {     [       1	                       1	        ]  }
C	      { sqrt[ -------------  +  ------------------------------- ]  }
C	      {     [ (PAR(1)*U)**2     pi*PAR(1)*PAR(2)*U*P*sqrt(T0/T) ]  }

	IF (IMOD.EQ.1) THEN
	  
	  KU  = ABS(PAR(1)) * U
	  PS  = PI * ABS(PAR(2)) * P * SQRT(T0/T)
	  TAU = KU * SQRT(PS / (KU + PS))
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  2: Goody model with Voigt lines, fitting only parameters 1,2.
C  3: Goody model with Voigt lines, fitting all 3 parameters.
C
C  These models are the same, except that PAR(3), and hence the line width, 
C  is held constant for IMOD=2.
C
C  The transmission is defined as:
C
C  	     {	       (        K * U * V(X,Y)             )                   }
C  TRAN = exp{ -2 * int(-----------------------------------) dX [limits 0-inf] }
C	     {         ( 1 + [K * U * d * V(X,Y) / q(v,T)] )                   }
C
C  where
C  K 	    = K(T0) = S0/d = absorption coefficient at T=T0
C  V(X,Y)   = Voigt line shape function
C  q(v,T)   = q0 * sqrt(T) = Doppler width 
C  a(v,T,P) = a0(self) * (P/P0) * sqrt(T0/T) = self-broadened Lorentz width 
C  S0       = mean line strength
C  d 	    = mean line spacing
C  X 	    = (v - v0) / q(v,T)
C  Y 	    = a(v,T,P) / q(v,T)
C
C  PAR(1) = K(T0) = S0/d
C  PAR(2) = d/q0 
C  PAR(3) = a0(self)/q0
C
C  so
C  	     {	      (        PAR(1) * U * V(X,Y)         )                   }
C  TRAN = exp{-2 * int(------------------------------------) dX [limits 0-inf] }
C	     {        ( 1 + PAR(1)*PAR(2)*U*V(X,Y)/sqrt(T) )                   }
C  with
C  Y = PAR(3) * (P/P0) * sqrt(T0)/T
C
C  V(X,Y) is the Voigt line shape function, called VOIGTF_EKS below. 
C  This function is Bob Wells' DRAYSN.FOR, corrected for the sqrt(pi) 
C  term by SBC.
C
C  Numerical Integration: ...................................................
C
C  The integral is evaluated using Simpson's Rule for 101 points.  
C  The code below actually applies Simpson's Rule twice, over two 
C  different ranges, producing terms TAU and TAU2.  
C
C  TAU is the summation for X = 0 to 30*Y*R1 in increments of 0.30*Y*R1. 
C  TAU2 is the summation for X2 = 30*Y*R1 to 3030*Y*R1 in increments of 30*Y*R1.
C  The choice of X/Y=30*R1 as the cross-over between the high and low
C  resolution calculations was chosen by EKS, because by this point
C  the integrand is less than 1% of its maximum. (Note: SBC found that the
C  Voigt function is < 1% when X/Y = 30Y.  
C
C  The integrand of the function may be written:
C       B(X,Y) =         V(X,Y)
C                   ---------------
C                     1 + KUDEL*V(X,Y)
C
C  Put C = 0.01*Bmax = 0.01*B(0,Y). When the function B is 1% of max
C  then:
C               V1 = V(X,Y) = C/(1-KUDEL*C)
C
C  In the wings, the Voigt line shape looks like a Lorentz line. Thus:
C               (nu - nu0) is proportional to 1/sqrt(k)
C  Now when X=30Y, V(X,Y) is always less than 1% of Vmax = V(0,Y). Thus
C  putting D = V(0,Y), when the factor KUDEL is large, the integration limit
C  must be scaled to:
C               X_int = 30*Y*SQRT(0.01*D/V1)
C
C  and the integration is then more appropriate. Further checks to ensure
C  numerical accuracy should be done.
C
C  Perform Simpson's Rule integration over the high and low resolution 
C  ranges as part of the same loop.
C
C  ..........................................................................

	ELSE IF (IMOD.EQ.2.OR.IMOD.EQ.3) THEN
	  KU 	= 2.0 * ABS(PAR(1)) * U
	  KUDEL = 0.5 * KU * ABS(PAR(2)) / SQRT(T)
	  Y 	= ABS(PAR(3)) * (P/P0) * SQRT(T0)/T

	  D1    = VOIGTF_EKS(0.,Y)
	  B1    = D1/(1.0 + KUDEL*D1)
	  C1    = 0.01*B1

	  V1    = C1 / (1. - C1*KUDEL)
	  R1    = SQRT(0.01*D1/V1)

	  NPTS 	= 100
	  DELX 	= 0.30 * Y * R1
	  DELX2 = 100.0 * DELX
	  TAU 	= 0.0
	  TAU2 	= 0.0

	  DO 20 I=0,NPTS
	    X 	= FLOAT(I) * DELX
	    VT 	= VOIGTF_EKS(X,Y)
	    X2 	= FLOAT(NPTS) * DELX + FLOAT(I+1) * DELX2
	    VT2 = VOIGTF_EKS(X2,Y)
	    TAU = TAU  + WT(I) * VT  / (1.0 + KUDEL * VT)
	    TAU2= TAU2 + WT(I) * VT2 / (1.0 + KUDEL * VT2)
20	  CONTINUE

	  TAU = (KU * (TAU*DELX + TAU2*DELX2)) / 3.0
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  4: Zachor-King model with Lorentz lineshape (3 parameters).
C
C  The transmission is defined as:
C
C             {		       -1	                }
C  TRAN = exp {-----------------------------------------}
C	      {     [    1	             1	     ]  }
C	      { sqrt[ --------  +  ----------------- ]  }
C	      {     [ (K*U)**2     (-ln(Tstrong))**2 ]  }
C
C  where
C  Tstrong is the strong limit for the transmission as defined by King (1964)
C          = 1 - P{n,x},  P(n,x) being the incomplete gamma function
C	   = 1 - P{n, [n * GAM(n) * sqrt(2*C*U*P/PI) ]**(1/n) }
C	   = 1 - P{n, [ n * GAM(n) * sqrt(4*K(T)*y*U*P*sqrt(T0/T)) ]**(1/n) }
C  
C  PAR(1) = K(T0) = S0/d
C  PAR(2) = y = a0(self) / P0*d 
C  PAR(3) = n = adjustable parameter defined in the model
C  so 
C  Tstrong = 1 - P{ PAR(3), 
C  [PAR(3) * GAM(PAR(3)) * sqrt(4*K(T)*PAR(2)*U*P*sqrt(T0/T)) ]**(1/PAR(3)) }
C
C  Subroutine STRONGLOR_EKS calculates the strong limit of the Zachor-King
C  model for this explicit Lorentz lineshape formulation using Numerical 
C  Recipes routines to evaluate the gamma and incomplete gamma functions.

	ELSE IF (IMOD.EQ.4) THEN
	  K = ABS(PAR(1))
	  Y = ABS(PAR(2))
C 	  CALL STRONGLOR_EKS(PAR(3),Y,K,T,U,P,ST,SQST)
	  TAU = 1.0 / SQRT( 1.0/(K*U)**2 + ST )
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  5: Gibson-Pierluissi (fudged Zachor-King) model with Lorentz lineshape
C     (5 parameters).
C
C  The transmission is defined as:
C
C          {				-1				    }
C  TRAN=exp{----------------------------------------------------------------}
C	   {     [    1	             PAR(4)                PAR(5)        ]  }
C	   { sqrt[ --------  +  -----------------  +  ------------------ ]  }
C	   {     [ (K*U)**2     (-ln(Tstrong))**2     K*U*(-ln(Tstrong)) ]  }
C
C  PAR(1) = K(T0) = S0/d
C  PAR(2) = y = a0(self) / P0*d 
C  PAR(3) = n = adjustable parameter defined in the model
C  PAR(4) = coefficient B_s > 0, in the Gibson-Pierluissi model 
C  PAR(5) = -K * coefficient B_ws, in the Gibson-Pierluissi model
C	  = M, in the Zachor-King formulation
C  and
C  Tstrong = 1 - P{ PAR(3), 
C  [PAR(3) * GAM(PAR(3)) * sqrt(4*K(T)*PAR(2)*U*P*sqrt(T0/T)) ]**(1/PAR(3)) }
C
C  Subroutine STRONGLOR_EKS calculates the strong limit of the Zachor-King
C  model for this explicit Lorentz lineshape formulation using Numerical 
C  Recipes routines to evaluate the gamma and incomplete gamma functions.
C
C  Parameters 4 and 5 provide complete flexibility in weighting
C  the three terms of the full Zachor-King model.  Note that in the full
C  Zachor-King model, PAR(4)=1 and PAR(5)=M.
C
C  Must have ABS(PAR(5)) < SQRT(PAR(4)).

	ELSE IF (IMOD.EQ.5) THEN
	  K = ABS(PAR(1))
	  Y = ABS(PAR(2))
C 	  CALL STRONGLOR_EKS(PAR(3),Y,K,T,U,P,ST,SQST)
	  KU    = K * U
	  WEAKT = (1.0/KU)**2.0
	  SQWT  = SQRT(WEAKT)
	  IF (ABS(PAR(5)).GE.SQRT(ABS(PAR(4)))) THEN
	    PAR(5) = (PAR(5)/ABS(PAR(5))) * 0.9999 * SQRT(ABS(PAR(4)))
	  ENDIF
	  TERM  = SQRT(WEAKT + ABS(PAR(4))*ST + PAR(5)*SQWT*SQST)
	  TAU   = 1.0/TERM
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  6: Smith model - simple empirical exponential
C
C  The transmission is defined as the simple empirical exponential:
C
C  TRAN = exp{ - K * U**a * P**b * T**c }
C  where
C  PAR(1) = K = absorption coefficient
C  PAR(2) = a = exponent of absorber amount
C  PAR(3) = b = exponent of pressure
C  PAR(4) = c = exponent of temperature
C  so
C  TRAN = exp{ - PAR(1) * U**PAR(2) * P**PAR(3) * T**PAR(4) }
C
C  The original code was
C	TAU=PAR(1)*U**PAR(2)*P**PAR(3)*T**PAR(4)
C  but this, and several variations of it, led to overflow problems.
C  Now using the form in Rodgers, 1976, p.28, and ABS(K).

	ELSE IF (IMOD.EQ.6) THEN
	  TERM = LOG(ABS(PAR(1))) + PAR(2)*LOG(U) + PAR(3)*LOG(P)
     #		 + PAR(4)*LOG(T)
	  IF (TERM.GT.70.0) TERM = 70.0
	  TAU = EXP(TERM)
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  7: The Strong limit.
C
C  The transmission is defined as the general expression for the strong limit 
C  of the random model (Rodgers, 1976, p.26), and is valid when U/P is large:
C
C  TRAN = exp{ - PAR(1) * sqrt(U*P) }
C  where
C  PAR(1) = a generalized absorption coefficient

	ELSE IF (IMOD.EQ.7) THEN
	  TAU = ABS(PAR(1)) * SQRT(U*P)
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  8: The Weak Limit
C
C  The transmission is defined the weak limit of the Goody-Lorentz and
C  Malkmus-Lorentz models, and is valid when U/P is small:
C
C  TRAN = exp{ - K(T0) * U }
C  where
C  PAR(1) = K(T0)

	ELSE IF (IMOD.EQ.8) THEN
	  TAU = ABS(PAR(1)) * U
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  9: The Goody model with Lorentz lines (as for IMOD=1) but with 
C     temperature dependence of K added as parameter 3.
C
C  The transmission is defined as:
C
C             {		          -1	                          }
C  TRAN = exp {---------------------------------------------------}
C	      {     [      1	             1		       ]  }
C	      { sqrt[ -----------  +  ------------------------ ]  }
C	      {     [ {K(T)*U}**2     pi*K(T)*Y*U*P*sqrt(T0/T) ]  }
C
C  where
C  K(T) = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C       = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(3) * {1/T0-1/T} )
C
C  PAR(1) = K(T0) = S0/d
C  PAR(2) = y = a0(self) / P0*d 
C  PAR(3) = E 

	ELSE IF (IMOD.EQ.9) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(3))*(T-T0)/(T0*T))
	  KU  = K * U 
	  PS  = PI * ABS(PAR(2)) * P * SQRT(T0/T)
	  TAU = KU * SQRT(PS / (KU + PS))
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  10: The Goody model with Voigt lines (as for IMOD=3) but with 
C      temperature dependence of K added as parameter 4.
C
C  The transmission is defined as:
C
C  	     {        (        K(T) * U * V(X,Y)           )                 }
C  TRAN = exp{-2 * int(------------------------------------)dX [limits 0-inf]}
C	     {        (1 + [K(T) * U * d * V(X,Y) / q(v,T)])                 }
C
C  where
C  K(T)   = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C  	  = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(4) * {1/T0-1/T} )
C  Y      = PAR(3) * (P/P0) * sqrt(T0) / T
C
C  PAR(1) = K(T0)
C  PAR(2) = d/q0 
C  PAR(3) = a0(self)/q0
C  PAR(4) = E 
C
C  See the notes above for IMOD=2 for full details of the numerical integration.

	ELSE IF (IMOD.EQ.10) THEN
	  NPTS 	= 100
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(4))*(T-T0)/(T0*T))
	  KU	= 2.0 * K * U
	  KUDEL = 0.5 * KU * ABS(PAR(2)) / SQRT(T)
	  Y	= ABS(PAR(3)) * (P/P0) * SQRT(T0)/T

	  D1    = VOIGTF_EKS(0.,Y)
	  B1    = D1/(1.0 + KUDEL*D1)
	  C1    = 0.01*B1

	  V1    = C1 / (1. - C1*KUDEL)
	  R1    = SQRT(0.01*D1/V1)

	  DELX 	= 0.30 * Y * R1
	  DELX2 = 100.0 * DELX
	  TAU 	= 0.0
	  TAU2 	= 0.0

	  DO 21 I=0,NPTS
	    X	= FLOAT(I) * DELX
	    VT	= VOIGTF_EKS(X,Y)
	    X2	= FLOAT(NPTS)*DELX + FLOAT(I+1)*DELX2
	    VT2	= VOIGTF_EKS(X2,Y)
	    TAU	= TAU  + WT(I) * VT  / (1.0 + KUDEL * VT)
	    TAU2= TAU2 + WT(I) * VT2 / (1.0 + KUDEL * VT2)
21	  CONTINUE

	  TAU = (KU * (TAU*DELX + TAU2*DELX2)) / 3.0
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  11: The Weak Limit (as for IMOD=8) but with temperature dependence 
C      of K added as parameter 2.
C
C  The transmission is defined as:
C  TRAN = exp{ -K(T) * U }
C
C  where
C  K(T) = K(T0) * (T0/T)**1.5 * exp( 1.439 * E * {1/T0-1/T} )
C
C  PAR(1) = K(T0)
C  PAR(2) = E 

	ELSE IF (IMOD.EQ.11) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(2))*(T-T0)/(T0*T))
	  TAU = K * U
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  12: The Goody model with Voigt lines and variable PAR(3) (as for IMOD=3) 
C      but a foreign-broadening parameter added as parameter 4.
C
C  The transmission is defined as:
C
C  	     {	       (        K * U * V(X,Y)             )                   }
C  TRAN = exp{ -2 * int(-----------------------------------) dX [limits 0-inf] }
C	     {         ( 1 + [K * U * d * V(X,Y) / q(v,T)] )                   }
C
C  where
C  a(v,T,P) = a0(self) * (P/P0) * sqrt(T0/T) * [C + (1-C)*a0(foreign)/a0(self)]
C	    = Lorentz width including both self and foreign broadening 
C  C 	    = mixing ratio
C  a0(self)/a0(foreign) = self-to-foreign broadening parameter ("SFB")
C
C  PAR(1)   = K(T0)
C  PAR(2)   = d/q0 
C  PAR(3)   = a0(self)/q0
C  PAR(4)   = a0(self)/a0(foreign) = self-to-foreign broadening parameter 
C
C  so
C  	     {	      (        PAR(1) * U * V(X,Y)         )                   }
C  TRAN = exp{-2 * int(------------------------------------) dX [limits 0-inf] }
C	     {        ( 1 + PAR(1)*PAR(2)*U*V(X,Y)/sqrt(T) )                   }
C  with
C  Y = a(v,T,P) / q(T,P)
C    = (P/P0) * sqrt(T0)/T * PAR(3) * [ C + (1-C)/PAR(4) ]
C
C  See the notes above for IMOD=2 for full details of the numerical integration.

	ELSE IF (IMOD.EQ.12) THEN
	  NPTS 	= 100
	  KU 	= 2.0 * ABS(PAR(1)) * U
	  KUDEL = 0.5 * KU * ABS(PAR(2)) / SQRT(T)
	  P3	= C + (1.0-C) / ABS(PAR(4))
	  Y	= (P/P0) * (SQRT(T0)/T) * ABS(PAR(3)) * P3

	  D1    = VOIGTF_EKS(0.,Y)
	  B1    = D1/(1.0 + KUDEL*D1)
	  C1    = 0.01*B1

	  V1    = C1 / (1. - C1*KUDEL)
	  R1    = SQRT(0.01*D1/V1)

	  DELX 	= 0.30 * Y * R1
	  DELX2 = 100.0 * DELX
	  TAU 	= 0.0
	  TAU2 	= 0.0

	  DO 420 I=0,NPTS
	    X 	= FLOAT(I) * DELX
	    VT 	= VOIGTF_EKS(X,Y)
	    X2 	= FLOAT(NPTS) * DELX + FLOAT(I+1) * DELX2
	    VT2 = VOIGTF_EKS(X2,Y)
	    TAU = TAU + WT(I) * VT  / (1.0 + KUDEL * VT)
	    TAU2= TAU2+ WT(I) * VT2 / (1.0 + KUDEL * VT2)
420	  CONTINUE

	  TAU = (KU * (TAU * DELX + TAU2 * DELX2)) / 3.0
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  13: The Goody model with Voigt lines (as for IMOD=3) but with the
C      temperature dependence of K added as parameter 4 (as in IMOD=10), and
C      a foreign-broadening parameter added as parameter 5 (as in IMOD=12).
C
C  The transmission is defined as:
C
C  	     {        (        K(T) * U * V(X,Y)           )                 }
C  TRAN = exp{-2 * int(------------------------------------)dX [limits 0-inf]}
C	     {        (1 + [K(T) * U * d * V(X,Y) / q(v,T)])                 }
C
C  where
C  K(T)     = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C           = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(4) * {1/T0-1/T} )
C  a(v,T,P) = a0(self) * (P/P0) * sqrt(T0/T) * [C + (1-C)*a0(foreign)/a0(self)]
C  Y = a(v,T,P) / q(T,P)
C    = (P/P0) * sqrt(T0)/T * PAR(3) * [ C + (1-C)/PAR(5) ]
C
C  PAR(1) = K(T0)
C  PAR(2) = d/q0 
C  PAR(3) = a0(self)/q0
C  PAR(4) = E 
C  PAR(5) = a0(self)/a0(foreign) self-to-foreign broadening parameter ("SFB")
C
C  See the notes above for IMOD=2 for full details of the numerical integration.

	ELSE IF (IMOD.EQ.13) THEN
	  NPTS	= 100
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(4))*(T-T0)/(T0*T))
	  KU	= 2.0 * K * U
	  KUDEL	= 0.5 * KU * ABS(PAR(2)) / SQRT(T)
	  P3	= C + (1.0-C) / ABS(PAR(5))
	  Y	= (P/P0) * (SQRT(T0)/T) * ABS(PAR(3)) * P3

	  D1    = VOIGTF_EKS(0.,Y)
	  B1    = D1/(1.0 + KUDEL*D1)
	  C1    = 0.01*B1

	  V1    = C1 / (1. - C1*KUDEL)
	  R1    = SQRT(0.01*D1/V1)

	  DELX	= 0.30 * Y * R1
	  DELX2	= 100.0 * DELX
	  TAU	= 0.0
	  TAU2	= 0.0

	  DO 421 I=0,NPTS
	    X	= FLOAT(I) * DELX
	    VT	= VOIGTF_EKS(X,Y)
	    X2	= FLOAT(NPTS) * DELX + FLOAT(I+1) * DELX2
	    VT2	= VOIGTF_EKS(X2,Y)
	    TAU	= TAU  + WT(I) * VT  / (1.0 + KUDEL * VT)
	    TAU2= TAU2 + WT(I) * VT2 / (1.0 + KUDEL * VT2)
421	  CONTINUE

	  TAU = (KU * (TAU * DELX + TAU2 * DELX2)) / 3.0
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  14: The Goody model with Voigt lines and K(T) (as for IMOD=10) but with 
C      a two-parameter far-wing continuum effect added.
C
C  The transmission is defined as:
C
C  	  {       	                 (     K(T) * U * V(X,Y)    )          }
C TRAN=exp{-C1*U*(P/P0)*(T/T0)**C2 -2*int(--------------------------)dX [0-inf]}
C	  {              	         (1+[K(T)*U*d*V(X,Y)/q(v,T)])          }
C  where
C
C  PAR(1) = K(T0)
C  PAR(2) = d/q0 
C  PAR(3) = a0(self)/q0
C  PAR(4) = E 
C  PAR(5) = C1 = coefficient of the far-wing continuum term
C  PAR(6) = C2 = exponent of temperature dependence in the far-wing term
C
C  so
C  	      {	                               
C  TRAN = exp { - PAR(5) * U * (P/P0) * (T/T0)**PAR(6) 
C	      {
C  	                 (        K(T) * U * V(X,Y)         )                  }
C  		- 2 * int(----------------------------------) dX limits[0-inf] }
C	                 ( 1 + K(Y)*PAR(2)*U*V(X,Y)/sqrt(T) )                  }
C  with
C  Y    = PAR(3) * (P/P0) * sqrt(T0) / T
C  K(T) = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(4) * {1/T0-1/T} )
C
C  See the notes above for IMOD=2 for full details of the numerical integration.

	ELSE IF (IMOD.EQ.14) THEN
	  NPTS 	= 100
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(4))*(T-T0)/(T0*T))
	  KU	= 2.0 * K * U
	  KUDEL	= 0.5 * KU * ABS(PAR(2)) / SQRT(T)
	  Y	= ABS(PAR(3)) * (P/P0) * SQRT(T0)/T

	  DELX	= 0.30 * Y
	  DELX2	= 100.0 * DELX

C          PRINT*,Y,DELX,DELX2

	  TAU	= 0.0
	  TAU2	= 0.0

	  DO 521 I=0,NPTS
	    X	= FLOAT(I) * DELX
	    VT	= SNGL(VOIGTF_EKS(DBLE(X),DBLE(Y)))
	    X2	= FLOAT(NPTS) * DELX + FLOAT(I+1) * DELX2
	    VT2	= SNGL(VOIGTF_EKS(DBLE(X2),DBLE(Y)))
	    TAU	= TAU  + WT(I) * VT  / (1.0 + KUDEL * VT)
	    TAU2= TAU2 + WT(I) * VT2 / (1.0 + KUDEL * VT2)
C            PRINT*,I,X,VT,TAU,X2,VT2,TAU2
521	  CONTINUE

C          PRINT*,KU,DELX,DELX2
	  TAU = (KU * (TAU * DELX + TAU2 * DELX2)) / 3.0
C          PRINT*,TAU
	  TAU = TAU + ABS(PAR(5)) * U * (P/P0) * (T/T0)**PAR(6)
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)


C          print*,PAR(1)*U,EXP(-TAU)

C---------------------------------------------------------------------------
C
C  15: Malkmus model with Lorentz lines.
C
C  The transmission is defined as:
C
C             { 	    - 2 * K * U               }
C  TRAN = exp {---------------------------------------}
C	      {     	[    	K * U * sqrt(T/T0)  ] }
C	      { 1 + sqrt[ 1 + --------------------- ] }
C	      {     	[ 	pi * y * P 	    ] }
C
C  where
C  PAR(1) = K(T0) 
C  PAR(2) = y 
C
C  so
C             { 	    - 2 * PAR(1) * U               }
C  TRAN = exp {--------------------------------------------}
C	      {     	[    	PAR(1) * U * sqrt(T/T0)  ] }
C	      { 1 + sqrt[ 1 + -------------------------- ] }
C	      {     	[ 	pi * PAR(2) * P     	 ] }

	ELSE IF (IMOD.EQ.15) THEN
	  KU   = ABS(PAR(1)) * U
	  TERM = 1.0 + KU * SQRT(T/T0) / (PI * ABS(PAR(2)) * P)
	  TAU  = (2.0 * KU) / (1.0 + SQRT(TERM))
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  16: Malkmus model with Voigt lines.
C
C  The transmission is defined as:
C
C  	     { -2	( 	     [     K * U * d	      )            }
C  TRAN = exp{ --- * int( q(v,T) * ln[ 1 + --------- * V(X,Y) ) dX [0-inf] }
C	     {  d       ( 	     [      q(v,T)            )            }
C
C  where
C  V(X,Y)   = Voigt line shape function
C  q(v,T)   = q0 * sqrt(T) = Doppler width 
C  a(v,T,P) = a0(self) * (P/P0) * sqrt(T0/T) = self-broadened Lorentz width 
C  S0       = mean line strength
C  d 	    = mean line spacing
C  X 	    = (v - v0) / q(v,T)
C  Y 	    = a(v,T,P) / q(v,T)
C
C  PAR(1)   = K(T0)
C  PAR(2)   = d/q0 
C  PAR(3)   = a0(self)/q0
C
C  so
C  	    {  -2	 ( 	      [    PAR(1)*PAR(2)*U*V(X,Y))             }
C TRAN = exp{------ * int(sqrt(T) * ln[1 + -----------------------) dX [0-inf] }
C	    {PAR(2)      ( 	      [    	   sqrt(T)        )            }
C  with
C  Y = PAR(3) * (P/P0) * sqrt(T0)/T
C
C  See the notes above for IMOD=2 for full details of the numerical integration.
C
C  Perform Simpson's Rule integration over the high and low resolution 
C  ranges as part of the same loop.

	ELSE IF (IMOD.EQ.16) THEN
	  TERM  = ABS(PAR(1)) * ABS(PAR(2)) * U / SQRT(T)
	  Y 	= ABS(PAR(3)) * (P/P0) * SQRT(T0)/T
	  NPTS 	= 100
	  DELX 	= 0.30 * Y
	  DELX2 = 100.0 * DELX
	  TAU 	= 0.0
	  TAU2 	= 0.0

	  DO 30 I=0,NPTS
	    X 	= FLOAT(I) * DELX
	    VT 	= VOIGTF_EKS(X,Y)
	    X2 	= FLOAT(NPTS) * DELX + FLOAT(I+1) * DELX2
	    VT2 = VOIGTF_EKS(X2,Y)
	    TAU = TAU  + WT(I) * LOG(1.0 + TERM * VT)
	    TAU2= TAU2 + WT(I) * LOG(1.0 + TERM * VT2)
30	  CONTINUE

	  TAU=((2.0*SQRT(T)/ABS(PAR(2)))*(TAU*DELX+TAU2*DELX2))/3.0
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C---------------------------------------------------------------------------
C
C  17: Malkmus model with Lorentz lines (as for IMOD=15) but with 
C      temperature dependence of K added as parameter 3.
C
C  The transmission is defined as:
C
C             { 	  - 2 * K(T) * U                }
C  TRAN = exp {-----------------------------------------}
C	      {     	[    	K(T) * U * sqrt(T/T0) ] }
C	      { 1 + sqrt[ 1 +   --------------------- ] }
C	      {     	[ 	  pi * y * P 	      ] }
C
C  where
C  K(T) = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C       = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(3) * {1/T0-1/T} )
C
C  PAR(1) = K(T0) = S0/d
C  PAR(2) = y = a0(self) / P0*d 
C  PAR(3) = E 

	ELSE IF (IMOD.EQ.17) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(3))*(T-T0)/(T0*T))
	  KU   = K * U
	  TERM = 1.0 + KU * SQRT(T/T0) / (PI * ABS(PAR(2)) * P)
	  TAU  = (2.0 * KU) / (1.0 + SQRT(TERM))
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  18: Malkmus model with Voigt lines (as for IMOD=16) but with 
C      temperature dependence of K added as parameter 4.
C
C  The transmission is defined as:
C
C  	     { -2	( 	     [     K(T) * U * d	         )            }
C  TRAN = exp{ --- * int( q(v,T) * ln[ 1 + ------------ * V(X,Y) ) dX [0-inf] }
C	     {  d       ( 	     [         q(v,T)            )            }
C
C  where
C  K(T)   = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C  	  = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(4) * {1/T0-1/T} )
C  Y      = PAR(3) * (P/P0) * sqrt(T0) / T
C
C  PAR(1) = K(T0)
C  PAR(2) = d/q0 
C  PAR(3) = a0(self)/q0
C  PAR(4) = E 
C
C  See the notes above for IMOD=2 for full details of the numerical integration.
C
C  Perform Simpson's Rule integration over the high and low resolution 
C  ranges as part of the same loop.

	ELSE IF (IMOD.EQ.18) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(4))*(T-T0)/(T0*T))
	  TERM  = K * ABS(PAR(2)) * U / SQRT(T)
	  Y 	= ABS(PAR(3)) * (P/P0) * SQRT(T0)/T
	  NPTS 	= 100
	  DELX 	= 0.30 * Y
	  DELX2 = 100.0 * DELX
	  TAU 	= 0.0
	  TAU2 	= 0.0

	  DO 130 I=0,NPTS
	    X 	= FLOAT(I) * DELX
	    VT 	= VOIGTF_EKS(X,Y)
	    X2 	= FLOAT(NPTS) * DELX + FLOAT(I+1) * DELX2
	    VT2 = VOIGTF_EKS(X2,Y)
	    TAU = TAU  + WT(I) * LOG(1.0 + TERM * VT)
	    TAU2= TAU2 + WT(I) * LOG(1.0 + TERM * VT2)
130	  CONTINUE

	  TAU=((2.0*SQRT(T)/ABS(PAR(2)))*(TAU*DELX+TAU2*DELX2))/3.0
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  19: The Malkmus model with Voigt lines and K(T) (as for IMOD=18) but with 
C      a two-parameter far-wing continuum effect added.
C
C  The transmission is defined as:
C
C  	  {                          2    ( 	    [  K(T)*U*d	      )     }
C TRAN=exp{-C1*U*(P/P0)*(T/T0)**C2 - -*int(q(v,T)*ln[1+--------*V(X,Y)) dX] }
C	  {                          d    ( 	    [  q(v,T)         )     }
C
C  where
C  K(T)   = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C  	  = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(4) * {1/T0-1/T} )
C  Y      = PAR(3) * (P/P0) * sqrt(T0) / T
C
C  PAR(1) = K(T0)
C  PAR(2) = d/q0 
C  PAR(3) = a0(self)/q0
C  PAR(4) = E 
C  PAR(5) = C1 = coefficient of the far-wing continuum term
C  PAR(6) = C2 = exponent of temperature dependence in the far-wing term
C
C  See the notes above for IMOD=2 for full details of the numerical integration.
C
C  Perform Simpson's Rule integration over the high and low resolution 
C  ranges as part of the same loop.

	ELSE IF (IMOD.EQ.19) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(4))*(T-T0)/(T0*T))
	  TERM  = K * ABS(PAR(2)) * U / SQRT(T)
	  Y 	= ABS(PAR(3)) * (P/P0) * SQRT(T0)/T
	  NPTS 	= 100
	  DELX 	= 0.30 * Y
	  DELX2 = 100.0 * DELX
	  TAU 	= 0.0
	  TAU2 	= 0.0

	  DO 131 I=0,NPTS
	    X 	= FLOAT(I) * DELX
	    VT 	= VOIGTF_EKS(X,Y)
	    X2 	= FLOAT(NPTS) * DELX + FLOAT(I+1) * DELX2
	    VT2 = VOIGTF_EKS(X2,Y)
	    TAU = TAU  + WT(I) * LOG(1.0 + TERM * VT)
	    TAU2= TAU2 + WT(I) * LOG(1.0 + TERM * VT2)
131	  CONTINUE

	  TAU=((2.0*SQRT(T)/ABS(PAR(2)))*(TAU*DELX+TAU2*DELX2))/3.0
	  TAU = TAU + ABS(PAR(5)) * U * (P/P0) * (T/T0)**PAR(6)
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  20: Smith model - simple empirical exponential (as for IMOD=6)
C      with temperature dependence of K added as parameter 5.
C
C  The transmission is defined as the simple empirical exponential:
C
C  TRAN = exp{ - K(T) * U**a * P**b * T**c }
C  where
C  PAR(1) = K = absorption coefficient
C  PAR(2) = a = exponent of absorber amount
C  PAR(3) = b = exponent of pressure
C  PAR(4) = c = exponent of temperature
C  PAR(5) = E = lower state energy
C  so
C  TRAN = exp{ - K(T) * U**PAR(2) * P**PAR(3) * T**PAR(4) }
C  where
C  K(T) = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C       = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(5) * {1/T0-1/T} )
C
	ELSE IF (IMOD.EQ.20) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(5))*(T-T0)/(T0*T))
	  TERM = LOG(K) + PAR(2)*LOG(U) + PAR(3)*LOG(P)
     #		 + PAR(4)*LOG(T)
	  IF (TERM.GT.70.0) TERM = 70.0
	  TAU = EXP(TERM)
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  21: Zachor-King model (as for IMOD=4) with Lorentz lineshape 
C      explicitly used and T dependence of K added as parameter 4.
C
C  The transmission is defined as:
C
C             {		          -1	                   }
C  TRAN = exp {--------------------------------------------}
C	      {     [      1	             1	        ]  }
C	      { sqrt[ -----------  +  ----------------- ]  }
C	      {     [ [K(T)*U]**2     [-ln(Tstrong)]**2 ]  }
C
C  where
C  Tstrong is the strong limit for the transmission as defined by King (1964)
C          = 1 - P{n,x},  P(n,x) being the incomplete gamma function
C	   = 1 - P{n, [n * GAM(n) * sqrt(2*C*U*P/PI) ]**(1/n) }
C	   = 1 - P{n, [ n * GAM(n) * sqrt(4*K(T)*y*U*P*sqrt(T0/T)) ]**(1/n) }
C  
C  PAR(1) = K(T0) = S0/d
C  PAR(2) = y = a0(self) / P0*d 
C  PAR(3) = n = adjustable parameter defined in the model
C  PAR(4) = E 
C  with
C  K(T) = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C       = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(4) * {1/T0-1/T} )
C
C  so 
C  Tstrong = 1 - P{ PAR(3), 
C  [PAR(3) * GAM(PAR(3)) * sqrt(4*K(T)*PAR(2)*U*P*sqrt(T0/T)) ]**(1/PAR(3)) }
C
C  Subroutine STRONGLOR_EKS calculates the strong limit of the Zachor-King
C  model for this explicit Lorentz lineshape formulation using Numerical 
C  Recipes routines to evaluate the gamma and incomplete gamma functions.

	ELSE IF (IMOD.EQ.21) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(4))*(T-T0)/(T0*T))
C 	  CALL STRONGLOR_EKS(PAR(3),PAR(2),K,T,U,P,ST,SQST)
	  TAU = 1.0 / SQRT( 1.0/(K*U)**2 + ST )
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  22: Gibson-Pierluissi (as for IMOD=5) with Lorentz lineshape 
C      explicitly used and T dependence of K added as parameter 6.
C
C  The transmission is defined as:
C
C          {				-1				    }
C  TRAN=exp{----------------------------------------------------------------}
C	   {    [     1	            PAR(4)                PAR(5)          ] }
C	   {sqrt[ ----------- + ----------------- + --------------------  ] }
C	   {    [ [K(T)*U]**2   (-ln(Tstrong))**2   K(T)*U*[-ln(Tstrong)] ] }
C  
C  PAR(1) = K(T0) = S0/d
C  PAR(2) = y = a0(self) / P0*d 
C  PAR(3) = n = adjustable parameter defined in the model
C  PAR(4) = coefficient B_s > 0, in the Gibson-Pierluissi model 
C  PAR(5) = -K * coefficient B_ws, in the Gibson-Pierluissi model
C	  = M, in the Zachor-King formulation
C  PAR(6) = E 
C  with
C  K(T) = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C       = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(6) * {1/T0-1/T} )
C  and
C  Tstrong = 1 - P{ PAR(3), 
C  [PAR(3) * GAM(PAR(3)) * sqrt(4*K(T)*PAR(2)*U*P*sqrt(T0/T)) ]**(1/PAR(3)) }
C
C  Subroutine STRONGLOR_EKS calculates the strong limit of the Zachor-King
C  model for this explicit Lorentz lineshape formulation using Numerical 
C  Recipes routines to evaluate the gamma and incomplete gamma functions.
C
C  Parameters 4 and 5 provide complete flexibility in weighting
C  the three terms of the full Zachor-King model.  Note that in the full
C  Zachor-King model, PAR(4)=1 and PAR(5)=M.
C
C  Must have ABS(PAR(5)) < SQRT(PAR(4)).

	ELSE IF (IMOD.EQ.22) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(6))*(T-T0)/(T0*T))
C 	  CALL STRONGLOR_EKS(PAR(3),PAR(2),K,T,U,P,ST,SQST)
	  KU    = K * U
	  WEAKT = (1.0/KU)**2.0
	  SQWT  = SQRT(WEAKT)
	  IF (ABS(PAR(5)).GE.SQRT(ABS(PAR(4)))) THEN
	    PAR(5) = (PAR(5)/ABS(PAR(5))) * 0.9999 * SQRT(ABS(PAR(4)))
	  ENDIF
	  TERM  = SQRT(WEAKT + ABS(PAR(4))*ST + PAR(5)*SQWT*SQST)
	  TAU   = 1.0/TERM
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  23: The Goody model with Lorentz lines and with temperature dependence 
C      of K (as for IMOD=9) but with foreign broadening added as parameter 4.
C
C  The transmission is defined as:
C
C            {		          -1	                                      }
C TRAN = exp {----------------------------------------------------------------}
C	     {     [      1	                       1                    ] }
C	     { sqrt[ -----------  +  -------------------------------------- ] }
C	     {     [ {K(T)*U}**2     pi*K(T)*Y*U*P*[C+(1-C)/SFB]*sqrt(T0/T) ] }
C
C  where
C  a(v,T,P) = a0(self) * (P/P0) * sqrt(T0/T) * [C + (1-C)*a0(foreign)/a0(self)]
C	    = Lorentz width including both self and foreign broadening 
C  C 	    = mixing ratio
C  a0(self)/a0(foreign) = self-to-foreign broadening parameter ("SFB")
C
C  K(T) = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C       = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(3) * {1/T0-1/T} )
C
C  PAR(1) = K(T0) = S0/d
C  PAR(2) = y = a0(self) / P0*d 
C  PAR(3) = E 
C  PAR(4) = a0(self)/a0(foreign) = self-to-foreign broadening parameter 

	ELSE IF (IMOD.EQ.23) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(3))*(T-T0)/(T0*T))
	  KU  = K * U 
	  P3  = C + (1.0-C) / ABS(PAR(4))
	  PS  = PI * ABS(PAR(2)) * P * SQRT(T0/T) * P3

	  TAU = KU * SQRT(PS / (KU + PS))
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C---------------------------------------------------------------------------
C
C  24: Malkmus model with Lorentz lines and k(T) (as for IMOD=17) but with 
C      foreign broadening added as parameter 4.
C
C  The transmission is defined as:
C
C             { 	  - 2 * K(T) * U                     }
C  TRAN = exp {----------------------------------------------}
C	      {     	[    	K(T) * U * sqrt(T/T0)      ] }
C	      { 1 + sqrt[ 1 +   -------------------------- ] }
C	      {     	[ 	pi * Y * P * [C+(1-C)/SFB] ] }
C
C  where
C  a(v,T,P) = a0(self) * (P/P0) * sqrt(T0/T) * [C + (1-C)*a0(foreign)/a0(self)]
C	    = Lorentz width including both self and foreign broadening 
C  C 	    = mixing ratio
C  a0(self)/a0(foreign) = self-to-foreign broadening parameter ("SFB")
C
C  K(T) = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C       = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(3) * {1/T0-1/T} )
C
C  PAR(1) = K(T0) = S0/d
C  PAR(2) = y = a0(self) / P0*d 
C  PAR(3) = E 
C  PAR(4) = a0(self)/a0(foreign) = self-to-foreign broadening parameter 

	ELSE IF (IMOD.EQ.24) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(3))*(T-T0)/(T0*T))
	  KU   = K * U
	  P3  = C + (1.0-C) / ABS(PAR(4))
	  TERM = 1.0 + KU * SQRT(T/T0) / (PI * ABS(PAR(2)) * P * P3)
	  TAU  = (2.0 * KU) / (1.0 + SQRT(TERM))
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C-------------------------------------------------------------------------
C
C  25: Malkmus model with Voigt lines and k(T) (as for IMOD=18) but with 
C      self-to-foreign broadening added as parameter 5.
C
C  The transmission is defined as:
C
C  	     { -2	( 	     [     K(T) * U * d	         )            }
C  TRAN = exp{ --- * int( q(v,T) * ln[ 1 + ------------ * V(X,Y) ) dX [0-inf] }
C	     {  d       ( 	     [         q(v,T)            )            }
C
C  where
C  K(T)   = K(T0)  * (T0/T)**1.5 * exp( 1.439 * E      * {1/T0-1/T} )
C  	  = PAR(1) * (T0/T)**1.5 * exp( 1.439 * PAR(4) * {1/T0-1/T} )
C
C  a(v,T,P) = a0(self) * (P/P0) * sqrt(T0/T) * [C + (1-C)*a0(foreign)/a0(self)]
C
C  Y = a(v,T,P) / q(T,P)
C    = (P/P0) * sqrt(T0)/T * PAR(3) * [ C + (1-C)/PAR(5) ]
C
C  PAR(1) = K(T0)
C  PAR(2) = d/q0 
C  PAR(3) = a0(self)/q0
C  PAR(4) = E 
C  PAR(5) = a0(self)/a0(foreign) self-to-foreign broadening parameter ("SFB")
C
C  See the notes above for IMOD=2 for full details of the numerical integration.
C
C  Perform Simpson's Rule integration over the high and low resolution 
C  ranges as part of the same loop.

	ELSE IF (IMOD.EQ.25) THEN
	  K = ABS(PAR(1))*(T0/T)**1.5*EXP(1.439*ABS(PAR(4))*(T-T0)/(T0*T))
	  TERM  = K * ABS(PAR(2)) * U / SQRT(T)
	  P3	= C + (1.0-C) / ABS(PAR(5))
	  Y	= (P/P0) * (SQRT(T0)/T) * ABS(PAR(3)) * P3
	  NPTS 	= 100
	  DELX 	= 0.30 * Y
	  DELX2 = 100.0 * DELX
	  TAU 	= 0.0
	  TAU2 	= 0.0

	  DO 630 I=0,NPTS
	    X 	= FLOAT(I) * DELX
	    VT 	= VOIGTF_EKS(X,Y)
	    X2 	= FLOAT(NPTS) * DELX + FLOAT(I+1) * DELX2
	    VT2 = VOIGTF_EKS(X2,Y)
	    TAU = TAU  + WT(I) * LOG(1.0 + TERM * VT)
	    TAU2= TAU2 + WT(I) * LOG(1.0 + TERM * VT2)
630	  CONTINUE

	  TAU=((2.0*SQRT(T)/ABS(PAR(2)))*(TAU*DELX+TAU2*DELX2))/3.0
	  IF (TAU.GT.70.0) TAU = 70.0
	  TRAN_EKS = EXP(-TAU)

C --------------------------------------------------------------------------

	ELSE
	  WRITE(5,1)
1	  FORMAT(' ERROR IN TRAN_EKS.FOR - MODEL NOT INCLUDED ')
	  STOP

	END IF

	RETURN
	END
