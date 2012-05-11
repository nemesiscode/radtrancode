C---------------------------------------------------------------------
	DOUBLE PRECISION FUNCTION VOIGTF_EKS(X,Y)
C--------------------------------------------------------------------
C  
C  This function evaluates the Voigt profile as
C
C  F(V) = Y/PI int{ exp(-T*T) / [Y*Y + (X-T)*(X-T)] dT
C	  integrated from - infinity to + infinity,
C
C  given the standard Voigt inputs X and Y.
C
C  Copied from [CALCUTT.LBL.SOURCE.GENLBL]DRAYSN.FOR by EKS, 21jan91.
C  Was known as FUNCTION SHAPE(X,Y) in SBC's code.
C  EKS renamed it VOIGTF_EKS to agree with the function name used in 
C  the spectral fitting programs.
C
C  Previously copied from Bob Wells by SBC, feb88.
C  Was known as FUNCTION DRAYSN(X,Y) in Bob's original code.
C
C  Since the Voigt profile actually has a factor of (PI)^(-1.5), not just 
C  (PI)^(-1) as defined above, SBC had to divide the outputs by sqrt(PI).
C  Comparing his version with Bob's original, each of the four 
C  outputs has been divided by a factor of sqrt(PI) in the code
C  below.  To do this, SBC added the following parameterisation:
C
	PARAMETER (PI=0.3183099,SQRTPI=0.5641896)
C
C  Original code, with SBC's modifications, follows ...
C
C---------------------------------------------------------------------
C
C THIS FUNCTION FOR THE EVALUATION OF THE VOIGT PROFILE
C IS TAKEN FROM DRAYSON JQSRT 1976
C
C *** ROUTINE COMPUTES THE VOIGT FUNCTION: Y/PI* INTEGRAL FROM ***
C *** - TO + INFINITY OF EXP (-T*T)/(Y*Y+(X-T)*(X-T))DT        ***

      DOUBLE PRECISION B(22),RI(15),XN(15),YN(15),D0(25),D1(25),
     & D2(25),D3(25),D4(25),HN(25),XX(3),HH(3),NBY2(19),C(21)

      DOUBLE PRECISION C0,DX,H,U,UU,V,VV,X,Y,Y2

      LOGICAL TRU

      DATA B(1),B(2)/0.0,0.7093602E-7/,XN/10.0,9.0,8.0,8.0,7.0,6.0,
     &     5.0,4.0,7*3.0/,H/0.201/,XX/.5246476,1.65068,.7071068/,
     &     HH/.2562121,.02588268,.2820948/,NBY2/9.5,9.0,8.5,8.0,7.5,
     &     7.0,6.5,6.0,5.5,5.0,4.5,4.0,3.5,3.0,2.5,2.0,1.5,1.0,0.5/
     &     C/.7093602E-7,-.2518434E-6,.8566874E-6,-.2787638E-5,
     &     .866074E-5,-.2565551E-4,.7228775E-4,-.1933631E-3,.489952E-3,
     &     -.1173267E-2,.2648762E-2,-.5623190E-2,.1119601E-1,
     &     -.2084976E-1,.3621573E-1,-.5851412E-1,.8770816E-1,-.121664,
     &     .15584,-.184,.2/,TRU/.FALSE./
     &      ,YN/3*0.6,0.5,2*0.4,4*0.3,1.0,0.9,0.8,2*0.7/

      IF(TRU)GOTO104

C *** REGION I. COMPUTE DAWSON'S FUNCTION AT MESH POINTS ***
      TRU=.TRUE.
      DO101I=1,15
  101 RI(I)=-I/2.
      DO103I=1,25
      HN(I)=H*(I-.5)
      C0=4.*HN(I)*HN(I)/25.-2.0
      DO102J=2,21
  102 B(J+1)=C0*B(J)-B(J-1)+C(J)
      D0(I)=HN(I)*(B(22)-B(21))/5.
      D1(I)=1.0-2.0*HN(I)*D0(I)
      D2(I)=(HN(I)*D1(I)+D0(I))/RI(2)
      D3(I)=(HN(I)*D2(I)+D1(I))/RI(3)
  103 D4(I)=(HN(I)*D3(I)+D2(I))/RI(4)
  104 IF(X-5.0)105,112,112
  105 IF(Y-1.0)110,110,106
  106 IF(X.GT.1.85*(3.6-Y))GOTO112

C *** REGION II: CONTINUED FRACTION. COMPUTE NUMBER OF TERMS NEEDED ***
      IF(Y.LT.1.45)GOTO107
      I=INT(Y*Y)
      GOTO108
  107 I=INT(11.0*Y)
  108 J=INT(X+X+1.85)
      MAX=INT(XN(J)*YN(I)+.46)
      MIN=MIN0(16,21-2*MAX)

C *** EVALUATE CONTINUED FRACTION ***
      UU=Y
      VV=X
      DO109J=MIN,19
      U=NBY2(J)/(UU*UU+VV*VV)
      UU=Y+U*UU
 109  VV=X-U*VV
      VOIGTF_EKS=UU/(UU*UU+VV*VV)*PI
      RETURN
  110 Y2=Y*Y
      IF(X+Y.GE.5.0)GOTO113

C *** REGION I. COMPUTE DAWSON'S FUNCTION AT X FROM TAYLOR SERIES ***
      N=INT(X/H)
      DX=X-HN(N+1)
      U=(((D4(N+1)*DX+D3(N+1))*DX+D2(N+1))*DX+D1(N+1))*DX+D0(N+1)
      V=1.-2.0*X*U

C *** TAYLOR SERIES EXPANSION ABOUT Y=0.0 ***
      VV=EXP(Y2-X*X)*COS(2.0*X*Y)/1.128379-Y*V
      UU=-Y
      MAX=INT(5.0+(12.5-X)*0.8*Y)
      DO111I=2,MAX,2
      U=(X*V+U)/RI(I)
      V=(X*U+V)/RI(I+1)
      UU=-UU*Y2
  111 VV=VV+V*UU
      VOIGTF_EKS=1.128379*VV*SQRTPI
      RETURN
  112 Y2=Y*Y
      IF(Y.LT.11.-.6875*X)GOTO113

C *** REGION IIIB: 2 POINT GAUSS-HERMITE QUADRATURE ***
      U=X-XX(3)
      V=X+XX(3)
      VOIGTF_EKS=SQRTPI*Y*(HH(3)/(Y2+U*U)+HH(3)/(Y2+V*V))
      RETURN

C *** REGION III A : 4 POINT GAUSS HERMITE QUADRATURE ***
  113 U=X-XX(1)
      V=X+XX(1)
      UU=X-XX(2)
      VV=X+XX(2)
      VOIGTF_EKS=SQRTPI*(Y*(HH(1)/(Y2+U*U)+HH(1)/(Y2+V*V)
     #  +HH(2)/(Y2+UU*UU)+HH(2)/(Y2+VV*VV)))

      RETURN
      END
