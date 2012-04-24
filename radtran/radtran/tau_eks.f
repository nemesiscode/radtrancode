      real function tau_eks(knu0,delad,y0,EL,SFB,qrot,P,T,U,q)
C     $Id:
C     *********************************************************************
C
C     Version 1.1
C     Calculate the transmission of a path using the the formulation of the
C     Goody-Voigt random band model given in Strong (1992) D.Phil.Thesis
C
C     Integration of the function over x is performed in two stages using 
C     Simpsons rule with 101 points over the ranges:
C	1)	x = 0 to 30*y
C 	2)	x = 30*y to 3030*y
C     This is the same method used in tran_eks.f by K.Strong.
C
C     ----------------------------------------------------------------------
C     Input Variables
C
C	knu0	real	Absorption coefficient at 296K (*1e20 (molec/cm2)**-1)
C	delad	real	line spacing / alpha_D0 (K**-0.5)
C	y0	real	alpha_L0 / alpha_D0 (K**-0.5)
C	EL	real	Lower state energy (cm-1)
C	SFB	real	Self-to-foreign broadening parameter
C	qrot	real	Rotational Partition Function Coefficient
C	P	real	Mean pressure of path (atm)
C	T	real	Mean temperature of path (K)
C	U	real	Absorber amount in path (*1e20 (molec/cm2))
C	q	real	Fractional abundance of active gas 
C
C     Output variable
C
C 	tau_eks	real	Transmission of path
C     ------------------------------------------------------------------------
C
C	Pat Irwin	26/9/94
C
C     *********************************************************************
      implicit none
      real U,delad,T,EL,knu0,y0,SFB,P,q,qrot
      integer IMOD,LPAR
      real T1,TRAN_EKS
      real CONDIT(6),PAR(6)

C      print*,'tau_eks knu0,delad,y0,EL,SFB,qrot,P,T,U,q'
C      print*,knu0,delad,y0,EL,SFB,qrot,P,T,U,q

      if(knu0.eq.0.or.U.eq.0.)then
       tau_eks=0.
      else

       CONDIT(1)=U
       CONDIT(2)=P
       CONDIT(3)=T
       CONDIT(4)=Q

       LPAR=6
       IMOD=14
       PAR(1)=KNU0
       PAR(2)=DELAD
       PAR(3)=Y0*(q + (1.0-q)/SFB)
       PAR(4)=EL
       PAR(5)=0.0
       PAR(6)=0.0

       T1=TRAN_EKS(IMOD,CONDIT,LPAR,PAR)


C       print*,'TAU_EKS = ',-LOG(T1) 

       tau_eks = -LOG(T1)

      endif

      RETURN

      END

