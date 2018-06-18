      real function benreuven(P,T,nu,nu0,LLQ,frac)
C     ****************************************************************
C     Function to calculate the Ben Reuven lineshape for the microwave
C	based on Devaraj et al. (2014).
C
C     Input parameters:
C	P	REAL	Pressure (atm)
C	T	REAL	Temperature (K)
C	FRAC	REAL	Gas vmr (usually NH3)
C	LLQ	CHARACTER*15	Lower state local quanta

C	IMPORTANT:  Experimental procedure
c	Currently hardcoded for partial pressures of H2, He for Jupiter.
c	The NH3 abundance comes from FRAC.
C	Unclear over whether conversion from GHz to cm-1 needed in all cases.
C
C	L.N.Fletcher 18/06/2018
C     ****************************************************************
      implicit none
      DOUBLE PRECISION nu,nu0
      CHARACTER*15 LLQ
      REAL P,T,c
      real ph2,phe,pnh3,fh2,fhe,fnh3,frac
      real pi,gamma,sum
      real delta,zeta
      real n1,n2,d1,d2
      real gamma0,gammah2,gammahe,gammanh3,zetah2,zetahe,zetanh3
      INTEGER J,K
      parameter (pi=3.1415927)
      parameter (c=299792458.)
      parameter (gammah2=1.7947, gammahe=0.75, gammanh3=0.7719)
      parameter (zetah2=1.2031, zetahe=0.3, zetanh3=0.5620)
      parameter (fh2=0.863, fhe=0.134)
      
      fnh3=frac
      
c     gamma0
c     Self-broadening linewidths of the inversion transitions of 
c	ammonia in MHz/Torr
      READ(LLQ(1:2),*)J
      READ(LLQ(3:4),*)K
      gamma0=25.923*K*(J*(J+1))**0.5
      
      ph2=fh2*P
      phe=fhe*P
      pnh3=fnh3*P

      
c     Linewidth from summing contributions from diff gases in GHz
c	[Exponents from Tab 3 of Deveraj et al. 2014)
      gamma=gammah2*ph2*(300/T)**0.8357 + 
     1   gammahe*phe*(300/T)**(2/3) +
     2	 gammanh3*gamma0*pnh3*(295/T)
     
c     Coupling parameter from summing contributions from diff gases in GHz
c	[Exponents from Tab 3 of Deveraj et al. 2014)
      zeta=zetah2*ph2*(300/T)**0.8610 + 
     1   zetahe*phe*(300/T)**(2/3) +
     2	 zetanh3*gamma0*pnh3*(295/T)**0.6206     
         
	 
c     Pressure shift parameter in GHz
      delta=-0.0404*gamma
      
      
c     Convert these parameters from GHz to cm-1 (LNF: do we need to do this??)
      gamma=1e4/(1e6*c/(gamma*1e9))
      zeta=1e4/(1e6*c/(zeta*1e9))
      delta=1e4/(1e6*c/(delta*1e9))
      
      

c     Calculate the Ben Reuven shape as numerator/denominator:

      n1=(gamma-zeta)*nu*nu
      n2=(gamma+zeta)*((nu0+delta)**2 + gamma**2 - zeta**2)
      d1=(nu**2 - (nu0+delta)**2 - gamma**2 + zeta**2)**2
      d2=(4*nu**2)*gamma**2
      
      sum=(n1+n2)/(d1+d2)

      benreuven = (2/pi)*sum*(nu/nu0)*(nu/nu0)
      

      return

      end
