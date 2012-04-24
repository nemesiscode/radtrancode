      real function rayleighj(v,P,T)
C     **************************************************************
C     Function to evaluate the Rayleigh scattering cross section for
C     Jovian air using data from Allen (1976) Astrophysical Quantities
C
C     Pat Irwin   8/6/98
C
C     Input variables
C	v	real	wavenumber (cm-1)
C	P	real	pressure (atm)
C	T	real	temperature (K)
C
C     Output variable
C	rayleighj	real	scattering cross section (cm2)
C
C     **************************************************************
      real P,T,k,nu2,Ns,P1,T1,xcorr,alpha,N0
      real AH2,AHe,BH2,BHe,x,nH2,nHe,fH2,nAir,P0,T0
      real emH2,emHe,v,lambda,delta,faniso,pi,m1,e1
      parameter (AH2=13.58E-5,BH2 = 7.52E-3)
      parameter (AHe= 3.48E-5,BHe = 2.30E-3)
      parameter (pi=3.1415927,fH2 = 0.864)
      parameter (k=1.37971e-23)
      parameter (P0=1.013e5,T0=273.15)

      lambda = 1e-2/v			! lambda in m
      x = 1.0/(lambda*1e6)

      nH2 = AH2*(1.0+BH2*x*x)
      nHe = AHe*(1.0+BHe*x*x)

C     calculate the Jupiter air's refractive index at STP (Actually n-1)
      nAir = fH2*nH2 + (1-fH2)*nHe
 
C     H2,He Seem pretty isotropic to me?...Hence delta = 0.
      delta = 0.0

      temp = 32*(pi**3)*nAir**2
 
      N0 = P0/(k*T0)

      x = N0*lambda*lambda

C     Factor 1e4 converts to cm2

      rayleighj = temp*1e4*faniso(delta)/(3*(x**2))

      return

      end

