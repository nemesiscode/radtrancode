      real function rayleighv(v,P,T)
C     **************************************************************
C     Function to evaluate the Rayleigh scattering cross section for
C     CO2 air using data from Ignatiev et al. (1997)
C
C     Con Tsang   13/5/08
C
C     Input variables
C	v	real	wavenumber (cm-1)
C	P	real	pressure (atm)
C	T	real	temperature (K)
C
C     Output variable
C	rayleighv	real	scattering cross section (cm2)
C
C     **************************************************************
     
      real :: v,lambda,lambda2,lambda4,x,sigma
      real, parameter :: C1=156.63,C2=189.94
      real, parameter :: nu1=0.965,nu2=0.035
      real, parameter :: A1=17.904,A2=11.041
      real, parameter :: kboltz=1.38064852e-23
 
      lambda = 10000/v			! lambda in microns
      lambda2 = lambda**2
      lambda4 = lambda2**2		
            
      x = A1*nu1/(C1-1.0/lambda2)^2
      x = x + A2*nu2/(C2-1.0/lambda2)^2
      
      rayleighv = kboltz*x/(10*lambda4)

      return

      end

