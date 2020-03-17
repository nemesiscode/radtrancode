      real function rayleighv(v,P,T)
C     **************************************************************
C     Function to evaluate the Rayleigh scattering cross section for
C     CO2 air using data Allen (1976) Astrophysical Quantities
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
     
      real :: v,lambda
      real, parameter :: C=8.8E-28		! provided by B. Bezard
 
      lambda = 10000/v			! lambda in microns
      lambda = lambda**4		
            
      rayleighv = C/lambda


      return

      end

