      real function rayleigha(v,P,T)
C     **************************************************************
C     Function to evaluate the Rayleigh scattering cross section using 
C     Atmospheric Radiation Liou (1980), p 79.
C     Modified with new data from Handbook of Chemistry and Physics 1996
C
C     Pat Irwin   8/6/98
C
C     Input variables
C	v	real	wavenumber (cm-1)
C	P	real	pressure (atm)
C	T	real	temperature (K)
C
C     Output variable
C	rayleigha	real	scattering cross section (cm2)
C
C     **************************************************************
      real P,T,k,nu2,Ns,P1,T1,xcorr
      real c1,c2,c3,v,lambda,delta,faniso,pi,m1
      parameter (c1=8342.12,c2=2406030.0,c3=15997.0)
      parameter (delta = 0.035,pi=3.1415927)
      parameter (k=1.37971e-23)

      lambda = 1e4/v			! lambda in um
      nu2 = 1.0/(lambda*lambda)		! convenient parameter

C     calculate the air's refractive index
      m1 = c1 + c2/(130.0 - nu2) + c3/(38.9 - nu2)
      m1 = m1*1e-8

      P1 = P*1.013e5			! Convert P to Pa
      T1 = T - 273.15			! Convert T to centigrade

C     Extra correction. From Handbook?
      xcorr = P1*(1 + P1*(61.3-T1)*1E-10)/(96095.4*(1+0.003661*T1))
      m1 = m1*xcorr

C     calculate number density of air molecules (m-3)
      Ns = P1/(k*T)
      lambda = lambda*1e-6		! Convert lambda to m
      x = Ns*lambda*lambda

      temp = 8*(pi**3)*(m1*(m1 + 2.0))**2
 
C     Calculate scattering x-section. faniso corrects for anisotropy.
C     Factor of 1e4 converts from m2 to cm2

      rayleigha = temp*1e4*faniso(delta)/(3*(x**2))

      return

      end


      real function faniso(delta)
C     *****************************************************************
C     Function to calculate anisotropic correction
C
C     *****************************************************************

      real delta

      faniso = (6.0+3.0*delta)/(6.0 - 7.0*delta)

      return
      end
