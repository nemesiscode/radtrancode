      real function hanning(FWHM,k)
C     ***********************************************************************************
C     Function to compute Hanning Function from parameters defined in
C     Wolfram MathWorld
C
C
C     Input parameters
C 	FWHM	real	FWHM of Hamming function
C	k	real	offset from centre of Hamming Function in same 
C			units as the FWHM
C    
C     Output parameter
C	hanning	real	Hanning Function weight
C
C     Pat Irwin		1/9/14
C
C     ***********************************************************************************
      implicit none
      real pi,FWHM,a,k,stest,arg,y
      parameter (pi=3.1415927)
      
      a = 1./fwhm
      stest = 1.0 - 4.0*a*a*k*k
      arg = 2*pi*k*a

      if(abs(arg).lt.1e-3)then
       y = 1 - (1.0/6.0)*arg**2 + (1.0/120.0)*arg**4
      else
       y = sin(arg)/arg
      endif

      if(abs(stest).le.1.0e-6)then
       hanning = 0.5
      else
       hanning = y/stest
      endif

      return

      end

