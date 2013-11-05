      real function hamming(FWHM,k)
C     ***********************************************************************************
C     Function to compute Hamming Function from parameters defined in
C     Wolfram MathWorld
C
C     Adapted from grid_hamm.pro, written by Nick Teanby
C
C     Input parameters
C 	FWHM	real	FWHM of Hamming function
C	k	real	offset from centre of Hamming Function in same 
C			units as the FWHM
C    
C     Output parameter
C	hamming	real	Hamming Function weight
C
C     Pat Irwin		5/11/13
C
C     ***********************************************************************************
      implicit none
      real pi,FWHM,a,k,stest
      parameter (pi=3.1415927)
      
      a = 0.9077/fwhm
      stest = 1.0 - 4.0*a*a*k*k

      if (abs(stest).gt.1.0e-6) then
        if (k.ne.0.0) then
         hamming=a*(1.08-0.64*a*a*k*k)*sin(2.0*pi*a*k)/
     &         ((1.0-4.0*a*a*k*k)*(2.0*pi*a*k))
        else
         hamming=a*1.08
        endif
      else
C       sin(..... converges to 0.5 close to (1.0-4.0*a*a*k*k)=0
        hamming=a*(1.08-0.64*a*a*k*k)*0.5
      endif

      return

      end

