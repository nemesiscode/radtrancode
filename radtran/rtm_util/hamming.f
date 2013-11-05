      real function hamming(FWHM,xoff)
C     ***********************************************************************************
C     Function to compute Hamming Function from parameters defined in Wikipedia
C
C     Input parameters
C 	FWHM	FWHM of Hamming function
C	xoff	offset from centre of Hamming Function in units of FWHM
C    
C     Output parameter
C	hamming	Hamming Function weight
C
C     Pat Irwin		5/11/13
C
C     ***********************************************************************************
      implicit none
      real pi,FWHM,alpha,beta,wfwhm,xx,xoff
      parameter (pi=3.1415927, alpha=0.53836, beta=0.46164)
      parameter (wfwhm = 0.5264)
      

C     Compute position of xoff in scaled window covering range 0 to 1.
      xx = 0.5 + 0.5*wfwhm*abs(xoff)/(0.5*FWHM)

      if(xx.le.1.0) then
       hamming = alpha-beta*cos(2.*pi*xx)
      else
       hamming=0.
      endif

      return

      end

