      real function calc_chi(ny,y,yn,se1,chisq)
C     $Id:
C     *******************************************************************
C     Function to calculate the retrieval cost function which combines
C     departure from a priori and closeness to spectrum.
C
C     Input variables
C	ny	integer	Number of points in measurement vector
C	y(my)	real	Measurement vector
C	yn(my)	real	Calculated measurement vector
C       se1(my) real    Measured radiance variances
C	
C     Output variables
C      calc_chi	real	Closeness of fit to measurement vector
C
C     Pat Irwin		Original	18/7/00
C     Pat Irwin		Tidied		17/10/03	
C     Pat Irwin		Modified from calc_phiret  11/6/13
C
C     *******************************************************************
      implicit none
      integer ny
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      real y(my),yn(my)
      real se1(my)
      real chisq
      integer i,j
     
      chisq = 0.0


      do i=1,ny
       chisq=chisq+(yn(i)-y(i)/se1(i))**2
      enddo

      calc_chi = chisq

      return

      end
