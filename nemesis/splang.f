      subroutine splang(ny,y,ntheta,theta,yy,yy2)
C     *********************************************************
C     Subroutine to work out the cubic spline setup for a series which
C     repeats a periodic vector, y(ny), three times as an approximation of 
C     a true periodic spline.
C     
C     The  input vector, y(ny), which is assumed to equally sample the angle 
C     range 0 to 2pi, with y(ny+1)=y(1).
C    
C     Input variables
C	y(ny)	real	Values to fourier analyse
C	ny	integer	Number of ordinates
C
C     Output variables
C       ntheta		integer	Number of spline points
C	theta(ntheta)	real	Spline angles (degrees)
C	yy(ntheta)	real	Spline points
C	yy2(ntheta)	real	Spline second derivatives
C
C     Pat Irwin  4/12/18
C
C     *********************************************************
      implicit none

      integer ny,i,j,mlon,mthet,ntheta
      real dtheta
      parameter(mlon=20,mthet=3*mlon+1)
      real y(mlon),yy(mthet),yy2(mthet),theta(mthet)

      dtheta = 360./float(ny)

      ntheta = 3*ny+1
      j=1
      do i=1,ntheta
       theta(i)=-360+(i-1)*dtheta   
       yy(i)=y(j)
       j=j+1
       if(j.gt.ny)j=1
      enddo

      call spline (theta,yy,ntheta,1.e30,1.e30,yy2)

      return

      end
       


