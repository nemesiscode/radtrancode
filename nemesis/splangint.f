      subroutine splangint(ny,ntheta,theta,yy,yy2,ang,result,grad)
C     *********************************************************
C     Subroutine to work out the cubic spline setup for a series which
C     repeats a periodic vector, y(ny), three times as an approximation of 
C     a true periodic spline.
C     
C     The  input vector, y(ny), which is assumed to equally sample the angle 
C     range 0 to 2pi, with y(ny+1)=y(1).
C    
C     Input variables
C	ny		integer	Number of ordinates
C       ntheta		integer	Number of spline points
C	theta(ntheta)	real	Spline angles (degrees)
C	yy(ntheta)	real	Spline points
C	yy2(ntheta)	real	Spline second derivatives
C	ang		real	Angle to interpolate to
C     Output variables
C	result		real	Intepolated value
C	grad(ny)	real	approximate ROC of result with ordinates
C				from linear interpolation
C
C     Pat Irwin  4/12/18
C
C     *********************************************************
      implicit none

      integer ny,i,j,mlon,mthet,ntheta
      real dtheta
      parameter(mlon=20,mthet=3*mlon+1)
      real grad(mlon),yy(mthet),yy2(mthet),theta(mthet)
      real result,ang,f

      call splint(theta,yy,yy2,ntheta,ang,result)

      dtheta = 360./float(ny)

      j=1+int((ang+360.)/dtheta)
      f = (ang-theta(j))/dtheta

      do i=1,ny
       grad(i)=0.
      enddo

      grad(j)=1.0-f
      if(j.lt.ny)then
       grad(j+1)=f
      else
       grad(1)=f
      endif

      return

      end
       


