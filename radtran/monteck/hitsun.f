      subroutine hitsun(solvec,dvec,scos,ihit)
C     ********************************************************
C     Subroutine to see if exiting photon comes within a certain
C     angular distance from the Sun.
C
C     Input variables
C   	solvec(3)	real	Direction vector of Sun
C	dvec(3)		real	Direction vector of exiting photon
C	scos		real	Cos(theta_min) - cosine of minimum
C				acceptable angle between Sun and photon  
C
C     Output variable
C	ihit		integer	0=if miss, 1 if hit
C
C     Pat Irwin		1/8/05
C
C     ********************************************************
      real solvec(3),dvec(3),scos,xx,yy,xy,cthet
      integer ihit,i

      xx=0.0
      yy=0.0
      xy=0.0
      do i=1,3
       xx=xx+solvec(i)**2
       yy=yy+dvec(i)**2
       xy=xy+solvec(i)*dvec(i)
      enddo

      cthet = xy/sqrt(xx*yy)

      if(cthet.ge.scos)then
       ihit=1
      else
       ihit=0
      endif

      return

      end



