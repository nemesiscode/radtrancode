      subroutine deflectray(vec,r,theta,vec1)
C     **************************************************************
C     Procedure to deflect a ray in the plane defined by the initial
C     vector and the normal to the surface.
C
C     Input variables
C       vec(3)	real  original vector
C	r(3)    real  position vector = normal
C	theta   real  Angle to deflect ray by. 
C		   If theta > 0 then angle of ray w.r.t. r increased 
C
C     Output variables
C	vec1(3)	real  refracted ray
C
C     Pat Irwin   13/3/15
C    
C     **************************************************************
      implicit none
      real vec(3),r(3),theta,vec1(3),dtr,sproduct,n(3)
      real rr,gamma,cos2thet,aa,bb,cc,alpha(2)
      real thet1,thet2(2)
      integer i,ierr,beta(2),j
      parameter(dtr=3.1415927/180.)

      rr=sqrt(sproduct(r,r))
      do i=1,3
       n(i)=r(i)/rr
      enddo

      gamma=sproduct(vec,n)

      cos2thet=cos(theta*dtr)**2
      aa = gamma**2-cos2thet 
      bb = 2.0*gamma*(1.0-cos2thet)
      cc = 1.0-cos2thet

      call quadsolve(aa,bb,cc,alpha,ierr)
      if(ierr.eq.1) then
        print*,'Error in deflectray - no solution'
        print*,vec
        print*,r
        print*,theta
        stop
      endif

      do i=1,2
       beta(i) = cos(theta*dtr)/(1.0+gamma*alpha(i))
      enddo
      thet1 = acos(sproduct(vec,n))/dtr

      do i=1,2
       do j=1,3
        vec1(j)=beta(i)*(vec(j)+alpha(i)*n(j))
       enddo
       thet2(i)=acos(sproduct(vec1,n))/dtr
      enddo

      if(theta.gt.0.0) then
       if(thet2(1).gt.thet2(2)) then 
         i=1 
       else
         i=2
       endif
      else
       if(thet2(1).lt.thet2(2)) then 
        i=1 
       else 
        i=2
       endif 
      endif

      do j=1,3
       vec1(j)=beta(i)*(vec(j)+alpha(i)*n(j))
      enddo

      return

      end
