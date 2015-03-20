      subroutine findnear(jlay,nlay,rad,r0,vec,ethick,ilay,
     1	r1,dist,theta)
C     *************************************************************
C     Procedure to find nearest intersection of ray with a spherical shell
C     above below or at the same level as the starting point. Intesection 
C     must be further away from start than ethick.
C
C     Input variables:
C 	jlay		integer	current radius ordinate
C 	nlay		integer	number of radius ordinates
C	rad(nlay+1)	real	radii
C 	r0(3)		real	Starting position vector
C 	vec(3) 		real	direction vector of ray
C	ethick 		real	Minimum distance away from valid solution
C
C     Output variables:
C 	ilay		integer	radius ordinate of intersection
C 	r1(3)		real	New position vector
C	dist		real	Distance between r0 and r1
C 	theta		real	Angle between ray and local normal at r1
C		 		> 90 and ray going down
C		 		< 90 and ray going up
C
C     Pat Irwin	Original	20/3/15
C
C     *************************************************************
      implicit none
      include '../includes/arrdef.f'
      integer jlay,nlay,ilay,i,j1,j2,j,ierrout,k
      real rad(maxlay+1),r0(3),vec(3),ethick,r1(3),dist,theta
      real dtr,sproduct,sum,rtmp,rsave(3)
      parameter (dtr=3.1415927/180.)


      j1=jlay-1
      if(j1.lt.1)j1=1
      j2=jlay+1
      if(j2.gt.nlay+1)j2=nlay+1

      sum=1e36
      dist=-1.
      ilay=-1
      do 10 j=j1,j2
       rtmp = rad(j)
       call sphintercept(r0,vec,rtmp,r1,ethick,ierrout)
       if(ierrout.eq.0) then
         dist=0.
         do i=1,3
          dist=dist+(r0(i)-r1(i))**2
         enddo
         dist=sqrt(dist)
         if(dist.lt.sum.and.dist.gt.ethick)then
          sum=dist
          ilay=j
          do k=1,3
           rsave(k)=r1(k)
          enddo
         endif
       endif 
10    continue

      if(ilay.lt.0) then 
       print*,'Error in findnear - solution not found'
      endif

      do k=1,3
           r1(k)=rsave(k)
      enddo
      dist=sum
      theta=acos(sproduct(r1,vec)/rad(ilay))/dtr

      return

      end

