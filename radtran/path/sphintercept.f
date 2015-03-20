      subroutine sphintercept(r0,dir,R,r1,thick,ierrout)
C     ************************************************************
C     Procedure to find point of intersection between a line and a sphere.
C
C     Input parameters
C	r0(3)	real	Starting position vector
C	dir(3)	real	Direction vector
C	R	real	Radius of sphere to intercept
C	thick	real	Minimum distance of intercept from start for
C			 solution to be considered significant
C     Output variables
C	r1(3)	real	Position vector of interception point
C	ierrout integer	Equals 1 if no solution, set to 0 otherwise
C
C     Pat Irwin	Original	20/3/15
C
C     ************************************************************
      implicit none
      real r0(3),dir(3),R,r1(3),thick
      real aa,bb,cc,soln(2),lambda,min,max,a
      integer ierrout,ierr,i
      aa=0.
      bb=0.
      cc=-R**2
      do i=1,3
       aa=aa+dir(i)**2
       bb=bb+2.0*r0(i)*dir(i)
       cc=cc+r0(i)**2
      enddo

      ierrout=0
      call quadsolve(aa,bb,cc,soln,ierr)

C      print*,soln(1),soln(2)
      if(ierr.eq.0) then
       if(soln(1).lt.-thick.and.soln(2).lt.-thick) then
C        print*,'Solution is in wrong direction'
        ierrout=1
        return
       else
        if(soln(1).lt.-thick)then
         lambda=soln(2)
        else
         if(soln(2).lt.-thick)then
          lambda=soln(1)
         else
          a=min(abs(soln(1)),abs(soln(2)))
C          print*,'a,thick',a,thick
          if(a.gt.thick)then
           lambda=a
          else
           lambda=max(soln(1),soln(2))
          endif
C          print*,'lambda',lambda
         endif
        endif
       endif
      else
       print*,'sphintercept: No real solution'
       ierrout=1
       return
      endif

      do i=1,3
       r1(i)=r0(i)+lambda*dir(i)
      enddo

      return

      end
