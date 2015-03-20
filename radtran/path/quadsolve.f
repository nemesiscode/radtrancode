      subroutine quadsolve(a,b,c,xx,ierr)
C     *********************************************************
C     Procedure to find real solutions of x satisfying the equation:
C	ax^2+bx+c=0.
C
C     Input variables
C	a,b,c	real	coefficients of quadratice equation
C
C     Output variables
C	x(2)	real	Solutions to x
C	ierr	integer	Set to 1 if no solution exists, Set to 0 otherwise.
C
C     Pat Irwin	20/3/15
C
C     *********************************************************
      implicit none
      integer ierr,i
      real a,b,c,xx(2)
      logical isnan,ntest1,ntest2
C      print*,'abc',a,b,c
      ierr=0
      do i=1,2
       xx(i)=0.
      enddo
C      test=b**2-4*a*c
C      print*,'test = ',test
      if(b**2.ge.4*a*c) then
        xx(1)=(-b+sqrt(b**2-4.0*a*c))/(2.*a)
        xx(2)=(-b-sqrt(b**2-4.0*a*c))/(2.*a)
      else 
        ierr=1
      endif

C      print*,'xx',xx(1),xx(2)
      ntest1=isnan(xx(1))
      ntest2=isnan(xx(2))
      if(ntest1.or.ntest2)then
       if(b.eq.0.0.and.c.eq.0.0)then
        xx(1)=-b
        xx(2)=-b
       else
        ierr=1
       endif
      endif
      return

      end
