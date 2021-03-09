      subroutine assess(nx,ny,kk_in,sx_in,se_in)
C     $Id:
C     ********************************************************************
C     This subroutine assesses the retrieval matrices to see
C     whether an exact retrieval may be expected.
C
C     One formulation of the gain matrix is
C        dd = sx*kk_T*(kk*sx*kk_T + se)^-1
C
C     If the retrieval is exact, the se will be very small. Since se is
C     diagonal all we need do is compare to  the diagonal elements of
C     
C     The matrices are defined following the convention
C				 A(a1,a2) where a1=number or rows
C				  	        a2=number of columns
C
C     The input variables are:
C	nx		integer	   number of elements of x-matrix
C	ny		integer	   number of elements of y-matrix
C       kk_in(my,mx)	real	   Rate of change of vector y with vector x
C				    calculated for xn
C       sx_in(mx,mx)	real	   Covariance matrix of a priori
C 	se_in(my,my)	real	   Diagonal elements of measurement vector
C                                  covariance matrix
C
C     All input variables are transferred to common dimension arrays of
C     dimension (my,my) for ease of internal manipulation. 
C     Internal arrays are double precision. 
C
C     The following subroutines are called:
C
C	dmult_mat	multiplies two matrices
C
C     Pat Irwin	9/9/04		Original (adapted from dretrieve)
C
C     ********************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
     
      integer nx,ny,a1,a2,m1,m2
      integer kk1,kk2,kt1,kt2,i,j
      integer icheck
      real kk_in(my,mx)
      real sx_in(mx,mx),se_in(my,my)

      double precision a(my,my)
      double precision m(my,my)
      double precision kk(my,my),kt(my,my)
      double precision sum1,sum2,sum3

      integer idiag,iquiet
      common/diagnostic/idiag,iquiet
     
C     load kk with kk_in
      do i=1,ny
       do j=1,nx
        kk(i,j) = dble(kk_in(i,j))
        kt(j,i) = dble(kk_in(i,j))
       end do
      end do
      kk1=ny
      kk2=nx
      kt1=nx
      kt2=ny

C     ******************************************************************
C     Calculate core of the Contribution Function
C     ******************************************************************

C     Load a with sx_in

      do i=1,nx
         do j=1,nx
          a(i,j)=dble(sx_in(i,j))
         end do
      end do

C     Multiply sx*kt, put answer in m for later use
      call dmult_mat(my,a,nx,nx,kt,kt1,kt2,m,m1,m2)

C     Multiply kk*m, put answer in a
      call dmult_mat(my,kk,kk1,kk2,m,m1,m2,a,a1,a2)

C     Add se to a
      sum1 = 0.0
      sum2 = 0.0
      sum3 = 0.0
      do i=1,ny
       do j=1,ny
         a(i,j) = a(i,j)+dble(se_in(i,j))
         if(i.eq.j)then
          sum1 = sum1+a(i,j)
          sum2 = sum2+dble(se_in(i,j))
          sum3 = sum3 + a(i,j)/dble(se_in(i,j))
         endif
       enddo
      end do
      sum1 = sum1/dble(ny)
      sum2 = sum2/dble(ny)
      sum3 = sum3/dble(ny) 
      if(idiag.gt.0)then
       print*,'Assess:'
       print*,'Average of diagonal elements of Kk*Sx*Kt : ',sum1 
       print*,'Average of diagonal elements of Se : ',sum2 
       print*,'Ratio = ',sum1/sum2
       print*,'Average of Kk*Sx*Kt/Se element ratio : ',sum3 
       if(sum3.gt.10.0)then
        print*,'******************* ASSESS WARNING *****************'
        print*,'Insufficient constraint. Solution likely to be exact'
        print*,'****************************************************'
       endif
      endif

      return

      end
