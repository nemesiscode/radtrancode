      subroutine calcnextxn(nx,ny,x0_in,xn_in,y_in,yn_in,
     1 dd_in,aa_in,x_out)
C     $Id:
C     ********************************************************************
C     This subroutine performs the optimal estimation retrieval of the
C     vector x from a set of measurements y and forward derivative matrix
C     kk. The equation solved is (re: p147 of Houghton, Taylor and Rodgers):

C   
C     xn+1 = x0 + dd*(y-yn) - aa*(x0 - xn)
C
C     The matrices are defined following the convention
C				 A(a1,a2) where a1=number or rows
C				  	        a2=number of columns
C
C     The input variables are:
C	nx		integer	   number of elements of x-matrix
C	ny		integer	   number of elements of y-matrix
C	x0_in(mx)	real	   a priori guess of x
C	xn_in(mx)	real       Retrieved vector x at last iteration
C	y_in(my)	real	   The measurement vector y that we are
C				    trying to fit
C	yn_in(my)	real	   Fitted value of y at last iteration
C	dd_in(mx,my)    double     Gain matrix
C	aa_in(mx,mx)	double     Averaging Kernels
C
C     All input variables are transferred to common dimension arrays of
C     dimension (my,my) for ease of internal manipulation. 
C     Internal arrays are double precision. 
C
C     The following subroutines are called:
C
C	dmult_mat	multiplies two matrices
C
C     Parameters output are:
C
C	x_out(nx)	real	   Newly retrieved state vector
C
C     Program developed from PMIRR water vapour retrieval study
C
C     Pat Irwin	30/1/94		Original
C     Pat Irwin	17/10/03	Tidied for Nemesis
C     Pat Irwin 9/9/04          Modified to use pre-calculated gain matrix
C
C     ********************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
     
      integer nx,ny,a1,a2,b1,b2,c1,c2,d1,d2
      integer kk1,kk2,n1,n2,i,j,x11,x12,x21,x22
      integer icheck
      real y_in(my),yn_in(my),x0_in(mx),x_out(mx)
      real xn_in(mx)

      double precision a(my,my),b(my,my),dd_in(mx,my),c(my,my)
      double precision x1(my,my),aa_in(mx,mx),d(my,my),x2(my,my)

C     load a with aa_in
      do i=1,nx
       do j=1,nx
        a(i,j) = aa_in(i,j)
       end do
      end do
      a1=nx
      a2=nx

C     Subtract x0 - xn, put answer in b
      do 20 i=1,nx
       b(i,1)=dble(x0_in(i) - xn_in(i))
20    continue
      b1=nx
      b2=1


C     Multiply a*b, put answer in x1
      call dmult_mat(my,a,a1,a2,b,b1,b2,x1,x11,x12) 


C     load d with dd_in
      do i=1,nx
       do j=1,ny
        d(i,j) = dd_in(i,j)
       end do
      end do
      d1=nx
      d2=ny

C     Put c = y-yn
      do 25 i=1,ny
       c(i,1)=dble(y_in(i) - yn_in(i))
25    continue
      c1=ny
      c2=1

C     Multiply d*c, put answer in x2
      call dmult_mat(my,d,d1,d2,c,c1,c2,x2,x21,x22) 

C     Subtract x1 from x2 and add to x0
      do 40 i=1,nx
       x_out(i)=x0_in(i)+sngl(x2(i,1)-x1(i,1))
40    continue

      return

      end
