      subroutine calc_gain_matrix(nx,ny,kk_in,sx_in,sxi,se_in,sei,
     1  dd_out,aa_out)
C     $Id:
C     ********************************************************************
C     This subroutine calculates the optimal estimation gain matrix
C     (sometimes known as contribution functions)
C
C     The gain matrix is either calculated as:   
C        dd = sx*kk_T*(kk*sx*kk_T + se)^-1
C
C     if nx >= ny or as
C        dd = ((sx^-1 + kk_T*se^-1*kk)^-1)*kk_T*se^-1
C     otherwise.
C
C     It is assumed that se^-1 and sx^-1 have already been
C     calculated.
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
C	sxi(mx,mx)	double     Inverse of sx
C 	se_in(my,my)	real	   Measurement vectorcovariance matrix
C 	sei(my,my)	double	   Inverse of se
C
C     All input variables are transferred to common dimension arrays of
C     dimension (my,my) for ease of internal manipulation. 
C     Internal arrays are double precision. 
C
C     The following subroutines are called:
C
C	dmult_mat	multiplies two matrices
C	dinvertm	inverts a matrix
C
C     Parameters output are:
C
C	dd_out(mx,my)	double	   Gain matrix
C	aa_out(mx,mx)	double	   Averaging kernels
C
C     Pat Irwin	9/9/04		Original (adapted from dretrieve)
C
C     ********************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
     
      integer nx,ny,a1,a2,b1,b2,m1,m2,dd1,dd2
      integer kk1,kk2,kt1,kt2,c1,c2,i,j
      integer icheck
      real kk_in(my,mx)
      real sx_in(mx,mx),se_in(my,my)

      double precision a(my,my),b(my,my),sxi(mx,mx)
      double precision m(my,my),c(my,my),dd(my,my)
      double precision kk(my,my),kt(my,my)
      double precision sum1,sum2,sum3,sei(my,my)
      double precision dd_out(mx,my),aa_out(mx,mx)


     
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
C     Calculate the Contribution Functions and Averaging Kernels
C     ******************************************************************

      if(nx.ge.ny)then
C       Load a with sx_in

        do i=1,nx
         do j=1,nx
          a(i,j)=dble(sx_in(i,j))
         end do
        end do

C       Multiply sx*kt, put answer in m for later use
        call dmult_mat(my,a,nx,nx,kt,kt1,kt2,m,m1,m2)

C       Multiply kk*m, put answer in a
        call dmult_mat(my,kk,kk1,kk2,m,m1,m2,a,a1,a2)

C       Add se to a
        do i=1,ny
         do j=1,ny
          a(i,j) = a(i,j)+dble(se_in(i,j))
         enddo
        end do

C       Invert a, put answer in b
C       jmod = 0 : LU decomposition
C       jmod = 1 : Gaussian elimination
C       jmod = 2 : Cholesky decomposition

c        jmod=2 
        print*,'Calc_gain_matrix: Inverting...'
        call dinvertm(jmod,icheck,a,ny,my,b)
        b1=ny
        b2=ny

C       Multiply m(= sx*kt from above)*b, put answer in dd
C       This the gain matrix 
        call dmult_mat(my,m,m1,m2,b,b1,b2,dd,dd1,dd2)

      else

       
C       Calculate kt*sei. Put answer in m and keep for later
        call dmult_mat(my,kt,kt1,kt2,sei,ny,ny,m,m1,m2)

C       Calculate m*kk, put answer in c  
        call dmult_mat(my,m,m1,m2,kk,kk1,kk2,c,c1,c2)
   
C       Add sxi to c
        do i=1,nx
         do j=1,nx
          c(j,i)=c(j,i)+sxi(j,i)
         enddo
        enddo

C       Invert c, put answer in a
c        jmod = 2
        print*,'Calc_gain_matrix: Inverting...'
        call dinvertm(jmod,icheck,c,nx,my,a)
        a1=nx
        a2=nx

C       Multiply a*(kt*sei) [=a*m, m kept from earlier]
C       This is the gain matrix
        call dmult_mat(my,a,a1,a2,m,m1,m2,dd,dd1,dd2)
      endif

C     multiply dd*kk, put answer in a
C     These are the averaging kernels
      call dmult_mat(my,dd,dd1,dd2,kk,kk1,kk2,a,a1,a2)


C     Transfer matrices to real*4 output matrices
      do i=1,nx
        do j=1,ny
         dd_out(i,j)=dd(i,j)
        enddo
        do j=1,nx
         aa_out(i,j)=a(i,j)
        enddo
      enddo

      return

      end
