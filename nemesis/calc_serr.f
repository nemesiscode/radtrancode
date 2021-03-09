      subroutine calc_serr(nx,ny,sx_in,se_in,aa_in,dd_in,
     &  st_out,sn_out,sm_out)
C     $Id:
C     ********************************************************************
C     This subroutine calculates the error covariance matrices after the
C     final iteration has been completed.
C
C     The subroutine calculates the MEASUREMENT covariance matrix according
C     to the equation (re: p130 of Houghton, Taylor and Rodgers) :
C
C     			sm = dd*se*dd_T
C
C     The subroutine calculates the SMOOTHING error covariance matrix
C     according to the equation:
C
C			sn = (aa-I)*sx*(aa-I)_T
C
C     The subroutine also calculates the TOTAL error matrix:
C
C			st=sn+sm
C
C     Input variables are:			
C
C	nx		integer    number of elements of vector x
C	ny		integer	   number of elements of vector y
C       sx_in(mx,mx)	real	   Covariance matrix of a priori x
C	se_in(my,my)	real 	   Covariance matrix of measurement vector y
C  	dd_in(mx,my)	double	   Gain matrix
C	aa_in(mx,mx)	double     Averaging kernels
C
C     The following subroutines are called:
C	mult_mat	multiplies two matrices
C	add_mat		adds two matrices
C	sub_mat		subtracts two matrices
C
C     The array dimensions mx and my are included in arraylen.f
C
C     Output variables are:
C
C	sm_out(mx,mx)	real	   Final measurement covariance matrix
C	sn_out(mx,mx)	real	   Final smoothing error covariance matrix
C	st_out(mx,mx)	real	   Final full covariance matrix
C
C     Pat Irwin	30/1/94		Original
C     Pat Irwin	17/10/03	Removed idim and replaced with my
C
C     ********************************************************************

      implicit none


      integer aa1,aa2,nx,ny,sn1,sn2,sm1,sm2,st1,st2,c1,c2,a1,a2
      integer dt1,dt2,dd1,dd2,se1,se2,i,j
      integer icheck
      integer m1,m2,b1,b2

C     ********************************************************************
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
C     ********************************************************************

      real sx_in(mx,mx),se_in(my,my)
      real st_out(mx,mx),sn_out(mx,mx),sm_out(mx,mx)
      double precision dd_in(mx,my),aa_in(mx,mx)
      
      double precision dd(my,my),dt(my,my)
      double precision a(my,my),b(my,my),c(my,my),m(my,my)
      double precision se(my,my),sx(my,my)
      double precision sm(my,my),sn(my,my),st(my,my)

      do i=1,nx
        do j=1,ny
         dd(i,j)=dd_in(i,j)
         dt(j,i)=dd_in(i,j)
        enddo
      enddo
      dd1 = nx
      dd2 = ny
      dt1 = ny
      dt2 = nx

      do i=1,ny
       do j=1,ny
        se(i,j)= se_in(i,j)
       end do
      end do
      se1=ny
      se2=ny

C     multiply dd*se, put answer in a
      call dmult_mat(my,dd,dd1,dd2,se,se1,se2,a,a1,a2)


C     multiply a*dt, put answer in sm
      call dmult_mat(my,a,a1,a2,dt,dt1,dt2,sm,sm1,sm2)

C     load b with sx_in
      do i=1,nx
       do j=1,nx
        b(i,j)=dble(sx_in(i,j))
       end do
      end do

C     load a with aa-ii

      do i=1,nx
       do j=1,nx
        a(i,j)=aa_in(i,j)
       end do
       a(i,i)=a(i,i)-1.0
      end do
      aa1 = nx
      aa2 = nx

C     multiply a*b, put answer in c
      call dmult_mat(my,a,aa1,aa2,b,nx,nx,c,c1,c2)


C     load a with (aa-ii)T

      do i=1,nx
       do j=1,nx
        a(i,j)=aa_in(j,i)
       end do
       a(i,i)=a(i,i)-1.0
      end do

C     multiply c*a, put answer in sn
      call dmult_mat(my,c,c1,c2,a,aa1,aa2,sn,sn1,sn2)


C     Add sn and sm together to give st

      call dadd_mat(my,sm,nx,nx,sn,nx,nx,st,st1,st2)

C     Transfer matrices to output
      do i=1,nx
       do j=1,nx
        sm_out(i,j)=sngl(sm(i,j))
        sn_out(i,j)=sngl(sn(i,j))
        st_out(i,j)=sngl(st(i,j))
       end do
      end do
   
      return

      end



