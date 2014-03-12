      real function calc_phiprior(nx,xn,xa,sai)
C     $Id:
C     *******************************************************************
C     Function to calculate the retrieval cost function in terms of 
C     departure from a priori.
C
C     Input variables
C	nx	integer	Number if elements of state vector
C	xn(mx) 	real	State vector
C	sai(mx,mx) double precision Inverse of a priori covariance 
C			            matrix
C	
C     Output variables
C	calc_phiprior real Departure from a priori
C
C     Pat Irwin		Original	18/7/00
C     Pat Irwin		Tidied		17/10/03
C     Pat Irwin		Computes departure of a priori only  11/6/11	
C
C     *******************************************************************
      implicit none
      integer nx,ny
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      real xn(mx),xa(mx)
      double precision sai(mx,mx),phi1,phi2,a(my,my),b(my,my)
      double precision bt(my,my),c(my,my)
      integer i,j,a1,a2,b1,b2,bt1,bt2,c1,c2
     

      do i=1,nx
       do j=1,nx
        b(j,i)=0.0
        bt(j,i)=0.0
       enddo
      enddo

C     Load xn-xa into b and bt
      do i=1,nx
       b(i,1)=dble(xn(i)-xa(i))
       bt(1,i)=b(i,1)
      enddo
      b1 = nx
      b2 = 1
      bt1 = 1
      bt2 = nx
C     Load sai into c
      do i=1,nx
       do j=1,nx
        c(i,j)=sai(i,j)
       enddo
      enddo

C     Multiply sai*b, put answer in a
      call dmult_mat(my,c,nx,nx,b,b1,b2,a,a1,a2)

C     Multiply bt*a, put answer in c
      call dmult_mat(my,bt,bt1,bt2,a,a1,a2,c,c1,c2)

      phi2 = c(1,1)

      print*,'calc_phiprior: phi2 = ',sngl(phi2)

      calc_phiprior = sngl(phi2)

      return

      end
