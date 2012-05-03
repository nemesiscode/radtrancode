      real function calc_phiret(ny,y,yn,sei,nx,xn,xa,sai,chisq)
C     $Id:
C     *******************************************************************
C     Function to calculate the retrieval cost function which combines
C     departure from a priori and closeness to spectrum.
C
C     Input variables
C	ny	integer	Number of points in measurement vector
C	y(my)	real	Measurement vector
C	yn(my)	real	Calculated measurement vector
C	sei(my,my) double Inverse of measurement covariance matrix
C	nx	integer	Number if elements of state vector
C	xn(mx) 	real	State vector
C	sai(mx,mx) double precision Inverse of a priori covariance 
C			            matrix
C	
C     Output variables
C	chisq	real	Closeness of fit to measurement vector
C	calc_phiret real Total cost function
C
C     Pat Irwin		Original	18/7/00
C     Pat Irwin		Tidied		17/10/03	
C
C     *******************************************************************
      implicit none
      integer nx,ny
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      real y(my),yn(my),xn(mx),xa(mx)
      double precision sei(my,my)
      real chisq
      double precision sai(mx,mx),phi1,phi2,a(my,my),b(my,my)
      double precision bt(my,my),c(my,my)
      integer i,j,a1,a2,b1,b2,bt1,bt2,c1,c2
     
      phi1 = 0.0


      do i=1,ny
       do j=1,ny
        b(j,i)=0.0
        bt(j,i)=0.0
       enddo
      enddo

C     Load yn-y into b and bt
      do i=1,ny
       b(i,1)=dble(yn(i)-y(i))
       bt(1,i)=b(i,1)
      enddo
      b1 = ny
      b2 = 1
      bt1 = 1
      bt2 = ny

C     Multiply sei*b, put answer in a
      call dmult_mat(my,sei,ny,ny,b,b1,b2,a,a1,a2)

C     Multiply bt*a, put answer in c
      call dmult_mat(my,bt,bt1,bt2,a,a1,a2,c,c1,c2)

      phi1 = c(1,1)

      chisq = sngl(phi1)

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

      print*,'calc_phiret: phi1,phi2 = ',sngl(phi1),sngl(phi2)

      calc_phiret = sngl(phi1+phi2)

      return

      end
