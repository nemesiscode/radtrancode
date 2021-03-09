      subroutine readnextiter(nx,ny,xn1,xa,y,se1,yn1,yn,
     & chisq,phi)
C     ***************************************************************
C     Subroutine to read in next iterated solution from the iter file
C
C     Input variables
C     	None   	It is assumed that .itr file is already open with file 
C		number = 37
C
C     Output variables
C	nx	integer	length of state vector
C	ny	integer	length of measurement vector
C	xn1(mx)	real	Current state vector
C	xa(mx)	real	A priori state vector
C	y(my)	real	Measurement vector
C	se1(my)	real	diagonal elements of measurement covariance matrix
C	yn1(my)	real	Current modelled measurement vector
C	yn(my)	real	Previous modelled best-fit measurement vector
C 	chisq	real	chisq parameter
C	phi	real	phi parameter
C
C     Pat Irwin		25/5/16
C
C     ***************************************************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      integer nx,ny,i,j
      real xn1(mx),xa(mx),y(my),se1(my),yn1(my),yn(my)
      real chisq,phi,kk(my,mx)
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      read(37,*)chisq,phi
      read(37,*)(xn1(i),i=1,nx)
      read(37,*)(xa(i),i=1,nx)
      read(37,*)(y(i),i=1,ny)
      read(37,*)(se1(i),i=1,ny)
      read(37,*)(yn1(i),i=1,ny)
      read(37,*)(yn(i),i=1,ny)
      do i=1,nx
          read(37,*)(kk(j,i),j=1,ny)
      enddo
      if(idiag.gt.0)print*,'readnextiter OK'
      return

      end
      
