      subroutine check_kk(nx,ny,xn,xa,sa,kk,nxf,xnf,xaf,saf,kkf,nset)
C     $Id:
C     *****************************************************************
C     Subroutine to check for any zero rows of kk and remove. This is done
C     to prevent the retrieval algorithm becoming unstable. The
C     routine also strips out the current and a priori state vectors 
C     accordingly.
C
C     Input variables
C	nx	integer	Number of elements in state vector
C	ny	integer	Number if elements in measurement vector
C	xn(mx)	real	Current state vector
C	xa(mx)	real	A Priori state vector
C	sa(mx,mx) real	A priori covariance matrix
C	kk(my,mx) real	Functional derivatives matrix
C	
C     Output variables
C 	nxf	integer	Number if state vector elements for which kk has
C				non-zero rows
C	xnf(mx)	real	Stripped current state vector
C	xaf(mx)	real	Stripped a priori state vector
C	saf(mx,mx) real	Stripped a priori covariance matrix
C	kkf(my,mx) real	Stripped functional derivatives matrix
C	nset(mx) integer Element numbers of xn which have non-zero kk rows
C
C     Pat Irwin	17/10/03 	Original
C	
C     *****************************************************************
      implicit none
C     Set measurement vector and source vector lengths here.
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

      real kk(my,mx),kkf(my,mx),xa(mx),xaf(mx),xn(mx),xnf(mx)
      real sa(mx,mx),saf(mx,mx),xk
      integer i,j,nx,ny,nxf,ix,nset(mx)
      logical ltest
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      if(idiag.gt.0)print*,'Check_kk, requested nx = ',nx
      ix=0
      DO j=1,nx
       ltest=.true.
       xk = kk(1,j)
       do i=2,ny
        if(kk(i,j).ne.xk)ltest=.false.
       enddo
       if(ltest)then
        if(idiag.gt.0)then
         print*,'Variable : ',j,' has flat weighting function. Discard.'
        endif
       else
        ix=ix+1
        nset(ix)=j
       endif
      ENDDO

      nxf = ix
      if(idiag.gt.0)print*,'Check_kk, allowed nxf = ',nxf

      if(nxf.gt.0)then
       do i=1,nxf
        xnf(i)=xn(nset(i))
        xaf(i)=xa(nset(i))
        do j=1,nxf
         saf(j,i)=sa(nset(j),nset(i))
        enddo
        do j=1,ny
         kkf(j,i)=kk(j,nset(i))
        enddo
       enddo
      endif

      return

      end
