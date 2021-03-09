      subroutine writerawx(lraw,ispec,xlat,xlon,npro,nvar,varident,
     1  varparam,nx,xn,st,nvarx,varidentx,varparamx,nxx,xnx,stx)
C     $Id:
C     ***********************************************************************
C     Output the results of retrieval code (as writeraw.f), but also writes the
C     results of previous retrievals in case LIN=4
C
C     Input variables
C	lraw		integer		Output unit number
C	ispec		integer		ID number of spectrum
C	xlat		real		Latitude
C	xlon		real		Longitude 
C	npro		integer		Number of levels in profiles
C	nvar		integer		Number of variable profiles
C	varident(mvar,3) integer	Identity of profiles and 
C						parameterisation
C	varparam(mvar,mparam) real 	Extra parameters as required
C	nx		integer		Number of elements in state vector
C	xn(mx)		real		Retrieved measurement vector
C	st(mx,mx)	real		A priori covariance matrix
C
C     Pat Irwin		21/2/05
C
C     ***********************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'

      integer nx,i,j,ispec,lraw,npro,nvartot,nxtot
      real xn(mx),st(mx,mx),xlat,xlon,varparam(mvar,mparam)
      integer varident(mvar,3),nvar,ivar
      real xnx(mx),stx(mx,mx),varparamx(mvar,mparam)
      integer varidentx(mvar,3),nvarx,ivarx,nxx
      real xntot(mx),sttot(mx,mx)
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      write(lraw,*) ispec,'   ! ispec'
      write(lraw,*)xlat,xlon,' ! Latitude, Longitude'

      nvartot = nvar + nvarx      
      write(lraw,*)npro,nvartot,'  ! npro,nvar'

      do ivar=1,nvar
        write(lraw,*)ivar,'  ! ivar'
        write(lraw,*)(varident(ivar,j),j=1,3)
        write(lraw,*)(varparam(ivar,j),j=1,mparam)
      enddo
      
      do ivarx=1,nvarx
        write(lraw,*)ivarx+nvar,'  ! ivar'
        write(lraw,*)(varidentx(ivarx,j),j=1,3)
        write(lraw,*)(varparamx(ivarx,j),j=1,mparam)
      enddo

      nxtot = nx + nxx

   
      do i=1,mx
       xntot(i)=0.0
       do j=1,mx
        sttot(i,j)=0.0
       enddo
      enddo


      do i=1,nx
       xntot(i)=xn(i)
       do j=1,nx
         sttot(i,j) = st(i,j)
       enddo
      enddo

      do i=nx+1,nxtot
       xntot(i)=xnx(i-nx)
       do j=nx+1,nxtot
         sttot(i,j) = stx(i-nx,j-nx)
       enddo
      enddo


      write(lraw,*)nxtot,'  ! nx'
      write(lraw,*)(xntot(i),i=1,nxtot)
      do i=1,nxtot
       write(lraw,*)(sttot(i,j),j=i,nxtot)
      enddo

      return

      end


