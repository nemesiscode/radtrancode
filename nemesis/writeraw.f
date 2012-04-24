      subroutine writeraw(lraw,ispec,xlat,xlon,npro,nvar,varident,
     1  varparam,nx,xn,st)
C     $Id:
C     ***********************************************************************
C     Output the results of retrieval code
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

      integer nx,i,j,ispec,lraw,npro
      real xn(mx),st(mx,mx),xlat,xlon,varparam(mvar,mparam)
      integer varident(mvar,3),nvar,ivar

      write(lraw,*) ispec,'   ! ispec'
      write(lraw,*)xlat,xlon,' ! Latitude, Longitude'
      
      write(lraw,*)npro,nvar,'  ! npro,nvar'

      do ivar=1,nvar
        write(lraw,*)ivar,'  ! ivar'
        write(lraw,*)(varident(ivar,j),j=1,3)
        write(lraw,*)(varparam(ivar,j),j=1,mparam)
      enddo
      write(lraw,*)nx,'  ! nx'
      write(lraw,*)(xn(i),i=1,nx)
      do i=1,nx
       write(lraw,*)(st(i,j),j=i,nx)
      enddo

      return

      end


