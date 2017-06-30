      subroutine readraw(lraw,xlat,xlon,npro,nvar,varident,
     1  varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,nx,xn,st)
C     $Id:
C     ***********************************************************************
C     Output the results of retrieval code
C
C     Output variables
C	lraw		integer		file unit number
C	xlat		real		Latitude
C	xlon		real		Longitude 
C	npro		integer		Number of levels in profiles
C	nvar		integer		Number of variable profiles
C	varident(mvar,3) integer	Identity of profiles and 
C						parameterisation
C	varparam(mvar,mparam) real 	Extra parameters as required
C	jsurf		integer		Position of surface temperature
C					element (if retrieved)
C	jalb		integer		Position of surface albedo spec
C	jxsc		integer		Position of x-section spec
C	jtan		integer		Position of tangent ht correction
C	jpre		integer		Position of tangent pressure
C	jrad		integer		Position of surface radius
C	jlogg		integer		Position of surface log_10(g)
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

      integer nx,i,j,ispec,lraw,npro,ivar1
      real xn(mx),st(mx,mx),xlat,xlon,varparam(mvar,mparam)
      integer varident(mvar,3),nvar,ivar,jsurf,ioff,np,jalb,jtan
      integer jpre, jrad, npvar, jlogg, jxsc

      read(lraw,*)ispec
      read(lraw,*)xlat,xlon
      
      read(lraw,*)npro,nvar

      jsurf = -1
      jalb = -1
      jxsc = -1
      jtan = -1
      jpre = -1
      jrad = -1
      jlogg = -1
      ioff=0
      do ivar=1,nvar
        read(lraw,*)ivar1
        read(lraw,*)(varident(ivar,j),j=1,3)
        read(lraw,*)(varparam(ivar,j),j=1,mparam)
        np = 1
        if(varident(ivar,1).le.100)then
          np = npvar(varident(ivar,3),npro)
        endif
        if(varident(ivar,1).eq.888)np = int(varparam(ivar,1)) 
        if(varident(ivar,1).eq.887)np = int(varparam(ivar,1)) 
        if(varident(ivar,1).eq.444)np = 2+int(varparam(ivar,1)) 
        if(varident(ivar,1).eq.445)np = 3+int(varparam(ivar,1)) 
        if(varident(ivar,1).eq.222)np = 8
        if(varident(ivar,1).eq.223)np = 9
        if(varident(ivar,1).eq.224)np = 9
        if(varident(ivar,1).eq.225)np = 11
        if(varident(ivar,1).eq.226)np = 8
        if(varident(ivar,1).eq.227)np = 7

        ioff=ioff+np
        if(varident(ivar,1).eq.999)jsurf=ioff
        if(varident(ivar,1).eq.887)jxsc=ioff
        if(varident(ivar,1).eq.888)jalb=ioff
        if(varident(ivar,1).eq.889)jalb=ioff
        if(varident(ivar,1).eq.777)jtan=ioff
        if(varident(ivar,1).eq.666)jpre=ioff
        if(varident(ivar,1).eq.555)jrad=ioff
        if(varident(ivar,1).eq.333)jlogg=ioff

      enddo
      read(lraw,*)nx
      read(lraw,*)(xn(i),i=1,nx)
      do i=1,nx
       read(lraw,*)(st(i,j),j=i,nx)
       do j=i,nx
        st(j,i)=st(i,j)
       enddo
      enddo

      

      return

      end


