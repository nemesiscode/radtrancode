      subroutine writeoutp(iform,runname,ispace,lout,ispec,xlat,xlon,
     1  npro,nvar,varident,varparam,nx,ny,y,yn,se,xa,sa,xn,err,ngeom,
     2  nconv,vconv,gasgiant,jpre,jrad,jlogg,jfrac,iscat,lin)
C     $Id:
C     ***********************************************************************
C     Output the results of retrieval code with the format required by the Matlab
C     visualization tools
C
C     Input variables
C	runname		character*100	run name
C	iform		integer		Output format: 0 = radiance
C						       1 = F_plan/F_star
C						       2 = 100*A_plan/A_star
C						       3 = planet spectral flux
C							    i.e. F_plan
C						       4 = Transmission*solar_flux
C						       5 = Transmission (solar occultation)
C	ispace		integer		0=cm-1,1=microns
C	lout		integer		Output unit number
C	ispec		integer		Spectrum ID
C	xlat		real		Latitude
C	xlon		real		Longitude 
C	npro		integer		Number of levels in .ref file
C	nvar		integer		Number of variable profiles
C	varident(mvar,3) integer	Identity of profiles and 
C						parameterisation
C	varparam(mvar,mparam) real 	Extra parameters as required
C	nx		integer		Number of elements in state vector
C	ny		integer		Number of elements in measurement
C						vector
C	y(my)		real		Measured spectrum
C	yn(my)		real		Best calculated spectrum
C	se(my)		real		Variances of measured spectrum
C	xa(mx)		real		A priori measurement vector
C	sa(mx,mx)	real		A priori covariance matrix
C	xn(mx)		real		Retrieved vector x
C	err(mx)		real		Retrieved errors
C	ngeom		integer		Number of observation geometries
C	nconv(mgeom)	integer		Number of convolution wavenumbers
C					 at each observation angle
C	vconv(mgeom,mconv) real		Convolution wavenumbers
C	gasgiant	logical		Gas giant flag
C	jpre		integer		Indicates if pressure retrieval 
C					performed.
C	jrad		integer		Indicates if radiusretrieval 
C					performed.
C	jlogg		integer		Indicates if surface gravity retrieval 
C					performed.
C	jfrac		integer		Indicates if profile fraction average 
C					performed.
C       iscat		integer		Flag to indicate scattering calc.
C	lin		integer		Previous retrieval flag
C
C     Pat Irwin		29/7/96
C     Steve Smith       10/12/96 Added extra space to xn op format
C     Steve Smith       17/3/97 New version for anneal_tng
C     Pat Irwin		17/10/03 Revised and recommented for Nemesis
C     Pat Irwin		11/5/12	 Updated for variable formats and new Nemesis
C     Pat Irwin		30/1/13	Updated to write out last retrieved .prf
C				 files.
C
C     ***********************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/planrad.f'
      include 'arraylen.f'

      integer ny,nx,i,ngeom,igeom,lout,ispec,nsubspec,iform
      real y(my),xa(mx),xn(mx),err(mx),se(my),yn(my)
      real err1,xerr1,sa(mx,mx),xfac,xlat,xlon
      integer nconv(mgeom),j,ioff,varident(mvar,3),nvar,npro
      integer nxtemp,ivar,np,ix,iflag,ispace,npvar
      integer logflag,xflag,jpara,flagh2p,jpre,ncont,npro1
      integer iscat,lin,jrad,jlogg,iplanet,jfrac
      real xa1,ea1,xn1,en1,iav,xdnu,RADIUS,Grav
      parameter (Grav=6.672E-11)
      real vconv(mgeom,mconv),varparam(mvar,mparam)
      real relerr
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      real xmapx(maxv,maxgas+2+maxcon,maxpro)
      real xlatx,varidentx(mvar,3),varparamx(mvar,mparam)
      real stx(mx,mx)
      integer nxx,xnx(mx),nvarx,nprox,jtanx,jprex,jradx
      integer icread,jloggx,ierr,ierrx,jfracx
      character*100 runname,aname,buffer,cellfile
      logical gasgiant,cellexist

      integer iwave,imode,nx1,nvmr,np1,nmode,nwave,max_wave
      integer max_mode,inorm
      parameter (max_wave = 1000,max_mode = 10)
      real r0,v0,clen,vm,nm,lambda0
      real wave(max_wave),xsec(max_mode,max_wave,2),nimag(max_wave)
      real v1(max_wave),k1(max_wave),vm1,n1(max_wave)
      real nreal(max_wave),minlam
      real srefind(max_wave,2),parm(3),rs(3)


1     FORMAT(A)

      write(lout,*)ispace
      write(lout,*)ngeom
      write(lout,*)ny
      write(lout,*)nx
      write(lout,*)xlat
      write(lout,*)xlon
C     We assume all observations have same number of points
      write(lout,*)nconv(1)
      
C     Spectra is in transmission units
      ioff = 0
      do igeom=1,ngeom
        do j=1,nconv(igeom)
          i = ioff+j
          err1 = sqrt(se(i))
          write(lout,1001)vconv(igeom,j)
          write(lout,1001)y(i)
          write(lout,1001)err1
          write(lout,1001)yn(i)
        enddo
        ioff = ioff+nconv(igeom)
      enddo          

      write(lout,*)nvar

      nxtemp=0
      do 299 ivar=1,nvar
        do j=1,3
          write(lout,*)varident(ivar,j)
        enddo
        if(varident(ivar,3).eq.25)then
          do j=1,mparam
            write(lout,*)varparam(ivar,j)
          enddo
        elseif(varident(ivar,1).eq.445)then
          do j=1,6
            write(lout,*)varparam(ivar,j)
          enddo
        else
          do j=1,5
            write(lout,*)varparam(ivar,j)
          enddo
        endif

        np=1
        if(varident(ivar,1).le.100)then
          np = npvar(varident(ivar,3),npro,varparam(ivar,1))
        endif
        if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))
        if(varident(ivar,1).eq.887)np = int(varparam(ivar,1))
        if(varident(ivar,1).eq.444)np = 2+int(varparam(ivar,1))
        if(varident(ivar,1).eq.445)np = 3+(2*int(varparam(ivar,1)))
        if(varident(ivar,1).eq.222)np = 8
        if(varident(ivar,1).eq.223)np = 9
        if(varident(ivar,1).eq.224)np = 9
        if(varident(ivar,1).eq.225)np = 11
        if(varident(ivar,1).eq.226)np = 8
        if(varident(ivar,1).eq.227)np = 7

        write(lout,*)np

        do i = 1,np
          ix = nxtemp+i
          xa1 = xa(ix)
          ea1 = sqrt(abs(sa(ix,ix)))
          xn1 = xn(ix)
          en1 = err(ix)

          iflag = logflag(varident(ivar,1),varident(ivar,3),
     &     varparam(ivar,1),i)
          print*,xa1,ea1,xn1,en1,iflag
          if(iflag.eq.1)then
            xa1 = exp(xa1)
            ea1 = xa1*ea1
            xn1 = exp(xn1)
            en1 = xn1*en1
          endif      

          write(lout,1002)xa1
          write(lout,1002)ea1
          write(lout,1002)xn1
          write(lout,1002)en1

        enddo
        nxtemp = nxtemp+np

299   continue

1001  format(f14.10)
1002  format(e12.5)

      return

      end


