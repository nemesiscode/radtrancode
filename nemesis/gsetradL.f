      subroutine gsetradL(runname,nconv,vconv,fwhm,ispace,iscat,
     1 gasgiant,layht,nlayer,laytyp,layint,xlat,lin,hcorrx,
     2 nvar,varident,varparam,nx,xn,jpre,tsurf,xmap)
C     $Id:
C     ************************************************************************
C     Subroutine to write out the .pat, .prf, .xsc and .sca and aerosol 
C     files needed for a CIRSradg run. Routine also returns xmap
C     (calculated by subprofretg) which relates the functional
C     derivatives calculated by CIRSRADG to the state vector elements
C
C     Input variables
C       runname         character*100    Root run name.
C       nconv           integer         Number of calculation wavelengths
C       vconv(mconv)   	real            Calculation wavelength array
C	fwhm		real		FWHM of convolved spectrum
C	ispace		integer		0=cm-1, 1= wavelengths
C	iscat		integer		Scattering ID
C	gasgiant	logical		Indicates if planet is a gas giant
C	layht		real		Altitude of base layer
C	nlayer		integer		Number of layers
C	laytyp		integer		How layers are separated
C	layint		integer		How layer amounts are integrated
C	xlat		real		latitude of observation
C	lin		integer		Unit number of previous retrieval (if any)
C       nvar    	integer 	Number of variable profiles 
C					  (e.g. gas,T,aerosol)
C       varident(mvar,3) integer 	identity of constituent to retrieved
C					 and parameterisation
C       varparam(mvar,mparam) real 	Additional parameters constraining
C					  profile.
C	nx		integer		Number of elements in state vector
C	xn(mx)		real		state vector
C       jpre            integer         Position of tangent pressure
C                                       in xn
C	tsurf		real		Surface temperature
C
C    Output variables      
C	xmap(maxv,maxgas+2+maxcon,maxpro) real Mapping array relate
C		functional derivatives calculated by CIRSRADG to the
C		state vector elements
C	tsurf		real		Updated surface temperature (if
C					previously retrieved)
C
C     Pat Irwin	17/10/03	Modified from oxcirsg for Nemesis
C
C     ************************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'

      integer nconv,lin,ispace,iscat
      real xlat,fwhm,xlatx,hcorrx,tsurf
      integer nlayer,laytyp,nx,nxx,ncont,jpre
      integer layint,jsurfx,jalbx,jtanx,jprex,nprox
      real layht
      real vconv(mconv)
      integer flagh2p
      double precision mu(maxmu),wtmu(maxmu)
      real xn(mx),xnx(mx),stx(mx,mx)
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      character*100 runname

      integer nvar,varident(mvar,3),i,j,ivar,ivarx
      real varparam(mvar,mparam)
      integer nvarx,varidentx(mvar,3)
      real varparamx(mvar,mparam)
      logical gasgiant


      call subprofretg(runname,ispace,iscat,gasgiant,xlat,nvar,
     1  varident,varparam,nx,xn,jpre,ncont,flagh2p,xmap)

      hcorrx=0.0

      if(lin.eq.1.or.lin.eq.3)then

       call readxtmp(runname,xlatx,nvarx,varidentx,varparamx,nprox,
     1 nxx,xnx,stx,jsurfx,jalbx,jtanx,jprex)

       call stripvar(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  nxx,xnx)

       print*,'gsetradL - variables to be updated from .pre'
       print*,'xlatx = ',xlatx
       print*,'nvarx = ',nvarx
       do i=1,nvarx
        print*,'varidentx',(varidentx(i,j),j=1,3)
        print*,'varparamx',(varparamx(i,j),j=1,mparam)
       enddo
       print*,'gsetradL: Calling subprofretgx'

       call subprofretgx(runname,ispace,iscat,gasgiant,nvarx,varidentx,
     1  varparamx,nxx,xnx,jprex,flagh2p)

       do ivarx=1,nvarx
        if(varidentx(ivarx,1).eq.777)hcorrx=xnx(jtanx)
        if(varidentx(ivarx,1).eq.999)tsurf=xnx(jsurfx)
       enddo

      endif

      call gwritepatL(runname,iscat,nconv,vconv,fwhm,layht,nlayer,
     2 laytyp,layint,flagh2p)

      return

      end
    
