      subroutine coreret(runname,ispace,iscat,ilbl,ica,kiter,phlimit,
     1  fwhm,xlat,xlon,ngeom,nav,nwave,vwave,nconv,vconv,angles,
     2  gasgiant,lin,lpre,nvar,varident,varparam,npro,jsurf,jalb,jxsc,
     3  jtan,jpre,jrad,jlogg,jfrac,wgeom,flat,flon,nx,lx,xa,sa,ny,y,se1,
     4  xn,sm,sn,st,yn,kk,aa,dd)
C     $Id:
C     ******************************************************************
C
C     Input variables
C	runname	character*100	Root name of associated run files
C       ispace           integer Indicates if wavelengths in vconv and 
C                               vwave are in wavenumbers(0) or 
C                               wavelengths (1)
C	iscat	integer	Set to 0 for thermal emission
C			Set to 1 to use plane-parallel scattering RTM
C			Set to 2 if internal radiation field to be 
C				    calculated first for limb/near-limb 
C				    observations.
C			Set to 3 for single-scattering calculations
C	ilbl	integer	Set to 0 for correlated-k caculation
C                       Set to 1 for lbl calculation
C                       Set to 2 for lbl-table calculation
C	ica	integer	1 if single retrieval, 0 otherwise.
C	kiter	integer	Maximum number of iterations
C	phlimit	real	Limiting % change in cost function to consider solution
C			converged.
C	fwhm	real	Required FWHM of final convoluted spectrum
C	xlat	real	Latitude of observed site.
C	xlon	real	Longitude of observed site.
C	ngeom	integer	Number of observation angles at which site is observed
C	nav(ngeom) integer  Number of synthetic spectra required
C                       to simulate each FOV-averaged measurement spectrum.
C	nwave(mgeom) integer Number of 'calculation wavelengths' (tabulated 
C			wavelengths in k-tables covering required
C			wavelength range.
C	vwave(mgeom,mwave) real	'Calculation' wavelengths
C	nconv(mgeom) integer Number of 'convolution wavelengths' (output
C			wavelengths where spectrum has been convolved 
C			with FWHM)
C	vconv(mgeom,mconv) real	'Convolution' wavelengths
C	angles(mgeom,mav,3) real Observation angles of each observation geometry
C			(solar,emission,azimuth)
C			if emission angle < 0, solar angle field holds the
C			tangent altitude (km)
C	gasgiant	logical Flag for gas giant planet
C	lin		integer Indicates if previous retrieval to be used
C			        to constrain profiles of to be used as the 
C				a priori for new retrieval
C	lpre		Previous retrieval unit number
C	nvar	integer	Number of variable profiles (gas,T,aerosol)
C	varident(mvar,3) integer Identity of constituent to retrieved and
C				 	method of parameterisation
C	varparam(mvar,mparam) real Additional parameters constraining
C					profile.
C	jsurf		integer	Position of surface temperature element in
C				xa (if included)
C	jalb		integer	Position of surface albedo spectrum in
C				xa (if included)
C	jxsc		integer	Position of x-section spectrum in
C				xa (if included)
C	jtan		integer	Position of tangent height correction in
C				xa (if included)
C	jpre		integer	Position of tangent pressure in
C				xa (if included)
C	wgeom(mgeom,mav) real	Integration weights 
C	flat(mgeom,mav)	real	Integration point latitudes 
C	flon(mgeom,mav)	real	Integration point longitudes 
C	nx		integer	Number of elements in measurement vector
C	lx(mx)		integer	1 if log, 0 otherwise	
C	xa(mx)		real	a priori state vector
C	sa(mx,mx)	real 	A priori covariance matrix
C	ny	integer	Number of elements in measured spectra array
C	y(my)	real	Measurement vector
C	se1(my)	real	Measured radiance variances
C
C     Output variables
C	xn(mx)	real	best fit state vector
C	st(mx,mx)	real	best fit covariance matrix
C	yn(my)	real	best fit calculated spectra
C	kk(my,mx)	real	Calculated dR/dx matrix
C	aa(mx,mx)	double	Averaging kernels
C	dd(mx,my)	double	Contribution functions
C
C     Pat Irwin		29/4/01
C	  Pat Irwin 10/10/03 conversion for Nemesis
C	  Mahmuda Afrin Badhan	05/08/14  Alternate convergence criteria added to address time-consuming retrievals 
C                                     caused by oscillating solutions
C
C     ************************ VARIABLES *******************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      integer iter,kiter,ica,iscat,i,j,icheck,j1,j2,j3,j4,jsurf
      integer jalb,jalbx,jtan,jtanx,jpre,jprex,ilbl,jrad,jradx
      integer iprfcheck,iplanet,lx(mx),jlogg,jloggx,jxsc,jxscx
      integer jfrac,jfracx
      real phlimit,alambda,xtry,tphi,abstphi
      integer xflag,ierr,ncont,flagh2p,npro1,jpara
      real xdnu,xmap(maxv,maxgas+2+maxcon,maxpro)
      CHARACTER*100 runname,itname,abort,aname,buffer,rdw

      
      real xn(mx),se1(my),se(my,my),calc_phiret,sf(my,my)
      real fwhm,xlat,xlatx,xdiff,xn1(mx),x_out(mx)
      real xlonx,RADIUS,xlon
      integer nprox,nvarx,varidentx(mvar,3),jsurfx,nxx,ix,np,npro
      integer n_alambda
      real st(mx,mx),varparamx(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xnx(mx),stx(mx,mx),ynx(my)
      integer ifix(mx),ifixx(mx)

      integer nvar,varident(mvar,3),lin,lin0,lpre,ispace,nav(mgeom),k
      real varparam(mvar,mparam)

      integer ngeom, nwave(mgeom), nconv(mgeom), nx, ny, nxf,ivar,ivarx
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real xa(mx),kk1(my,mx),sa(mx,mx),y(my),yn(my),kkx(my,mx)
      real yn1(my),s1(mx,mx),kk(my,mx),kkold(my,mx),ynold(my)
      real wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real vwaveT(mwave),vconvT(mconv)
      integer nwaveT,nconvT,npvar,jj,nyfull,tmpvarint,iterred
      logical gasgiant,abexist,ntest,isnan,redwavbool,finalbool

      double precision s1d(mx,mx),sai(mx,mx)
      double precision s1e(my,my),sei(my,my)
      double precision dd(mx,my),aa(mx,mx)

      integer nconvfull(mgeom),nwavefull(mgeom),igeom,igeom2,irank
      integer nrdw(mgeom),rdwindices(mgeom,mconv),rdwi1,rdwi2,maxirank
      real vwavefull(mgeom,mwave),vconvfull(mgeom,mconv),yfull(my)
      real vconvi(mgeom,mconv),vwavei(mgeom,mwave),seold(my,my),phfac
      integer nconvi(mgeom),nwavei(mgeom),rdwindicesi(mgeom,mconv)
      integer rank(mgeom,mconv),nconvold(mgeom),rankdiff

      real phi,ophi,chisq,xchi,oxchi,tmpvar
C     **************************** CODE ********************************

C     ++++++++++++++++++ Read in extra parameters to test vmr profile +++
C     Look to see if the CIA file refined has variable para-H2 or not.
      call file(runname,runname,'cia')
      open(12,FILE=runname,STATUS='OLD')
       read(12,1) aname
       read(12,*) xdnu
       read(12,*) jpara
      close(12)

      flagh2p=0
      ierr=0
      if(jpara.ne.0)then
       flagh2p=1
      endif

C     Read in number of aerosol types from the aerosol.ref file
      OPEN(UNIT=1,FILE='aerosol.ref',STATUS='OLD')
C     First skip header
55     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 55
       READ(BUFFER,*)NPRO1,NCONT
      CLOSE(1)

1     FORMAT(A)

C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C     If a reduced wavelength grid exists, do as many iterations of the retrieval as possible
C     using the reduced grid to save time, and then do the final few iterations using the full grid.
      finalbool = .false.
      call file(runname,rdw,'rdw')
      inquire(file=rdw,exist=redwavbool)
      if(redwavbool)then
C      Store full wavelength grid parameters
       nconvfull=nconv
       vconvfull=vconv
       nwavefull=nwave
       vwavefull=vwave
       nyfull=ny
       yfull=y

C      Read in reduced wavelength grid
       irank=1
       rankdiff = 1
       call readrdw(rdw,nconvfull,vconvfull,ngeom,nconv,vconv,
     1  rdwindices,rank,irank,nconvi,vconvi,rdwindicesi,maxirank,
     2  rankdiff)
       nwave=nconv
       ny=0
       do igeom=1,ngeom
        ny=ny+nconv(igeom)
        vwave(igeom,1:nwave(igeom))=vconv(igeom,1:nconv(igeom))
       enddo
       do i=1,my
        y(i)=0.0
       enddo
       j1=0
       j2=0
       do igeom=1,ngeom
        do i=1,nwave(igeom)
         tmpvarint=rdwindices(igeom,i)+j1
         y(i+j2)=yfull(tmpvarint)
        enddo
        j1=j1+nwavefull(igeom)
        j2=j2+nwave(igeom)
       enddo
      endif

      iterred=1		!iteration at which to commence retrieval

      if(ilbl.eq.0.or.ilbl.eq.2)then
C      Find all the calculation and convolution wavelengths and rank
C      in order

       call rankwave(ngeom,nwave,vwave,nconv,vconv,nwaveT,vwaveT,
     1  nconvT,vconvT)

      endif

C     Initialise s1d
      do i=1,mx
       do j=1,mx
        s1d(i,j)=0.0
       enddo
      enddo

      do i=1,nx
       do j=1,nx
        s1d(i,j)=dble(sa(i,j))
       enddo
      enddo
      
C     Calculate inverse of sa.
      jmod = 2
      icheck=0
      print*,'coreret: Checking can invert Sa'
      call dinvertm(jmod,icheck,s1d,nx,mx,sai)
      if(icheck.eq.1)then
       print*,'************* WARNING *************'
       print*,'Coreret, sa does not invert cleanly'
       do i=1,nx
        print*,i,xa(i),(sa(i,j),j=1,nx)
       enddo
       print*,'Aborting...'
       print*,'***********************************'
       stop
      endif

C     Find if any of the variables have such small error that we can 
C     fix them
      call setifix(xa,sa,nvar,varident,varparam,npro,ifix)

C     Initialise s1e and se
      do i=1,my
       do j=1,my
        sei(i,j)=0.0
        se(i,j)=0.0
       enddo
      enddo

      if(redwavbool.eqv..true.)then
       j1=0
       j2=0
       do igeom=1,ngeom
        do i=1,nconv(igeom)
         rdwi1=rdwindices(igeom,i)+j1
         se(j2+i,j2+i)=se1(rdwi1)
         sei(j2+i,j2+i)=1.0/dble(se1(rdwi1))
        enddo
        j1=j1+nwavefull(igeom)
        j2=j2+nconv(igeom)
       enddo
      else
       do i=1,ny
        se(i,i)=se1(i)
        sei(i,i)=1.0/dble(se1(i))
       enddo
      endif


      CALL readrefiplan(runname,iplanet,xlat,RADIUS)


C     Calculate first spectrum and k-matrix

C     Load state vector with a priori
      do i=1,nx
       xn(i)=xa(i)
      enddo

      if(lin.eq.1.or.lin.eq.3)then

       if(lin.eq.1) then
C        Just substituting parameters from .pre file
         call readraw(lpre,xlatx,xlonx,nprox,nvarx,varidentx,varparamx,
     1  jsurfx,jalbx,jxscx,jtanx,jprex,jradx,jloggx,jfracx,nxx,xnx,stx)
       
        xdiff = abs(xlat-xlatx)
        if(xdiff.gt.lat_tolerance)then
          print*,'Coreret: Aborting - latitudes inconsistent'
          print*,xlatx,xlat
          stop
        endif
         
        do ivarx=1,nvarx
         do ivar=1,nvar
          if(varidentx(ivarx,1).eq.varident(ivar,1))then
           if(varidentx(ivarx,2).eq.varident(ivar,2))then
             print*,'Coreret: Can not use previous retrieval to add'
             print*,'Radiance error, since identity of variable is'
             print*,'identical to one of those being retrieved in this'
             print*,'retrieval'
          print*,'ivar,varident : ',ivar,(varident(ivar,j),j=1,3)
          print*,'ivarx,varidentx : ',ivarx,(varidentx(ivarx,j),j=1,3)
             stop
           endif
          endif
         enddo
        enddo

C       Write out x-data to temporary .str file for later routines.
        call writextmp(runname,xlatx,xlonx,nvarx,varidentx,varparamx,
     1   nprox,nxx,xnx,stx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,jloggx,
     2   jfracx)

       else
C       substituting and retrieving parameters from .pre file. 
C       Current record from .pre file already read in by
C       readapriori.f. Hence just read in from temporary .str file
        call readxtmp(runname,xlatx,xlonx,nvarx,varidentx,varparamx,
     1   nprox,nxx,xnx,stx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,
     2   jloggx,jfracx)

       endif
 
C      Calc. gradient of all elements of xnx matrix.
       do i=1,nxx
        ifixx(i)=1
       enddo

       lin0 = 0

       if(ilbl.eq.1)then
         print*,'Calling forwardnoglbl - A'
         CALL forwardnoglbl(runname,ispace,iscat,fwhm,ngeom,nav,
     1    wgeom,flat,flon,nconv,vconv,angles,gasgiant,lin0,
     2    nvarx,varidentx,varparamx,jsurfx,jalbx,jxscx,jtanx,jprex,
     3    jradx,jloggx,jfracx,RADIUS,nxx,xnx,ifixx,ny,ynx,kkx,kiter)
       else
 
       if(iscat.eq.0)then
         if(jalbx.lt.0.and.jxscx.lt.0)then
          print*,'Calling forwardavfovX - A'
          CALL forwardavfovX(runname,ispace,iscat,fwhm,ngeom,
     1     nav,wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,
     2     lin0,nvarx,varidentx,varparamx,jsurfx,jalbx,jxscx,jtanx,
     3     jprex,jradx,jloggx,jfracx,RADIUS,nxx,xnx,ny,ynx,kkx)
         else
          print*,'Calling forwardnogX - A'
          CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1     wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin0,
     2     nvarx,varidentx,varparamx,jsurfx,jalbx,jxscx,jtanx,jprex,
     3     jradx,jloggx,jfracx,RADIUS,nxx,xnx,ifixx,ny,ynx,kkx,kiter,
     4	   iprfcheck)
         endif
       elseif(iscat.eq.1.or.iscat.eq.3.or.iscat.eq.4)then
         print*,'Calling forwardnogX - A1. iscat = ',iscat
         CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1    wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin0,
     2    nvarx,varidentx,varparamx,jsurfx,jalbx,jxscx,jtanx,jprex,
     3    jradx,jloggx,jfracx,RADIUS,nxx,xnx,ifixx,ny,ynx,kkx,kiter,
     4	  iprfcheck)
         print*,'forwardnogX - A1 - called OK'
       elseif(iscat.eq.2)then
         print*,'Calling intradfield - A1'
         CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1    vconvT,gasgiant,lin0,nvarx,varidentx,
     2    varparamx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,jloggx,
     3    jfracx,RADIUS,nxx,xnx)
         print*,'intradfield called OK'
         CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1    wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin0,
     2    nvarx,varidentx,varparamx,jsurfx,jalbx,jxscx,jtanx,jprex,
     3    jradx,jloggx,jfracx,RADIUS,nxx,xnx,ifixx,ny,ynx,kkx,kiter,
     4	  iprfcheck)
       else
         print*,'Coreret: iscat invalid',iscat
         stop
        endif
       endif

       if(lin.eq.3) then

C        strip out variables from kkx and stx that will be retrieved in 
C        this run.
         call scankkx(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  varparam,kkx,stx,nxx)
       endif

       print*,'Calc forward model error'
       call calcfwderr(nxx,ny,kkx,stx,sf)

C      Add effect of previous retrieval errors to measurement covariance
C      matrix

       print*,'Writing sef file'
       call file(runname,runname,'sef')
       open(35,file=runname,status='unknown')
       write(35,*)'Additional measurement covariance matrix due to'
       write(35,*)'uncertainty in previous retrieval'
       write(35,*)ny,'   ! ny'
       do i=1,ny
        write(35,*)(sf(i,j),j=i,ny)
        do j=1,ny       
         se(i,j)=se(i,j)+sf(i,j)
        enddo
       enddo

       close(35)

C      Recalculate inverse of se
       do i=1,ny
        do j=1,ny
         s1e(i,j)=dble(se(i,j))
        enddo
       enddo

       print*,'Calculating inverse'
C      Calculate inverse of se
       jmod = 2
       icheck=0
       call dinvertm(jmod,icheck,s1e,ny,my,sei)
       if(icheck.eq.1)then
        print*,'************* WARNING *************'
        print*,'Coreret, se does not invert cleanly'
        print*,'Aborting...'
        print*,'***********************************'
        stop
       endif
       print*,'Inverse OK'
      endif 


198   if(ilbl.eq.1)then
         print*,'Calling forwardnoglbl - B'
         CALL forwardnoglbl(runname,ispace,iscat,fwhm,ngeom,nav,
     1    wgeom,flat,flon,nconv,vconv,angles,gasgiant,lin,
     2    nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,
     3    jrad,jlogg,jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter)
         print*,'call OK'
      else

       if(iscat.eq.0)then
        if(jalb.lt.0.and.jxsc.lt.0)then
         print*,'Calling forwardavfovX - B'
         CALL forwardavfovX(runname,ispace,iscat,fwhm,ngeom,nav,
     1    wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2    nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3    jlogg,jfrac,RADIUS,nx,xn,ny,yn,kk)
        else
         print*,'Calling forwardnogX - B' 
         CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1    wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2    nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3    jlogg,jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter,iprfcheck)
        endif

C        print*,'forwardavfovX OK, jpre = ',jpre

       elseif(iscat.eq.1.or.iscat.eq.3.or.iscat.eq.4)then

        print*,'Calling forwardnogX - B1, iscat = ',iscat

        CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3   jlogg,jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter,iprfcheck)

C        print*,'forwardnogX OK, jpre = ',jpre


       elseif(iscat.eq.2)then

        print*,'Calling intradfield - B',gasgiant
        CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1   vconvT,gasgiant,lin,nvar,varident,varparam,jsurf,jalb,
     2   jxsc,jtan,jpre,jrad,jlogg,jfrac,RADIUS,nx,xn)
        print*,'intradfield - B called OK'
        CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3   jlogg,jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter,iprfcheck)
        print*,'B - OK'
       endif
      endif

199   open(12,file='kk.dat',status='unknown')
      write(12,*)nx,ny
      do i=1,ny
       write(12,*)(kk(i,j),j=1,nx)
      enddo
      close(12)


C     Now calculate the gain matrix and averaging kernels
      call calc_gain_matrix(nx,ny,kk,sa,sai,se,sei,dd,aa)

C     Calculate initial value of cost function phi.
      phi = calc_phiret(ny,y,yn,sei,nx,xn,xa,sai,chisq)
      ophi = phi
      oxchi = chisq/float(ny)
c      write(*,*)ophi,oxchi

C     Assess whether retrieval is likely to be OK
      call assess(nx,ny,kk,sa,se)
 
      if(ica.eq.1)then       ! Open and write only for single spec ret.
       inquire(file=itname, number=tmpvarint)!check that .itr file hasn't already been opened before opening it
       if(tmpvarint.ne.37)then
        call file(runname,itname,'itr')
        open(37,file=itname,status='unknown')
        write(37,*)'' 
       else
        write(37,*)'###'
       endif
      endif

C     alambda is a Marquardt-Levenberg-type 'braking parameter'
      
      alambda = 1.0

C     Set the trial vectors xn1, and yn1 to be the same as the initial
C     vectors xn, yn
      do i=1,nx
       xn1(i)=xn(i)
      enddo
      do i=1,ny
       yn1(i)=yn(i)
      enddo

      print*,'iterred,kiter',iterred,kiter

      do 401 iter = iterred, kiter

        if(ica.eq.1)then
         write(37,*)''
         write(37,*)nx,ny,kiter
         write(37,*)chisq,phi
         write(37,*)(xn1(i),i=1,nx)
         write(37,*)(xa(i),i=1,nx)
         write(37,*)(y(i),i=1,ny)
         write(37,*)(se1(i),i=1,ny)
         write(37,*)(yn1(i),i=1,ny)
         write(37,*)(yn(i),i=1,ny)
         do i=1,nx
          write(37,*)(kk(j,i),j=1,ny)
         enddo
        endif

C       Now calculate next iterated xn1
        call calcnextxn(nx,ny,xa,xn,y,yn,dd,aa,x_out)

        do i=1,nx
         xn1(i) = xn(i) + (x_out(i)-xn(i))/(1.0+alambda)

C        Add additional brake for model 102 to stop silly fractions.
         if(jfrac.gt.0)then
          if(xn1(jfrac).lt.0.01.or.xn1(jfrac).gt.0.99)then
           alambda=alambda*10
           if(alambda.gt.1e10)alambda=1e10
           goto 401
          endif
         endif

         print*,'i, x_old, x_next',i,xn(i),xn1(i)
         ntest=isnan(xn1(i))
         if(ntest)then
          print*,'NAN detected. Increase brake and try again'
          alambda=alambda*10.
          if(alambda.gt.1e10)alambda=1e10
          goto 401
         endif
C        Check to see if log numbers have gone out of range
         if(lx(i).eq.1)then
          if(xn1(i).gt.85.or.xn1(i).lt.-85)then
           print*,'Coreret - log(number gone out of range)'
           print*,'Increasing brake'
           alambda=alambda*10
           if(alambda.gt.1e10)alambda=1e10
           goto 401
          else
           print*,'exp(x_old),exp(x_next)',exp(xn(i)),exp(xn1(i))
          endif
         endif
        enddo

C       Test to see if any vmrs have gone negative.
        xflag=0
        call subprofretg(xflag,runname,ispace,iscat,gasgiant,xlat,
     1    xlon,nvar,varident,varparam,nx,xn1,jpre,ncont,flagh2p,
     2    xmap,ierr)
        if (ierr.eq.1)then
          alambda = alambda*10.0             ! increase Marquardt brake
          if(alambda.gt.1e10)alambda=1e10
          goto 401
        endif


        ix=1
        do ivar = 1,nvar
         np=1
         if(varident(ivar,1).le.100)then
           np=npvar(varident(ivar,3),npro,varparam(ivar,1))
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

         do j=ix,ix+np-1
          if(varident(ivar,1).eq.0)then
           if(varident(ivar,3).eq.0) then
            if(xn1(j).lt.1.0) then
             print*,'Temperature has gone negative, Increase alambda'
             alambda = alambda*10.0		! increase Marquardt brake
             if(alambda.gt.1e10)alambda=1e10
             goto 401
            endif
           endif

           if(varident(ivar,3).eq.16.and.j.eq.ix) then
            if(xn1(j).lt.1.0) then
             print*,'Temperature has gone negative, Increase alambda'
             alambda = alambda*10.0		! increase Marquardt brake
             if(alambda.gt.1e10)alambda=1e10
             goto 401
            endif
           endif
          endif

          if(varident(ivar,1).eq.444.and.j.eq.ix+1)then
           if(exp(xn1(j)).lt.0.01)then
             print*,'Variance of size distribution gone too small',
     1		exp(xn1(j))
             print*,'Increase alambda',alambda
             alambda = alambda*10.0		! increase Marquardt brake
             if(alambda.gt.1e10)alambda=1e10
             goto 401
           endif
          endif

          if(varident(ivar,1).eq.444)then
           if(j.ge.ix+2.and.j.le.ix+np-1)then
            if(exp(xn1(j)).gt.1)then
             print*,'Imaginary refractive index too large',
     1		exp(xn1(j))
             print*,'Increase alambda',alambda
             alambda = alambda*10.0		! increase Marquardt brake
             if(alambda.gt.1e10)alambda=1e10
             goto 401
            endif
           endif
          endif

          if(varident(ivar,1).eq.445.and.j.eq.ix+1)then
           if(exp(xn1(j)).lt.0.01)then
             print*,'Variance of size distribution gone too small',
     1		exp(xn1(j))
             print*,'Increase alambda',alambda
             alambda = alambda*10.0		! increase Marquardt brake
             if(alambda.gt.1e10)alambda=1e10
             goto 401
           endif
          endif

          if(varident(ivar,1).eq.445.and.j.eq.ix+2)then
           if(xn1(j).ge.0)then
             print*,'Particle shell:core ratio greater than 1',
     1		exp(xn1(j))
             print*,'Increase alambda',alambda
             alambda = alambda*10.0		! increase Marquardt brake
             if(alambda.gt.1e10)alambda=1e10
             goto 401
           endif
          endif

          if(varident(ivar,1).eq.227.and.j.eq.ix)then!if TC altitude is too low (MIGHT NEED TO CHANGE THIS FOR OTHER WAVELENGTHS)
           if(exp(xn1(j)).gt.5)then
              print*,'TC altitude too low:',exp(xn1(j))
              print*,'Increase alambda',alambda
              alambda = alambda*10.0		! increase Marquardt brake
              if(alambda.gt.1e10)alambda=1e10
              goto 401
           endif
          endif

          if(varident(ivar,1).eq.227.and.j.eq.ix+1)then!if TC opacity becomes too large
           if(exp(xn1(j)).ge.1000)then
              print*,'TC opacity too high:',exp(xn1(j))
              print*,'Increase alambda',alambda
              alambda = alambda*10.0		! increase Marquardt brake
              if(alambda.gt.1e10)alambda=1e10
              goto 401
           endif
          endif

          if(varident(ivar,1).eq.227.and.j.eq.ix+3)then!if CB altitude goes above SH altitude or below TC base altitude, or if TC becomes too narrow
           if(varparam(ivar,1).eq.0)then
            tmpvar = 0.9
           else
            tmpvar = varparam(ivar,1)
           endif
           if(exp(xn1(j-3)).lt.1.5*exp(xn1(j)).or.
     1       (tmpvar*exp(xn1(j)))-0.02.lt.exp(xn1(j+2)))then
              print*,'Breakdown of cloud structure:'
              print*,'TC base pressure = ',exp(xn1(j-3))
              print*,'CB pressure = ',exp(xn1(j))
              print*,'SH pressure = ',exp(xn1(j+2))
              print*,'Increase alambda',alambda
              alambda = alambda*10.0		! increase Marquardt brake
              if(alambda.gt.1e10)alambda=1e10
              goto 401
           endif
          endif


         enddo
         ix=ix+np
        enddo


C       Calculate test spectrum using trial state vector xn1. 
C       Put output spectrum into temporary spectrum yn1 with
C       temporary kernel matrix kk1. Does it improve the fit? 

        if(ilbl.eq.1)then
         print*,'Calling forwardnoglbl - C'
         CALL forwardnoglbl(runname,ispace,iscat,fwhm,ngeom,nav,     
     1    wgeom,flat,flon,nconv,vconv,angles,gasgiant,lin,
     2    nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3    jlogg,jfrac,RADIUS,nx,xn1,ifix,ny,yn1,kk1,kiter)
        else
         if(iscat.eq.0)then
          if(jalb.lt.0.and.jxsc.lt.0)then
           print*,'Calling forwardavfovX - C'
           CALL forwardavfovX(runname,ispace,iscat,fwhm,ngeom,nav,
     1      wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,
     2      lin,nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,
     3      jrad,jlogg,jfrac,RADIUS,nx,xn1,ny,yn1,kk1)
          else
           print*,'Calling forwardnogX - C'
           CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1      wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,
     2      lin,nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3      jlogg,jfrac,RADIUS,nx,xn1,ifix,ny,yn1,kk1,kiter,iprfcheck)
          endif

         elseif(iscat.eq.1.or.iscat.eq.3.or.iscat.eq.4)then

         print*,'Calling forwardnogX - C'
          CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1     wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2     nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3     jlogg,jfrac,RADIUS,nx,xn1,ifix,ny,yn1,kk1,kiter,iprfcheck)
         elseif(iscat.eq.2)then
         print*,'Calling intradfield - C'
          CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1     vconvT,gasgiant,lin,nvar,varident,varparam,jsurf,jalb,
     2     jxsc,jtan,jpre,jrad,jlogg,jfrac,RADIUS,nx,xn1)
          CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1     wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2     nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3     jlogg,jfrac,RADIUS,nx,xn1,ifix,ny,yn1,kk1,kiter,iprfcheck)
         endif
        endif

        if(iprfcheck.eq.1)then
C        iteration has led to negative temperatures, vmrs or dust amounts
C        Increase brake and try again.
         print*,'Profile gone wobbly, increase brake and try again'
         alambda = alambda*10.0              ! increase Marquardt brake
         if(alambda.gt.1e10)alambda=1e10
         goto 401
        endif

C       Calculate the cost function for this trial solution.
        phi = calc_phiret(ny,y,yn1,sei,nx,xn1,xa,sai,chisq)

        xchi = chisq/float(ny)
        print*,'chisq/ny = ',xchi
        print*,'it.,al.,ophi.,phi.',
     1   iter,alambda,ophi,phi

C       What's %phi between last and this iteration?    
        tphi = 100.0*(ophi-phi)/ophi
        abstphi = abs(tphi)
        print*,'%phi, abs(%phi) : ',tphi,abstphi

        if(redwavbool.eqv..true.)then!reduce phlimit in first few iterations of reduced wavelength scheme using a fudge factor to waste less time
         phfac=float(nconv(1))/float(nconvfull(1))
        else
         phfac=1.0
        endif

C       Does trial solution fit the data better?
        if(phi.le.ophi)then
            print*,'Successful iteration. Updating xn,yn and kk'
            do i=1,nx
             xn(i)=xn1(i)         		! update xn to new value
            enddo
            do i=1,ny
             yn(i)=yn1(i)				! update yn and kk
             do j=1,nx
              kk(i,j)=kk1(i,j)
             enddo
            enddo

            print*,'Calculating new gain matrix and averaging kernels'
C           Now calculate the gain matrix and averaging kernels
            call calc_gain_matrix(nx,ny,kk,sa,sai,se,sei,dd,aa)

            print*,'calc_gain_matrix OK'
          
C           Has solution converged?
          
            if(tphi.ge.0.0.and.tphi.le.phlimit/phfac.and.
     1       alambda.lt.1.0)then
              print*,'%phi, phlimit : ',tphi,phlimit/phfac
              print*,'Phi has converged'
              print*,'Terminating retrieval'
              oxchi=xchi
              if((redwavbool.eqv..true.).and.(irank.lt.maxirank))then
               GOTO 201
              else
               GOTO 202     
              endif              
            else
              ophi=phi
              oxchi = xchi
              alambda = alambda*0.3		! reduce Marquardt brake
            endif
          
        elseif (iter.ge.5.and.abstphi.le.phlimit/phfac)then
C       If phi > ophi, accept new xn and kk only if current solution 
C       would converge under one of the alternate criterions:

C           If alambda is small enough, increase it to decrease abs(tphi) value. 						
            if (alambda.lt.0.1) then        ! don't allow lambda to increase beyond 1.0
	       alambda = alambda*10.0		! increase Marquardt brake further
            else
C              If lambda is close to 1.0 or greater when condition met, accept that iteration.
               print*,'Accepting iteration. Updating xn,yn and kk'
               do i=1,nx
	        xn(i)=xn1(i)         		! update xn to new value
	       enddo

	       do i=1,ny
		 yn(i)=yn1(i)				! update yn and kk
		 do j=1,nx
		  kk(i,j)=kk1(i,j)
		 enddo
	       enddo

	       print*,'%phi, phlimit, alambda : ',tphi,phlimit,alambda
	       print*,'Phi has converged under the alternate criteria'
	       if (alambda.ge.1.0) then
		 print*,'In addition, alambda is >= 1.0'
	       endif	
	       print*,'Terminating retrieval under alternate conditions'
               oxchi=xchi
               if((redwavbool.eqv..true.).and.(irank.lt.maxirank))then
                GOTO 201
               else
                GOTO 202     
               endif   
	    endif														
							
C	    If alternate criterions aren't met either, leave xn and kk alone and try again with more braking
	else
            alambda = alambda*10.0		! increase Marquardt brake
            if(alambda.gt.1e10)alambda=1e10
        endif


C       If a reduced wavelength grid is being used but the retrieval still hasn't converged to a solution, 
C       perform the last 2 iterations using the highest-ranked grid
        if((redwavbool.eqv..true.).and.(iter.ge.kiter-2).and.
     1   (finalbool.eqv..false.))then
         if((iter.ge.2).and.(irank.lt.maxirank)) then 
          rankdiff = maxirank - irank
          irank = maxirank - 1
          finalbool = .true.
          print*, 'Replacing current wavelength grid with'
          print*, 'highest-ranked grid for last few iterations'
          goto 201
         else
          goto 202
         endif
        endif


        call file(runname,runname,'abo')
        inquire(file=runname,exist=abexist)
        if(abexist)then
         open(83,file=runname,status='old')
         read(83,'(A)')abort
         close(83)
         if(abort.eq.'stop'.or.abort.eq.'STOP')then
           print*,'Terminating retrieval'
           GOTO 202
         endif
        endif
       
401   continue       

C     If using the reduced wavelength scheme
201   if(redwavbool.eqv..true.)then

       irank = irank + 1
      
C      Extend wavelength grid to next-ranked wavelengths, re-initialise measurement vector/covariance matrix
C      and run initial forward model again using extended grid
       if(irank.le.maxirank)then
       nconvold = nconv
       call readrdw(rdw,nconvfull,vconvfull,ngeom,nconv,vconv,
     1  rdwindices,rank,irank,nconvi,vconvi,rdwindicesi,maxirank,
     2  rankdiff)
        nwave=nconv
        nwavei=nconvi
        ny=0
        do igeom=1,ngeom
         ny=ny+nconv(igeom)
         vwave(igeom,1:nwave(igeom))=vconv(igeom,1:nconv(igeom))
         vwavei(igeom,1:nwavei(igeom))=vconvi(igeom,1:nconvi(igeom))
        enddo
        do i=1,my
         y(i)=0.0
        enddo
        j1=0
        j2=0
        do igeom=1,ngeom
         do i=1,nwave(igeom)
          tmpvarint=rdwindices(igeom,i)+j1
          y(i+j2)=yfull(tmpvarint)
         enddo
         j1=j1+nwavefull(igeom)
         j2=j2+nwave(igeom)
        enddo

C       Re-initialise s1e and se
        do i=1,my
         do j=1,my
          sei(i,j)=0.0
          se(i,j)=0.0
         enddo
        enddo

        j1=0
        j2=0
        do igeom=1,ngeom
         do i=1,nconv(igeom)
          rdwi1=rdwindices(igeom,i)+j1
          se(j2+i,j2+i)=se1(rdwi1)
          sei(j2+i,j2+i)=1.0/dble(se1(rdwi1))
         enddo
         j1=j1+nwavefull(igeom)
         j2=j2+nconv(igeom)
        enddo
        

        iterred=iter+1

        if((iscat.eq.1.or.iscat.ge.3).and.(ilbl.ne.1))then

         kkold = kk
         ynold = yn

         print*,'Performing forward model on new wavelengths'
         print*,'Calling forwardnogX - D'
         CALL forwardnogXrdw(runname,ispace,iscat,fwhm,ngeom,nav,
     1     wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2     nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3     jlogg,jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter,iprfcheck,
     4     rdwindicesi,kkold,nconvi,vconvi,ynold)

           if(ilbl.eq.0.or.ilbl.eq.2)then
C          Find all the calculation and convolution wavelengths and rank
C          in order

            call rankwave(ngeom,nwave,vwave,nconv,vconv,nwaveT,vwaveT,
     1       nconvT,vconvT)

           endif
           goto 199
        else
         print*,'Warning coreret.f: '
         print*,'forwardnoglbl.f, forwardavfovX.f and intradfield.f'
         print*,'not yet completely optimised for reduced wavelength'
         print*,'scheme. Refer to forwardnogXrdw.f for reference.'

         if(ilbl.eq.0.or.ilbl.eq.2)then
C        Find all the calculation and convolution wavelengths and rank
C        in order

          call rankwave(ngeom,nwave,vwave,nconv,vconv,nwaveT,vwaveT,
     1     nconvT,vconvT)

         endif
         goto 198
        endif
C      
       else
        
        print*, 'Highest-ranked wavelength grid reached. Ending run.'
        goto 202         

       endif

      endif

202   if(ica.eq.1)close(37)

      if(ica.eq.1)then
C      Write out k-matrix for reference
       OPEN(52,file='kk.out',form='unformatted',status='unknown')
       write(52)y,yn
       write(52)kk
       CLOSE(52)

       close(37)

      endif

      if(redwavbool.eqv..true.)then
C	Make sure input and output spectral errors are consistent when ending run
        j2=0
        do igeom=1,ngeom
         do i=1,nconv(igeom)
          se1(j2+i) = se(j2+i,j2+i)
         enddo
         j2=j2+nconv(igeom)
        enddo

	do i=1,nconv(1)
	 print*, 'Se=', vconv(1,i), se(i,i)
	 print*, 'Se1=', vconv(1,i), se1(i)
	enddo
      endif

C      print*,'chisq/ny is equal to : ',chisq/float(ny)
      print*,'chisq/ny is equal to : ',oxchi
      if(chisq.gt.ny)then
       print*,'Coreret: WARNING'
       print*,'chisq/ny should be less than 1 if correctly retrieved'
      endif

      print*,'Calculating final covariance matrix'
      CALL calc_serr(nx,ny,sa,se,aa,dd,st,sn,sm)
      print*,'Matrix calculated'

C     Make sure errors stay as a priori for kiter < 0
      if(kiter.lt.0)then
       do i=1,nx
        do j=1,nx
         st(i,j)=sa(i,j)
        enddo
       enddo
      endif
      
      return

      end


