      subroutine coreretdisc(runname,ispace,iscat,ica,kiter,phlimit,
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
C       lx(mx)          integer 1 if log, 0 otherwise
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
C			10/10/03 conversion for Nemesis
C
C     ************************ VARIABLES *******************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      integer iter,kiter,ica,iscat,i,j,icheck,j1,j2,jsurf
      integer jalb,jalbx,jtan,jtanx,jpre,jprex,jrad,jradx
      integer lx(mx),jlogg,jloggx,jxsc,jxscx,jfrac,jfracx
      real phlimit,alambda,xtry,tphi
      integer xflag,ierr,ncont,flagh2p,npro1,jpara
      real xdnu,xmap(maxv,maxgas+2+maxcon,maxpro)
      CHARACTER*100 runname,itname,abort,aname,buffer

      real xn(mx),se1(my),se(my,my),calc_phiret,sf(my,my)
      real fwhm,xlat,xlatx,xdiff,xn1(mx),x_out(mx)
      real xlonx,xlon
      integer nprox,nvarx,varidentx(mvar,3),jsurfx,nxx,ix,np,npro
      real st(mx,mx),varparamx(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xnx(mx),stx(mx,mx),ynx(my)

      integer nvar,varident(mvar,3),lin,lin0,lpre,ispace,nav(mgeom),k
      real varparam(mvar,mparam),RADIUS
      integer iplanet

      integer ngeom, nwave(mgeom), nconv(mgeom), nx, ny, nxf,ivar,ivarx
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real xa(mx),kk1(my,mx),sa(mx,mx),y(my),yn(my),kkx(my,mx)
      real yn1(my),s1(mx,mx),kk(my,mx)
      real wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real vwaveT(mwave),vconvT(mconv)
      integer nwaveT,nconvT,npvar
      logical gasgiant,abexist

      double precision s1d(mx,mx),sai(mx,mx)
      double precision s1e(my,my),sei(my,my)
      double precision dd(mx,my),aa(mx,mx)

      real phi,ophi,chisq,xchi,oxchi
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet
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

C     Find all the calculation and convolution wavelengths and rank
C     in order

      call rankwave(ngeom,nwave,vwave,nconv,vconv,nwaveT,vwaveT,
     1 nconvT,vconvT)

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
      icheck=0
      jmod=2
      call dinvertm(jmod,icheck,s1d,nx,mx,sai)
      if(icheck.eq.1)then
       print*,'************* WARNING *************'
       print*,'Coreret, sa does not invert cleanly'
       print*,'Aborting...'
       print*,'***********************************'
       stop
      endif


C     Initialise s1e and se
      do i=1,my
       do j=1,my
        sei(i,j)=0.0
        se(i,j)=0.0
       enddo
      enddo

      do i=1,ny
        se(i,i)=se1(i)
        sei(i,i)=1.0/dble(se1(i))
      enddo

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
             if(idiag.gt.0)then
              print*,'Coreret: Can not use previous retrieval to add'
              print*,'Radiance error, since identity of variable is'
              print*,'identical to one of those being retrieved in this'
              print*,'retrieval'
          print*,'ivar,varident : ',ivar,(varident(ivar,j),j=1,3)
          print*,'ivarx,varidentx : ',ivarx,(varidentx(ivarx,j),j=1,3)
             endif
             stop
           endif
          endif
         enddo
        enddo

C       Write out x-data to temporary .str file for later routines.
        call writextmp(runname,xlatx,xlonx,nvarx,varidentx,varparamx,
     1   nprox,nxx,xnx,stx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,
     2   jloggx,jfracx)

       else
C       substituting and retrieving parameters from .pre file. 
C       Current record frrom .pre file already read in by
C       readapriori.f. Hence just read in from temporary .str file
 
        call readxtmp(runname,xlatx,xlonx,nvarx,varidentx,varparamx,
     1   nprox,nxx,xnx,stx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,
     2   jloggx,jfracx)

       endif
 
       lin0 = 0

       if(iscat.eq.0)then
        if(idiag.gt.0)print*,'Calling forwarddisc - A'
        CALL forwarddisc(runname,ispace,iscat,fwhm,ngeom,
     1   nav,wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,
     2   lin0,nvarx,varidentx,varparamx,jsurfx,jalbx,jxscx,jtanx,
     3   jprex,jradx,jloggx,jfracx,RADIUS,nxx,xnx,ny,ynx,kkx)
       else
        if(idiag.gt.0)print*,'Option not supported for disc-averaging'
       endif

       if(lin.eq.3) then
C        strip out variables from kkx and stx that will be retrieved in 
C        this run.
         call scankkx(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  varparam,kkx,stx,nxx)
       endif

       call calcfwderr(nxx,ny,kkx,stx,sf)

C      Add effect of previous retrieval errors to measurement covariance
C      matrix

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

      endif

      if(iscat.eq.0)then
        if(idiag.gt.0)print*,'Calling forwarddisc - B'
        CALL forwarddisc(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3   jlogg,jfrac,RADIUS,nx,xn,ny,yn,kk)

      else

         if(idiag.gt.0)print*,'Option not supported'

      endif

      if(idiag.gt.0)then
       open(12,file='kk.dat',status='unknown')
       write(12,*)nx,ny
       do i=1,ny
        write(12,*)(kk(i,j),j=1,nx)
       enddo
       close(12)
      endif

C     Now calculate the gain matrix and averaging kernels
      call calc_gain_matrix(nx,ny,kk,sa,sai,se,sei,dd,aa)

C     Calculate initial value of cost function phi.
      phi = calc_phiret(ny,y,yn,sei,nx,xn,xa,sai,chisq)
      ophi = phi
      oxchi = chisq/float(ny)

C     Assess whether retrieval is likely to be OK
      call assess(nx,ny,kk,sa,se)
 
      if(ica.eq.1)then		! Open and write only for single spec ret.
       call file(runname,itname,'itr')
       open(37,file=itname,status='unknown')
       write(37,*)nx,ny,kiter
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

      do 401 iter = 1, kiter

        if(ica.eq.1)then
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

C       x_out(nx) is the next iterated value of xn using classical N-L
C       optimal estimation. However, we want to apply a braking parameter
C       alambda to stop the new trial vector xn1 being too far from the
C       last 'best-fit' value xn

145     continue
        do i=1,nx
         xn1(i) = xn(i) + (x_out(i)-xn(i))/(1.0+alambda)


C        Add additional brake for model 102 to stop silly fractions.
         if(jfrac.gt.0)then
          if(xn1(jfrac).lt.0.01.or.xn1(jfrac).gt.0.99)then
           alambda=alambda*10
           if(alambda.gt.1e10)alambda=1e10
           goto 145
          endif
         endif

C        Check to see if log numbers have gone out of range
         if(lx(i).eq.1)then
          if(xn1(i).gt.85.or.xn1(i).lt.-85)then
           if(idiag.gt.0)then
            print*,'Coreretdisc - log(number gone out of range)'
            print*,'Increasing brake'
           endif
           alambda = alambda*10.0               ! increase Marquardt brake
           if(alambda.gt.1e30)then
            print*,'Death spiral - stopping'
            stop
           endif
           goto 145
          endif
         endif

        enddo

C       Test to see if any vmrs have gone negative.
        xflag=0
        call subprofretg(xflag,runname,ispace,iscat,gasgiant,xlat,xlon,
     1    nvar,varident,varparam,nx,xn1,jpre,ncont,flagh2p,xmap,ierr)
        if (ierr.eq.1)then
          alambda = alambda*10.0             ! increase Marquardt brake
          if(alambda.gt.1e10)alambda=1e10
          goto 145
        endif

        ix=1
        do ivar = 1,nvar
         np=1
         if(varident(ivar,1).le.100)then
           np=npvar(varident(ivar,3),npro,varparam(ivar,1))
         endif
         if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))
         if(varident(ivar,1).eq.887)np = int(varparam(ivar,1))
         if(varident(ivar,1).eq.444)then
          if(varparam(ivar,2).gt.0.0)then
           np = 2+int(varparam(ivar,1))
          else
           np = 3
          endif
         endif
         if(varident(ivar,1).eq.446)then
          if(varparam(ivar,2).gt.0.0)then
           np = 3+int(varparam(ivar,1))
          else
           np = 5
          endif
         endif

         if(varident(ivar,1).eq.445)np = 3+(2*int(varparam(ivar,1)))

         do j=ix,ix+np-1
          if(varident(ivar,1).eq.0)then
           if(varident(ivar,3).eq.0) then
            if(xn1(j).lt.1.0) then
             if(idiag.gt.0)then
              print*,'Temperature has gone negative, Increase alambda'
             endif
             alambda = alambda*10.0             ! increase Marquardt brake
             if(alambda.gt.1e10)alambda=1e10
             goto 145
            endif
           endif
           if(varident(ivar,3).eq.16.and.j.eq.ix) then
            if(xn1(j).lt.1.0) then
             if(idiag.gt.0)then
              print*,'Temperature has gone negative, Increase alambda'
             endif
             alambda = alambda*10.0             ! increase Marquardt brake
             if(alambda.gt.1e10)alambda=1e10
             goto 145
            endif
           endif
          endif
         enddo
         ix=ix+np
        enddo
C       Calculate test spectrum using trial state vector xn1. 
C       Put output spectrum into temporary spectrum yn1 with
C       temporary kernel matrix kk1. Does it improve the fit? 

        if(iscat.eq.0)then
          if(idiag.gt.0)print*,'Calling forwarddisc - C'
          CALL forwarddisc(runname,ispace,iscat,fwhm,ngeom,nav,
     1     wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,
     2     lin,nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,
     3     jrad,jlogg,jfrac,RADIUS,nx,xn1,ny,yn1,kk1)
        else
          if(idiag.gt.0)print*,'Option not supported'
        endif

C       Calculate the cost function for this trial solution.
        phi = calc_phiret(ny,y,yn1,sei,nx,xn1,xa,sai,chisq)

        xchi = chisq/float(ny)
        if(idiag.gt.0)print*,'chisq/ny = ',xchi
        if(idiag.gt.0)print*,'it.,al.,ophi.,phi.',
     1   iter,alambda,ophi,phi

C       Does trial solution fit the data better?
        if(phi.le.ophi)then
          if(idiag.gt.0)then
           print*,'Successful iteration. Updating xn,yn and kk'
          endif
          do i=1,nx
           xn(i)=xn1(i)         		! update xn to new value
          enddo
          do i=1,ny
           yn(i)=yn1(i)				! update yn and kk
           do j=1,nx
            kk(i,j)=kk1(i,j)
           enddo
          enddo

C         Now calculate the gain matrix and averaging kernels
          call calc_gain_matrix(nx,ny,kk,sa,sai,se,sei,dd,aa)

C         Has solution converged?
          tphi = 100.0*(ophi-phi)/ophi
          if(tphi.ge.0.0.and.tphi.le.phlimit.and.alambda.lt.1.0)then
            if(idiag.gt.0)print*,'%phi, phlimit : ',tphi,phlimit
            if(idiag.gt.0)print*,'Phi has converged'
            if(idiag.gt.0)print*,'Terminating retrieval'
            GOTO 202
          else
            ophi=phi
            oxchi = xchi
            alambda = alambda*0.3		! reduce Marquardt brake
          endif
        else
C	  Leave xn and kk alone and try again with more braking
          alambda = alambda*10.0		! increase Marquardt brake
          if(alambda.gt.1e10)alambda=1e10
        endif

        call file(runname,runname,'abo')
        inquire(file=runname,exist=abexist)
        if(abexist)then
         open(83,file=runname,status='old')
         read(83,'(A)')abort
         close(83)
         if(abort.eq.'stop'.or.abort.eq.'STOP')then
           if(idiag.gt.0)print*,'Terminating retrieval'
           GOTO 202
         endif
        endif
               
401   continue       

202   if(ica.eq.1)close(37)

      if(ica.eq.1)then
C      Write out k-matrix for reference
       OPEN(52,file='kk.out',form='unformatted',status='unknown')
       write(52)y,yn
       write(52)kk
       CLOSE(52)

       close(37)

      endif

      if(idiag.gt.0)print*,'chisq/ny is equal to : ',chisq/float(ny)
      if(chisq.gt.ny)then
       if(idiag.gt.0)then
        print*,'Coreret: WARNING'
        print*,'chisq/ny should be less than 1 if correctly retrieved'
       endif
      endif

      CALL calc_serr(nx,ny,sa,se,aa,dd,st,sn,sm)

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


