      subroutine coreretMCMC(runname,ispace,iscat,ilbl,ica,miter,
     1  niter,fwhm,xlat,ngeom,nav,nwave,vwave,nconv,vconv,
     2  angles,gasgiant,nvar,varident,varparam,npro,jsurf,jalb,jtan,
     3  jpre,jrad,jlogg,wgeom,flat,nx,lx,xa,sa,ny,y,se1,idum)
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
C	ica	integer	1 if single retrieval, 0 otherwise.
C	miter	integer	Total number of report steps
C	niter	integer Number of iterations between reporting steps
C	fwhm	real	Required FWHM of final convoluted spectrum
C	xlat	real	Latitude of observed site.
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
C	nvar	integer	Number of variable profiles (gas,T,aerosol)
C	varident(mvar,3) integer Identity of constituent to retrieved and
C				 	method of parameterisation
C	varparam(mvar,mparam) real Additional parameters constraining
C					profile.
C	jsurf		integer	Position of surface temperature element in
C				xa (if included)
C	jalb		integer	Position of surface albedo spectrum in
C				xa (if included)
C	jtan		integer	Position of tangent height correction in
C				xa (if included)
C	jpre		integer	Position of tangent pressure in
C				xa (if included)
C	wgeom(mgeom,mav) real	Integration weights 
C	flat(mgeom,mav)	real	Integration point latitudes 
C	nx		integer	Number of elements in measurement vector
C       lx(mx)          integer 1 if log, 0 otherwise
C	xa(mx)		real	a priori state vector
C	sa(mx,mx)	real 	A priori covariance matrix
C	ny	integer	Number of elements in measured spectra array
C	y(my)	real	Measurement vector
C	se1(my)	real	Measured radiance variances
C	idum	integer	negative seed number of random number generator
C
C     Pat Irwin		29/4/01
C			10/10/03 conversion for Nemesis
c			10/6/13  converted from coreret.f
C
C     ************************ VARIABLES *******************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      integer iter,kiter,ica,iscat,i,j,icheck,j1,j2,jsurf,lin
      integer jalb,jtan,jpre,ilbl,jrad,jlogg,miter,niter,idum,itry
      CHARACTER*100 runname,itname,abort

      real xn(mx),se1(my),calc_chi,kk(my,mx),calc_phiprior
      real fwhm,xlat,xn1(mx),ran11,test,RADIUS
      integer ix,np,npro,lout,iplanet,lx(mx)
      real alpha,alpha_chi,alpha_phi
      logical accept
      integer nvar,varident(mvar,3),ispace,nav(mgeom),k
      real varparam(mvar,mparam)

      integer ngeom, nwave(mgeom), nconv(mgeom), nx, ny, nxf,ivar,ivarx
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real xa(mx),sa(mx,mx),y(my),yn(my)
      real yn1(my),sx(mx,mx)
      real wgeom(mgeom,mav),flat(mgeom,mav)
      real vwaveT(mwave),vconvT(mconv)
      integer nwaveT,nconvT,npvar,iprfcheck
      logical gasgiant
      double precision s1d(mx,mx),sai(mx,mx)

      real chisq,xchi,ochisq,oxchi,phi,ophi,xphi
C     **************************** CODE ********************************

      open(41,file='testxtry.dat',status='unknown')

      if(ilbl.eq.0)then
C      Find all the calculation and convolution wavelengths and rank
C      in order

       call rankwave(ngeom,nwave,vwave,nconv,vconv,nwaveT,vwaveT,
     1  nconvT,vconvT)

      endif

      CALL readrefiplan(runname,iplanet,xlat,RADIUS)

C     Calculate first spectrum

C     Load state vector with a priori
      do i=1,nx
       xn(i)=xa(i)
      enddo


C     Invert covariance matrix
      do i=1,nx
       do j=1,nx
        s1d(i,j)=dble(sa(i,j))
       enddo
      enddo

C     Calculate inverse of sa.   
      jmod = 2
      icheck=0
      call dinvertm(jmod,icheck,s1d,nx,mx,sai)
      if(icheck.eq.1)then
       print*,'************* WARNING *************'
       print*,'CoreretMCMC, sa does not invert cleanly'
       print*,'Aborting...'
       print*,'***********************************'
       stop
      endif

C     Set kiter to -1 to stop the code calculating any gradients
      kiter = -1
C     Set lin=0 to prevent code looking for previous retrievals
      lin=0

      if(ilbl.eq.1)then
         print*,'Calling forwardnoglbl - B'
         CALL forwardnoglbl(runname,ispace,iscat,fwhm,ngeom,nav,
     1    wgeom,flat,nconv,vconv,angles,gasgiant,lin,
     2    nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,
     3    jlogg,RADIUS,nx,xn,ny,yn,kk,kiter)
      else

       if(iscat.ne.2)then
 
        CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,jlogg,
     3   RADIUS,nx,xn,ny,yn,kk,kiter,iprfcheck)
       else

        CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1   vconvT,gasgiant,lin,nvar,varident,varparam,jsurf,jalb,
     2   jtan,jpre,jrad,jlogg,RADIUS,nx,xn)

        CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,jlogg,
     3   RADIUS,nx,xn,ny,yn,kk,kiter,iprfcheck)

       endif
      endif

C     Calculate initial value of cost function chisq.
      open(12,file='test.dat',status='unknown')
       write(12,*)nx,ny
       write(12,*)(xn(i),i=1,nx)
       do i=1,ny
        write(12,*)y(i),se1(i),yn(i)
       enddo
      close(12)

      chisq = calc_chi(ny,y,yn,se1)
      ochisq = chisq
      oxchi = chisq/float(ny)
      
      print*,'Initial chisq/ny = ',oxchi
      phi = calc_phiprior(nx,xn,xa,sai)
      ophi = phi

      lout=38

C     Set the trial vectors xn1, and yn1 to be the same as the initial
C     vectors xn, yn
      do i=1,nx
       xn1(i)=xn(i)
      enddo
      do i=1,ny
       yn1(i)=yn(i)
      enddo
C     Set proposal covariance of proposal distribution to be proportional
C     to the apriori distribution
      do i=1,nx
       do j=1,nx
        sx(j,i)=sa(j,i)
       enddo
      enddo


C     Set proposal distribution to be same as a priori distribution (for now)

      do 401 iter = 1, miter

       do 402 itry = 1,niter

        print*,'iter,itry',iter,itry

C       Now calculate next iterated xn
        call modxvecMCMCA(idum,npro,nvar,varident,varparam,
     1   nx,xn,sx,xn1)


C       Calculate test spectrum using trial state vector xn1. 
C       Put output spectrum into temporary spectrum yn1.

        if(ilbl.eq.1)then
         CALL forwardnoglbl(runname,ispace,iscat,fwhm,ngeom,nav,     
     1    wgeom,flat,nconv,vconv,angles,gasgiant,lin,
     2    nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,jlogg,
     3    RADIUS,nx,xn1,ny,yn1,kk,kiter)
        else
         if(iscat.ne.2)then
          CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1     wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2     nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,jlogg,
     3     RADIUS,nx,xn1,ny,yn1,kk,kiter,iprfcheck)
         else
          CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1     vconvT,gasgiant,lin,nvar,varident,varparam,jsurf,jalb,
     2     jtan,jpre,jrad,jlogg,RADIUS,nx,xn1)
          CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1     wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2     nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,jlogg,
     3     RADIUS,nx,xn1,ny,yn1,kk,kiter,iprfcheck)
         endif
        endif

C       Calculate the cost function for this trial solution.
        chisq = calc_chi(ny,y,yn1,se1)

C        open(12,file='test1.dat',status='unknown')
C         write(12,*)nx,ny
C         write(12,*)(xn1(i),i=1,nx)
C         do i=1,ny
C          write(12,*)y(i),se1(i),yn1(i)
C         enddo
C        close(12)

        phi = calc_phiprior(nx,xn,xa,sai)

        xchi = chisq/float(ny)
        print*,'chisq/ny = ',xchi
        xphi = phi/float(nx)
        print*,'phi,phi/nx = ',xphi

        print*,'XXX',nx
        do i=1,nx
         print*,xa(i),xn(i),xn(i)-xa(i),sqrt(sa(i,i)),
     1  (xn(i)-xa(i))/sqrt(sa(i,i))
        enddo

C       Does trial solution fit the data better?

C       Calculate probability of acceptance
        alpha_chi = exp(-0.5*(chisq-ochisq))
        print*,'AAA',chisq,ochisq,alpha_chi
        alpha_phi = exp(-0.5*(phi-ophi))
        print*,'BBB',phi,ophi,alpha_phi
C       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C       Benneke and Seager include the probability of the solution as
C       determined from apriori information. I'm not yet sure how to do
C       this and my first attempt calc_phiprior seems to screw things up
C       so here I force the code not to worry about phi!
C        alpha_phi = exp(-0.5*(phi-ophi)/float(nx))
C        print*,'CCC',phi,ophi,alpha_phi
        alpha_phi = 1.
        ophi=1.
        phi=1.
C       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        alpha = alpha_chi*alpha_phi
        print*,'Alpha',alpha
        accept=.false.
C       If fit is better and the solution is closer to the priori then accept
        if(chisq.le.ochisq.and.phi.le.ophi)then
         print*,'XX'
         accept=.true.
        else
C        If not automatically accepted then choose depending on combined
C        probability
         test=ran11(idum)
         print*,'YY',alpha,test
         if(alpha.ge.test)then
          accept=.true.
          print*,'YYX'
         endif
        endif

        print*,'Accept = ',accept

        if(accept)then
         do i=1,nx
          xn(i)=xn1(i)         		! update xn to new value
         enddo
         do i=1,ny
          yn(i)=yn1(i)				! update yn
         enddo
         ochisq=chisq
         oxchi=xchi
         ophi=phi
        endif

        if(ica.eq.1)then
         write(39,*)(xn(i),i=1,nx)
         write(39,*)(yn(i),i=1,ny)
         write(39,*)chisq
        endif

402    continue

       write(lout,*)iter
       write(lout,*)(xn(i),i=1,nx)
       write(lout,*)(yn(i),i=1,ny)
       write(lout,*)chisq

401   continue       

      if(chisq.ge.ochisq)then
       xchi=oxchi
      endif
      print*,'chisq/ny is equal to : ',xchi

      if(xchi.gt.1.)then
       print*,'Coreret: WARNING'
       print*,'chisq/ny should be less than 1 if correctly retrieved'
      endif

      close(41)

      return

      end


