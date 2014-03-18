      subroutine coreretMC(runname,ispace,iscat,ica,kiter,phlimit,
     1  fwhm,xlat,ngeom,nav,nwave,vwave,nconv,vconv,angles,
     2  gasgiant,lin,lpre,nvar,varident,varparam,npro,jsurf,jalb,jtan,
     3  jpre,jrad,wgeom,flat,nx,lx,xa,sa,ny,y,se1,xn,sm,sn,st,yn,kk,
     4  aa,dd)
C     $Id:
C     ******************************************************************
C
C     Input variables
C	runname	character*100	Root name of associated run files
C       ispace           integer Indicates if wavelengths in vconv and 
C                               vwave are in wavenumbers(0) or 
C                               wavelengths (1)
C	iscat	integer	Set to 1 to use scattering RTM, 0 for thermal emission
C			Set iscat=2 if internal radiation field to be
C			calculated first
C	ica	integer	1 if single retrieval, 0 otherwise.
C	kiter	integer	Maximum number of iterations
C	phlimit	real	Limiting % change in cost function to consider solution
C			converged.
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
      integer jalb,jalbx,jtan,jtanx,jpre,jprex,iscat1
      integer jrad,jradx,npvar,iplanet,lx(mx)
      real phlimit,alambda,xtry,tphi, RADIUS
      CHARACTER*100 runname,itname,abort

      real xn(mx),se1(my),se(my,my),calc_phiret,sf(my,my)
      real fwhm,xlat,xlatx,xdiff,xn1(mx),x_out(mx)
      real xlonx
      integer nprox,nvarx,varidentx(mvar,3),jsurfx,nxx,ix,np,npro
      real st(mx,mx),varparamx(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xnx(mx),stx(mx,mx),ynx(my)

      integer nvar,varident(mvar,3),lin,lin0,lpre,ispace,nav(mgeom),k
      real varparam(mvar,mparam)

      integer ngeom, nwave(mgeom), nconv(mgeom), nx, ny, nxf,ivar,ivarx
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real xa(mx),kk1(my,mx),sa(mx,mx),y(my),yn(my),kkx(my,mx)
      real yn1(my),s1(mx,mx),kk(my,mx)
      real wgeom(mgeom,mav),flat(mgeom,mav)
      real vwaveT(mwave),vconvT(mconv)
      integer nwaveT,nconvT
      logical gasgiant

      double precision s1d(mx,mx),sai(mx,mx)
      double precision s1e(my,my),sei(my,my)
      double precision dd(mx,my),aa(mx,mx)

      real phi,ophi,chisq,xchi,oxchi
C     **************************** CODE ********************************

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
C       print*,(sa(i,j),j=1,nx)
      enddo
      
C     Calculate inverse of sa.
      jmod = 2
      icheck=0
      call dinvertm(jmod,icheck,s1d,nx,mx,sai)
      if(icheck.eq.1)then
       print*,'************* WARNING *************'
       print*,'CoreretMC, sa does not invert cleanly'
C       print*,'Setting to simple diagonal matrix'
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

      CALL readrefiplan(runname,iplanet,RADIUS)


C     Calculate first spectrum and k-matrix

C     Load state vector with a priori
      do i=1,nx
       xn(i)=xa(i)
      enddo

C      print*,'coreretMC: lin = ',lin
      if(lin.eq.1.or.lin.eq.3)then

       if(lin.eq.1) then
C        Just substituting parameters from .pre file
         call readraw(lpre,xlatx,xlonx,nprox,nvarx,varidentx,varparamx,
     1  jsurfx,jalbx,jtanx,jprex,jradx,nxx,xnx,stx)
       
        xdiff = abs(xlat-xlatx)
        if(xdiff.gt.lat_tolerance)then
          print*,'CoreretMC: Aborting - latitudes inconsistent'
          print*,xlatx,xlat
          stop
        endif
         
        do ivarx=1,nvarx
         do ivar=1,nvar
          if(varidentx(ivarx,1).eq.varident(ivar,1))then
           if(varidentx(ivarx,2).eq.varident(ivar,2))then
             print*,'CoreretMC: Can not use previous retrieval to add'
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
        call writextmp(runname,xlatx,nvarx,varidentx,varparamx,nprox,
     1   nxx,xnx,stx,jsurfx,jalbx,jtanx,jprex,jradx)

       else
C       substituting and retrieving parameters from .pre file. 
C       Current record frrom .pre file already read in by
C       readapriori.f. Hence just read in from temporary .str file
 
        call readxtmp(runname,xlatx,nvarx,varidentx,varparamx,nprox,
     1   nxx,xnx,stx,jsurfx,jalbx,jtanx,jprex,jradx)

       endif
 
       lin0 = 0

       if(iscat.eq.0)then
        print*,'Calling forwardavfovX'
        CALL forwardavfovX(runname,ispace,iscat,fwhm,ngeom,
     1   nav,wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,
     2   lin0,nvarx,varidentx,varparamx,jsurfx,jalbx,jtanx,jprex,
     3   jradx,RADIUS,nxx,xnx,ny,ynx,kkx)
       elseif(iscat.eq.1)then
        print*,'Calling forwardnogMC'
        CALL forwardnogMC(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin0,
     2   nvarx,varidentx,varparamx,jsurfx,jalbx,jtanx,jprex,
     3   nxx,xnx,ny,ynx,kkx,kiter)
       elseif(iscat.eq.2)then
        print*,'Calling intradfield'
        CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1   vconvT,gasgiant,lin0,nvarx,varidentx,varparamx,
     2   jsurfx,jalbx,jtanx,jprex,nxx,xnx)

        iscat1=1
        CALL forwardnogX(runname,ispace,iscat1,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin0,
     2   nvarx,varidentx,varparamx,jsurfx,jalbx,jtanx,jprex,jradx,
     3   RADIUS,nxx,xnx,ny,ynx,kkx,kiter)
       else
        print*,'CoreretMC: iscat invalid',iscat
        stop
       endif

       if(lin.eq.3) then
C        strip out variables from kkx and stx that will be retrieved in 
C        this run.
         call scankkx(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  kkx,stx,nxx)
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
        print*,'CoreretMC, se does not invert cleanly'
        print*,'Aborting...'
        print*,'***********************************'
        stop
       endif

      endif

C      print*,'OK, iscat = ',iscat

      if(iscat.eq.0)then
        print*,'Calling forwardavfovX - A'
C        print*,runname,ispace,iscat,fwhm
C        print*,ngeom
C        do i=1,ngeom
C         print*,nav(i)
C         print*,(wgeom(i,j),j=1,nav(i))
C         print*,(flat(i,j),j=1,nav(i))
C         do j=1,nav(i)
C          print*,(angles(i,j,k),k=1,3)
C         enddo
C         print*,nwave(i),nconv(i)
C        enddo
C        print*,gasgiant,lin
C        print*,nvar
C        do i=1,nvar
C         print*,(varident(i,j),j=1,3)
C         print*,(varparam(i,j),j=1,5)
C        enddo
C        print*,'jsurf,jalb,jtan,jpre',jsurf,jalb,jtan,jpre
C        print*,nx,(xn(i),i=1,nx)
C        print*,ny

        CALL forwardavfovX(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,RADIUS,
     3   nx,xn,ny,yn,kk)

C        print*,'forwardavfovX OK, jpre = ',jpre

      elseif(iscat.eq.1)then
 
       CALL forwardnogMC(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,nx,xn,ny,
     3   yn,kk,kiter)

C        print*,'forwardnogMC OK, jpre = ',jpre


      elseif(iscat.eq.2)then

       print*,'testing'
       print*,'nx = ',nx
       CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1   vconvT,gasgiant,lin,nvar,varident,varparam,jsurf,jalb,
     2   jtan,jpre,nx,xn)

       iscat1=1
C       print*,'Calling forwardnogX'
       CALL forwardnogX(runname,ispace,iscat1,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,
     3   RADIUS,nx,xn,ny,yn,kk,kiter)

      else
       iscat1=1
       CALL forwardnogX(runname,ispace,iscat1,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,RADIUS,
     3   nx,xn,ny,yn,kk,kiter)

C        print*,'forwardnogX OK, jpre = ',jpre
      endif

      open(12,file='kk.dat',status='unknown')
      write(12,*)nx,ny
      do i=1,ny
       write(12,*)(kk(i,j),j=1,nx)
      enddo
      close(12)


C      print*,'Calling calc_gain_matrix'
C     Now calculate the gain matrix and averaging kernels
      call calc_gain_matrix(nx,ny,kk,sa,sai,se,sei,dd,aa)

C     Calculate initial value of cost function phi.
      phi = calc_phiret(ny,y,yn,sei,nx,xn,xa,sai,chisq)
      ophi = phi
      oxchi = chisq/float(ny)

C      print*,'Calling assess'
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

C        print*,'Calling calcnextxn'
C       Now calculate next iterated xn1
        call calcnextxn(nx,ny,xa,xn,y,yn,dd,aa,x_out)

C       x_out(nx) is the next iterated value of xn using classical N-L
C       optimal estimation. However, we want to apply a braking parameter
C       alambda to stop the new trial vector xn1 being too far from the
C       last 'best-fit' value xn
C        print*,'Next xn'

145     continue
        do i=1,nx
         xn1(i) = xn(i) + (x_out(i)-xn(i))/(1.0+alambda)
C         print*,i,xn(i),xn1(i)
C         if(xn1(i).lt.0.0)then
C          print*,'*********** Warning ************'
C          print*,'Possible error in coreretMC. Trial vector element'
C          print*,'has gone negative.'
C          print*,i,xa(i),xn1(i),x_out(i)
C         endif

C        Check to see if log numbers have gone out of range
         if(lx(i).eq.1)then
          if(xn1(i).gt.85.or.xn1(i).lt.-85)then
           print*,'Coreret - log(number gone out of range)'
           print*,'Increasing brake'
           alambda = alambda*10.0               ! increase Marquardt brake
           if(alambda.gt.1e10)alambda=1e10
           goto 145
          endif
         endif

        enddo

        ix=1
        do ivar = 1,nvar
         np=1
         if(varident(ivar,1).le.100)then
           np=npvar(varident(ivar,3),npro)
         endif
         if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))
         if(varident(ivar,1).eq.444)np = 2+int(varparam(ivar,1))

         do j=ix,ix+np-1
          if(varident(ivar,1).eq.0)then
           if(varident(ivar,3).eq.0) then
            if(xn1(j).lt.1.0) then
             print*,'Temperature has gone negative, Increase alambda'
             alambda = alambda*10.0             ! increase Marquardt brake
             if(alambda.gt.1e10)alambda=1e10
             goto 145
            endif
           endif
           if(varident(ivar,3).eq.16.and.j.eq.ix) then
            if(xn1(j).lt.1.0) then
             print*,'Temperature has gone negative, Increase alambda'
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
          CALL forwardavfovX(runname,ispace,iscat,fwhm,ngeom,nav,
     1     wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,
     2     lin,nvar,varident,varparam,jsurf,jalb,jtan,jpre,
     3     jrad,RADIUS,nx,xn1,ny,yn1,kk1)
        elseif(iscat.eq.1)then
          CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1     wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2     nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,RADIUS,
     3     nx,xn1,ny,yn1,kk1,kiter)
        elseif(iscat.eq.2)then
          CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1     vconvT,gasgiant,lin,nvar,varident,varparam,jsurf,jalb,
     2     jtan,jpre,nx,xn1)
          iscat1=1
          CALL forwardnogX(runname,ispace,iscat1,fwhm,ngeom,nav,
     1     wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2     nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,RADIUS,
     3     nx,xn1,ny,yn1,kk1,kiter)
        else
          iscat1=1
          CALL forwardnogX(runname,ispace,iscat1,fwhm,ngeom,nav,
     1     wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2     nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,RADIUS,
     3     nx,xn1,ny,yn1,kk1,kiter)
        endif

C       Calculate the cost function for this trial solution.
        phi = calc_phiret(ny,y,yn1,sei,nx,xn1,xa,sai,chisq)

        xchi = chisq/float(ny)
        print*,'chisq/ny = ',xchi
        print*,'it.,al.,ophi.,phi.',
     1   iter,alambda,ophi,phi

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

C          print*,'Calculating new gain matrix and averaging kernels'
C         Now calculate the gain matrix and averaging kernels
          call calc_gain_matrix(nx,ny,kk,sa,sai,se,sei,dd,aa)

C          print*,'calc_gain_matrix OK'

C         Has solution converged?
          tphi = 100.0*(ophi-phi)/ophi
          if(tphi.ge.0.0.and.tphi.le.phlimit.and.alambda.lt.1.0)then
            print*,'%phi, phlimit : ',tphi,phlimit
            print*,'Phi has converged'
            print*,'Terminating retrieval'
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
        open(83,file=runname,status='old')
        read(83,'(A)')abort
        close(83)
        if(abort.eq.'stop'.or.abort.eq.'STOP')then
          print*,'Terminating retrieval'
          GOTO 202
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

      print*,'chisq/ny is equal to : ',chisq/float(ny)
      if(chisq.gt.ny)then
       print*,'CoreretMC: WARNING'
       print*,'chisq/ny should be less than 1 if correctly retrieved'
      endif

C      print*,'Calculating final covariance matrix'
      CALL calc_serr(nx,ny,sa,se,aa,dd,st,sn,sm)
C      print*,'Matrix calculated'

      return

      end


