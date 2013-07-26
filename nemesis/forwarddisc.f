      subroutine forwarddisc(runname,ispace,iscat,fwhm,ngeom,
     1 nav,wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,
     2 lin,nvar,varident,varparam,jsurf,jalb,jtan,jpre,nx,xn,ny,
     3 yn,kk)
C     $Id:
C     **************************************************************
C     Subroutine to calculate an FOV-averaged spectrum and
C     KK-matrix using CIRSRADG gradient code. Uses numerical integration
C     for cases where the limb component is also important. 
C     
C     Input variables:
C       runname(60)   character Name of run.
C       ispace           integer Indicates if wavelengths in vconv and
C                               vwave are in wavenumbers(0) or
C                               wavelengths (1)
C	iscat		integer 0=non-scattering, 1=scattering
C       fwhm            real    Desired FWHM of final spectrum
C       ngeom           integer Number of observation geometries to average
C       nav(mgeom)      integer         Number of synthetic spectra required
C                                       to simulate each FOV-averaged
C                                       measurement spectrum.
C	wgeom(mgeom,mav)real	Integration weights to use
C	flat(mgeom,mav)	real	Integration point latitudes
C       nwave(mgeom) 	integer Number of calculation wavelengths
C       vwave(mgeom,mwave) real Calculation wavelengths
C       nconv(mgeom)    integer Number of convolution wavelengths
C       vconv(mgeom,mconv) real    Convolution wavelengths
C       angles(mgeom,mav,3) real    Observation angles
C	gasgiant	logical Indicates if planet is a gas giant
C	lin		integer integer to indicate role of previous
C                               retrieval (if any)
C       nvar    integer Number of variable profiles (gas,T,aerosol)
C       varident(nvar,3) integer identity of constituent to retrieved and
C					parameterisation method
C       varparam(nvar,mparam) real Additional arameters constraining profile.
C	jsurf		integer	Position of surface temperature element in
C				xn (if included)
C	jalb		integer	Position of surface albedo spectrum in
C				xn (if included)
C	jtan		integer	Position of tangent height correction in
C				xn (if included)
C	jpre		integer	Position of tangent pressure in
C				xn (if included)
C       nx              integer Number of elements in state vector
C       xn(mx)          real	State vector
C       ny      	integer Number of elements in measured spectra array
C
C     Output variables
C       yn(my)          real    Synthetic radiances
C       kk(my,mx)       real    dR/dx matrix
C
C     Pat Irwin	4/4/01		Original
C     Pat Irwin 17/10/03	Tidied for Nemesis
C     Pat Irwin 8/2/04		Modified from forward.f to integrate
C				several observations
C     Pat Irwin 20/9/10		Updated from forwardavfovX for disc-average
C				calculations.
C     Pat Irwin 1/3/12		Updated for Nemesis2.0
C
C     **************************************************************

      implicit none
      integer i,j,lin,ispace,iav,ispace1
      integer ngeom,ioff,igeom
      real interpem
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/gascom.f'
      include 'arraylen.f'
      INCLUDE '../radtran/includes/planrad.f'
      real xlat,planck_wave,planckg_wave,Bg
      real wgeom(mgeom,mav),flat(mgeom,mav),pressR,delp
      integer layint,inormal,iray,itype,nlayer,laytyp,iscat
      integer nwave(mgeom),jsurf,nem,nav(mgeom),nwave1
      integer jalb,jtan,jpre,k,iptf
      real vwave(mgeom,mwave),angles(mgeom,mav,3),vwave1(mwave)
      real pi
      parameter(pi=3.1415927)
      real calcout(maxout3),fwhm
      real gradients(maxout4),vv
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      real vconv(mgeom,mconv),vconv1(mconv),xfac
      real layht,tsurf,esurf,gtsurf
      real xn(mx),yn(my),kk(my,mx),yn1(my)
      integer ny,iconv
      integer nphi,ipath
      integer nmu,isol,lowbc,nf
      real dist,xdist,star,Rstar,galb,sol_ang,emiss_ang,aphi
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,solfile,solname
      real xmap(maxv,maxgas+2+maxcon,maxpro)

      integer nvar,varident(mvar,3)
      real varparam(mvar,mparam),RADIUS1
      logical gasgiant,fexist
      real vem(maxsec),emissivity(maxsec)
      
             
c  ** variables for solar reflected cloud **
      real solar,xsolar
      real refl_cloud_albedo
      logical reflecting_atmos
      common /refl_cloud_params/refl_cloud_albedo,reflecting_atmos

C     Initialise arrays
      do i=1,my
       yn(i)=0.0
       yn1(i)=0.0
       do j=1,mx
        kk(i,j)=0.0
       enddo
      enddo

      print*,'forwarddisc: jsurf,jalb,jtan,jpre',jsurf,jalb,jtan,jpre
      if(jalb.gt.0.or.jtan.gt.0.or.jpre.gt.0)then
       print*,'Warning from forwarddisc'
       print*,'jtan, jalb or jpre > 0'
       stop
      endif

      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf,nphi,layht,tsurf,nlayer,laytyp,layint)

      if(jsurf.gt.0)then
       tsurf = xn(jsurf)
       print*,'tsurf = ',tsurf
      endif

      ioff = 0
      do 100 igeom=1,ngeom
       print*,'Forwarddisc. Spectrum ',igeom,' of ',ngeom

       nconv1 = nconv(igeom)
       nwave1 = nwave(igeom)
 
       do 105 i=1,nconv1
        vconv1(i)=vconv(igeom,i)
105    continue
       do 106 i=1,nwave1
        vwave1(i)=vwave(igeom,i)
106    continue

       if(nwave1.gt.1)call sort(nwave1,vwave1)
       if(nconv1.gt.1)call sort(nconv1,vconv1)

       do 110 iav = 1,nav(igeom)
         sol_ang = angles(igeom,iav,1)
         emiss_ang = angles(igeom,iav,2)
         aphi = angles(igeom,iav,3)

         xlat = flat(igeom,iav)   

C        Set up parameters for non-scattering cirsrad run.
         CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF)
         itype=11			! scloud11. not used here


C        Set up all files for a direct cirsrad run
         call gsetraddisc(runname,iscat,nmu,mu,wtmu,isol,dist,
     1    lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,
     2    layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,lin,
     3    nvar,varident,varparam,nx,xn,jalb,jtan,jpre,tsurf,xmap)

C        If planet is not a gas giant and observation is not at limb then
C        we need to read in the surface emissivity spectrum
         if(.not.gasgiant.and.emiss_ang.ge.0)then
          call readsurfem(runname,nem,vem,emissivity)
         else
           nem=2
           vem(1)=-100.0
           vem(2)=1e7
           emissivity(1)=1.0
           emissivity(2)=1.0
         endif


         call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1    vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gtsurf,
     2    nx, xmap, vconv1, nconv1, npath, calcout, gradients)


         ipath = 1
         xfac = 1.0
         do j=1,nconv1
          iconv=-1
          do k = 1,nconv1
           if(vconv(igeom,j).eq.vconv1(k))iconv=k
          enddo
          if(iconv.lt.0)then
           print*,'Error in forwarddisc iconv<1'
           stop
          endif

	  ioff1=nconv1*(ipath-1)+iconv
          yn(ioff+j)=yn(ioff+j)+wgeom(igeom,iav)*calcout(ioff1)
         enddo
    
         do i=1,nx
           do j=1,nconv1
            iconv=-1
            do k = 1,nconv1
             if(vconv(igeom,j).eq.vconv1(k))iconv=k
            enddo
            ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
            if(i.ne.jsurf)then 
             kk(ioff+j,i)=kk(ioff+j,i)+wgeom(igeom,iav)*
     1		gradients(ioff2)
            else
             kk(ioff+j,i)=kk(ioff+j,i)+wgeom(igeom,iav)*
     1		gtsurf
            endif
           enddo

         enddo

         if (reflecting_atmos) then
          print*,'ADDING IN REFLECTING ATMOSPHERE CONTRIBUTION'
          print*,'cloud albedo=',refl_cloud_albedo

          ipath = 2

c         initialise solar spectrum
          CALL FILE(runname,solfile,'sol')
          inquire(file=solfile,exist=fexist)
          if(fexist) then
           call opensol(solfile,solname)
          else
           print*,'Error in forwarddisc. solar flux file not defined'
           print*,'for reflectivity calculation'
           stop
          endif

          CALL init_solar_wave(ispace,solname)

          do j=1,nconv1
           vv = vconv(igeom,j)
           iconv=-1
           do k = 1,nconv1
            if(vv.eq.vconv1(k))iconv=k
           enddo

             ioff1=nconv1*(ipath-1)+iconv
C            Get star flux at planet's radius
             CALL get_solar_wave(vconv1(j),dist,solar)

C            Get overall star spectral power
             xdist=-1.
             CALL get_solar_wave(vconv1(j),xdist,xsolar)

           Bg = solar*refl_cloud_albedo/pi
           
           radius1=radius2

           print*,'forwarddisc: radius = ',radius1
           xfac=(2.*pi)*4.*pi*RADIUS1**2/xsolar

           yn(ioff+j)=yn(ioff+j) + xfac*wgeom(igeom,iav)*
     1					calcout(ioff1)*Bg

           do i=1,nx
            ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
            kk(ioff+j,i)=kk(ioff+j,i) +
     1         xfac*wgeom(igeom,iav)*gradients(ioff2)*Bg
           enddo
     
          enddo

	 endif

110    continue

       ioff = ioff + nconv1

100   continue

      return

      end
