      subroutine calc_input_files(runname,ispace,iscat,fwhm,flat,
     1 flon,nconv,vconv,angles,gasgiant,lin,nvar,varident,varparam,
     2 jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,RADIUS,nx,xn)
C     $Id:
C     **************************************************************
C     Subroutine to write out input files for last successful iteration of
C     xn.
C
C     Input variables:
C       runname(60)   character Name of run.
C       ispace           integer Indicates if wavelengths in vconv and
C                               vwave are in wavenumbers(0) or
C                               wavelengths (1)
C       fwhm            real    Desired FWHM of final spectrum
C       flat(mgeom,mav)  real    Integration point latitudes
C       flon(mgeom,mav)  real    Integration point longitudes
C       nconv(mgeom)    integer Number of convolution wavelengths
C       vconv(mgeom,mconv) real    Convolution wavelengths
C       angles(mgeom,mav,3) real    Observation angles
C       gasgiant        logical Indicates if planet is a gas giant
C       lin             integer indicates role of previous retrieval (if any)
C       nvar    integer Number of variable profiles (gas,T,aerosol)
C       varident(nvar,3) integer identity of constituent to retrieved and
C					parameterisation method
C       varparam(nvar,mparam) real Additional arameters constraining profile.
C       jsurf           integer Position of surface temperature element in
C                               xn (if included)
C	jalb		integer position of first surface albedo element in
C				xn (if included)
C	jxsc		integer position of first x-section element in
C				xn (if included)
C	jtan		integer position of tangent ht. correction element in
C				xn (if included)
C	jpre		integer position of tangent pressure element in
C				xn (if included)
C       jrad		integer position radius element in
C                               xn (if included)
C       jlogg		integer position surface gravity (log(g)) element in
C                               xn (if included)
C       RADIUS		real    Planetary radius at 0km altitude
C       nx              integer Number of elements in state vector
C       xn(mx)          real	State vector
C
C     Output variables
C       None.
C
C     Pat Irwin	22/9/16		Original
C
C     **************************************************************

      implicit none
      integer i,j
      integer igeom,lin
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/gascom.f'
      include '../radtran/includes/planrad.f'
      include 'arraylen.f'
      real xlat,Grav,fwhm,xlon
      parameter (Grav=6.672E-11)
      integer layint,inormal,iray,itype,nlayer,laytyp,iscat
      integer iptf,jrad,j1,iav,ispace
      real RADIUS
      integer imie,imie1,jlogg
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      real vconv(mgeom,mconv),flat(mgeom,mav)
      real layht,tsurf,esurf,angles(mgeom,mav,3)
      real xn(mx),flon(mgeom,mav)
      real vconv1(mconv)
      integer ny,jsurf,jalb,jxsc,jtan,jpre,nem
      integer nphi,ipath,iconv,k
      integer nmu,isol,lowbc,nf,nf1,nx2,kiter
      real dist,galb,sol_ang,emiss_ang,z_ang,aphi,vv
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,logname
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      common /imiescat/imie1

      integer nvar,varident(mvar,3)
      real varparam(mvar,mparam)
      logical gasgiant
C     jradf and jloggf are passed via the planrad common block   
      jradf=jrad
      jloggf=jlogg


      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf1,nphi,layht,tsurf,nlayer,laytyp,layint)

      if(jsurf.gt.0)then
       tsurf = xn(jsurf)
      endif

      igeom=1
      iav=1
      sol_ang = angles(igeom,iav,1)
      emiss_ang = angles(igeom,iav,2)
      aphi = angles(igeom,iav,3)
      nconv1=nconv(igeom)   
      do 105 i=1,nconv1
        vconv1(i)=vconv(igeom,i)
105   continue

      if(sol_ang.lt.emiss_ang) then
           z_ang = sol_ang
        else
           z_ang = emiss_ang
      endif

C     New bit to increase number of Fourier components depending on
C     miniumum zenith angle
         if(iscat.eq.1.and.z_ang.ge.0.0)then
           nf = int(30*z_ang/90.0)
      else
           nf=nf1
      endif

      print*,'Angles : ',sol_ang,emiss_ang,aphi
      print*,'nf = ',nf
      xlat = flat(igeom,iav)
      xlon = flon(igeom,iav)


C     If we're retrieving planet radius then add correction to reference
C     radius
C     N.B.radius2 is passed via the planrad common block.
      if(jrad.gt.0)then
          radius2 = xn(jrad) + radius
      else      
          radius2 = radius
      endif

C     If we're retrieving surface gravity then modify the planet mass
C     N.B. mass2 is passed via the planrad common block. Assume xn(jlogg)
C     holds log_10(surface gravity in cm/s^2). Need factor of 1e-20 to convert
C     mass to units of 1e24 kg.       
      if(jlogg.gt.0)then
          mass2 = 1e-20*10**(xn(jlogg))*(radius2**2)/Grav
      endif

         
      CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,INH3,
     1    IPTF,IMIE, iuvscat)
      IMIE1=IMIE
      itype=11			! scloud11wave


C     Set up all files for a direct cirsrad run
      print*,'calling gsetrad'
      call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,
     1 lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,
     2 layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,xlon,
     2 lin,nvar,varident,varparam,nx,xn,jalb,jxsc,jtan,jpre,tsurf,xmap)
      print*,'gsetrad called OK'

      return

      end
