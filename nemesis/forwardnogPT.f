      subroutine forwardnogPT(runname,ispace,fwhm,ngeom,nav,
     1 wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,nvar,
     2 varident,varparam,jrad,RADIUS,nx,xn,ny,yn,kk,kiter)
C     $Id:
C     **************************************************************
C     Subroutine to calculate a primary transit spectrum of an exoplanet.
C  
C     Adapted from earlier limb observation Nemesis subroutine forwardL.f
C
C     Calculates the ratio of the planet area to star area and outputs
C     result as a percentage. i.e. Result = 100*planet_area/star_area
C     
C     Input variables:
C       runname(60)   character Name of run.
C       ispace           integer Indicates if wavelengths in vconv and
C                               vwave are in wavenumbers(0) or
C                               wavelengths (1)
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
C       vconv(mgeom,mconv) real Convolution wavelengths
C       angles(mgeom,mav,3) real    Observation angles
C	gasgiant	logical Indicates if planet is a gas giant
C       lin             integer integer to indicate role of previous
C                               retrieval (if any)
C       nvar    integer Number of variable profiles (gas,T,aerosol)
C       varident(nvar,3) integer identity of constituent to retrieved and
C					parameterisation method
C       varparam(nvar,mparam) real Additional arameters constraining profile.
C       RADIUS		real    Radius of planet at 0km tangent altitude
C       nx              integer Number of elements in state vector
C       xn(mx)          real	State vector
C       ny      	integer Number of elements in measured spectra array
C       kiter		integer	Number of iterations of NemesisPT
C
C     Output variables
C       yn(my)          real    Synthetic radiances
C       kk(my,mx)       real    dR/dx matrix
C
C     Pat Irwin	4/4/01		Original
C     Pat Irwin 17/10/03	Tidied for Nemesis
C     Pat Irwin 8/2/04		Modified from forward.f to integrate
C				several observations
C     Pat Irwin */*/11		Modified from forwardl to model primary transits
C     Pat Irwin 25/4/12         Overhauled and tidied.
C
C     **************************************************************

      implicit none
      integer i,j,k,lin,ispace,iav,jpath
      integer ngeom,ioff,igeom,iread
      real interpem,dv
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/gascom.f'
      include 'arraylen.f'
      include '../radtran/includes/planrad.f'
      real xlat,planck_wave,planckg_wave,Bg,height(100),htan
      real wgeom(mgeom,mav),flat(mgeom,mav),fh,thetrot
      integer layint,inormal,iray,iptf,itype,nlayer,laytyp
      integer nwave(mgeom),jsurf,jrad,nem,nav(mgeom),nwave1
      real vwave(mgeom,mwave),angles(mgeom,mav,3),vwave1(mwave)
      real calcout(maxout3),fwhm,output(maxout3)
      real vv,xref,dx
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      integer ipixA,ipixB,ichan,ix
      real vconv(mgeom,mconv),vconv1(mconv)
      real layht,tsurf,esurf,pressR,delp,altbore,thbore,gtsurf
      real xn(mx),yn(my),kk(my,mx),yn1(my),caltbore
      double precision y1(mconv),y2(mconv,mx)
      integer ny,iscat,kiter,nx2,ix1
      integer nphi,ipath
      integer nmu,isol,lowbc,nf
      real dist,galb,sol_ang,emiss_ang,aphi
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      real hcorr,hcorrx

      integer nconvtmp,nwavetmp,iflag,i1,iswitch
      integer nview,iconv,nfov
      real vem(maxsec),emissivity(maxsec),PI,h1,dh,trans
      parameter (PI=3.1415927)

      real vtmp,tmp,delh,dtr,radius,radius1
      double precision area,area0,area1
      integer nvar,varident(mvar,3),ivar
      real varparam(mvar,mparam)
      logical gasgiant

      real stelrad,solwave(maxbin),solrad(maxbin)
      integer solnpt

      common /solardat/iread, stelrad, solwave, solrad,  solnpt
 
C     jradf is passed via the planrad common block
      jradf=jrad

C     Initialise arrays
      do i=1,my
       yn(i)=0.0
       yn1(i)=0.0
       do j=1,mx
        kk(i,j)=0.0
       enddo
      enddo


C     Get all the wavelengths for which calculations are required 
C     at this location, regardless of which geometry.

      nconv1=0
      do igeom=1,ngeom
       if(igeom.eq.1)then
        nconv1=nconv(igeom)
        nwave1=nwave(igeom)
        do i=1,nconv1
         vconv1(i)=vconv(igeom,i)
        enddo
        do i=1,nwave1
         vwave1(i)=vwave(igeom,i)
        enddo
       else
        nconvtmp = nconv(igeom)
        nwavetmp = nwave(igeom)
        do i=1,nconvtmp
         vtmp=vconv(igeom,i)
         iflag=0
         do i1=1,nconv1
          if(vtmp.eq.vconv1(i1))iflag=1
         enddo
         if(iflag.eq.0)then
          nconv1=nconv1+1
          vconv1(nconv1)=vtmp
         endif
        enddo
        do i=1,nwavetmp
         vtmp=vwave(igeom,i)
         iflag=0
         do i1=1,nwave1
          if(vtmp.eq.vwave1(i1))iflag=1
         enddo
         if(iflag.eq.0)then
          nwave1=nwave1+1
          vwave1(nwave1)=vtmp
         endif
        enddo
       endif
      enddo

C     Now sort wavelength arrays
31    continue
      iswitch=0
      do i=1,nconv1-1
       if(vconv1(i).gt.vconv1(i+1))then
        tmp = vconv1(i+1)
        vconv1(i+1)=vconv1(i)
        vconv1(i)=tmp
        iswitch=1
       endif
      enddo      
      if(iswitch.eq.1)goto 31

41    continue
      iswitch=0
      do i=1,nwave1-1
       if(vwave1(i).gt.vwave1(i+1))then
        tmp = vwave1(i+1)
        vwave1(i+1)=vwave1(i)
        vwave1(i)=tmp
        iswitch=1
       endif
      enddo      
      if(iswitch.eq.1)goto 41

      print*,'ForwardnogPT vconv:'
      do i=1,nconv1
       print*,i,vconv1(i)
      enddo

      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf,nphi,layht,tsurf,nlayer,laytyp,layint)

      if(jsurf.gt.0)then
       tsurf = xn(jsurf)
      endif


      iscat=0
      igeom=1
      iav=1

      xlat=flat(igeom,iav)

      if(kiter.ge.0)then
           nx2 = nx+1
      else
           nx2 = 1
      endif


      print*,'Nav = ',nav(igeom)

      nconv1 = nconv(igeom)
      nwave1 = nwave(igeom)

      if(nav(igeom).gt.1)then
         print*,'You should not be doing explicit FOV averaging with'
         print*,'forwardnogPT. All averaging is done implicitly.'
         print*,'Aborting...'
         stop
      endif

      do 111 ix1=1,nx2

         ix=ix1-1

         print*,'forwardnogPT, ix,nx = ',ix,nx

         if(ix.gt.0)then
          xref = xn(ix)
          dx = 0.05*xref
          if(dx.eq.0)dx = 0.1
          if(ix.eq.jrad)dx=10.
          xn(ix)=xn(ix)+dx
         endif


C        If we're retrieving planet radius then add correction to reference
C        radius
         if(jrad.gt.0)then
          radius1 = xn(jrad) + radius
         else
          radius1 = radius
         endif
         print*,'radius=',radius1

C        radius2 is passed via the planrad common block
         radius2=radius1


C        Set up all files for a direct cirsrad run of limb spectra and
C        near-limb observations
         call gsetradPT(runname,nconv1,vconv1,fwhm,ispace,iscat,
     1    gasgiant,layht,nlayer,laytyp,layint,xlat,lin,
     3    nvar,varident,varparam,nx,xn,xmap)

C        If planet is not a gas giant then we need to read in the surface
C         emissivity spectrum
         if(.not.gasgiant)then
          call readsurfem(runname,nem,vem,emissivity)
         else
           nem=2
           vem(1)=-100.0
           vem(2)=1e7
           emissivity(1)=1.0
           emissivity(2)=1.0
         endif

         CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF)


C        Set up parameters for non-scattering cirsrad run.
         itype=12                  ! scloud12. not used here

         call CIRSrtf_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1    vwave1,nwave1,npath,output,vconv1,nconv1, itype, 
     2    nem, vem, emissivity, tsurf, calcout)

         if(ix.eq.0)then
C         Read in base heights from '.drv' file
          call readdrvh(runname,height)        
          area0 = pi*stelrad**2
         endif

         area1 = pi*(radius1+height(1))**2
         print*,area0,area1
   
         do 206 iconv=1,nconv1
 
          area = 0.0

          do 205 ipath=1,npath        

           h1 = radius1 + height(ipath)
           ioff1=ipath+(iconv-1)*npath
           trans = calcout(ioff1)
           y1(ipath) = 2.*pi*h1*(1.-trans)
205       continue        

          do 207 ipath=1,npath-1
           dh = height(ipath+1)-height(ipath)
           area = area + 0.5*(y1(ipath)+y1(ipath+1))*dh
207       continue

        
          yn1(iconv)=sngl(100.0*(area1+area)/area0)
          if(ix.eq.0)then
           yn(iconv)=yn1(iconv)
          else
           kk(iconv,ix)=(yn1(iconv)-yn(iconv))/dx
          endif

206      continue

C       Set a priori element back to original value
        if(ix.gt.0) xn(ix)=xref
        
111   continue

      return

      end
