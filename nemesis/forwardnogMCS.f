      subroutine forwardnogMCS(runname,ispace,iscat,fwhm,xlat,ngeom,nav,
     1 wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,
     2 lin,nvar,varident,varparam,jsurf,jalb,jtan,jpre,RADIUS,
     3 SATRAD,thetrot,altbore,nx,xn,ny,yn,kk,kiter)
C     $Id:
C     **************************************************************
C     Subroutine to calculate a synthetic spectrum and KK-matrix using
C     FINITE DIFFERENCES. The routine is identical in operation to 
C     forwardavfovX.f but calculates the K-matrix using old-fashioned 
C     finite differences. Routine exists to check that forwardavfovX.f, 
C     using the internal gradients is operating correctly, and also for
C     scattering calculations.
C
C     Input variables:
C       runname(60)   character Name of run.
C       ispace           integer Indicates if wavelengths in vconv and
C                               vwave are in wavenumbers(0) or
C                               wavelengths (1)
C	iscat		integer 0=non-scattering, 1=scattering
C       fwhm            real    Desired FWHM of final spectrum
C       ngeom           integer Number of observation geometries included.
C       nav(mgeom)      integer         Number of synthetic spectra required
C                                       to simulate each FOV-averaged
C                                       measurement spectrum.
C       wgeom(mgeom,mav)real     Integration weights to use
C       flat(mgeom,mav)  real    Integration point latitudes
C       nwave(mgeom) integer Number of calculation wavelengths
C       vwave(mgeom,mwave) real    Calculation wavelengths
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
C	jtan		integer position of tangent ht. correction element in
C				xn (if included)
C	jpre		integer position of tangent pressure element in
C				xn (if included)
C       RADIUS          real    Radius of planet at 0km tangent altitude
C       SATRAD          real    Distance of satellite from centre
C                               of planet
C       thetrot         real    Angle of rotation of pixel array from 
C                               vertical
C	altbore		real	Altitude of tangent boresight above 
C				surface
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
C     Pat Irwin 28/10/03	Modified from forward.f for testing 
C				   purposes
C
C     **************************************************************

      implicit none
      integer i,j,ispace,ulog
      parameter (ulog=17)
      integer ngeom,ioff,igeom,lin
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/gascom.f'
      include 'arraylen.f'
      real xlat,xref,dx
      integer layint,inormal,iray,itype,nlayer,laytyp,iscat,iscat1
      integer nwave(mgeom),ix,ix1,iav,nwave1,lin0,nlayer1,fmod
      real vwave(mgeom,mwave),interpem,height(100),vv
      real calcout(maxout3),fwhm,planck_wave,output(maxout3)
      real calcouts(maxout3)
      integer iptf
      real htan,thetrot,altbore,caltbore
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      integer ipixA,ipixB,ichan
      real vconv(mgeom,mconv),wgeom(mgeom,mav),flat(mgeom,mav)
      real layht,tsurf,esurf,angles(mgeom,mav,3)
      real xn(mx),yn(my),kk(my,mx),ytmp(my),yn1(my)
      real vconv1(mconv),vwave1(mwave)
      integer ny,jsurf,jalb,jtan,jpre,nem,nav(mgeom)
      integer nphi,ipath
      integer nmu,isol,lowbc,nf,nf1,nx2,kiter
      real dist,galb,sol_ang,emiss_ang,z_ang,aphi
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,logname,intname
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      real hcorr,hcorrx

      integer nconvtmp,nwavetmp,iflag,i1,iswitch
      integer nview,mview,iconv,nfov,mfov,iread
      parameter (mfov=200,mview=100)

      real vtmp,tmp,topht,delh,dtr,thcentre,RADIUS
      real SATRAD,hview(mview),thview(mview),thbore
      real thfov(mfov),rfov(mfov),radmean

      integer nvar,varident(mvar,3),ivar
      real varparam(mvar,mparam)
      logical gasgiant
      real vem(maxsec),emissivity(maxsec)


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

      print*,'ForwardnogMCS vconv:'
      do i=1,nconv1
       print*,i,vconv1(i)
      enddo

      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf1,nphi,layht,tsurf,nlayer,laytyp,layint)

      if(jsurf.gt.0)then
       tsurf = xn(jsurf)
      endif


C     Read in base heights from '.drv' file
      call readdrvh(runname,height)

      topht=100.0

C     Define tangent height array
      nview = nlayer + 1 + int(0.5*nlayer)
      do i=1,nlayer
        hview(i)=height(i)
      enddo
      hview(nlayer+1)=layht-0.001
      delh = (topht-layht)/nlayer
      do i=1,int(0.5*nlayer)
        hview(nlayer+1+i)=layht-(0.001+i*delh)       
      enddo
C     Now calculate viewing angle w.r.t spacecraft
      print*,'Pre-calc grid'
      print*,'Radius, satrad, thetrot = ',radius,satrad,thetrot
      do i=1,nview
        thview(i) = asin((RADIUS+hview(i))/SATRAD)
        print*,i,hview(i),thview(i)
      enddo


C     If planet is not a gas giant then we need to read in the surface 
C       emissivity spectrum
      if(.not.gasgiant)then
          call readsurfem(runname,nem,vem,emissivity)
      else
           nem=2
           vem(1)=-100.0
           vem(2)=1e7
           emissivity(1)=1.0
           emissivity(2)=1.0
      endif

      if(kiter.ge.0)then
        nx2 = nx+1
      else
        nx2 = 1
      endif

      iread=1

      do 111 ix1=1,nx2

       ix = ix1-1

       print*,'forwardnogMCS, ix,nx,nx2 = ',ix,nx,nx2

       if(ix.gt.0)then
            xref = xn(ix)
            dx = 0.05*xref
            if(dx.eq.0)dx = 0.1

            if(ix.eq.jtan)then
             if(emiss_ang.lt.0)then
               dx=1.0
             else
               goto 111
             endif
            endif                             
            xn(ix)=xn(ix)+dx
            if(ix.eq.jsurf)tsurf=xn(ix)
       endif

       hcorr = 0.0
       print*,'Calling intradfield'
       CALL intradfield(runname,ispace,xlat,nwave1,vwave1,nconv1,
     1   vconv1,gasgiant,lin,nvar,varident,varparam,jsurf,jalb,
     2   jtan,jpre,nx,xn)
  

C      Modifying internal radiation field.
       fmod=2
C       print*,'Calling mradfield'
C       call mradfield(fmod)

       print*,'Setting up for interpolation calculations'       
C      Set up parameters for scattering cirsrad run.

       CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF)


       itype=11			! scloud11wave

       iscat=0
C      Set up all files for a direct cirsrad run of limb spectra and
C      near-limb observations
       call gsetradMCS(runname,nconv1,vconv1,fwhm,ispace,iscat,
     1     gasgiant,layht,topht,nlayer,laytyp,layint,xlat,lin,hcorrx,
     3     nvar,varident,varparam,nx,xn,jpre,tsurf,xmap)

       print*,'Calculating grid of radiances'
C      Compute suite of limb/near-limb radiances 
          call CIRSrtf_waveS(runname, dist, inormal, iray, fwhm, ispace,
     1     vwave1,nwave1,npath, output, vconv1, nconv1, itype, nem, 
     2     vem, emissivity,tsurf, calcout)

       print*,'Raw cirsrtf_waveS output : '
       do ipath=1,npath
        print*,ipath,(calcout(ipath+(j-1)*npath),j=1,nconv1)
       enddo

       ioff=0

       do 100 igeom=1,ngeom
        print*,'ForwardnogMCS. Spectrum ',igeom,' of ',ngeom
       
        nwavetmp = nwave(igeom)
        nconvtmp = nconv(igeom)

C       Search for B channels to get pixel ID of tangent point. If no B-channels
C       listed, assume that channels are in usual order

        ipixB=igeom

        do i=1,nconvtmp
         vv = vconv(igeom,i)
         if (vv.lt.350.)then
C         Get B pixel numbers
          ipixB = int(100*(0.001+vv*10-int(vv*10)))
          goto 646
         endif
        enddo
646     continue
C       Compute A pixel numbers
        ipixA = 1 + 21-ipixB

        print*,'ipixA,ipixB : ',ipixA,ipixB

        if(nav(igeom).gt.1)then
         print*,'You should not be doing explicit FOV averaging with'
         print*,'forwardnogMCS. All averaging is done implicitly.'
         print*,'Aborting...'
         stop
        endif
        iav=1

        
        sol_ang = angles(igeom,iav,1)
        emiss_ang = angles(igeom,iav,2)
        aphi = angles(igeom,iav,3)
        print*,'Iav = ',iav
        print*,'Angles : ',sol_ang,emiss_ang,aphi

        xlat = flat(igeom,iav)
        dtr = 3.1415927/180.0
 

        if(emiss_ang.lt.0.or.emiss_ang.gt.80)then
           print*,'Interpolating pre-calculated spectra'
           if(emiss_ang.lt.0)then
             htan = sol_ang+hcorr+hcorrx
           else
             htan = -RADIUS*(1.0-SIN(emiss_ang*dtr))+hcorr+hcorrx
           endif
           caltbore = altbore+hcorr+hcorrx
           print*,'forwardnogMCS: sol_ang, hcorr, hcorrx',
     1       sol_ang, hcorr,hcorrx
           print*,'htan = ',htan
           print*,'caltbore = ',caltbore
           thcentre = asin((RADIUS+htan)/SATRAD)
           thbore = asin((RADIUS+caltbore)/SATRAD)
C          Now run through wavelengths and find right wavelength in
C          pre-calculated array

           do 206 i=1,nconvtmp
            iconv=-1
            do j=1,nconv1
             if(vconv(igeom,i).eq.vconv1(j))iconv=j
            enddo
            if(iconv.lt.0)then
             print*,'forwardnogMCS iconv < 0'
             stop
            endif
            vv = vconv(igeom,i)
            esurf = interpem(nem,vem,emissivity,vv)
            ichan = int(10*(0.001+vv-int(vv)))

            print*,'vv,esurf,ichan,thetrot,thcentre,thbore',
     1        vv,esurf,ichan,thetrot,thcentre,thbore

            if(abs(thetrot).le.1.0) then
              call computeFOVA(iread,ichan,ipixA,ipixB,thcentre,
     1         thbore,thetrot,nfov,thfov,rfov)
              iread=0
            else
              call computeFOVB(ichan,ipixA,ipixB,thcentre,thbore,
     1         thetrot,nfov,thfov,rfov)
            endif

C            print*,'nfov,thcentre',nfov,thcentre

C            do j=1,nfov
C             print*,thfov(j),rfov(j)
C            enddo

C            stop

C           Interpolate radiances

            call interpolateMCSnog(nview,npath,thview,
     1 nconv1,iconv,vconv1,nlayer,nx,esurf,jsurf,ispace,tsurf,
     2 calcout,nfov,thfov,rfov,radmean)     

            if(ix.eq.0)then
             print*,iconv,vconv1(iconv),ioff,i,radmean
             yn(ioff+i) = radmean
            else
             yn1(ioff+i) = radmean
             kk(ioff+i,ix) = (yn1(ioff+i)-yn(ioff+i))/dx
            endif

206        continue

        else

C         Set up parameters for scattering cirsrad run.

          intname='intrad'
          call setup(intname,gasgiant,nmu,mu,wtmu,isol,
     1 dist,lowbc,galb,nf,nphi,layht,tsurf,nlayer1,laytyp,layint)


          CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF)

          itype=11                      ! scloud11wave
          iscat1=1

C         Set up all files for a direct cirsrad run
          call gsetrad(intname,iscat1,nmu,mu,wtmu,isol,dist,   
     1     lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,
     2     layht,nlayer1,laytyp,layint,sol_ang,emiss_ang,aphi,
     3     xlat,lin,nvar,varident,varparam,nx,xn,jalb,jtan,
     4     jpre,tsurf,xmap)

C          print*,'Calling scattering calc : sol_ang, emiss_ang',
C     1     sol_ang,emiss_ang

          call CIRSrtf_wave(intname, dist, inormal, iray, fwhm, ispace,
     1     vwave1,nwave1,npath, output, vconv1, nconv1, itype, nem, 
     2     vem, emissivity,tsurf, calcouts)

C           print*,'scattering nadir calculation napth = ',npath
C           do i=1,nconv1
C            print*,calcouts(i)
C           enddo

           do 207 i=1,nconvtmp
            iconv=-1
            do j=1,nconv1
             if(vconv(igeom,i).eq.vconv1(j))iconv=j
            enddo
            if(iconv.lt.0)then
             print*,'forwardnogMCS iconv < 0'
             stop
            endif
           
            ioff1=iconv

            if(ix.eq.0)then
             yn(ioff+i) = calcouts(ioff1)
            else
             yn1(ioff+i) = calcouts(ioff1)
             kk(ioff+i,ix) = (yn1(ioff+i)-yn(ioff+i))/dx
            endif

207        continue

        endif

        ioff = ioff + nconvtmp

100    continue

       if(ix.gt.0)then
         xn(ix)=xref
         if(ix.eq.jsurf)tsurf=xn(ix)
       endif
 
111   continue

      return

      end
