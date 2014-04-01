      subroutine forwardavfovMCS(runname,ispace,fwhm,xlat,ngeom,nav,
     1 wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,nvar,
     2 varident,varparam,jsurf,jalb,jtan,jpre,RADIUS,SATRAD,thetrot,
     3 altbore,nx,xn,ny,yn,kk)
C     $Id:
C     **************************************************************
C     Subroutine to calculate an FOV-averaged limb spectra and
C     KK-matrix using CIRSRADG gradient code. Uses numerical integration
C     for cases where the limb component is also important. 
C     
C     Input variables:
C       runname(100)   character Name of run.
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
C	jsurf		integer	Position of surface temperature element in
C				xn (if included)
C       jalb            integer Position of surface albedo spectrum in
C                               xn (if included)
C       jtan            integer Position of tangent height correction in
C                               xn (if included)
C       jpre            integer Position of tangent pressure in
C                               xn (if included)
C       RADIUS		real    Radius of planet at 0km tangent altitude
C       SATRAD          real    Distance of satellite from centre
C				of planet
C       thetrot         real    Angle of rotation of pixel array from 
C                               vertical
C       altbore         real    Altitude of tangent boresight above
C                               surface
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
C
C     **************************************************************

      implicit none
      integer i,j,k,lin,ispace,iav,jpath
      integer ngeom,ioff,igeom,iread
      real interpem,dv
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/gascom.f'
      include 'arraylen.f'
      real xlat,planck_wave,planckg_wave,Bg,height(100),htan
      real wgeom(mgeom,mav),flat(mgeom,mav),fh,thetrot
      integer layint,inormal,iray,itype,nlayer,laytyp
      integer nwave(mgeom),jsurf,nem,nav(mgeom),nwave1
      real vwave(mgeom,mwave),angles(mgeom,mav,3),vwave1(mwave)
      real calcout(maxout3),fwhm,calcoutL(maxout3)
      real calcout1(maxout3),gradients1(maxout4)
      real calcout2(maxout3),gradients2(maxout4)
      real gradients(maxout4),vv,gradientsL(maxout4)
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      integer ipixA,ipixB,ichan,iptf,imie,imie1
      real vconv(mgeom,mconv),vconv1(mconv)
      real layht,tsurf,esurf,gradtsurf(maxout3),pressR
      real delp,altbore,thbore
      real xn(mx),yn(my),kk(my,mx),yn1(my),caltbore
      integer ny,iscat,jalb,jtan,jpre
      integer nphi,ipath
      integer nmu,isol,lowbc,nf
      real dist,galb,sol_ang,emiss_ang,aphi
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      real hcorr,hcorrx
      common /imiescat/imie1

      integer nconvtmp,nwavetmp,iflag,i1,iswitch
      integer nview,mview,iconv,nfov,mfov
      real vem(maxsec),emissivity(maxsec)
      parameter (mfov=200,mview=100)

      real vtmp,tmp,topht,delh,dtr,thcentre,RADIUS
      real SATRAD,hview(mview),thview(mview)
      real thfov(mfov),rfov(mfov),radmean,gradmean(mx)
      integer nvar,varident(mvar,3),ivar
      real varparam(mvar,mparam)
      logical gasgiant


C      print*,'ForwardavfovMCS'
C      print*,runname
C      print*,ispace,fwhm,ngeom
C      do i=1,ngeom
C       print*,nav(i)
C       do j=1,nav(i)
C        print*,wgeom(i,j),flat(i,j),(angles(i,j,k),k=1,3)
C       enddo
C       print*,nwave(i)
C       print*,(vwave(i,j),j=1,nwave(i))
C       print*,nconv(i)
C       print*,(vconv(i,j),j=1,nconv(i))
C       print*,gasgiant,lin,nvar
C      enddo
C      do i=1,nvar
C       print*,(varident(i,j),j=1,3)     
C       print*,(varparam(i,j),j=1,5)     
C      enddo
C      print*,jsurf,jalb,jtan,jpre,nx
C      print*,(xn(i),i=1,nx)
C      print*,ny


C     Initialise arrays
      do i=1,my
       yn(i)=0.0
       yn1(i)=0.0
       do j=1,mx
        kk(i,j)=0.0
       enddo
      enddo


      if(jalb.gt.1)then
       print*,'Warning from forwardavfovMCS'
       print*,'Can not currently do surface albedo retrievals with'
       print*,'gradient code'   
       stop
      endif

      print*,'forwardavfovMCS. jpre = ',jpre

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

      print*,'ForwardavfovMCS vconv:'
      do i=1,nconv1
       print*,i,vconv1(i)
      enddo

      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf,nphi,layht,tsurf,nlayer,laytyp,layint)

      if(jsurf.gt.0)then
       tsurf = xn(jsurf)
      endif

      iscat=0
      topht = 100.0


C     Set up all files for a direct cirsrad run of limb spectra and
C     near-limb observations
      call gsetradMCS(runname,nconv1,vconv1,fwhm,ispace,iscat,
     1    gasgiant,layht,topht,nlayer,laytyp,layint,xlat,lin,hcorrx,
     3    nvar,varident,varparam,nx,xn,jpre,tsurf,xmap)

C     If planet is not a gas giant then we need to read in the surface 
C      emissivity spectrum
      if(.not.gasgiant)then
          call readsurfem(runname,nem,vem,emissivity)
      else
           nem=2
           vem(1)=-100.0
           vem(2)=1e7
           emissivity(1)=1.0
           emissivity(2)=1.0
      endif

C     Set up parameters for non-scattering cirsrad run.

      CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF,IMIE)
      IMIE1=IMIE

      itype=12                  ! scloud12. not used here

      print*,'hcorrx = ',hcorrx

      call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1  vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf, 
     2  nx, xmap, vconv1, nconv1, npath,calcoutL, gradientsL)


      if(jpre.gt.0)then

        print*,'Setting up arrays for calculating  RoC'
        print*,'with tangent pressure'
        pressR = xn(jpre)
        delp = pressR*0.01
        xn(jpre)=pressR+delp

C       Set up all files to recalculate limb spectra
        call gsetradMCS(runname,nconv1,vconv1,fwhm,ispace,iscat,
     1    gasgiant,layht,topht,nlayer,laytyp,layint,xlat,lin,hcorrx,
     3    nvar,varident,varparam,nx,xn,jpre,tsurf,xmap)

        call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace,
     1   vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf, 
     2   nx, xmap, vconv1, nconv1, npath, calcout1, gradients1)

         xn(jpre)=pressR

      endif

      print*,'hcorrx = ',hcorrx

C     Read in base heights from '.drv' file
      call readdrvh(runname,height)

      print*,'hcorrx = ',hcorrx

      ioff = 0

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
      do i=1,nview
        thview(i) = asin((RADIUS+hview(i))/SATRAD)
        print*,i,hview(i),thview(i)
      enddo

C      print*,'nconv1 = ',nconv1
C      do ipath=1,nlayer
C       print*,ipath,(calcoutL(nconv1*(ipath-1)+iconv),iconv=1,nconv1)
C      enddo
C      do i = nlayer+1,nview  
C       ipath = nlayer+1 + (i-nlayer-1)*2
C       print*,ipath,(calcoutL(nconv1*(ipath-1)+iconv),iconv=1,nconv1)
C       ipath=ipath+1
C       print*,ipath,(calcoutL(nconv1*(ipath-1)+iconv),iconv=1,nconv1)       
C      enddo

      iread=1

      do 100 igeom=1,ngeom
       print*,'ForwardavfovMCS. Obs ',igeom,' of ',ngeom

       print*,'Nav = ',nav(igeom)

       nconvtmp = nconv(igeom)
       nwavetmp = nwave(igeom)
C       Search for B channels to get pixel ID of tangent point. If no B-channels
C       listed, assume that channels are in usual order

       ipixB = igeom

       do i=1,nconvtmp
        vv = vconv(igeom,i)
        if (vv.lt.350.)then
C        Get B pixel numbers 
         ipixB = int(100*(0.001+vv*10-int(vv*10)))
         goto 646 
        endif
       enddo
646    continue
C      Compute A pixel numbers
       ipixA = 1 + 21-ipixB

       print*,'vv, ipixA, ipixB',vv,ipixA,ipixB

       hcorr = 0.0
       do ivar = 1,nvar
        if(varident(ivar,1).eq.777)hcorr = xn(jtan)
       enddo

       if(nav(igeom).gt.1)then
        print*,'You should not be doing explicit FOV averaging with'
        print*,'forwardavfovMCS. All averaging is done implicitly.'
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
           print*,'forwardavfovMCS: sol_ang, hcorr, hcorrx',
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
             print*,'forwardavfovMCS iconv < 0'
             stop
            endif           
            vv = vconv(igeom,i)
            esurf = interpem(nem,vem,emissivity,vv)
            ichan = int(10*(0.001+vv-int(vv)))

            print*,'vv,esurf,ichan,thetrot',vv,esurf,ichan,thetrot

            if(abs(thetrot).le.1.0) then
              call computeFOVA(iread,ichan,ipixA,ipixB,thcentre,thbore,
     1        thetrot,nfov,thfov,rfov)
              iread=0
            else
              call computeFOVB(ichan,ipixA,ipixB,thcentre,thbore,
     1        thetrot,nfov,thfov,rfov)
            endif

            call interpolateMCS(nview,thview,
     1 nconv1,iconv,vconv1,nlayer,nx,esurf,jsurf,ispace,tsurf,
     2 calcoutL,gradientsL,nfov,thfov,rfov,radmean,gradmean)


            yn(ioff+i) = radmean

            do j=1,nx
             kk(ioff+i,j) = gradmean(j)
            enddo

206        continue

       else

           print*,'Calculating new nadir-spectra'
           iscat = 0
           call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,
     1      galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,layht,
     2      nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,lin,
     3      nvar,varident,varparam,nx,xn,jalb,jtan,jpre,tsurf,xmap)
      

           call CIRSrtfg_wave(runname,dist,inormal,iray,fwhm,ispace,
     1      vwave1,nwave1,itype, nem, vem, emissivity, tsurf, 
     2      gradtsurf, nx, xmap, vconv1, nconv1, npath, calcout, 
     3      gradients)

C          First path is here assumed to be thermal emission
     
      
           ipath = 1
           do j=1,nconvtmp
            vv = vconv(igeom,j)
            iconv=-1
            do k=1,nconv1
             if(vv.eq.vconv1(k))iconv=k
            enddo
            ioff1=nconv1*(ipath-1)+iconv
            print*,j,ioff1,calcout(ioff1)
            yn(ioff+j)=yn(ioff+j)+wgeom(igeom,iav)*calcout(ioff1)
           enddo

           do i=1,nx
            if(i.ne.jtan.and.i.ne.jpre)then
             do j=1,nconvtmp
              vv = vconv(igeom,j)
              iconv=-1
              do k=1,nconv1
               if(vv.eq.vconv1(k))iconv=k
              enddo
              ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
              kk(ioff+j,i)=kk(ioff+j,i)+
     1		wgeom(igeom,iav)*gradients(ioff2)
             enddo
            endif   
           enddo
          
       endif


       if(jpre.gt.0)then

        print*,'Calculating RoC with tangent pressure'

        iav=1
        sol_ang = angles(igeom,iav,1)
        emiss_ang = angles(igeom,iav,2)
        aphi = angles(igeom,iav,3)
         
        print*,'Iav = ',iav
        print*,'Angles : ',sol_ang,emiss_ang,aphi

        xlat = flat(igeom,iav)   

        if(emiss_ang.lt.0.or.emiss_ang.gt.80)then
           print*,'Interpolating limb-calculated spectra'
           if(emiss_ang.lt.0)then
             htan = sol_ang+hcorr+hcorrx
           else
             htan = -RADIUS*(1.0-SIN(emiss_ang*dtr))+hcorr+hcorrx
           endif

           htan = sol_ang + hcorr
           print*,'forwardfovMCS: sol_ang,hcorr,htan',sol_ang,
     1      hcorr,htan
           thcentre = asin((RADIUS+htan)/SATRAD)
 
C          Now run through wavelengths and find right wavelength in 
C          pre-calculated array
           do 207 i=1,nconvtmp
            iconv=-1
            do j=1,nconv1
             if(vconv(igeom,i).eq.vconv1(j))iconv=j
            enddo
            print*,i,vconv(igeom,i),iconv
            if(iconv.lt.0)then
             print*,'forwardavfovMCS iconv < 0'
             stop
            endif           
            esurf = interpem(nem,vem,emissivity,vconv1(iconv))
            vv = vconv(igeom,i)
            ichan = int(10*(0.001+vv-int(vv)))

            print*,'PRESS vv,ichan,thetrot',vv,ichan,thetrot

            if(abs(thetrot).le.1.0) then
              call computeFOVA(iread,ichan,ipixA,ipixB,thcentre,thbore,
     1        thetrot,nfov,thfov,rfov)
              iread=0
            else
              call computeFOVB(ichan,ipixA,ipixB,thcentre,thbore,
     1        thetrot,nfov,thfov,rfov)
            endif

            call interpolateMCS(nview,thview,
     1 nconv1,iconv,vconv1,nlayer,nx,esurf,jsurf,ispace,tsurf,
     2 calcout1,gradients1,nfov,thfov,rfov,radmean,gradmean)

            yn1(ioff+i) = radmean

            kk(ioff+i,jpre) = kk(ioff+i,jpre) +
     1            (yn1(ioff+i)-yn(ioff+i))/delp

207        continue
        else

           print*,'Calculating new nadir-spectra'
           iscat = 0
           xn(jpre)=xn(jpre)+delp
           call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,
     1      galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,layht,
     2      nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,lin,
     3      nvar,varident,varparam,nx,xn,jalb,jtan,jpre,tsurf,xmap)


           call CIRSrtfg_wave(runname,dist,inormal,iray,fwhm,ispace,
     1      vwave1,nwave1,itype, nem, vem, emissivity, tsurf, 
     2      gradtsurf, nx, xmap, vconv1, nconv1, npath, calcout2, 
     3      gradients2)

C          Need to assume order of paths. First path is assumed to be
C          thermal emission

           ipath=1
           do j=1,nconvtmp
            vv = vconv(igeom,j)
            iconv=-1
            do k=1,nconv1
              if(vv.eq.vconv1(k))iconv=k
            enddo
            ioff1=nconvtmp*(ipath-1)+iconv

            yn1(ioff+j)=calcout2(ioff1)

            kk(ioff+j,jpre) = (yn1(ioff+j)-yn(ioff+j))/
     1       (pressR*delp)

           enddo
        
           xn(jpre)=xn(jpre)-delp

        endif

       endif

       if(jtan.gt.0)then

        iav=1
        sol_ang = angles(igeom,iav,1)
        emiss_ang = angles(igeom,iav,2)
        aphi = angles(igeom,iav,3)
         
        print*,'Iav = ',iav
        print*,'Angles : ',sol_ang,emiss_ang,aphi

        xlat = flat(igeom,iav)   

        if(emiss_ang.lt.0.or.emiss_ang.gt.80)then
           print*,'Interpolating limb-calculated spectra'
           htan = sol_ang + hcorr + 1.0
           caltbore = altbore + hcorr + 1.0
           print*,'forwardfovMCS: sol_ang,hcorr,htan',sol_ang,
     1      hcorr,htan

           thcentre = asin((RADIUS+htan)/SATRAD)
           thbore = asin((RADIUS+caltbore)/SATRAD)
 
C          Now run through wavelengths and find right wavelength in 
C          pre-calculated array
           do 208 i=1,nconvtmp
            iconv=-1
            do j=1,nconv1
             if(vconv(igeom,i).eq.vconv1(j))iconv=j
            enddo
            if(iconv.lt.0)then
             print*,'forwardavfovMCS iconv < 0'
             stop
            endif           
            esurf = interpem(nem,vem,emissivity,vconv1(iconv))
            vv = vconv(igeom,i)
            ichan = int(10*(0.001+vv-int(vv)))

            if(abs(thetrot).le.1.0) then
              call computeFOVA(iread,ichan,ipixA,ipixB,thcentre,thbore,
     1        thetrot,nfov,thfov,rfov)
              iread=0
            else
              call computeFOVB(ichan,ipixA,ipixB,thcentre,thbore,
     1        thetrot,nfov,thfov,rfov)
            endif

            call interpolateMCS(nview,thview,
     1 nconv1,iconv,vconv1,nlayer,nx,esurf,jsurf,ispace,tsurf,
     2 calcoutL,gradientsL,nfov,thfov,rfov,radmean,gradmean)

            yn1(ioff+i) = radmean

            kk(ioff+i,jtan) = kk(ioff+i,jtan) +
     1            (yn1(ioff+i)-yn(ioff+i))

208        continue

        endif

       endif

       ioff = ioff + nconvtmp

100   continue

      return

      end
