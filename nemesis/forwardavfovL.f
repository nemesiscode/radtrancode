      subroutine forwardavfovL(runname,ispace,fwhm,ngeom,nav,
     1 wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,nvar,
     2 varident,varparam,jsurf,jalb,jtan,jpre,nx,xn,ny,yn,kk)
C     $Id:
C     **************************************************************
C     Subroutine to calculate an FOV-averaged limb spectra and
C     KK-matrix using CIRSRADG gradient code. Uses numerical integration
C     for cases where the limb component is also important. 
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
C	jsurf		integer	Position of surface temperature element in
C				xn (if included)
C       jalb            integer Position of surface albedo spectrum in
C                               xn (if included)
C       jtan            integer Position of tangent height correction in
C                               xn (if included)
C       jpre            integer Position of tangent pressure in
C                               xn (if included)
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
      integer ngeom,ioff,igeom
      real interpem,dv
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/gascom.f'
      include 'arraylen.f'
      real xlat,planck_wave,planckg_wave,Bg,height(100),htan
      real wgeom(mgeom,mav),flat(mgeom,mav),fh
      integer layint,inormal,iray,itype,nlayer,laytyp
      integer iptf,imie,imie1
      integer nwave(mgeom),jsurf,nem,nav(mgeom),nwave1
      real vwave(mgeom,mwave),angles(mgeom,mav,3),vwave1(mwave)
      real calcout(maxout3),fwhm,calcoutL(maxout3)
      real calcout1(maxout3),gradients1(maxout3)
      real gradients(maxout4),vv,gradientsL(maxout3)
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      real vconv(mgeom,mconv),vconv1(mconv)
      real layht,tsurf,esurf,gradtsurf(maxout3),pressR,delp
      real xn(mx),yn(my),kk(my,mx),yn1(my)
      integer ny,iscat,jalb,jtan,jpre
      integer nphi,ipath
      integer nmu,isol,lowbc,nf
      real dist,galb,sol_ang,emiss_ang,aphi
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      real hcorr,hcorrx
      common /imiescat/ imie1
      integer nvar,varident(mvar,3),ivar
      real vem(maxsec),emissivity(maxsec)
      real varparam(mvar,mparam)
      logical gasgiant


C      print*,'ForwardavfovL'
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
       print*,'Warning from forwardavfovL'
       print*,'Can not currently do surface albedo retrievals with'
       print*,'gradient code'   
       stop
      endif

      print*,'forwardavfovL. jpre = ',jpre

C     Check that all observations have same wavelengths - necessary for
C     NemesisL

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
        if(nconv(igeom).ne.nconv1)then
         print*,'forwardavfovL: All spectra must be of same length'
         print*,'igeom,nconv(1),nconv(igeom)',igeom,nconv1,nconv(igeom)
         stop
        endif
        do i=1,nconv1
         dv = 100.0*abs(vconv(igeom,i)-vconv1(i))/vconv1(i)
         if(dv.gt.1.0)then
          print*,'forwardavfovL: wavenumbers do not match'
          print*,igeom,vconv1(i),vconv(igeom,i),dv
          stop
         endif
        enddo
       endif
      enddo

      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf,nphi,layht,tsurf,nlayer,laytyp,layint)

      if(jsurf.gt.0)then
       tsurf = xn(jsurf)
      endif

      iscat=0
C     Set up all files for a direct cirsrad run of limb spectra
      call gsetradL(runname,nconv1,vconv1,fwhm,ispace,iscat,
     1    gasgiant,layht,nlayer,laytyp,layint,xlat,lin,hcorrx,
     3    nvar,varident,varparam,nx,xn,jpre,tsurf,xmap)


C     If planet is not a gas giant then we need to read in the 
C      surface emissivity spectrum
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

      itype=12			! scloud12. not used here

      print*,'hcorrx = ',hcorrx

      call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1  vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf, 
     2  nx, xmap, vconv1, nconv1, npath, calcoutL, gradientsL)

      print*,'hcorrx = ',hcorrx

C     Read in base heights from '.drv' file
      call readdrvh(runname,height)

      print*,'hcorrx = ',hcorrx

      ioff = 0

      do 100 igeom=1,ngeom
       print*,'ForwardavfovL. Spectrum ',igeom,' of ',ngeom
       print*,'Nav = ',nav(igeom)

       nconv1 = nconv(igeom)
       nwave1 = nwave(igeom)

       do 105 i=1,nconv1
        vconv1(i)=vconv(igeom,i)
105    continue
       do 106 i=1,nwave1
        vwave1(i)=vwave(igeom,i)
106    continue

       hcorr = 0.0
       do ivar = 1,nvar
        if(varident(ivar,1).eq.777)hcorr = xn(jtan)
       enddo

       do 110 iav = 1,nav(igeom)
         sol_ang = angles(igeom,iav,1)
         emiss_ang = angles(igeom,iav,2)
         aphi = angles(igeom,iav,3)
         
         print*,'Iav = ',iav
         print*,'Angles : ',sol_ang,emiss_ang,aphi

         xlat = flat(igeom,iav)   

         if(emiss_ang.lt.0)then
           print*,'Interpolating limb-calculated spectra'
           htan = sol_ang+hcorr+hcorrx
           print*,'forwardavfovL: sol_ang, hcorr,hcorrx',
     1       sol_ang, hcorr,hcorrx
           print*,'htan = ',htan
           jpath=-1
           do i=1,nlayer-1
            if(height(i).le.htan.and.height(i+1).gt.htan)jpath=i
           enddo

           if(jpath.lt.1)then
            print*,'Error in forwardavfovL: tangent height is'
            print*,'out of range'
            print*,nlayer,height(1),height(nlayer),htan
           else
            fh = (htan-height(jpath))/(height(jpath+1)-height(jpath))
           endif

           if(htan.lt.height(1))then
            jpath=1
            fh=0.0
           endif
           if(htan.ge.height(nlayer))then
            jpath=nlayer-1
            fh=1.0
           endif
         
           do ipath=jpath,jpath+1
            fh=1.0-fh
            do j=1,nconv1
 	     ioff1=nconv1*(ipath-1)+j
             yn(ioff+j)=yn(ioff+j)+wgeom(igeom,iav)*fh*calcoutL(ioff1)
            enddo
    
            do i=1,nx
             if(i.ne.jtan.and.i.ne.jpre)then
              do j=1,nconv1 
              ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + j
              kk(ioff+j,i)=kk(ioff+j,i)+
     1	  	wgeom(igeom,iav)*fh*gradientsL(ioff2)
              enddo
             endif
            enddo
           enddo
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
     3 	    gradients)

C          Need to assume order of paths. First path is assumed to be
C          thermal emission, 2nd path is transmission to ground (if planet
C          is not a gas giant)
     
      
           ipath = 1
           do j=1,nconv1
            ioff1=nconv1*(ipath-1)+j
            yn(ioff+j)=yn(ioff+j)+wgeom(igeom,iav)*calcout(ioff1)
           enddo

           do i=1,nx
            if(i.ne.jtan.and.i.ne.jpre)then
             do j=1,nconv1
             ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + j
             kk(ioff+j,i)=kk(ioff+j,i)+
     1		wgeom(igeom,iav)*gradients(ioff2)
             enddo
            endif   
           enddo
          
C          If planet is not a gas giant and observation is not at limb then
C          we need to add the radiation from the ground
           if(.not.gasgiant.and.emiss_ang.ge.0)then
            ipath = 2
          
            do j=1,nconv1
             vv = vconv1(j)
             esurf = interpem(nem,vem,emissivity,vv)
             ioff1=nconv1*(ipath-1)+j

             Bg = planck_wave(ispace,vconv1(j),tsurf)*esurf

             yn(ioff+j)=yn(ioff+j) + wgeom(igeom,iav)*calcout(ioff1)*Bg

             do i=1,nx
              ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + j
              if(i.eq.jsurf)then
                kk(ioff+j,i)=kk(ioff+j,i)+
     1    wgeom(igeom,iav)*calcout(ioff1)*esurf*planckg_wave(ispace,
     2    vconv1(j),tsurf)
              else
                kk(ioff+j,i)=kk(ioff+j,i) +
     1               wgeom(igeom,iav)*gradients(ioff2)*Bg
              endif
             enddo
          
            enddo

           endif

         endif

110    continue

       if(jpre.gt.0)then

        print*,'Calculating RoC with tangent pressure'
        pressR = xn(jpre)
        delp = pressR*0.01
        xn(jpre)=pressR+delp

C       Set up all files to recalculate limb spectra
        call gsetradL(runname,nconv1,vconv1,fwhm,ispace,iscat,
     1    gasgiant,layht,
     2    nlayer,laytyp,layint,xlat,lin,hcorrx,
     3    nvar,varident,varparam,nx,xn,jpre,tsurf,xmap)

        call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1   vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf, 
     2   nx, xmap, vconv1, nconv1, npath, calcout1, gradients1)


        do 112 iav = 1,nav(igeom)
         sol_ang = angles(igeom,iav,1)
         emiss_ang = angles(igeom,iav,2)
         aphi = angles(igeom,iav,3)
         
         print*,'Iav = ',iav
         print*,'Angles : ',sol_ang,emiss_ang,aphi

         xlat = flat(igeom,iav)   

         if(emiss_ang.lt.0)then
           print*,'Interpolating limb-calculated spectra'
           htan = sol_ang + hcorr
           print*,'forwardfovL1: sol_ang,hcorr,htan',sol_ang,
     1      hcorr,htan
           jpath=-1
           do i=1,nlayer-1
            if(height(i).le.htan.and.height(i+1).gt.htan)jpath=i
           enddo
           if(jpath.lt.1)then
            print*,'Error in forwardavfovL: tangent height is'
            print*,'out of range'
            print*,nlayer,height(1),height(nlayer),htan
           endif
           fh = (htan-height(jpath))/(height(jpath+1)-height(jpath))
         
           print*,jpath,fh
           do ipath=jpath,jpath+1
            fh=1.0-fh
            do j=1,nconv1
 	     ioff1=nconv1*(ipath-1)+j
             yn1(ioff+j)=yn1(ioff+j)+wgeom(igeom,iav)*fh*calcout1(ioff1)
C            tangent pressure taken as logs so need to adjust gradient
             kk(ioff+j,jpre) = kk(ioff+j,jpre) +
     1            (yn1(ioff+j)-yn(ioff+j))/(pressR*delp)
             print*,ipath,ioff+j,yn1(ioff+j),yn(ioff+j),delp,pressR
            enddo
           enddo
         endif

         xn(jpre)=pressR

112     continue



       endif

       if(jtan.gt.0)then

        do 111 iav = 1,nav(igeom)
         sol_ang = angles(igeom,iav,1)
         emiss_ang = angles(igeom,iav,2)
         aphi = angles(igeom,iav,3)
         
         print*,'Iav = ',iav
         print*,'Angles : ',sol_ang,emiss_ang,aphi

         xlat = flat(igeom,iav)   

         if(emiss_ang.lt.0)then
           print*,'Interpolating limb-calculated spectra'
           htan = sol_ang + hcorr + 1.0
           print*,'forwardfovL1: sol_ang,hcorr,htan',sol_ang,
     1      hcorr,htan
           jpath=-1
           do i=1,nlayer-1
            if(height(i).le.htan.and.height(i+1).gt.htan)jpath=i
           enddo
           if(jpath.lt.1)then
            print*,'Error in forwardavfovL: tangent height is'
            print*,'out of range'
            print*,nlayer,height(1),height(nlayer),htan
           endif
           fh = (htan-height(jpath))/(height(jpath+1)-height(jpath))
         
           do ipath=jpath,jpath+1
            fh=1.0-fh
            do j=1,nconv1
 	     ioff1=nconv1*(ipath-1)+j
             yn1(ioff+j)=yn1(ioff+j)+wgeom(igeom,iav)*fh*calcoutL(ioff1)
             kk(ioff+j,jtan) = kk(ioff+j,jtan) + yn1(ioff+j)-yn(ioff+j)
            enddo
           enddo
         endif

111     continue

       endif

       ioff = ioff + nconv1

100   continue

      return

      end
