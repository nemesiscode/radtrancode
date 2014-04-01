      subroutine forwardavfovX(runname,ispace,iscat,fwhm,ngeom,
     1 nav,wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,
     2 lin,nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,RADIUS,
     3 nx,xn,ny,yn,kk)
C     $Id:
C     **************************************************************
C     Subroutine to calculate an FOV-averaged spectrum and
C     KK-matrix using CIRSRADG gradient code. Uses numerical integration
C     for cases where the limb component is also important. 
C     
C     Input variables:
C       runname(100)   character Name of run.
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
C       jrad            integer position radius element in
C                               xn (if included)
C       RADIUS          real    Planetary radius at 0km altitude
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
C     Pat Irwin 1/3/12		Updated for Nemesis2.0
C
C     **************************************************************

      implicit none
      integer i,j,lin,ispace,iav
      integer ngeom,ioff,igeom
      real interpem
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/gascom.f'
      include '../radtran/includes/planrad.f'
      include 'arraylen.f'
      real xlat,planck_wave,planckg_wave,Bg
      real wgeom(mgeom,mav),flat(mgeom,mav),pressR,delp
      integer layint,inormal,iray,itype,nlayer,laytyp,iscat
      integer nwave(mgeom),jsurf,nem,nav(mgeom),nwave1
      integer jalb,jtan,jpre,k,iptf,jrad,imie,imie1
      real vwave(mgeom,mwave),angles(mgeom,mav,3),vwave1(mwave)
      real calcout(maxout3),fwhm,RADIUS
      real gradients(maxout4),vv
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      real vconv(mgeom,mconv),vconv1(mconv)
      real layht,tsurf,esurf,gradtsurf(maxout3)
      real xn(mx),yn(my),kk(my,mx),yn1(my)
      integer ny,iconv
      integer nphi,ipath
      integer nmu,isol,lowbc,nf
      real dist,galb,sol_ang,emiss_ang,aphi
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,solfile,solname
      real xmap(maxv,maxgas+2+maxcon,maxpro),xfac,pi,xdist
      parameter (pi=3.1415927)

      integer nvar,varident(mvar,3)
      real varparam(mvar,mparam)
      logical gasgiant,fexist
      real vem(maxsec),emissivity(maxsec)
      common /imiescat/imie1

      integer cellngas,cellid(maxgas),celliso(maxgas),icread
      real cellength,cellpress,celltemp,cellvmr(maxgas)
      common/celldat/icread,cellngas,cellid,celliso,cellvmr,cellength,
     1  cellpress,celltemp


c  ** variables for solar refelcted cloud **
	real solar
      real refl_cloud_albedo
      logical reflecting_atmos
      common /refl_cloud_params/refl_cloud_albedo,reflecting_atmos

      real stelrad,solwave(maxbin),solrad(maxbin)
      integer solnpt,iform,iread

      common /solardat/iread, iform, stelrad, solwave, solrad,  solnpt

C     jradf is passed via tha planrad common block
      jradf=jrad



C      print*,runname
C      print*,ispace,iscat,fwhm,ngeom
C      print*,(nav(i),i=1,ngeom)
C      do j=1,ngeom
C       print*,(wgeom(j,i),i=1,nav(j))
C       print*,(flat(j,i),i=1,nav(j))
C       do i=1,nav(j)
C        print*,j,i,angles(j,i,k),k=1,3)  
C       enddo
C      enddo
C      do j=1,ngeom
C       print*,nwave(j),(vwave(j,i),i=1,nwave(j))
C       print*,nconv(j),(vconv(j,i),i=1,nconv(j))
C      enddo
C      print*,gasgiant,lin,nvar
C      do i=1,nvar
C       print*,(varident(i,j),j=1,3)
C       print*,(varparam(i,j),j=1,mparam)
C      enddo
C      print*,jsurf,jalb,jtan,jpre
C      print*,nx,ny
C      print*,(xn(i),i=1,nx)

             
C     Initialise arrays
      do i=1,my
       yn(i)=0.0
       yn1(i)=0.0
       do j=1,mx
        kk(i,j)=0.0
       enddo
      enddo

      print*,'forwardavfovX: jsurf,jalb,jtan,jpre',jsurf,jalb,jtan,jpre
      if(jalb.gt.1)then
       print*,'Warning from forwardavfovX'
       print*,'Can not currently do surface albedo retrievals with'
       print*,'gradient code'
       stop
      endif

      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf,nphi,layht,tsurf,nlayer,laytyp,layint)

      if(jsurf.gt.0)then
       tsurf = xn(jsurf)
       print*,'tsurf = ',tsurf
      endif


C     If we're retrieving planet radius then add correction to reference
C     radius
C     N.B.radius2 is passed via the planrad common block.
      if(jrad.gt.0)then
        radius2 = xn(jrad) + radius
      else
        radius2 = radius
      endif

      ioff = 0
      do 100 igeom=1,ngeom
       print*,'ForwardavfovX. Spectrum ',igeom,' of ',ngeom

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

C        Set up parameters for non-scattering cirsrad run.
         CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF,IMIE)
         IMIE1=IMIE
         itype=11			! scloud11. not used here


C        Set up all files for a direct cirsrad run
         call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,
     1    lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,
     2    layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,lin,
     3    nvar,varident,varparam,nx,xn,jalb,jtan,jpre,tsurf,xmap)


         call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1    vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf, 
     2    nx, xmap, vconv1, nconv1, npath, calcout, gradients)

C        Need to assume order of paths. First path is assumed to be
C        thermal emission


         print*,'ForwardavfovX: Npath = ',npath

         if(icread.ne.1)then
C         Not an SCR calculation. Assume 1st path is the thermal emission
          ipath = 1
          do j=1,nconv1
           iconv=-1
           do k = 1,nconv1
            if(vconv(igeom,j).eq.vconv1(k))iconv=k
           enddo
           if(iconv.lt.0)then
            print*,'Error in forwardavfovX iconv<1'
            stop
           endif
 	   ioff1=nconv1*(ipath-1)+iconv
           yn(ioff+j)=yn(ioff+j)+wgeom(igeom,iav)*calcout(ioff1)
          enddo
    
C         Calculate gradients
          do i=1,nx

           if(i.ne.jtan.and.i.ne.jpre.and.i.ne.jrad)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
             ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
             kk(ioff+j,i)=kk(ioff+j,i)+wgeom(igeom,iav)*gradients(ioff2)
            enddo
           endif

           if(i.eq.jrad)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+j,i)=kk(ioff+j,i)+wgeom(igeom,iav)*
     1          2.*calcout(ioff1)/RADIUS
            enddo
           endif

           if(i.eq.jsurf)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+j,i)=kk(ioff+j,i)+wgeom(igeom,iav)*
     1          gradtsurf(ioff1)
            enddo
           endif

          enddo

  
          if (gasgiant.and.reflecting_atmos) then
           print*,'ADDING IN REFLECTING ATMOSPHERE CONTRIBUTION'
           print*,'cloud albedo=',refl_cloud_albedo
           
C          Assume this calculation is in path 2.
           ipath = 2
 
c	   initialise solar spectrum
           CALL FILE(runname,solfile,'sol')
           inquire(file=solfile,exist=fexist)
           if(fexist) then
            call opensol(solfile,solname)
           else
            print*,'Error in forwardavfovX. solar flux file not defined'
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
         
 	    CALL get_solar_wave(vconv1(j),dist,solar)

            Bg = solar*refl_cloud_albedo/3.141592654
           
            yn(ioff+j)=yn(ioff+j) + wgeom(igeom,iav)*calcout(ioff1)*Bg
                      
            do i=1,nx
             ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
             kk(ioff+j,i)=kk(ioff+j,i) + 
     1                wgeom(igeom,iav)*gradients(ioff2)*Bg
            enddo

           enddo
          
          endif

         else
C         We have a cell defined, which means we have two outputs, wideband and 
C         sideband. Output the sideband first (path=4) and then 
C         the wideband (path=5)        

          do j=1,nconv1
           iconv=-1
           do k = 1,nconv1
            if(vconv(igeom,j).eq.vconv1(k))iconv=k
           enddo
           if(iconv.lt.0)then
            print*,'Error in forwardavfovX iconv<1'
            stop
           endif
           ipath=4
 	   ioff1=nconv1*(ipath-1)+iconv
           yn(ioff+j)=yn(ioff+j)+wgeom(igeom,iav)*calcout(ioff1)
           ipath=5
 	   ioff1=nconv1*(ipath-1)+iconv
           yn(ioff+nconv1+j)=yn(ioff+nconv1+j)+wgeom(igeom,iav)*
     1		calcout(ioff1)
          enddo
    
          do i=1,nx

C          Now the gradients
           if(i.ne.jtan.and.i.ne.jpre.and.i.ne.jrad)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
             ipath=4
             ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
             kk(ioff+j,i)=kk(ioff+j,i)+wgeom(igeom,iav)*gradients(ioff2)
             ipath=5
             ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
             kk(ioff+nconv1+j,i)=kk(ioff+nconv1+j,i)+wgeom(igeom,iav)*
     1		gradients(ioff2)
            enddo
           endif

           if(i.eq.jrad)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
             ipath=4
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+j,i)=kk(ioff+j,i)+wgeom(igeom,iav)*
     1          2.*calcout(ioff1)/RADIUS
             ipath=5
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+nconv1+j,i)=kk(ioff+nconv1+j,i)+wgeom(igeom,iav)*
     1          2.*calcout(ioff1)/RADIUS
            enddo
           endif

           if(i.eq.jsurf)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
             ipath=4
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+j,i)=kk(ioff+j,i)+wgeom(igeom,iav)*
     1          gradtsurf(ioff1)
             ipath=5
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+nconv1+j,i)=kk(ioff+nconv1+j,i)+wgeom(igeom,iav)*
     1          gradtsurf(ioff1)
            enddo
           endif

          enddo


         endif


110    continue

       if(jtan.gt.0.or.jpre.gt.0)then

        if(jpre.gt.0)then
          print*,'Calculating RoC with tangent pressure'
          pressR = xn(jpre)
          delp = pressR*0.01
          xn(jpre)=pressR+delp
        endif

        do 111 iav = 1,nav(igeom)
         sol_ang = angles(igeom,iav,1)
         emiss_ang = angles(igeom,iav,2)
         aphi = angles(igeom,iav,3)

         if(jtan.gt.0)then
          if(emiss_ang.lt.0)sol_ang = sol_ang+1.0
         endif

         xlat = flat(igeom,iav)   

C        Set up parameters for non-scattering cirsrad run.

         CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF,IMIE)
         IMIE1=IMIE

         itype=11			! scloud11. not used here


C        Set up all files for a direct cirsrad run
         call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,
     1    lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,
     2    layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,lin,
     3    nvar,varident,varparam,nx,xn,jalb,jtan,jpre,tsurf,xmap)


         call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1    vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf, 
     2    nx, xmap, vconv1, nconv1, npath, calcout, gradients)


         if(icread.ne.1)then
C         First path is assumed to be thermal emission if not SCR calculation
        
          ipath = 1
          do j=1,nconv1
           iconv=-1
           do k=1,nconv1
            if(vconv(igeom,j).eq.vconv1(k))iconv=k
           enddo 
 	   ioff1=nconv1*(ipath-1)+iconv
           yn1(ioff+j)=yn1(ioff+j)+wgeom(igeom,iav)*calcout(ioff1)
          enddo

         else

          do j=1,nconv1
           iconv=-1
           do k=1,nconv1
            if(vconv(igeom,j).eq.vconv1(k))iconv=k
           enddo 
           ipath=4
 	   ioff1=nconv1*(ipath-1)+iconv
           yn1(ioff+j)=yn1(ioff+j)+wgeom(igeom,iav)*calcout(ioff1)
           ipath=5
 	   ioff1=nconv1*(ipath-1)+iconv
           yn1(ioff+nconv1+j)=yn1(ioff+nconv1+j)+
     1		wgeom(igeom,iav)*calcout(ioff1)
          enddo

         endif


111      continue

         do j=1,nconv1

          if(jtan.gt.0)then
C          Assume change in tangent height pressure of 1km.
           kk(ioff+j,jtan) = kk(ioff+j,jtan) + yn1(ioff+j)-yn(ioff+j)
           if(icread.eq.1)then
            kk(ioff+nconv1+j,jtan) = kk(ioff+nconv1+j,jtan) + 
     1		yn1(ioff+nconv1+j)-yn(ioff+nconv1+j)
           endif
          elseif(jpre.gt.0)then
           kk(ioff+j,jpre) = kk(ioff+j,jpre) +
     1        (yn1(ioff+j)-yn(ioff+j))/delp  
           if(icread.eq.1)then
            kk(ioff+nconv1+j,jpre) = kk(ioff+nconv1+j,jpre) +
     1        (yn1(ioff+nconv1+j)-yn(ioff+nconv1+j))/delp  
           endif
          endif

         enddo
    

         if(jpre.gt.0)then
          xn(jpre)=pressR
         endif
 
       endif

       if(icread.ne.1)then
        ioff = ioff + nconv1
       else
        ioff = ioff + 2*nconv1
       endif

100   continue

      return

      end
