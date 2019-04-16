      subroutine forwardavfovX(runname,ispace,iscat,fwhm,ngeom,
     1 nav,wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,
     2 lin,nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3 jlogg,jfrac,RADIUS,nx,xn,ny,yn,kk)
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
C	flon(mgeom,mav)	real	Integration point longitudes
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
C	jxsc		integer	Position of x-section spectrum in
C				xn (if included)
C	jtan		integer	Position of tangent height correction in
C				xn (if included)
C	jpre		integer	Position of tangent pressure in
C				xn (if included)
C       jrad            integer position radius element in
C                               xn (if included)
C	jlogg		integer	position of surface log_10(g) in
C                               xn (if included)
C	jfrac		integer	position of profile fraction in
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
      real xlat,planck_wave,planckg_wave,Bg,Grav
      parameter (Grav=6.672E-11)
      real wgeom(mgeom,mav),flat(mgeom,mav),pressR,delp
      real loggR,dellg,flon(mgeom,mav),xlon,xgeom
      integer layint,inormal,iray,itype,nlayer,laytyp,iscat
      integer nwave(mgeom),jsurf,nem,nav(mgeom),nwave1
      integer jalb,jxsc,jtan,jpre,k,iptf,jrad,imie,imie1,jlogg
      integer jfrac
      real vwave(mgeom,mwave),angles(mgeom,mav,3),vwave1(mwave)
      real calcout(maxout3),fwhm,RADIUS
      real gradients(maxout4),vv
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      real vconv(mgeom,mconv),vconv1(mconv)
      real layht,tsurf,esurf,gradtsurf(maxout3)
      real xn(mx),yn(my),kk(my,mx),yn1(my),ytmp(my)
      integer ny,iconv,iextra
      integer nphi,ipath,jx
      integer nmu,isol,lowbc,nf
      real dist,galb,sol_ang,emiss_ang,aphi
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,solfile,solname
      real xmap(maxv,maxgas+2+maxcon,maxpro),xfac,pi,xdist
      real v1(3),v0(3),xx,xp,xlat1,xlon1
      integer ix,ivar,npvar,i1,ifov
      parameter (pi=3.1415927)
      character*1 ans

      integer nvar,varident(mvar,3),npro,nvmr
      real varparam(mvar,mparam)
      logical gasgiant,fexist,gasgiant1,ipfov
      real vem(maxsec),emissivity(maxsec)
      common /imiescat/imie1

      integer cellngas,cellid(maxgas),celliso(maxgas),icread
      real cellength,cellpress,celltemp,cellvmr(maxgas)
      common/celldat/icread,cellngas,cellid,celliso,cellvmr,cellength,
     1  cellpress,celltemp


c  ** variables for solar reflected cloud **
	real solar
      real refl_cloud_albedo
      logical reflecting_atmos
      common /refl_cloud_params/refl_cloud_albedo,reflecting_atmos

      real stelrad,solwave(maxbin),solrad(maxbin)
      integer solnpt,iform,iread

      common /solardat/iread, iform, stelrad, solwave, solrad,  solnpt


1     format(a)

C     jradf and jloggf are passed via the planrad common block
      jradf=jrad
      jloggf=jlogg


      call readrefhead(runname,npro,nvmr,gasgiant1)

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

      if(jxsc.gt.1)then
       print*,'Warning from forwardavfovX'
       print*,'Can not currently do x-section retrievals with'
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

C     If we're retrieving surface gravity then modify the planet mass
C     N.B. mass2 is passed via the planrad common block. Assume xn(jlogg)
C     holds log_10(surface gravity in cm/s^2). Need factor of 1e-20 to convert
C     mass to units of 1e24 kg.
      if(jlogg.gt.0)then
        mass2 = 1e-20*10**(xn(jlogg))*(radius2**2)/Grav
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

       ipfov=.false.
       if(nav(igeom).gt.100)then
        ipfov=.true.
        ifov=76
        open(ifov,file='fov.txt',status='unknown')
        write(ifov,*)nav(igeom),nconv1
       endif

       do 110 iav = 1,nav(igeom)
         sol_ang = angles(igeom,iav,1)
         emiss_ang = angles(igeom,iav,2)
         aphi = angles(igeom,iav,3)

         xlat = flat(igeom,iav)   
         xlon = flon(igeom,iav)   
         xgeom = wgeom(igeom,iav)

C         print*,'t1',jfrac,xgeom
         if(jfrac.gt.0)then
C          print*,'test1'
          if(nav(igeom).ne.2)then
           print*,'Error in forwardavfovX'
           print*,'Model 102 only suitable for NAV=2'
           stop
          endif
          if(iav.eq.1)then
           xgeom=xn(jfrac)
          else
           xgeom=1.0 - xn(jfrac)
          endif
          print*,'t2',iav,xgeom
       
         endif

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
         CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,INH3,
     1    IPTF,IMIE, iuvscat)
         IMIE1=IMIE
         itype=11			! scloud11. not used here


C        Set up all files for a direct cirsrad run
         call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,
     1    lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,
     2    layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,xlon,
     3    lin,nvar,varident,varparam,nx,xn,jalb,jxsc,jtan,jpre,tsurf,
     4    xmap)


         call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1    vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf, 
     2    nx, xmap, vconv1, nconv1, npath, calcout, gradients,iscat)

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
           yn(ioff+j)=yn(ioff+j)+xgeom*calcout(ioff1)
           if(ipfov)then
             ytmp(j)=calcout(ioff1)
           endif
          enddo
    
          if(ipfov)then
           write(ifov,*)iav,xlat,xlon,xgeom,
     1       (ytmp(j),j=1,nconv1)
          endif

C          print*,'A'
C         Calculate gradients
          do i=1,nx

           if(i.ne.jtan.and.i.ne.jpre.and.i.ne.jrad.and.i.ne.jlogg.
     1 and.i.ne.jfrac)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
             ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
             kk(ioff+j,i)=kk(ioff+j,i)+xgeom*gradients(ioff2)
C             print*,i,j,ioff,ioff2,xgeom,gradients(ioff2),kk(ioff+j,i)
            enddo
           endif           

C          Model 102 - weighted average of two profiles.
           if(i.eq.jfrac)then

            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
             ioff1=nconv1*(ipath-1)+iconv
             if(iav.eq.1)then
              kk(ioff+j,i)=kk(ioff+j,i)+calcout(ioff1)
             else
              kk(ioff+j,i)=kk(ioff+j,i)-calcout(ioff1)
             endif
            enddo
           endif

           if(i.eq.jrad)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+j,i)=kk(ioff+j,i)+xgeom*
     1          2.*calcout(ioff1)/RADIUS2
            enddo
           endif

           if(i.eq.jsurf)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+j,i)=kk(ioff+j,i)+xgeom*
     1          gradtsurf(ioff1)
            enddo
           endif

          enddo


C          open(12,file='grad2.txt',status='unknown')
C          write(12,*)sol_ang,emiss_ang,xgeom,nconv1,nx
C          do j=1,nconv1
C           write(12,*)(kk(ioff+j,i),i=1,nx)
C          enddo
C          close(12)
C          read(5,1)ans



C          print*,'B'
  
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
           
            yn(ioff+j)=yn(ioff+j) + xgeom*calcout(ioff1)*Bg
                      
            do i=1,nx
             ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
             kk(ioff+j,i)=kk(ioff+j,i) + 
     1                xgeom*gradients(ioff2)*Bg
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
           yn(ioff+j)=yn(ioff+j)+xgeom*calcout(ioff1)
           ipath=5
 	   ioff1=nconv1*(ipath-1)+iconv
           yn(ioff+nconv1+j)=yn(ioff+nconv1+j)+xgeom*
     1		calcout(ioff1)
          enddo
    
          do i=1,nx

C          Now the gradients
           if(i.ne.jtan.and.i.ne.jpre.and.i.ne.jrad.and.i.ne.jlogg)then
            do j=1,nconv1
             iconv=-1
             do k = 1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
             ipath=4
             ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
             kk(ioff+j,i)=kk(ioff+j,i)+xgeom*gradients(ioff2)
             ipath=5
             ioff2 = nconv1*nx*(ipath-1)+(i-1)*nconv1 + iconv
             kk(ioff+nconv1+j,i)=kk(ioff+nconv1+j,i)+xgeom*
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
             kk(ioff+j,i)=kk(ioff+j,i)+xgeom*
     1          2.*calcout(ioff1)/RADIUS2
             ipath=5
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+nconv1+j,i)=kk(ioff+nconv1+j,i)+xgeom*
     1          2.*calcout(ioff1)/RADIUS2
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
             kk(ioff+j,i)=kk(ioff+j,i)+xgeom*
     1          gradtsurf(ioff1)
             ipath=5
  	     ioff1=nconv1*(ipath-1)+iconv
             kk(ioff+nconv1+j,i)=kk(ioff+nconv1+j,i)+xgeom*
     1          gradtsurf(ioff1)
            enddo
           endif

          enddo


         endif


110    continue

       if(ipfov)then
        close(ifov)
        stop
       endif

C       print*,'C'
       if(jtan.gt.0.or.jpre.gt.0.or.jlogg.gt.0)then

        do 113 iextra=1,3

         if(iextra.eq.1.and.jpre.lt.1)goto 113
         if(iextra.eq.2.and.jtan.lt.1)goto 113
         if(iextra.eq.3.and.jlogg.lt.1)goto 113

         if(iextra.eq.1.and.jpre.gt.0)then
          print*,'Calculating RoC with tangent pressure'
          pressR = xn(jpre)
          delp = pressR*0.01
          xn(jpre)=pressR+delp
         endif

         if(iextra.eq.3.and.jlogg.gt.0)then
          print*,'Calculating RoC with log(g)'
          loggR = xn(jlogg)
          dellg = loggR*0.01
          xn(jlogg)=loggR+dellg
          mass2 = 1e-20*10**(xn(jlogg))*(radius2**2)/Grav
         endif


         do 111 iav = 1,nav(igeom)
          sol_ang = angles(igeom,iav,1)
          emiss_ang = angles(igeom,iav,2)
          aphi = angles(igeom,iav,3)

          if(iextra.eq.2.and.jtan.gt.0)then
           if(emiss_ang.lt.0)sol_ang = sol_ang+1.0
          endif

          xlat = flat(igeom,iav)   
          xlon = flon(igeom,iav)   
          xgeom = wgeom(igeom,iav)   

C         Set up parameters for non-scattering cirsrad run.
 
          CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,INH3,
     1     IPTF,IMIE, iuvscat)
          IMIE1=IMIE

          itype=11			! scloud11. not used here


C         Set up all files for a direct cirsrad run
          call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,
     1     lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,
     2     layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,xlon,
     3     lin,nvar,varident,varparam,nx,xn,jalb,jxsc,jtan,jpre,tsurf,
     4     xmap)


          call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1     vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf, 
     2     nx, xmap, vconv1, nconv1, npath, calcout, gradients,iscat) 


          if(icread.ne.1)then
C          First path is assumed to be thermal emission if not SCR calculation
        
           ipath = 1
           do j=1,nconv1
            iconv=-1
            do k=1,nconv1
             if(vconv(igeom,j).eq.vconv1(k))iconv=k
            enddo 
 	    ioff1=nconv1*(ipath-1)+iconv
            yn1(ioff+j)=yn1(ioff+j)+xgeom*calcout(ioff1)
           enddo

          else

           do j=1,nconv1
            iconv=-1
            do k=1,nconv1
             if(vconv(igeom,j).eq.vconv1(k))iconv=k
            enddo 
            ipath=4
   	    ioff1=nconv1*(ipath-1)+iconv
            yn1(ioff+j)=yn1(ioff+j)+xgeom*calcout(ioff1)
            ipath=5
 	    ioff1=nconv1*(ipath-1)+iconv
            yn1(ioff+nconv1+j)=yn1(ioff+nconv1+j)+
     1		xgeom*calcout(ioff1)
           enddo

          endif


111      continue

         do j=1,nconv1

           if(iextra.eq.1.and.jpre.gt.0)then
            kk(ioff+j,jpre) = kk(ioff+j,jpre) +
     1         (yn1(ioff+j)-yn(ioff+j))/delp  
            if(icread.eq.1)then
             kk(ioff+nconv1+j,jpre) = kk(ioff+nconv1+j,jpre) +
     1         (yn1(ioff+nconv1+j)-yn(ioff+nconv1+j))/delp  
            endif
           endif

           if(iextra.eq.2.and.jtan.gt.0)then
C           Assume change in tangent height pressure of 1km.
            kk(ioff+j,jtan) = kk(ioff+j,jtan) + yn1(ioff+j)-yn(ioff+j)
            if(icread.eq.1)then
             kk(ioff+nconv1+j,jtan) = kk(ioff+nconv1+j,jtan) + 
     1		yn1(ioff+nconv1+j)-yn(ioff+nconv1+j)
            endif
           endif

           if(iextra.eq.3.and.jlogg.gt.0)then
            kk(ioff+j,jlogg) = kk(ioff+j,jlogg) +
     1         (yn1(ioff+j)-yn(ioff+j))/dellg  
            if(icread.eq.1)then
             kk(ioff+nconv1+j,jlogg) = kk(ioff+nconv1+j,jlogg) +
     1         (yn1(ioff+nconv1+j)-yn(ioff+nconv1+j))/dellg  
            endif
           endif

         enddo
    

         if(iextra.eq.1.and.jpre.gt.0)then
          xn(jpre)=pressR
         endif

         if(iextra.eq.3.and.jlogg.gt.0)then
          xn(jlogg)=loggR
          mass2 = 1e-20*10**(xn(jlogg))*(radius2**2)/Grav
         endif

113     enddo

       endif

       if(icread.ne.1)then
        ioff = ioff + nconv1
       else
        ioff = ioff + 2*nconv1
       endif

100   continue

C      print*,'Done'

      return

      end
