      subroutine forwardnogL(runname,ispace,fwhm,ngeom,nav,
     1 wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2 nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,occult,
     3 ionpeel,nx,xn,ny,yn,kk)
C     $Id:
C     **************************************************************
C     Subroutine to calculate an FOV-averaged limb spectra using CIRSRADG,
C     but the gradients KK are calculated using numerical integration.
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
C	flon(mgeom,mav)	real	Integration point longitudes
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
C       jxsc            integer Position of x-section spectrum in
C                               xn (if included)
C       jtan            integer Position of tangent height correction in
C                               xn (if included)
C       jpre            integer Position of tangent pressure in
C                               xn (if included)
C       occult          integer Solar occultation flag
C       ionpeel         integer Onion-peeling method flag (just one tangent height)
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
      integer i,j,k,lin,ispace,iav,jpath,ix
      integer ngeom,ioff,igeom,nlay1
      real interpem,dv
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/gascom.f'
      include 'arraylen.f'
      real flon(mgeom,mav),xlon
      real xlat,planck_wave,planckg_wave,Bg,height(maxlay),htan
      real wgeom(mgeom,mav),flat(mgeom,mav),fh,baseh1(maxlay)
      integer layint,inormal,iray,itype,nlayer,laytyp
      integer iptf,imie,imie1,occult,iconv,idummy,ionpeel,jlevlo,jlevhi
      integer nwave(mgeom),jsurf,nem,nav(mgeom),nwave1
      real vwave(mgeom,mwave),angles(mgeom,mav,3),vwave1(mwave)
      real calcout(maxout3),fwhm,calcoutL(maxout3),htanlo,htanhi
      real calcout1(maxout3),gradients1(maxout4),xn1(mx)
      real gradients(maxout4),vv,gradientsL(maxout4),pi,AU
      parameter (pi=3.1415927, AU=1.49597870e8)
      real solradius,solomeg,solar,xsol
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      real vconv(mgeom,mconv),vconv1(mconv)
      real layht,tsurf,esurf,gradtsurf(maxout3),pressR,delp
      real xn(mx),yn(my),kk(my,mx),yn1(my),xref,dx
      integer ny,iscat,jalb,jxsc,jtan,jpre
      integer nphi,ipath,lhay
      integer nmu,isol,lowbc,nf
      real dist,galb,sol_ang,emiss_ang,aphi
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,header
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      real hcorr,hcorrx
      common /imiescat/ imie1
      integer nvar,varident(mvar,3),ivar
      real vem(maxsec),emissivity(maxsec)
      real varparam(mvar,mparam)
      logical gasgiant,layexist

      real stelrad,solwave(maxbin),solrad(maxbin)
      integer solnpt,iform,iread

      common /solardat/iread, iform, stelrad, solwave, solrad,  solnpt


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
       print*,'Warning from forwardnogL'
       print*,'Can not currently do surface albedo retrievals with'
       print*,'gradient code'   
       stop
      endif

      print*,'forwardnogL. jpre = ',jpre

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
         print*,'forwardnogL: All spectra must be of same length'
         print*,'igeom,nconv(1),nconv(igeom)',igeom,nconv1,nconv(igeom)
         stop
        endif
        do i=1,nconv1
         dv = 100.0*abs(vconv(igeom,i)-vconv1(i))/vconv1(i)
         if(dv.gt.1.0)then
          print*,'forwardnogL: wavenumbers do not match'
          print*,igeom,vconv1(i),vconv(igeom,i),dv
          stop
         endif
        enddo
       endif
      enddo

      xlat = flat(1,1)   
      xlon = flon(1,1)   

      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf,nphi,layht,tsurf,nlayer,laytyp,layint)



C     Calculating the atmospheric path to calculate for the onion-peeling 
      jlevlo=1
      jlevhi=1
      if(ionpeel.eq.1)then
C     It is assumed that the base level of the atmospheric layers computed
C     by RADTRAN is read from the height.lay file

       if(laytyp.ne.5)then
         print*,'error in forwardavfovL - when using the onion-peeling'
         print*,'scheme it is assumed that the base level of each layer'
         print*,'is read from the height.lay file (laytyp = 5)'
       endif

       inquire(file='height.lay',exist=layexist)
       if(layexist)then
         lhay=12
1        FORMAT(A)
         open(lhay,file='height.lay',status='old')
C        Read in one line of header
         read(lhay,1)HEADER
         read(lhay,*)nlay1
         DO 107 i=1,nlay1
          READ(12,*)baseh1(i)
107      CONTINUE
         close(lhay)
       else
         print*,'error in forwardavfovL - when using the onion-peeling'
         print*,'scheme it is assumed that the base level of each layer'
         print*,'is read from the height.lay file, but it doesnt exist'
       endif

       do i=1,nlay1
        print*,baseh1(i)
       enddo

       htanlo = angles(1,1,1)
       htanhi = angles(ngeom,1,1)

       do i=1,nlay1
        if(baseh1(i).le.htanlo.and.baseh1(i+1).gt.htanlo)jlevlo=i
        if(baseh1(i).le.htanhi.and.baseh1(i+1).gt.htanhi)jlevhi=i
       enddo
      endif



      if(jsurf.gt.0)then
       tsurf = xn(jsurf)
      endif

      iscat=0
C     Set up all files for a direct cirsrad run of limb spectra
      call gsetradL(runname,nconv1,vconv1,fwhm,ispace,iscat,
     1    gasgiant,layht,nlayer,laytyp,layint,xlat,xlon,lin,hcorrx,
     2    nvar,varident,varparam,nx,xn,jpre,tsurf,occult,ionpeel,
     3    jlevlo,jlevhi,xmap)


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

      CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,INH3,
     1 IPTF,IMIE, iuvscat)
      IMIE1=IMIE

      itype=12			! scloud12. not used here

      print*,'hcorrx = ',hcorrx
 
      call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace, 
     1  vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf, 
     2  nx, xmap, vconv1, nconv1, npath, calcoutL, gradientsL,iscat)

C     Read in base heights from '.drv' file
      call readdrvh(runname,height)

      print*,'hcorrx = ',hcorrx

      ioff = 0


c     MAIN LOOP - CALCULATING SPECTRUM GIVEN AT EACH GEOMETRY

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
         

C          In the onion peeling method it is assumed that the tangent
C          heights coincide with the base altitude of each layer
           if(ionpeel.eq.1)then
            jpath=igeom
           endif


           if(occult.eq.2)then
C           just calculate transmission of path. Assumes no thermal emission 
C           into beam and no other nadir paths.
            do ipath=jpath,jpath+1
             fh=1.0-fh
             do j=1,nconv1
              ioff1=nconv1*(ipath-1)+j
              yn(ioff+j)=yn(ioff+j)+wgeom(igeom,iav)*fh*calcoutL(ioff1)
             enddo             
            enddo

           else

C           if OCCULT > 0, but <> 2, then first nlayer paths are
C           thermal emission and 2nd nlayer paths are transmission.
C           if OCCULT = 0, then there are nlayer thermal emission paths only

C           Add in thermal emission paths.
            do ipath=jpath,jpath+1
             fh=1.0-fh
             do j=1,nconv1
 	      ioff1=nconv1*(ipath-1)+j
              yn(ioff+j)=yn(ioff+j)+wgeom(igeom,iav)*fh*calcoutL(ioff1)
             enddo
            enddo

            if(occult.eq.1.or.occult.eq.3)then
C            Add solar occultation radiance to thermal emission

             solradius = stelrad/(dist*AU)
             solomeg = pi*solradius**2

             do ipath=jpath+nlayer,jpath+nlayer+1
              fh=1.0-fh
              do j=1,nconv1
C              Get Solar irradiance
               CALL get_solar_wave(vconv1(j),dist,solar)
               xsol =  solar/solomeg
  	       ioff1=nconv1*(ipath-1)+j
               yn(ioff+j)=yn(ioff+j)+
     1          wgeom(igeom,iav)*fh*calcoutL(ioff1)*xsol             
              enddo
             enddo

             if(occult.eq.3)then
C             divide thermal emission + solar occultation by solar radiance at top of atmosphere
              do j=1,nconv1
                CALL get_solar_wave(vconv1(j),dist,solar)
                xsol = solar/solomeg
                yn(ioff+j)=yn(ioff+j)/xsol
              enddo
             endif

            endif

           endif 

         else

           print*,'Calculating new nadir-spectra'
           iscat = 0
           call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,lowbc,
     1      galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,layht,
     2      nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,xlon,lin,
     3      nvar,varident,varparam,nx,xn,jalb,jxsc,jtan,jpre,tsurf,
     4      xmap)
      

           call CIRSrtfg_wave(runname,dist,inormal,iray,fwhm,ispace,
     1      vwave1,nwave1,itype, nem, vem, emissivity, tsurf, 
     2      gradtsurf, nx, xmap, vconv1, nconv1, npath, calcout, 
     3 	    gradients,iscat)

C          First path is assumed to be thermal emission
           
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
         endif

110    continue

       ioff = ioff + nconv1

100   continue



c     MAIN LOOP - CALCULATING NEW SPECTRUM AFTER PERTURBATION IN STATE VECTOR
 
      do i=1,nx
       xn1(i)=xn(i)
      enddo

      do 115 ix=1,nx

c       Initialise yn1 array
        do i=1,my
         yn1(i)=0.0
        enddo

        print*,'ForwardnogL :: Calculating gradients numerically'
        print*,'ix,xn(ix) = ',ix,xn(ix)

        xref = xn(ix)
        dx = 0.05*xref
        if(dx.eq.0)dx = 0.1

        if(ix.eq.jpre)then
          if(emiss_ang.lt.0)then
            dx=0.01*xref
          else
            goto 115
          endif
        endif

        if(ix.eq.jtan)then
          if(emiss_ang.lt.0)then
            dx=1.0
          else
            goto 115
          endif
        endif

        xn(ix)=xn(ix)+dx
        do i = 1,nx
         print*,i,xn1(i),xn(i),xn(i)-xn1(i)
        enddo

C       Set up all files to recalculate limb spectra after perturbation
        call gsetradL(runname,nconv1,vconv1,fwhm,ispace,iscat,
     1    gasgiant,layht,nlayer,laytyp,layint,xlat,xlon,lin,hcorrx,
     2    nvar,varident,varparam,nx,xn,jpre,tsurf,occult,ionpeel,
     3    jlevlo,jlevhi,xmap)

        call CIRSrtfg_wave(runname, dist, inormal, iray, fwhm, ispace,
     1   vwave1,nwave1,itype, nem, vem, emissivity, tsurf, gradtsurf,
     2   nx, xmap, vconv1, nconv1, npath, calcout1,gradients1,iscat) 

        ioff=0
        do 200 igeom=1,ngeom
         print*,'ForwardnogL. Spectrum ',
     1    igeom,' of ',ngeom
         print*,'Nav = ',nav(igeom)

         do 113 iav = 1,nav(igeom)
          sol_ang = angles(igeom,iav,1)
          emiss_ang = angles(igeom,iav,2)
          aphi = angles(igeom,iav,3)

          print*,'Iav = ',iav
          print*,'Angles : ',sol_ang,emiss_ang,aphi

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
            print*,'Error in forwardnogL: tangent height is'
            print*,'out of range'
            print*,nlayer,height(1),height(nlayer),htan
           endif
           fh = (htan-height(jpath))/(height(jpath+1)-height(jpath))
       

C          In the onion peeling method it is assumed that the tangent
C          heights coincide with the base altitude of each layer
           if(ionpeel.eq.1)then
            jpath=igeom
           endif

 
c          Now we calculate the actual spectra at the given tangent height

           if(occult.eq.2)then
C           just calculate transmission of path. Assumes no thermal emission 
C           into line of sight. First nlayers are transmission.

            print*,'occult = ',occult
            print*,'jpath,fh,nconv1 = ',jpath,fh,nconv1
            do ipath=jpath,jpath+1
             fh=1.0-fh
             do j=1,nconv1
              ioff1=nconv1*(ipath-1)+j
              yn1(ioff+j)=yn1(ioff+j)+
     1          wgeom(igeom,iav)*fh*calcout1(ioff1)
             enddo
            enddo

           else

C           First nlayer paths are thermal emission
            do ipath=jpath,jpath+1
             fh=1.0-fh
             do j=1,nconv1
              ioff1=nconv1*(ipath-1)+j
              yn1(ioff+j)=yn1(ioff+j)+
     1          wgeom(igeom,iav)*fh*calcout1(ioff1)
             enddo
            enddo

            if(occult.eq.1.or.occult.eq.3)then

             solradius = stelrad/(dist*AU)
             solomeg = pi*solradius**2

             do ipath=jpath+nlayer,jpath+nlayer+1
              fh=1.0-fh
              do j=1,nconv1
C              Get Solar irradiance
               CALL get_solar_wave(vconv1(j),dist,solar)
               xsol = solar/solomeg

               ioff1=nconv1*(ipath-1)+j
               yn1(ioff+j)=yn1(ioff+j)+
     1          wgeom(igeom,iav)*fh*calcout1(ioff1)*xsol
              enddo
             enddo

             if(occult.eq.3)then
              do j=1,nconv1
                CALL get_solar_wave(vconv1(j),dist,solar)
                xsol = solar/solomeg
                yn1(ioff+j)=yn1(ioff+j)/xsol
              enddo
             endif

            endif

           endif
          endif

113      continue     !loop in nav

 
         do j=1,nconv1
           kk(ioff+j,ix) = (yn1(ioff+j)-yn(ioff+j))/dx
         enddo

         ioff = ioff + nconv1

200     continue !loop in ngeom

        xn(ix) = xref

115   continue   !loop in state vector ix

      return

      end

