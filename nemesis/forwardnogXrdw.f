      subroutine forwardnogXrdw(runname,ispace,iscat,fwhm,ngeom,nav,
     1 wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,
     2 lin,nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3 jlogg,jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter,icheck,rdwindicesi,
     4 kkold,nconvi,vconvi,ynold)
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
C	iscat		integer 0=non-scattering
C				1=plane-parallel scattering
C				2=non-plane limb/near-limb scattering
C				3=single-scattering
C       fwhm            real    Desired FWHM of final spectrum
C       ngeom           integer Number of observation geometries included.
C       nav(mgeom)      integer         Number of synthetic spectra required
C                                       to simulate each FOV-averaged
C                                       measurement spectrum.
C       wgeom(mgeom,mav)real     Integration weights to use
C       flat(mgeom,mav)  real    Integration point latitudes
C       flon(mgeom,mav)  real    Integration point longitudes
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
C       jfrac		integer position of profile fraction element in
C                               xn (if included)
C       RADIUS		real    Planetary radius at 0km altitude
C       nx              integer Number of elements in state vector
C       xn(mx)          real	State vector
C	ifix(mx)	integer Vector showing which elements we need
C				 gradients for
C       ny      	integer Number of elements in measured spectra array
C	kiter		integer Number of iterations of Nemesis
C
C     Output variables
C       yn(my)          real    Synthetic radiances
C       kk(my,mx)       real    dR/dx matrix
C	icheck		integer	Check to see if temperature, vmr or dust has
C				gone negative
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
      include '../radtran/includes/planrad.f'
      include 'arraylen.f'
      real xlat,xref,dx,Grav
      parameter (Grav=6.672E-11)
      integer layint,inormal,iray,itype,nlayer,laytyp,iscat
      integer nwave(mgeom),ix,ix1,iav,nwave1,iptf,jrad,j1
      real vwave(mgeom,mwave),interpem,RADIUS
      real calcout(maxout3),fwhm,planck_wave,output(maxout3)
      real gradients(maxout4),pi
      parameter (pi=3.1415927)
      integer check_profile,icheck,imie,imie1,jlogg,ifix(mx)
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1,jfrac
      real vconv(mgeom,mconv),wgeom(mgeom,mav),flat(mgeom,mav)
      real layht,tsurf,esurf,angles(mgeom,mav,3),flon(mgeom,mav)
      real xn(mx),yn(my),kk(my,mx),ytmp(my),ystore(my),ynold(my)
      real vconv1(mconv),vwave1(mwave),xlon
      integer ny,jsurf,jalb,jtan,jpre,nem,nav(mgeom)
      integer nphi,ipath,iconv,k,jxsc
      integer nmu,isol,lowbc,nf,nf1,nx2,kiter
      real dist,galb,sol_ang,emiss_ang,z_ang,aphi,vv
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,logname
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      common /imiescat/imie1

      integer nvar,varident(mvar,3),rdwindicesi(mgeom,mconv)
      real varparam(mvar,mparam)
      logical gasgiant
      real vem(maxsec),emissivity(maxsec)

      real stelrad,solwave(maxbin),solrad(maxbin)
      integer solnpt,iform,iread

      integer nwavei1,nconvi1,nconvc1,nwavei(mwave)
      integer rdwindicesc(mgeom,mconv),nconvi(mgeom),nconvc(mgeom)
      real vconvi(mgeom,mconv),vconvi1(mconv),kkold(my,mx)
      real vwavei(mgeom,mwave),vwavei1(mwave)

      common /solardat/iread, iform, stelrad, solwave, solrad,  solnpt

      integer cellngas,cellid(maxgas),celliso(maxgas),icread
      real cellength,cellpress,celltemp,cellvmr(maxgas)
      common/celldat/icread,cellngas,cellid,celliso,cellvmr,cellength,
     1  cellpress,celltemp
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet



C     jradf and jloggf are passed via the planrad common block   
      jradf=jrad
      jloggf=jlogg

C     Calculate complement of rdwindicesi (ie subscripts of wavelengths NOT to be put through forward model)
      do igeom=1,ngeom
       j=1
       k=1
       do i=1,nwave(igeom)
        if(i.ne.rdwindicesi(igeom,j))then
         rdwindicesc(igeom,k) = i
         k=k+1
        else
         j=j+1
        endif
       enddo
       nconvc(igeom)=k-1
      enddo

C      print*,'-----------------'
C      print*,runname
C      print*,ispace,iscat,fwhm,ngeom
C      print*,(nav(i),i=1,ngeom)
C      do j=1,ngeom
C       print*,(wgeom(j,i),i=1,nav(j))
C       print*,(flat(j,i),i=1,nav(j))
C       do i=1,nav(j)
C        print*,j,i,(angles(j,i,k),k=1,3)  
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



      if(iform.eq.2) then
       print*,'forwardnogxrdw.f cannot be used '
       print*, 'to calculate A_plan/A_star'
       stop
      endif

      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf1,nphi,layht,tsurf,nlayer,laytyp,layint)



      call file(runname,logname,'log')

      open(ulog,file=logname,status='unknown')

C     Initialise arrays
      do i=1,my
       yn(i)=0.0
       do j=1,mx
        kk(i,j)=0.0
       enddo
      enddo

      call setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf1,nphi,layht,tsurf,nlayer,laytyp,layint)

      if(jsurf.gt.0)then
       tsurf = xn(jsurf)
      endif

      ioff = 0
      ioff2 = 0

      do 100 igeom=1,ngeom
       if(idiag.gt.0)print*,'ForwardnogX. Spectrum ',igeom,' of ',ngeom
       
       nwave1 = nwave(igeom)
       nconv1 = nconv(igeom)
       nwavei(igeom) = nconvi(igeom)
       nwavei1 = nwavei(igeom)
       nconvi1 = nconvi(igeom)
       nconvc1 = nconvc(igeom)
       vwavei(igeom,1:nwave(igeom))=vconvi(igeom,1:nconv(igeom))
       do 105 i=1,nconv1
        vconv1(i)=vconv(igeom,i)
        vconvi1(i)=vconvi(igeom,i)
105    continue
       do 106 i=1,nwave1
        vwave1(i)=vwave(igeom,i)
        vwavei1(i)=vwavei(igeom,i)
106    continue

       if(nwave1.gt.1)call sort(nwave1,vwave1)
       if(nconv1.gt.1)call sort(nconv1,vconv1)
       if(nwavei1.gt.1)call sort(nwavei1,vwavei1)
       if(nconvi1.gt.1)call sort(nconvi1,vconvi1)


        do 110 iav=1,nav(igeom)

         sol_ang = angles(igeom,iav,1)
         emiss_ang = angles(igeom,iav,2)
         aphi = angles(igeom,iav,3)
         
         if(sol_ang.lt.emiss_ang) then
           z_ang = sol_ang
         else
           z_ang = emiss_ang
         endif

C        New bit to increase number of Fourier components depending on
C        miniumum zenith angle
         if(iscat.eq.1.and.z_ang.ge.0.0)then
           nf = int(30*z_ang/90.0)
C            nf=9
C            nf=0
C            nf=20
         else
           nf=nf1
         endif

         if(idiag.gt.0)print*,'Angles : ',sol_ang,emiss_ang,aphi
         if(idiag.gt.0)print*,'nf = ',nf
         xlat = flat(igeom,iav)
         xlon = flon(igeom,iav)

         if(kiter.ge.0)then
           nx2 = nx+1
         else
           nx2 = 1
         endif


         do 111 ix1=1,nx2

          ix = ix1-1

          if(idiag.gt.0)print*,'forwardnogX, ix,nx = ',ix,nx
          if(ix.gt.0)then
            xref = xn(ix)
            dx = 0.05*xref
            if(dx.eq.0)dx = 0.1
            if(ix.eq.jrad)dx=10.
            if(ix.eq.jtan)then
             if(emiss_ang.lt.0)then
               dx=1.0
             else
               goto 111
             endif
            endif                             
            xn(ix)=xn(ix)+dx
          endif
C          if(idiag.gt.0)print*,'ix,xref,dx,xn(ix)',ix,xref,dx,xn(ix)
          if(jsurf.gt.0)then
           tsurf = xn(jsurf)
          endif

C         Check to see if this variable is unconstrained enough to bother
C         calculating its gradient.
          if(ix.gt.0.and.ifix(ix).eq.1)then
           xn(ix)=xref
           goto 111
          endif

C        If we're retrieving planet radius then add correction to reference
C        radius
C        N.B.radius2 is passed via the planrad common block.
         if(jrad.gt.0)then
          radius2 = xn(jrad) + radius
         else      
          radius2 = radius
         endif

C        If we're retrieving surface gravity then modify the planet mass
C        N.B. mass2 is passed via the planrad common block. Assume xn(jlogg)
C        holds log_10(surface gravity in cm/s^2). Need factor of 1e-20 to convert
C        mass to units of 1e24 kg.       
         if(jlogg.gt.0)then
          mass2 = 1e-20*10**(xn(jlogg))*(radius2**2)/Grav
         endif

         
C        Set up parameters for scattering cirsrad run.

         CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,INH3,
     1    IPTF,IMIE, iuvscat)
         IMIE1=IMIE
          itype=11			! scloud11wave

          if(idiag.gt.0)print*,'************** FORWARDNOGX ***********'
          if(idiag.gt.0)print*,'******** INORMAL = ',INORMAL
          if(idiag.gt.0)print*,'******** ITYPE = ',ITYPE


C         Set up all files for a direct cirsrad run
          if(idiag.gt.0)print*,'calling gsetrad'
          call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,
     1     lowbc,galb,nf,nconvi1,vconvi1,fwhm,ispace,gasgiant,
     2     layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,xlon,
     3     lin,nvar,varident,varparam,nx,xn,jalb,jxsc,jtan,jpre,tsurf,
     4     xmap)
          if(idiag.gt.0)print*,'gsetrad called OK'

C         If planet is not a gas giant and observation is not at limb then
C         we need to read in the surface emissivity spectrum.
C
C         Read in emissivity file  
          if(.not.gasgiant.and.emiss_ang.ge.0)then
           call readsurfem(runname,nem,vem,emissivity)
          else
           nem=2
           vem(1)=-100.0
           vem(2)=1e7
           emissivity(1)=1.0
           emissivity(2)=1.0
          endif

C         Check to see if any temperatures or vmrs have gone
C         negative and if so abort
          icheck = check_profile(runname)
          if(idiag.gt.0)print*,'forwardnogX, icheck = ',icheck
C         Check also to see if surface temperature has gone negative
          if(tsurf.lt.0.0)icheck=1   
          if(idiag.gt.0)print*,tsurf,icheck

          if(icheck.eq.1)then
           print*,'Profiles have gone awry. Abort, increase brakes and'
           print*,'try again'
           return
          endif
          
          if(idiag.gt.0)print*,'forwardnogX',runname,runname

          call CIRSrtf_wave(runname, dist, inormal, iray,fwhm, ispace,
     1     vwavei1,nwavei1,npath, output, vconvi1, nconvi1, itype,
     2     nem,vem,emissivity,tsurf, calcout)

          if(idiag.gt.0)print*,'Npath = ',npath

C         Unless an SCR calculation, first path is assumed to be thermal emission

          if(icread.ne.1)then       
           ipath=1
           do j=1,nconvi1
             iconv=-1
             do k=1,nconvi1
              if(vconvi(igeom,j).eq.vconvi1(k))iconv=k
             enddo
             if(iconv.lt.0)then
              print*,'Error in forwardnogXrdw iconv < 0'
              stop
             endif
             ioff1 = ipath + (iconv-1)*npath
             ytmp(ioff+rdwindicesi(igeom,j))=calcout(ioff1)
           enddo

           if(ix.eq.0)then
            do j=1,nconvi1
             k = rdwindicesi(igeom,j)
             yn(ioff+k)=yn(ioff+k)+wgeom(igeom,iav)*ytmp(ioff+k)
             ystore(ioff+k)=ytmp(ioff+k)
            enddo
C           Read old wavelengths (ones modelled in forwardnogX - C) back into yn
            do j=1,nconvc1
             k = rdwindicesc(igeom,j)
             yn(ioff+k)=ynold(ioff2+j)
            enddo
           else
            do j=1,nconvi1
             k = rdwindicesi(igeom,j)
             kk(ioff+k,ix)=kk(ioff+k,ix)+wgeom(igeom,iav)*
     1                      (ytmp(ioff+k) - ystore(ioff+k))/dx  
            enddo
            do j=1,nconvc1
C           Read old gradients (ones modelled in forwardnogX - C) back into Jacobian
             k = rdwindicesc(igeom,j)
             kk(ioff+k,ix)=kkold(ioff2+j,ix)
            enddo 
            xn(ix)=xref
            if(ix.eq.jlogg)then
             mass2 = 1e-20*10**(xn(jlogg))*(radius2**2)/Grav
            endif
           endif

          else

           print*,'Error forwardnogXrdw: '
           print*,'Reduced wavelength scheme not yet set up for '
           print*,'SCR cell retrievals. Please replace this message'
           print*,'in forwardnogXrdw.f with corresponding code for '
           print*,'SCR cell retrievals. Refer also to forwardnogX.f'
           print*,'for reference.'
           stop

          endif

111      continue

110     continue

       if(icread.ne.1)then
        ioff = ioff + nconv1
        ioff2 = ioff2 + nconvc1
       else
        ioff = ioff + 2*nconv1
       endif

100   continue

      close(ulog)


      return

      end
