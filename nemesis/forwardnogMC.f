      subroutine forwardnogMC(runname,ispace,iscat,fwhm,ngeom,nav,
     1 wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,
     2 lin,nvar,varident,varparam,jsurf,jalb,jtan,jpre,
     3 nx,xn,ny,yn,kk,kiter)
C     $Id:
C     **************************************************************
C     Subroutine to calculate a synthetic spectrum and KK-matrix using
C     FINITE DIFFERENCES. The routine is identical in operation to 
C     forwardavfovX.f but calculates the K-matrix using old-fashioned 
C     finite differences. Routine exists to check that forwardavfovX.f, 
C     using the internal gradients is operating correctly, and also for
C     scattering calculations.
C
C     Based on forwardnogX.f
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
      integer layint,inormal,iray,itype,nlayer,laytyp,iscat
      integer nwave(mgeom),ix,ix1,iav,nwave1,iptf
      real vwave(mgeom,mwave),interpem
      real calcout(maxout3),fwhm,planck_wave,output(maxout3)
      real gradients(maxout4)
      integer nx,nconv(mgeom),npath,ioff1,ioff2,nconv1
      real vconv(mgeom,mconv),wgeom(mgeom,mav),flat(mgeom,mav)
      real layht,tsurf,esurf,angles(mgeom,mav,3)
      real xn(mx),yn(my),kk(my,mx),ytmp(my),ystore(my)
      real vconv1(mconv),vwave1(mwave)
      integer ny,jsurf,jalb,jtan,jpre,nem,nav(mgeom)
      integer nphi,ipath,iconv,k
      integer nmu,isol,lowbc,nf,nf1,nx2,kiter
      real dist,galb,sol_ang,emiss_ang,z_ang,aphi,vv
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,logname
      real xmap(maxv,maxgas+2+maxcon,maxpro)

      integer nvar,varident(mvar,3)
      real varparam(mvar,mparam)
      logical gasgiant
      real vem(maxsec),emissivity(maxsec)

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

      do 100 igeom=1,ngeom
       print*,'ForwardnogMC. Spectrum ',igeom,' of ',ngeom
       
       nwave1 = nwave(igeom)
       nconv1 = nconv(igeom)
       do 105 i=1,nconv1
        vconv1(i)=vconv(igeom,i)
105    continue
       do 106 i=1,nwave1
        vwave1(i)=vwave(igeom,i)
106    continue

       if(nwave1.gt.1)call sort(nwave1,vwave1)
       if(nconv1.gt.1)call sort(nconv1,vconv1)


        do 110 iav=1,nav(igeom)

         sol_ang = angles(igeom,iav,1)
         emiss_ang = angles(igeom,iav,2)
         aphi = angles(igeom,iav,3)
 
         if(sol_ang.lt.emiss_ang) then
           z_ang = sol_ang
         else
           z_ang = emiss_ang
         endif

         print*,'Angles : ',sol_ang,emiss_ang,aphi
         xlat = flat(igeom,iav)

         if(kiter.ge.0)then
           nx2 = nx+1
         else
           nx2 = 1
         endif

         do 111 ix1=1,nx2

          ix = ix1-1

          print*,'forwardnocgMC, ix,nx = ',ix,nx

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
          endif

         
C         Set up parameters for scattering cirsrad run.

          CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF)


          print*,'************** FORWARDNOGMC ***********'
          print*,'******** INORMAL = ',INORMAL


C         Set up all files for a direct cirsrad run
          call gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,
     1     lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,
     2     layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,lin,
     3     nvar,varident,varparam,nx,xn,jalb,jtan,jpre,tsurf,xmap)

C         MonteCarlo code needs surface emissivity to be read in here.
          call readsurfem(runname,nem,vem,emissivity)


          print*,'runname = ',runname
          
          call CIRSrtf_waveMC(runname, dist, inormal, iray, fwhm, 
     1     ispace, vwave1,nwave1,sol_ang, emiss_ang, aphi, npath,
     2     output,vconv1, nconv1, nem,vem,emissivity,tsurf,
     3     calcout)


          ipath=1
          do j=1,nconv1
             iconv=-1
             do k=1,nconv1
              if(vconv(igeom,j).eq.vconv1(k))iconv=k
             enddo
             if(iconv.lt.0)then
              print*,'Error in forwardnogMC iconv < 0'
              stop
             endif
             ioff1=nconv1*(ipath-1)+iconv
             ytmp(ioff+j)=calcout(ioff1)
          enddo

          if(ix.eq.0)then
           do j=1,nconv1 
            yn(ioff+j)=yn(ioff+j)+wgeom(igeom,iav)*ytmp(ioff+j)
            ystore(ioff+j)=ytmp(ioff+j)
           enddo
          else
           do j=1,nconv1
            kk(ioff+j,ix)=kk(ioff+j,ix)+wgeom(igeom,iav)*
     1                      (ytmp(ioff+j) - ystore(ioff+j))/dx  
           enddo 
           xn(ix)=xref
          endif

111      continue
110     continue
     
       ioff = ioff + nconv1

100   continue

      close(ulog)

      return

      end
