      subroutine forwardnogXVenus(runname,ispace,iscat,fwhm,ngeom,nav,
     1 wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,
     2 lin,nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,
     3 jlogg,jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter,icheck)
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
C				3=single-scattering plane parallel
C				4=single-scattering spherical
C				5=internal flux calculations
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
      real xlat,xref,dx,Grav,xgeom
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
      real xn(mx),yn(my),kk(my,mx),ytmp(my),ystore(my)
      real yy(my),y1
      real vconv1(mconv),vwave1(mwave),xlon,yout(my)
      integer ny,jsurf,jalb,jtan,jpre,nem,nav(mgeom)
      integer nphi,ipath,iconv,k,jxsc
      integer nmu,isol,lowbc,nf,nf1,nx2,kiter
      real dist,galb,sol_ang,emiss_ang,z_ang,aphi,vv
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,logname
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      common /imiescat/imie1

      integer nvar,varident(mvar,3)
      real varparam(mvar,mparam)
      logical gasgiant
      real vem(maxsec),emissivity(maxsec)

      real stelrad,solwave(maxbin),solrad(maxbin)
      integer solnpt,iform,iread

      common /solardat/iread, iform, stelrad, solwave, solrad,  solnpt

      integer cellngas,cellid(maxgas),celliso(maxgas),icread
      real cellength,cellpress,celltemp,cellvmr(maxgas)
      common/celldat/icread,cellngas,cellid,celliso,cellvmr,cellength,
     1  cellpress,celltemp

      integer nlayv,nwave1v
      real vwavev(mconv),basepv(maxlay),basehv(maxlay)
      real solarv(mconv),radgroundv(mconv),galbv(mconv)
      real venera(mconv,maxlay),xx,trad(maxlay),hnow
      integer mfwhm,nfwhm
      parameter (mfwhm=1000)
      logical :: fwhmexist

C     Arrays
      real :: vfwhm(mfwhm),xfwhm(mfwhm)

C     jradf and jloggf are passed via the planrad common block   
      jradf=jrad
      jloggf=jlogg

      FWHMEXIST=.FALSE.


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
       print*,'forwardnogXVenus.f cannot be used to calculate'
       print*,'A_plan/A_star'
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

      nwave1 = nwave(1)
      nconv1 = nconv(1)
      do 105 i=1,nconv1
        vconv1(i)=vconv(1,i)
105   continue
      do 106 i=1,nwave1
        vwave1(i)=vwave(1,i)
106   continue

      if(nwave1.gt.1)call sort(nwave1,vwave1)
      if(nconv1.gt.1)call sort(nconv1,vconv1)



      sol_ang = angles(1,1,1)
      emiss_ang = 0.
      aphi = angles(1,1,3)         
      nf = 0

C      print*,'Angles : ',sol_ang,emiss_ang,aphi
C      print*,'nf = ',nf
      xlat = flat(1,1)
      xlon = flon(1,1)
      xgeom = wgeom(1,1)


      if(kiter.ge.0)then
           nx2 = nx+1
      else
           nx2 = 1
      endif


      do 110 ix1=1,nx2

         ix = ix1-1

         print*,'forwardnogXVenus, ix,nx = ',ix,nx
         if(ix.gt.0)then
            xref = xn(ix)
            dx = 0.05*xref
            if(dx.eq.0)dx = 0.1
            xn(ix)=xn(ix)+dx
         endif
         print*,'ix,xref,dx,xn(ix)',ix,xref,dx,xn(ix)
         if(jsurf.gt.0)then
           tsurf = xn(jsurf)
         endif

C        Check to see if this variable is unconstrained enough to bother
C        calculating its gradient.
         if(ix.gt.0.and.ifix(ix).eq.1)then
           print*,'Fix ',ix,xref
           xn(ix)=xref
           goto 110
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

          print*,'************** FORWARDNOGXVenus ***********'
          print*,'******** INORMAL = ',INORMAL
          print*,'******** ITYPE = ',ITYPE


C         Set up all files for a direct cirsrad run
          print*,'calling gsetradV'
          call gsetradV(runname,iscat,nmu,mu,wtmu,isol,dist,
     1     lowbc,galb,nf,nconv1,vconv1,fwhm,ispace,gasgiant,
     2     layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,xlon,
     3     lin,nvar,varident,varparam,nx,xn,jalb,jxsc,jtan,jpre,tsurf,
     4     xmap)
          print*,'gsetradV called OK'

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
          print*,'forwardnogXVenus, icheck = ',icheck
C         Check also to see if surface temperature has gone negative
          if(tsurf.lt.0.0)icheck=1   
          print*,tsurf,icheck

          if(icheck.eq.1)then
           print*,'Profiles have gone awry. Abort, increase brakes and'
           print*,'try again'
           return
          endif
          
          print*,'forwardnogXVenus',runname,runname

          call CIRSrtf_wave(runname, dist, inormal, iray,fwhm, ispace,
     1     vwave1,nwave1,npath, output, vconv1, nconv1, itype,
     2     nem,vem,emissivity,tsurf, calcout)


C         Now need to read in internal radiances from the venera.dat file

          open(12,file='venera.dat',status='old',err=222)
           read(12,*)nlayv
           read(12,*)nwave1v
C           print*,'nlayv = ',nlayv
           if(nwave1v.ne.nwave1)then
            print*,'nwave1v <> nwave1'
            stop
           endif
           read(12,*)(vwavev(i),i=1,nwave1v)
           read(12,*)(basepv(i),i=1,nlayv)
           read(12,*)(basehv(i),i=1,nlayv)
           do 100 j=1,nwave1v
            read(12,*)xx,solarv(j),radgroundv(j),galbv(j)
            read(12,*)(venera(j,i),i=1,nlayv)
100        continue
222       close(12)


C         Now need to interpolate internal radiances in the venera.dat file 
C         to the requested heights in the .spx file.

          ipath=1
          ioff=0

          do 111 igeom=1,ngeom
           hnow = angles(igeom,1,2)
C           print*,'igeom, hnow = ',igeom,hnow
           do k=1,nwave1v
            do i=1,nlayv
             trad(i)=venera(k,i)
            enddo
            call verint(basehv,trad,nlayv,y1,hnow)
            yy(k)=y1
           enddo
C          Smooth as required
C           print*,nwave1,(yy(i),i=1,nwave1)
           call cirsconv(runname,fwhm,nwave1,vwave1,yy,nconv1,
     & vconv1,yout,FWHMEXIST,NFWHM,VFWHM,XFWHM)
C           print*,nconv1,(yout(i),i=1,nconv1)
           do j=1,nconv1
            ytmp(ioff+j)=yout(j)
           enddo
 
           if(ix.eq.0)then
            do j=1,nconv1
             yn(ioff+j)=ytmp(ioff+j)
            enddo
           else
            do j=1,nconv1
             kk(ioff+j,ix)= (ytmp(ioff+j)-yn(ioff+j))/dx
            enddo
           endif
           ioff=ioff+nconv1
111       continue

          if(ix.gt.0)xn(ix)=xref
         
110   continue

      close(ulog)


      return

      end
