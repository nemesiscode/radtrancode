      PROGRAM generatespx
C     $Id:
C     ******************************************************************
C
C     CIRS retrieval code utilising correlated-k, thermal emission 
C     fast gradient radiative transfer model CIRSRADG. Extension of
C     Nemesis to retrieve from a number of locations in a row.  
C
C     CIRSRADG cannot currently deal with scattering calculations so this
C     gas to be done with CIRSRAD if a scattering calculation is required.
C     It is intended to upgrade CIRSRADG later.  
C
C     Minimisation is achieved using a modified non-linear estimation
C     which uses a Marquardt-Levenburg type brake.
C
C     Code can simultaneously retrieve to several measurements of the same
C     area at different viewing geometries.
C
C     Code can also average spectra over range of viewing angles.
C
C     Pat Irwin	        Modified from NIMS retrieval code 21/3/00
C			Updated	4/4/01
C			Updated for continuous vmr profiles 7/10/03
C			Updated for FOV-averaging 9/2/04
C
C     ******************************************************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

C     New compiler time
      real tot_time
      real rate
      integer c1,c2,cr,time1,time2,cm
     
      
C     TIME: Temporary variable returned by GETTIME containing the system time.
C     TIME1: System time at the beginning of program execution.
C     TIME2: System time at the end of program execution.

      character*100 buffer,ename,solfile,solname
      integer i,j,iscat,ica,k,lspec,lout,ispec,nspec,nspecx,ioff
      real xn(mx),se(my),err1(mx),woff,xdiff
      real fwhm,xlat,xlon,st(mx,mx),varparam(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xlatx,varparamx(mvar,mparam)
      real stx(mx,mx),xlonx,RADIUS
      integer varident(mvar,3),varidentx(mvar,3),iscat1,iplanet
      integer npro,nvmr,ispace,nav(mgeom),lraw,nprox,lpre
      integer lx(mx),iprfcheck,ifix(mx)
      character*100 runname
      integer ngeom, nwave(mgeom), nconv(mgeom), nx, ny, jsurf, jsurfx
      integer np,lin1,ioffx,ivarx,npx
      integer nwaveT(mgeom), nconvT(mgeom),jxsc
      integer ngas,ncont,nvar,nvarx,lin,nxx,igeom,nconv1,nwave1,jalb
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real vwaveT(mgeom,mwave),vconvT(mgeom,mconv)
      real xa(mx),rerr(mgeom,mconv),sa(mx,mx),y(my),yn(my)
      real xnx(mx),kk(my,mx),xerr
      real wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real vwave1(mwave),vconv1(mconv)
      double precision aa(mx,mx),dd(mx,my)
      real vkstart,vkend,vkstep
      integer idump,kiter,jtan,jtanx,jalbx,jpre,jprex,idum,lvec
      integer ivar,npvar,jrad,jlogg,iform,jradx,jloggx,jxscx
      integer jfrac,jfracx
C     ********** Scattering variables **********************
      real xwave(maxsec),xf(maxcon,maxsec),xg1(maxcon,maxsec)
      real xg2(maxcon,maxsec)
      real tnco,twave,frac,tico
      real phlimit,kkcor(mx,mx)
      logical gasgiant,solexist,percbool
      COMMON /hgphas/xwave,xf,xg1,xg2,tnco,twave,frac,tico
      COMMON /scatdump/ idump

      include '../radtran/includes/ciacom.f'

      CHARACTER*100 ANAME
      REAL DNU
      INTEGER IPARA
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C     ******************************************************


C     *******************************************************
C     ****************   CODE *******************************
C     *******************************************************
     
      idiag=1

C     Read in reference gas information data
      CALL RESERVEGAS

C     ----------- Scattering phase function initialisation --------------
      xwave(1)=-1                       ! Reset to force read of hgphase*
C                                         files.
C     ------------ Scattering phase function initialisation -------------

C     New compiler time
      CALL system_clock(count_rate=cr)
      CALL system_clock(count_max=cm)
      rate = REAL(cr)
      call system_clock(time1)

      CALL prompt('Enter run name : ')
      READ(5,1)buffer
1     FORMAT(a)
      runname = buffer(1:36)

      if(idiag.gt.0)print*,'checking files'
C     Make sure input files are OK
      CALL checkfiles(runname)

      CALL readrefhead(runname,npro,nvmr,gasgiant)
      if(npro.gt.maxpro)then
       print*,'Error in Nemesis. npro > maxpro : ',npro,maxpro
       stop
      endif

C     Read start, end and step of tables
      CALL readkkhead(runname,vkstart,vkend,vkstep)

      CALL file(runname,runname,'inp')
      OPEN(32,file=runname,status='old')

C     Read in whether to calculate with wavenumbers(0) or wavelength(1)
C     Also read in whether scattering is required (iscat)
      READ(32,*)ispace,iscat
      CALL readrefiplan(runname,iplanet,xlat,RADIUS)

      kiter=-1
      woff=0.0

C     Read in total number of spectra to simulate
      READ(32,*)nspec

C     Read in random -ve seed number
      READ(32,*)idum
C     Read in name of forward modelling error file
      READ(32,1)ename

C     Read in lin identifier in case want to use retrieved T
      READ(32,*)lin

      READ(32,*)percbool

      CLOSE(32)

      iform=0
     


C     See if there is a solar or stellar reference spectrum and read in
C     if present.
      call file(runname,solfile,'sol')
      inquire(file=solfile,exist=solexist)
      if(solexist)then
         call opensol(solfile,solname)
         CALL init_solar_wave(ispace,solname)
      else
         if(iform.eq.1)then
          print*,'Error in Nemesis. Flux-ratio calculation defined'
          print*,'but no solar file exists'
          stop
         endif
      endif

C     Open output files
      lout=38
      open(lout,file='generatespx.spx',status='unknown')
      lvec = 36
      open(lvec,file='modvector.dat',status='unknown')

      if(lin.gt.0)then
       lpre=39
C      if previous retrieval to be considered,
C      open previous raw retrieval file (copied to .pre)
       CALL file(runname,runname,'pre')
       open(lpre,file=runname,status='old')
       read(lpre,*)nspecx
       if(nspec.gt.nspecx)then
        if(idiag.gt.0)print*,'.pre file does not contain enough'
        if(idiag.gt.0)print*,'retrievals'
        stop
       endif
      endif

C     Open spectra file
      lspec=37
      CALL file(runname,runname,'spx')
      if(idiag.gt.0)write(6,1)runname
      open(lspec,file=runname,status='old')
C     Read in sample measurement vector, obs. geometry and covariances
       call readnextspavX(lspec,iform,woff,xlat,xlon,ngeom,nav,ny,y,se,
     1  fwhm,nconv,vconv,angles,wgeom,flat,flon)
      close(lspec)

C     Calculate the tabulated wavelengths of c-k look up tables
      do igeom=1,ngeom
       nconv1=nconv(igeom)
       do j=1,nconv1
        vconv1(j)=vconv(igeom,j)
       enddo
       CALL wavesetb(runname,vkstart,vkend,vkstep,nconv1,vconv1,
     1  fwhm,nwave1,vwave1)
       do j=1,nwave1
        vwave(igeom,j)=vwave1(j)
       enddo
       nwave(igeom)=nwave1
      enddo

C     set up a priori of x and its covariance
      lin1=0
      CALL readapriori(runname,lin1,lpre,xlat,npro,nvar,varident,
     1  varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,jfrac,nx,xa,
     2  sa,lx)
	
      write(lvec,*)nvar,'   ! nvar'     ! Number of variable profiles

      do 10 ivar=1,nvar
          write(lvec,*)(varident(ivar,j),j=1,3)
          write(lvec,*)(varparam(ivar,j),j=1,mparam)
10    continue


C     Find all the calculation and convolution wavelengths and rank
C     in order

      call rankwave(ngeom,nwave,vwave,nconv,vconv,nwaveT,vwaveT,
     1 nconvT,vconvT)

C     Read in forward modelling errors
      call forwarderr(ename,ngeom,nconv,vconv,woff,rerr)

C     Add forward errors to measurement covariances
      k=0
      DO i=1,ngeom
       do j=1,nconv(i)
        k = k+1
        xerr=rerr(i,j)
        if(percbool.eqv..true.)xerr=rerr(i,j)*(y(j)/100) 
        if(iform.eq.3)xerr=xerr*1e-18
        se(k)=se(k)+xerr**2
       enddo
      ENDDO

      CALL FILE(runname,runname,'cia')

      OPEN(12,FILE=runname,STATUS='OLD')
       READ(12,1)ANAME
       READ(12,*) DNU
       READ(12,*) IPARA
      CLOSE(12)
      IREAD1=1
      IREAD2=1
      IF(IPARA.EQ.0)THEN
       ANAME1=ANAME
       DNU1=DNU
      ELSE
       ANAME2=ANAME
       DNU2=DNU
       IPARA2=IPARA
      ENDIF

      call setifix(xa,sa,nvar,varident,varparam,npro,ifix)

      do 2999 ispec=1,nspec

C      Pseudo-random number generator seems to repeat cycle every 20 iterations, 
C      therefore need to reseed it to continue to get pseudo-random numbers
       if((mod(ispec-1,20).eq.0).and.(ispec.gt.1))then
        idum = idum*2
       endif

       call rmodapriori(idum,npro,nvar,varident,varparam,
     1  nx,xa,sa,xn)

       if(lin.gt.0) then
        call readraw(lpre,xlatx,xlonx,nprox,nvarx,varidentx,
     1   varparamx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,jloggx,
     2   jfracx,nxx,xnx,stx)

        ioff=0
        do ivar=1,nvar
         np=1
         if(varident(ivar,1).le.100)then
           np=npvar(varident(ivar,3),npro,varparam(ivar,1))
         endif
         if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))
         if(varident(ivar,1).eq.887)np = int(varparam(ivar,1))
         if(varident(ivar,1).eq.444)then
          if(varparam(ivar,2).gt.0.0)then
           np = 2+int(varparam(ivar,1))
          else
           np = 3
          endif
         endif
         if(varident(ivar,1).eq.446)then
          if(varparam(ivar,2).gt.0.0)then
           np = 3+2*int(varparam(ivar,1))
          else
           np = 5
          endif
         endif
         if(varident(ivar,1).eq.445)np = 3+int(varparam(ivar,1))
         if(varident(ivar,1).eq.222)np = 8
         if(varident(ivar,1).eq.223)np = 9
         if(varident(ivar,1).eq.224)np = 9
         if(varident(ivar,1).eq.225)np = 11
         if(varident(ivar,1).eq.226)np = 8
         if(varident(ivar,1).eq.227)np = 7

         if(varident(ivar,1).eq.0.and.np.eq.npro)then
          ioffx=0
          do ivarx=1,nvarx
           npx=1
           if(varidentx(ivarx,1).le.100)then
             npx=npvar(varidentx(ivarx,3),npro,varparamx(ivarx,1))
           endif
           if(varidentx(ivarx,1).eq.888)npx = int(varparamx(ivarx,1))
           if(varidentx(ivarx,1).eq.887)npx = int(varparamx(ivarx,1))
           if(varidentx(ivarx,1).eq.444)then
            if(varparamx(ivarx,2).gt.0.0)then
             npx = 2+int(varparamx(ivarx,1))
            else
             npx = 3
            endif
           endif
           if(varidentx(ivarx,1).eq.446)then
            if(varparamx(ivarx,2).gt.0.0)then
             npx = 3+2*int(varparamx(ivarx,1))
            else
             npx = 5
            endif
           endif
           if(varidentx(ivarx,1).eq.445)npx = 3+int(varparamx(ivarx,1))

           if(varidentx(ivarx,1).eq.0.and.npx.eq.npro)then
            do i=1,npro
             xn(ioff+i)=xnx(ioffx+i)
            enddo
           endif
           ioffx=ioffx+npx
          enddo
         endif
         ioff=ioff+np
        enddo
       endif

       idump=0	! flag for diagnostic print dumps

       ica=0
	
       if(iscat.eq.0)then
        CALL forwardavfovX(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,
     3   jfrac,RADIUS,nx,xn,ny,yn,kk)
      elseif(iscat.eq.1)then 
       if(idiag.gt.0)print*,'Calling forwardnogX'
       CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,
     3   jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter,iprfcheck)
      elseif(iscat.eq.2)then
       CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1   vconvT,gasgiant,lin,nvar,varident,varparam,jsurf,jalb,
     2   jxsc,jtan,jpre,jrad,jlogg,jfrac,RADIUS,nx,xn)
       iscat1=1
       CALL forwardnogX(runname,ispace,iscat1,fwhm,ngeom,nav,
     1   wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,
     3   jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter,iprfcheck)
      else
       iscat1=1
       CALL forwardnogX(runname,ispace,iscat1,fwhm,ngeom,nav,
     1   wgeom,flat,flon,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,
     3   jfrac,RADIUS,nx,xn,ifix,ny,yn,kk,kiter,iprfcheck)
      endif

      call writenextspavX(lout,iform,idum,woff,xlat,xlon,ngeom,nav,
     1  ny,yn,se,fwhm,nconv,vconv,angles,wgeom,flat,flon)

      call writemvec(lvec,npro,nvar,varident,varparam,jsurf,
     1  nx,xa,xn)

2999  continue

      close(lout)
      close(lvec)

      if(lin.gt.0)close(lpre)

C     New compiler time
      call system_clock(time2)
      tot_time=(time2-time1)/rate

      if(idiag.gt.0)write(6,*)'Model run OK'
      if(idiag.gt.0)write(6,244)tot_time
244   FORMAT(' Elapsed time (s) = ',F8.1)


      END

