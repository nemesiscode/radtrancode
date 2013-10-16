      PROGRAM generateMCSspx
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
      double precision time,time1,time2
C     TIME: Temporary variable returned by GETTIME containing the system time.
C     TIME1: System time at the beginning of program execution.
C     TIME2: System time at the end of program execution.

      character*200 buffer,ename
      integer i,j,iscat,ica,k,lspec,lout,ispec,nspec,nspecx,ioff
      integer linfo,npvar,iplanet
      real altbore,marsradius,satrad,thetrot,RADIUS
      real xn(mx),se(my),err1(mx),woff,xdiff
      real fwhm,xlat,xlon,st(mx,mx),varparam(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xlatx,varparamx(mvar,mparam)
      real stx(mx,mx),xlonx
      integer varident(mvar,3),varidentx(mvar,3),iscat1
      integer npro,nvmr,ispace,nav(mgeom),lraw,nprox,lpre
      character*100 runname    
      character*150 dummy	
C dummy: Character variable, used for reading-in header on info file   
      integer ngeom, nwave(mgeom), nconv(mgeom), nx, ny, jsurf, jsurfx
      integer np,lin1,ioffx,ivarx,npx,iform
      integer nwaveT(mgeom), nconvT(mgeom)
      integer ngas,ncont,nvar,nvarx,lin,nxx,igeom,nconv1,nwave1,jalb
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real vwaveT(mgeom,mwave),vconvT(mgeom,mconv)
      real xa(mx),rerr(mgeom,mconv),sa(mx,mx),y(my),yn(my)
      real xnx(mx),kk(my,mx)
      real wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real vwave1(mwave),vconv1(mconv)
      double precision aa(mx,mx),dd(mx,my)
      real vkstart,vkend,vkstep
      integer idump,kiter,jtan,jtanx,jalbx,jpre,jprex,idum,lvec
      integer ivar,jrad
C     ********** Scattering variables **********************
      real xwave(maxsec),xf(maxcon,maxsec),xg1(maxcon,maxsec)
      real xg2(maxcon,maxsec)
      real tnco,twave,frac,tico
      real phlimit,kkcor(mx,mx)
      logical gasgiant
      COMMON /hgphas/xwave,xf,xg1,xg2,tnco,twave,frac,tico
      COMMON /scatdump/ idump

      include '../radtran/includes/ciacom.f'

      CHARACTER*100 ANAME
      REAL DNU
      INTEGER IPARA


C     ******************************************************

C     ******* xlat Mod Variables (see code for info) *******
      integer inav, iav
      integer xlat_lmb_cnt
      real emiss_ang
      real xlat_lmb, xlat_lmb_store
C     ******************************************************


C     *******************************************************
C     ****************   CODE *******************************
C     *******************************************************

C     ----------- Scattering phase function initialisation --------------
      xwave(1)=-1                       ! Reset to force read of hgphase*
C                                         files.
C     ------------ Scattering phase function initialisation -------------

C     New compiler time
      call gettime(time)
      time1=time

      CALL prompt('Enter run name : ')
      READ(5,1)buffer
1     FORMAT(a)
      runname = buffer(1:36)

      print*,'checking files'
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

      CALL readrefiplan(runname,iplanet,RADIUS)

      kiter=0
      woff=0.0

C     Read in total number of spectra to simulate
      READ(32,*)nspec


C     Read in random -ve seed number
      READ(32,*)idum

C     Read in name of forward modelling error file
      READ(32,1)ename

C     Read in lin identifier in case want to use retrieved T
      READ(32,*)lin

      CLOSE(32)

      iform=0
     
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
        print*,'.pre file does not contain enough'
        print*,'retrievals'
        stop
       endif
      endif

C     Open spectra file
      lspec=37
      CALL file(runname,runname,'spx')
      write(6,1)runname
      open(lspec,file=runname,status='old')
      linfo=40
      CALL file(runname,runname,'info')
      open(linfo,file=runname,status='old')

C     Skip header lines of info file
182   continue
       read(linfo,'(a150)') dummy
       if (dummy(1:1).eq.'#') then
          goto 182
       endif
      backspace(linfo)

C     Read in sample measurement vector, obs. geometry and covariances
       call readnextspavX(lspec,iform,woff,xlat,xlon,ngeom,nav,ny,y,se,
     1  fwhm,nconv,vconv,angles,wgeom,flat,flon)

       call readnextinfo(linfo,altbore,marsradius,satrad,thetrot)

      close(lspec)
      close(linfo)

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
     1  varparam,jsurf,jalb,jtan,jpre,jrad,nx,xa,sa)
	
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
      print*, 'ename: (', ename, ')' 
      call forwarderr(ename,ngeom,nconv,vconv,woff,rerr)

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

      do 2999 ispec=1,nspec

       call rmodapriori(idum,npro,nvar,varident,varparam,
     1  nx,xa,sa,xn)

       if(lin.gt.0) then
        call readraw(lpre,xlatx,xlonx,nprox,nvarx,varidentx,
     1   varparamx,jsurfx,jalbx,jtanx,jprex,nxx,xnx,stx)

        ioff=0
        do ivar=1,nvar
         np=1
         if(varident(ivar,1).le.100)then
          np=npvar(varident(ivar,3).eq.0,npro)
         endif

         if(varident(ivar,1).eq.0.and.np.eq.npro)then
          ioffx=0
          do ivarx=1,nvarx
           npx=1
           if(varidentx(ivarx,1).le.100)then
            npx=npvar(varidentx(ivarx,3),npro)
           endif
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
        CALL forwardavfovMCS(runname,ispace,fwhm,xlat,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,marsradius,
     3   satrad,thetrot,altbore,nx,xn,ny,yn,kk)
      elseif(iscat.eq.1)then 
       CALL forwardnogX(runname,ispace,iscat,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,RADIUS,
     3   nx,xn,ny,yn,kk,kiter)
      elseif(iscat.eq.2)then
       CALL intradfield(runname,ispace,xlat,nwaveT,vwaveT,nconvT,
     1   vconvT,gasgiant,lin,nvar,varident,varparam,jsurf,jalb,
     2   jtan,jpre,nx,xn)
       iscat1=1
       CALL forwardnogX(runname,ispace,iscat1,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,RADIUS,
     3   nx,xn,ny,yn,kk,kiter)
      else
       iscat1=1
       CALL forwardnogX(runname,ispace,iscat1,fwhm,ngeom,nav,
     1   wgeom,flat,nwave,vwave,nconv,vconv,angles,gasgiant,lin,
     2   nvar,varident,varparam,jsurf,jalb,jtan,jpre,jrad,RADIUS,
     3   nx,xn,ny,yn,kk,kiter)
      endif

C     Mod to make sure XLAT is set to LMB value (also checks all LMB values the same)	
      xlat_lmb_cnt=0
      xlat_lmb_store=0.
      do igeom=1,ngeom
       do inav=1,nav(igeom)
        iav=nav(inav)
        emiss_ang = angles(igeom,iav,2)
        if (emiss_ang.lt.0.) then 
         xlat_lmb=flat(igeom,iav)
         if(xlat_lmb.ne.xlat_lmb_store.and.xlat_lmb_cnt.gt.0) then
          print*, 'xlat_lmb NE xlat_lmb_store (Limb xlat do not match)'
          STOP
         endif
         xlat_lmb_store=xlat_lmb
         xlat_lmb_cnt=xlat_lmb_cnt+1
         print*,'iav, emiss_ang = ', iav, emiss_ang
         print*,'xlat_lmb_cnt = ',xlat_lmb_cnt
         print*,'xlat_lmb, xlat_lmb_store = ',xlat_lmb,xlat_lmb_store
        endif
       enddo
      enddo
      xlat=xlat_lmb
      
C      Write out k-matrix for reference
       open(52,file='kk.out',form='unformatted',status='unknown')
       write(52)y,yn
       write(52)kk
       close(52)

      call writenextspavX(lout,iform,idum,woff,xlat,xlon,ngeom,nav,ny,
     1  yn,se,fwhm,nconv,vconv,angles,wgeom,flat,flon)

      call writemvec(lvec,npro,nvar,varident,varparam,jsurf,
     1  nx,xa,xn)

2999  continue

      close(lout)
      close(lvec)

      if(lin.gt.0)close(lpre)

C     New compiler time
      call gettime(time)
      time2=time
      tot_time=sngl(time2-time1)

      write(6,*)'Model run OK'
      WRITE(6,244)tot_time
244   FORMAT(' Elapsed time (s) = ',F8.1)


      END




