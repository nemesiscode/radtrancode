      subroutine nemesisEMCEE(runname, specsize, mcntemp, 
     1 mcnvmr, ith, intemp, invmr, vflag, gflag,
     2 inmass, inrad, MCMCspec)
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
      include '../radtran/includes/planrad.f'
      include '../radtran/includes/emcee.f'
      INCLUDE 'arraylen.f'

C     New compiler time
      real tot_time
      double precision time,time1,time2
C     TIME: Temporary variable returned by GETTIME containing the system time.
C     TIME1: System time at the beginning of program execution.
C     TIME2: System time at the end of program execution.

      character*100 buffer,ename
      integer i,j,iscat,ilbl,ica,k,lspec,lout,ispec,nspec,nspecx,ioff
      real xn(mx),se(my),err1(mx),woff,xdiff
      real fwhm,xlat,xlon,st(mx,mx),varparam(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xlatx,varparamx(mvar,mparam)
      real stx(mx,mx),xlonx
      integer varident(mvar,3),varidentx(mvar,3),iform
      integer npro,ispace,nav(mgeom),lraw,nprox,lpre, nvmr
      integer ngeom, nwave(mgeom),nconv(mgeom), nx, ny, jsurf, jsurfx
      integer ngas,ncont,nvar,nvarx,lin,nxx,igeom,nconv1,nwave1,jalb
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real xa(mx),rerr(mgeom,mconv),sa(mx,mx),y(my),yn(my)
      real xnx(mx),kk(my,mx),xerr
      real wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real vwave1(mwave),vconv1(mconv)
      double precision aa(mx,mx),dd(mx,my)
      real vkstart,vkend,vkstep
      integer idump,kiter,jtan,jtanx,jalbx,jpre,jprex
      integer jrad,jradx,lx(mx),jlogg,jloggx
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

C     Solar spectrum variables and flags
      integer iform1,iread,solnpt
      real solwave(maxbin),solrad(maxbin),solradius
      character*100 solfile,solname,cellfile
      integer cellngas,cellid(maxgas),celliso(maxgas),icread
      real cellength,cellpress,celltemp,cellvmr(maxgas)
      logical solexist,cellexist
      common/solardat/iread, iform, solradius, solwave, solrad,  solnpt
      common/celldat/icread,cellngas,cellid,celliso,cellvmr,cellength,
     1  cellpress,celltemp

c     RG EMCEE INPUTS
c 999 = nconv number of wavelengths, must be exact for python to not read in zeros
c and change the shape of the array
      character*100 runname
      integer arrsize, specsize,ith, vflag, gflag
      real MCMCspec(specsize), xfac
      real intemp(mcntemp), invmr(mcnvmr)
      real inrad, inmass
      integer mcntemp, mcnvmr
      character*3 sith
      character*255 path, path1, oldpath, path2
      logical countexist

cf2py intent(in) runname, specsize,mcntemp,mcnvmr
cf2py intent(in) intemp,invmr,ith
cf2py intent(in) inmass, inrad, vflag, gflag
cf2py intent(out) MCMCspec
cf2py depend(mcntemp) intemp
cf2py depend(mcnvmr) invmr
cf2py depend(specsize) MCMCspec

C     ******************************************************


C     *******************************************************
C     ****************   CODE *******************************
C     *******************************************************

C RG find free directory to do calculation

      CALL getcwd(path1)

      if (ith.lt.10)then
       write( sith, '(i1)' ) ith
      elseif (ith.ge.10.and.ith.lt.100)then
       write( sith, '(i2)' ) ith
      else
       write( sith, '(i3)' ) ith
      endif
      print*, path, ith, sith
      CALL chdir(TRIM(path1)//'/'//TRIM(sith))
      CALL getcwd(path)
      print*, ith, sith, path, path1

c RG assign array elements. This intermediate way of passing is due to
c being unable to pass the allocatable arrays (necessary for the F77-Py
c interface so that dimensions match during passing) in a common block. 

      VMRflag = vflag
      GRAVflag = gflag
      MCMCmass = inmass
      MCMCrad = inrad


      DO i=1, mcnvmr
       MCMCvmr(i) = 10**(invmr(i))
      ENDDO

      DO i=1, mcntemp
       MCMCtemp(i) = intemp(i)
      ENDDO

C     Read in reference gas information data
      CALL RESERVEGAS

C     ----------- Scattering phase function initialisation --------------
      xwave(1)=-1                       ! Reset to force read of hgphase*
C                                         files.
C     ------------ Scattering phase function initialisation -------------
      jradf=-1
      jloggf=-1


C     New compiler time
      call gettime(time)
      time1=time

c      CALL prompt('Enter run name : ')
c      READ(5,1)buffer
1     FORMAT(a)
c      runname = buffer(1:36)
c      runname = 'cirstest'

      print*,'checking files'
C     Make sure input files are OK
      CALL checkfiles(runname)

      CALL readrefhead(runname,npro,nvmr,gasgiant)
      if(npro.gt.maxpro)then
       print*,'Error in Nemesis. npro > maxpro : ',npro,maxpro
       stop
      endif

      print*,'Files OK',gasgiant

      CALL file(runname,runname,'inp')
      OPEN(32,file=runname,status='old')

C     Read in whether to calculate with wavenumbers(0) or wavelength(1)
C     Also read in whether scattering is required (iscat)
C     Also read in whether lbl calculation is required (ilbl)
      READ(32,*)ispace,iscat,ilbl

      if(ilbl.gt.0) then 
       print*,'Nemesis - LBL calculation'
      endif
      if(ilbl.eq.0)CALL readkkhead(runname,vkstart,vkend,vkstep)

C     Read any wavenumber offset to add to measured spectra
      READ(32,*)woff   

C     Read in name of forward modelling error file
      READ(32,1)ename     

C     Read in number of iterations
      READ(32,*)kiter 
      
C     Read limiting % change in phi
      READ(32,*)phlimit

C     Read in total number of spectra to fit and starting offset
      READ(32,*)nspec,ioff
C     Read in integer indicating if previous retrieval is to be
C     used to set some elements of the prf file (e.g. T-profile may
C     have already been retrieved from other wavelengths)'
C     lin = 0  indicates no previous retrievals
C     lin = 1  indicates that previous retrieval should be considered
C              and effect of retrieval errors accounted for
C     lin = 2  indicates that previous retrieval should be considered 
C              and used as a priori for current retrieval.
C     lin = 3  indicates that previous retrieval should be considered
C              and used as a priori for all parameters that match, and
C              used to fix all other parameters (including effect of 
C              propagation of retrieval errors).
      READ(32,*)lin
      iform1=0
      READ(32,*,END=999)iform1
999   continue
      CLOSE(32)


      print*,'iform1 = ',iform1

      if(iform1.eq.2)then
       print*,'Error in input file. Iform can be 0, 1 or 3 for Nemesis'
       stop
      endif
      
C     See if there is a solar or stellar reference spectrum and read in
C     if present.
      call file(runname,solfile,'sol')
      inquire(file=solfile,exist=solexist)
      if(solexist)then
         call opensol(solfile,solname)
         CALL init_solar_wave(ispace,solname)
      else
         if(iform1.eq.1)then
          print*,'Error in Nemesis. Flux-ratio calculation defined'
          print*,'but no solar file exists'
          stop
         endif
      endif
   
      iform=iform1


C     See if there is a SCR cell file present and read in
C     if present.
      icread=0
      call file(runname,cellfile,'cel')
      inquire(file=cellfile,exist=cellexist)
      if(cellexist)then
         icread=1
         call readcell(cellfile)
      endif

C     Open spectra file
      lspec=37
      CALL file(runname,runname,'spx')
      open(lspec,file=runname,status='old')
     
C     Open output file
      lout=38
      lraw=36
      CALL file(runname,runname,'mre')
      open(lout,file=runname,status='unknown')
      CALL file(runname,runname,'raw')
      open(lraw,file=runname,status='unknown')
      write(lout,*)nspec,' ! Total number of retrievals'
      write(lraw,*)nspec,' ! Total number of retrievals'
      print*,'lin = ',lin

      if(lin.gt.0)then
C      if previous retrieval to be considered, 
C      open previous raw retrieval file (copied to .pre)
       lpre=39
       CALL file(runname,runname,'pre')
       print*,'Nemesis: reading previous retrieval : ',runname
       open(lpre,file=runname,status='old')
       read(lpre,*)nspecx
       if(nspec+ioff-1.gt.nspecx)then
        print*,'.pre file does not contain enough'
        print*,'retrievals'
        stop
       endif
      endif


C     skip first ioff-1 spectra
      do ispec=1,ioff-1
       call readnextspavX(lspec,iform,woff,xlat,xlon,ngeom,nav,ny,y,se,
     & fwhm,nconv,vconv,angles,wgeom,flat,flon)

C      Look to see if previously retrieved information is to be used
C      and if so, read in
       if(lin.gt.0)then
      
        call readraw(lpre,xlatx,xlonx,nprox,nvarx,varidentx,
     1   varparamx,jsurfx,jalbx,jtanx,jprex,jradx,jloggx,nxx,
     2   xnx,stx)
      
       endif

      enddo

      do 2999 ispec=ioff,ioff-1+nspec

C     Read in measurement vector, obs. geometry and covariances
      call readnextspavX(lspec,iform,woff,xlat,xlon,ngeom,nav,ny,y,se,
     1  fwhm,nconv,vconv,angles,wgeom,flat,flon)
     
C		reset jradf and jloggf to -1 for start of each run

		jradf=-1
		jloggf=-1

      if(icread.eq.1)then
       do igeom=1,ngeom
        nconv(igeom)=nconv(igeom)/2
       enddo
      endif
C     Read in forward modelling errors
      call forwarderr(ename,ngeom,nconv,vconv,woff,rerr)

C     Add forward errors to measurement covariances
      k=0
      DO i=1,ngeom
       do j=1,nconv(i)
        k = k+1
        xerr=rerr(i,j)
        if(iform.eq.3)xerr=xerr*1e-18
        se(k)=se(k)+xerr**2
       enddo
      ENDDO

C     If we are doing an SCR cell, then we only need half the 
C     wavelengths in the spx file as we're calculating both WB and SB
C     outputs

      if(ilbl.eq.0)then
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
      endif
C     set up a priori of x and its covariance
      CALL readaprioriMCMC(runname,lin,lpre,xlat,npro,nvar,varident,
     1  varparam,jsurf,jalb,jtan,jpre,jrad,jlogg,nx,xa,sa,lx)
	
      DO i = 1, nx
        xn(i)=xa(i)
      ENDDO 

      idump=0	! flag for diagnostic print dumps

      if(nspec.eq.1)then
        ica = 1		! 1 = single retrieval
      else
        ica = 0		! 0 = multiple retrievals
      endif


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
      call coreret(runname,ispace,iscat,ilbl,ica,kiter,phlimit,
     1  fwhm,xlat,ngeom,nav,nwave,vwave,nconv,vconv,angles,
     2  gasgiant,lin,lpre,nvar,varident,varparam,npro,jsurf,jalb,jtan,
     3  jpre,jrad,jlogg,wgeom,flat,nx,lx,xa,sa,ny,y,se,xn,sm,sn,st,yn,
     4  kk,aa,dd)
C     Calculate retrieval errors.
C     Simple errors, set to sqrt of diagonal of ST
      do i=1,nx
       err1(i)=sqrt(abs(st(i,i)))
      enddo

C     write output
      CALL writeout(iform,runname,ispace,lout,ispec,xlat,xlon,npro,
     1 nvar,varident,varparam,nx,ny,y,yn,se,xa,sa,xn,err1,ngeom,
     2 nconv,vconv,gasgiant,jpre,jrad,jlogg,iscat,lin)

      CALL writeraw(lraw,ispec,xlat,xlon,npro,nvar,varident,
     1 varparam,nx,xn,st)

      if(ica.eq.1)then
C       Write out all the error matrices if only one case retrieved
        call write_covariance(runname,npro,nvar,varident,varparam,
     1    nx,ny,sa,sm,sn,st,se,aa,dd,kk)
      endif
2999  continue

      close(lspec)
      close(lout)
      close(lraw)
      if(lin.gt.0)close(lpre)

C     New compiler time
      call gettime(time)
      time2=time
      tot_time=sngl(time2-time1)

      write(6,*)'Model run OK'
      WRITE(6,244)tot_time
244   FORMAT(' Elapsed time (s) = ',F8.1)
      if(ispace.eq.0) then
C      Wavenumber space

C      Default format
       if(iform.eq.0)then
        write(*,*)'Radiances expressed as nW cm-2 sr-1 cm'
        xfac=1e9

C      F_plan/F_star format
       elseif(iform.eq.1)then
         write(*,*)'F_plan/F_star Ratio of planet'
         xfac=1.0

C      Spectral power format
       elseif(iform.eq.3)then
      write(*,*)'Spectral Radiation of planet: W (cm-1)-1'
         xfac=1e18

C      NemesisPT format
       elseif(iform.eq.2) then
        write(*,*)'Transit depth: 100*Planet_area/Stellar_area'
        xfac=1.

C      Default
       else
        print*,'Error in writeout - iform not defined. Default=0'
        write(*,*)'Radiances expressed as nW cm-2 sr-1 cm'
        xfac = 1e9
       endif

      else
C      Wavelength space

C      Default format
       if(iform.eq.0)then
        write(*,*)'Radiances expressed as uW cm-2 sr-1 um-1'
        xfac = 1e6

C      F_plan/F_star format
       elseif(iform.eq.1)then
         write(*,*)'F_plan/F_star Ratio of planet'
         xfac=1.0

C      Spectral irradiance format
       elseif(iform.eq.3)then
      write(*,*)'Spectral Radiation of planet: W um-1'
         xfac=1e18

C      NemesisPT format
       elseif(iform.eq.2)then
        write(*,*)'Transit depth: 100*Planet_area/Stellar_area'
        xfac=1.

C      Default format
       else
        print*,'Error in writeout - iform not defined. Default=0'
        write(*,*)'Radiances expressed as uW cm-2 sr-1 um-1'
        xfac=1e6
       endif

      endif

      do i=1,ngeom
       do j=1,nconv(i)
         MCMCspec(j) = yn(j)*xfac
        enddo
      enddo

c      open(unit=200, file='out.dat', status='append')
c      write(200,*), 'FINISHED:',ith, sith, TRIM(path), MCMCspec
c      close(200)

      CALL chdir(TRIM(path1))
      CALL getcwd(path1)
      print*, path1
      END



