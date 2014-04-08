      PROGRAM nemesisdisc
C     $Id:
C     ******************************************************************
C
C     CIRS retrieval code utilising correlated-k, thermal emission 
C     fast gradient radiative transfer model CIRSRADG. Extension of
C     Nemesis to calculated disc-averaged spectra either for standalone
C     planets (like brown dwarfs), or exoplanets in orbit about a star.  
C
C     This is a special version of the code for computing disc-averaged 
C     spectra (for thermal emission, non scattering cases only) using the
C     E3 hemispheric integration formalism. 
C
C     Minimisation is achieved using a modified non-linear estimation
C     which uses a Marquardt-Levenburg type brake.
C
C     Code looks to see if there is a <runname.sol> file here and in several
C     subsequent subroutines
C
C     If a *.sol file exists then it is assumed that we're calculating the
C     secondary transit of an exoplanet, in which case the code calculates
C     the flux ratio of the integrated star and planet fluxes. 
C
C     If a *.sol file does not exist, then the code calculates the surface 
C     spectral irradiace of the planet/brown-dwarf in units of W cm-2 um-1 or
C     W cm-2 (cm-1)-1

C     Pat Irwin	        Modified from NIMS retrieval code 21/3/00
C			Updated	4/4/01
C			Updated for continuous vmr profiles 7/10/03
C			Updated for FOV-averaging 9/2/04
C			disc-averaged version   27/9/10
C			Overhauled              25/4/12
C
C     ******************************************************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/planrad.f'
      INCLUDE 'arraylen.f'

C     New compiler time
      real tot_time
      double precision time,time1,time2
C     TIME: Temporary variable returned by GETTIME containing the system time.
C     TIME1: System time at the beginning of program execution.
C     TIME2: System time at the end of program execution.

      character*100 buffer,ename
      integer i,j,iscat,ica,k,lspec,lout,ispec,nspec,nspecx,ioff
      real xn(mx),se(my),err1(mx),woff,xdiff
      real fwhm,xlat,xlon,st(mx,mx),varparam(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xlatx,varparamx(mvar,mparam)
      real stx(mx,mx),xlonx
      integer varident(mvar,3),varidentx(mvar,3),lx(mx)
      integer npro,nvmr,ispace,nav(mgeom),lraw,nprox,lpre
      character*100 runname
      integer ngeom, nwave(mgeom), nconv(mgeom), nx, ny, jsurf, jsurfx
      integer ngas,ncont,nvar,nvarx,lin,nxx,igeom,nconv1,nwave1,jalb
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real xa(mx),rerr(mgeom,mconv),sa(mx,mx),y(my),yn(my)
      real xnx(mx),kk(my,mx),xerr
      real wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real vwave1(mwave),vconv1(mconv)
      integer iform,jlogg,jloggx
      double precision aa(mx,mx),dd(mx,my)
      real vkstart,vkend,vkstep
      integer idump,kiter,jtan,jtanx,jalbx,jpre,jrad,jradx,jprex
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
      character*100 solfile,solname
      logical solexist
      common/solardat/iread, iform, solradius, solwave, solrad,  solnpt


C     ******************************************************


C     *******************************************************
C     ****************   CODE *******************************
C     *******************************************************

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

      print*,'Files OK',gasgiant

C     Read start, end and step of tables
      CALL readkkhead(runname,vkstart,vkend,vkstep)

      CALL file(runname,runname,'inp')
      OPEN(32,file=runname,status='old')

C     Read in whether to calculate with wavenumbers(0) or wavelength(1)
C     Also read in whether scattering is required (iscat)
      READ(32,*)ispace,iscat

      if(iscat.ne.0)then
       print*,'Nemesisdisc does not work for iscat <> 0'
       stop
      endif

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


      iform1=1
      READ(32,*,END=999)iform1
999   continue
      CLOSE(32)

      print*,'iform1 = ',iform1
      if(iform1.ne.1.and.iform1.ne.3)then
       print*,'Error in input file. Iform can be 1 or 3 for Nemesisdisc'
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
        print*,'Error in Nemesisdisc. Flux-ratio calculation defined'
        print*,'or assumed, but no solar file present'
        stop
       endif
      endif

      iform=iform1


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
      CALL readapriori(runname,lin,lpre,xlat,npro,nvar,varident,
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
       READ(12,1) ANAME
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


	
      call coreretdisc(runname,ispace,iscat,ica,kiter,phlimit,
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

      CALL writeout(iform,runname,ispace,lout,ispec,xlat,xlon,
     1 npro,nvar,varident,varparam,nx,ny,y,yn,se,xa,sa,xn,err1,
     2 ngeom,nconv,vconv,gasgiant,jpre,jrad,jlogg,iscat,lin)

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


      END




