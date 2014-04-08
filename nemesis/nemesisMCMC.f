      PROGRAM nemesisMCMC
C     $Id:
C     ******************************************************************
C
C     Retrieval code utilising correlated-k and MCMC chain.
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
C			Modified for MCMC 11/6/13
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

      character*100 buffer,ename
      integer i,j,iscat,ilbl,ica,k,lspec,lout,ispec,nspec,nspecx,ioff
      real xn(mx),se(my),err1(mx),woff,xdiff
      real fwhm,xlat,xlon,st(mx,mx),varparam(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xlatx,varparamx(mvar,mparam)
      real stx(mx,mx),xlonx
      integer varident(mvar,3),varidentx(mvar,3),iform
      integer npro,nvmr,ispace,nav(mgeom),lraw,nprox,lpre
      character*100 runname
      integer ngeom, nwave(mgeom),nconv(mgeom), nx, ny, jsurf, jsurfx
      integer ngas,ncont,nvar,nvarx,lin,nxx,igeom,nconv1,nwave1,jalb
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real xa(mx),rerr(mgeom,mconv),sa(mx,mx),y(my),yn(my)
      real xnx(mx),kk(my,mx)
      real wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real vwave1(mwave),vconv1(mconv)
      double precision aa(mx,mx),dd(mx,my)
      real vkstart,vkend,vkstep
      integer idump,kiter,jtan,jtanx,jalbx,jpre,jprex
      integer jrad,jradx,lx(mx),jlogg,jloggx
      integer miter,niter,idum
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


C     *******************************************************
C     ****************   CODE *******************************
C     *******************************************************

C     Read in reference gas information data
      CALL RESERVEGAS

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

C     Read in number or reporting steps, miter, and number of iterations
C     between steps, niter 
      READ(32,*)miter,niter

C     Read in random negative integer as seed for random number 
C     generator
      READ(32,*)idum

C     Read in total number of spectra to fit and starting offset
      READ(32,*)nspec,ioff


      CLOSE(32)

C     Open spectra file
      lspec=37
      CALL file(runname,runname,'spx')
      open(lspec,file=runname,status='old')
     
C     Open output file
      lout=38
      CALL file(runname,runname,'mcm')
      open(lout,file=runname,status='unknown')
      write(lout,*)nspec,' ! Total number of fitted spectra'
      write(lout,*)miter,' ! Total number of reporting steps'
      write(lout,*)niter,' ! Number of iterations between reports'

      iform=0

C     skip first ioff-1 spectra
      do ispec=1,ioff-1

       call readnextspavX(lspec,iform,woff,xlat,xlon,ngeom,nav,ny,y,se,
     & fwhm,nconv,vconv,angles,wgeom,flat,flon)

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
        se(k)=se(k)+(rerr(i,j))**2
       enddo
      ENDDO

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
      CALL readapriori(runname,lin,lpre,xlat,npro,nvar,varident,
     1  varparam,jsurf,jalb,jtan,jpre,jrad,jlogg,nx,xa,sa,lx)
	

      write(lout,*)nx,' ! nx'
      write(lout,*)ny,' ! ny'

      DO i = 1, nx
        xn(i)=xa(i)
      ENDDO 

      idump=0	! flag for diagnostic print dumps

      ica=1     ! flag for intermediate iterations reports
      if(ica.eq.1)then
       CALL FILE(runname,runname,'itm')
       open(39,file=runname,status='unknown')
       write(39,*)nx,ny,niter,miter
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

	
      call coreretMCMC(runname,ispace,iscat,ilbl,ica,miter,niter,
     1  fwhm,xlat,ngeom,nav,nwave,vwave,nconv,vconv,angles,
     2  gasgiant,nvar,varident,varparam,npro,jsurf,jalb,jtan,
     3  jpre,jrad,jlogg,wgeom,flat,nx,lx,xa,sa,ny,y,se,idum)

      if(ica.eq.1)close(39)
 
2999  continue

      close(lspec)
      close(lout)
      

C     New compiler time
      call gettime(time)
      time2=time
      tot_time=sngl(time2-time1)

      write(6,*)'Model run OK'
      WRITE(6,244)tot_time
244   FORMAT(' Elapsed time (s) = ',F8.1)


      END




