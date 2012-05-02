      PROGRAM nemesisPT
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
C TIME: Temporary variable returned by GETTIME containing the system time.
C TIME1: System time at the beginning of program execution.
C TIME2: System time at the end of program execution.

      character*100 buffer,ename
      integer i,j,iscat,ica,k,lspec,lout,ispec,nspec,nspecx,ioff
      real xn(mx),se(my),err1(mx),woff,err2(mx),xdiff
      real fwhm,xlat,xlon,st(mx,mx),varparam(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xlatx,varparamx(mvar,mparam)
      real stx(mx,mx),xlonx
      integer varident(mvar,3),varidentx(mvar,3),igeom
      integer npro,nvmr,ispace,nav(mgeom),lraw,nprox,lpre
      character*100 runname
      character*150 dummy
      integer ngeom, nwave(mgeom), nconv(mgeom), nx, ny, jsurf
      integer ngas,ncont,nvar,nvarx,lin,nxx,jsurfx,nconv1,nwave1
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real kk(my,mx),xa(mx),rerr(mgeom,mconv),sa(mx,mx)
      real y(my),yn(my),xnx(mx),radius
      real wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real vconv1(mconv),vwave1(mwave)
      double precision aa(mx,mx),dd(mx,my)
      real vkstart,vkend,vkstep
      integer idump,kiter,jtan,jalb,jalbx,jpre,jtanx,jprex,iplanet
      integer jrad,jradx
C     ********** Scattering variables **********************
      real xwave(maxsec),xf(maxcon,maxsec),xg1(maxcon,maxsec)
      real xg2(maxcon,maxsec)
      real tnco,twave,frac,tico
      real phlimit,kkcor(mx,mx)
      logical gasgiant
      COMMON /hgphas/xwave,xf,xg1,xg2,tnco,twave,frac,tico
      COMMON /scatdump/ idump

      INCLUDE '../radtran/includes/ciacom.f'
      INCLUDE '../radtran/includes/gascom.f'

      CHARACTER*100 ANAME
      REAL DNU
      INTEGER IPARA

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


      CALL readrefiplan(runname,iplanet,radius)

      print*,'Files OK - plan radius = ',radius

C     Read start, end and step of tables
      CALL readkkhead(runname,vkstart,vkend,vkstep)

      CALL file(runname,runname,'inp')
      OPEN(32,file=runname,status='old')

C     Read in whether to calculate with wavenumbers(0) or wavelength(1)
C     Also read in whether scattering is required (iscat)
      READ(32,*)ispace,iscat

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
C              and used a apriori for next current retrieval.
C     lin = 3  indicates that previous retrieval should be considered
C              and used as a priori for all parameters that match, and
C              used to fix all other parameters (including effect of
C              propagation of retrieval errors).     
      READ(32,*)lin

      CLOSE(32)

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

       call readnextspavX(lspec,woff,xlat,xlon,ngeom,nav,ny,y,se,
     & fwhm,nconv,vconv,angles,wgeom,flat,flon)

C      Look to see if previously retrieved information is to be used
C      and if so, skipped
       if(lin.gt.0)then
      
       print*,'Reading in previously retrieved parameters'
        call readraw(lpre,xlatx,xlonx,nprox,nvarx,varidentx,
     1    varparamx,jsurfx,jalbx,jtanx,jprex,jradx,nxx,xnx,stx)
      
       endif

      enddo

      do 2999 ispec=ioff,ioff-1+nspec

C     Read in measurement vector, obs. geometry and covariances
      call readnextspavX(lspec,woff,xlat,xlon,ngeom,nav,ny,y,se,
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

C     Calculate the tabulated wavelengths of c-k look up tables
      do igeom=1,ngeom
       nconv1 = nconv(igeom)
       if(igeom.eq.1)print*,nconv1
       do j=1,nconv1
        vconv1(j)=vconv(igeom,j)
        if(igeom.eq.1)print*,vconv1(j)
       enddo
       CALL wavesetb(runname,vkstart,vkend,vkstep,nconv1,vconv1,fwhm,
     1  nwave1,vwave1)
       if(igeom.eq.1)print*,nwave1
       do j=1,nwave1
        vwave(igeom,j)=vwave1(j)
        if(igeom.eq.1)print*,vwave1(j)
       enddo
       nwave(igeom)=nwave1
      enddo

C     set up a priori of x and its covariance
      CALL readapriori(runname,lin,lpre,xlat,npro,nvar,varident,
     1  varparam,jsurf,jalb,jtan,jpre,jrad,nx,xa,sa)

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

      jpre=-1
      jalb=-1
      jtan=-1
      jsurf=-1
      print*,'Radius - 0 = ',jrad
      jrad=-1

      call coreretPT(runname,ispace,iscat,ica,kiter,phlimit,
     1  fwhm,xlat,ngeom,nav,nwave,vwave,nconv,vconv,angles,
     2  gasgiant,lin,lpre,nvar,varident,varparam,npro,jsurf,jalb,jtan,
     3  jpre,jrad,radius,wgeom,flat,nx,xa,sa,ny,y,se,xn,sm,sn,st,
     4  yn,kk,aa,dd)

C     Calculate retrieval errors.
C     Simple errors, set to sqrt of diagonal of ST
      do i=1,nx
       err1(i)=sqrt(abs(st(i,i)))
      enddo

C     Analyse final covariance matrix into its eigenvectors and hence get
C     more representative estimate of the final errors
C      call eigenerr(nx,st,err2)

C     write output
      print*,'Nemesis: Writing output'

      CALL writeoutPT(ispace,lout,ispec,xlat,xlon,npro,nvar,
     1 varident,varparam,nx,ny,y,yn,se,xa,sa,xn,err1,ngeom,
     2 nconv,vconv)

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




