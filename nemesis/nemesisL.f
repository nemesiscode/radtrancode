      PROGRAM nemesisL
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
      INCLUDE 'arraylen.f'

C     New compiler time
      real tot_time
      real rate
      integer c1,c2,cr,time1,time2,cm
C     TIME: Temporary variable returned by GETTIME containing the system time.
C     TIME1: System time at the beginning of program execution.
C     TIME2: System time at the end of program execution.

      character*100 buffer,ename
      integer i,j,iscat,ica,k,lspec,lout,loutp,ispec,nspec,nspecx,ioff
      real xn(mx),se(my),err1(mx),woff,xdiff
      real fwhm,xlat,xlon,st(mx,mx),varparam(mvar,mparam)
      real sn(mx,mx),sm(mx,mx),xlatx,varparamx(mvar,mparam)
      real stx(mx,mx),xlonx
      integer varident(mvar,3),varidentx(mvar,3),igeom,iform,iform1
      integer npro,nvmr,ispace,nav(mgeom),lraw,nprox,lpre,lprx
      integer ilbl,ilbl1,ishape,inum,ionpeel
      character*100 runname,solfile,solname,sfile,plotname
      logical solexist,percbool
      integer ngeom, nwave(mgeom), nconv(mgeom), nx, ny, jsurf
      integer ngas,ncont,nvar,nvarx,lin,nxx,jsurfx,nconv1,nwave1
      integer lx(mx)
      real vwave(mgeom,mwave),vconv(mgeom,mconv),angles(mgeom,mav,3)
      real kk(my,mx),xa(mx),rerr(mgeom,mconv),sa(mx,mx),xerr
      real y(my),yn(my),xnx(mx)
      real wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real vconv1(mconv),vwave1(mwave)
      double precision aa(mx,mx),dd(mx,my)
      real vkstart,vkend,vkstep
      integer idump,kiter,jtan,jalb,jalbx,jpre,jtanx,jprex
      integer jrad,jradx,jlogg,jloggx,jxsc,jxscx,occult
      integer jloggx1,jxscx1,jtanx1,jprex1,jradx1,jalbx1,jsurfx1
      integer nprox1,nvarx1,varidentx1(mvar,3),nxx1
      integer jfrac,jfracx,jfracx1
      real xlatx1,xlonx1,varparamx1(mvar,mparam),xnx1(mx),stx1(mx,mx)
C     ********** Scattering variables **********************
      real xwave(maxsec),xf(maxcon,maxsec),xg1(maxcon,maxsec)
      real xg2(maxcon,maxsec)
      real tnco,twave,frac,tico
      real phlimit,kkcor(mx,mx)
      logical gasgiant
      COMMON /hgphas/xwave,xf,xg1,xg2,tnco,twave,frac,tico
      COMMON /scatdump/ idump
      COMMON /lbltable/ ilbl1

      INCLUDE '../radtran/includes/ciacom.f'

      CHARACTER*100 ANAME
      REAL DNU
      INTEGER IPARA
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C     ******************************************************


C     *******************************************************
C     ****************   CODE *******************************
C     *******************************************************

C     Read in reference gas information data
      CALL RESERVEGAS

      idiag=1
C     ----------- Scattering phase function initialisation --------------
      xwave(1)=-1                       ! Reset to force read of hgphase*
C                                         files.
C     ------------ Scattering phase function initialisation -------------
      jradf=-1
      jloggf=-1


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
       print*,'Error in NemesisL. npro > maxpro : ',npro,maxpro
       stop
      endif

      if(idiag.gt.0)print*,'Files OK'

      CALL file(runname,runname,'inp')
      OPEN(32,file=runname,status='old')

C     Read in whether to calculate with wavenumbers(0) or wavelength(1)
C     Also read in whether scattering is required (iscat)
C     Also read in whether solar occultation is required (occult)
C     Also read in whether lbl calculation is required (ilbl)
C     Also read in whether gradients are calculated analytically or numerically (inum)
C     Also read in whether onion peeling method is to be used ionpeel

      READ(32,*)ispace,iscat,occult,ilbl,inum,ionpeel

      if(ilbl.eq.1) then
       if(idiag.gt.0)then
        print*,'NemesisL - LBL calculation. Not yet implemented'
       endif
      endif
      if(ilbl.eq.0)then
       if(idiag.gt.0)print*,'NemesisL - corr-k calculation'
       CALL readkkhead(runname,vkstart,vkend,vkstep)
      endif
      if(ilbl.eq.2)then
       if(idiag.gt.0)print*,'Nemesis - lbl-table calculation'
       CALL readkklblhead(runname,vkstart,vkend,vkstep)
      endif

C     Read any wavenumber offset to add to measured spectra
      READ(32,*)woff   

C     Read in name of forward modelling error file
      READ(32,1)ename     

C     Read in number of iterations
      READ(32,*)kiter 

      if(kiter.lt.0)then
        inum = 0   !If just using forward model, use fast analytical derivatives
      endif

      if(inum.eq.1.and.idiag.gt.0)then
       print*,'NemesisL - gradients from numerical differentiation'
      endif
      if(inum.eq.0.and.idiag.gt.0)then
       print*,'NemesisL - gradients calculated analytically'
      endif
      
      if(ionpeel.eq.1.and.idiag.gt.0)then
      print*,'NemesisL - Onion-peeling method (just one tangent height)'
      endif
 

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
C     lin = 4  indicates the same as lin = 3. However, in this case it
C              also indicates that current retrieved variables and the
C              previous ones must be saved in the output .raw file.
C              Useful for sequential retrieval, such as onion-peeling method

      READ(32,*)lin
      iform1=0
      percbool = .false.
      READ(32,*,END=999)iform1
      READ(32,*,END=999)percbool
999   continue
      CLOSE(32)

      if(idiag.gt.0)print*,'iform1 = ',iform1
      if(idiag.gt.0)print*,'percbool = ', percbool

      iform=iform1


C     See if there is a solar or stellar reference spectrum and read in
C     if present.
      call file(runname,solfile,'sol')
      inquire(file=solfile,exist=solexist)
      if(solexist)then
         call opensol(solfile,solname)
         CALL init_solar_wave(ispace,solname)
      else
         if(occult.eq.1)then
          print*,'Error in NemesisL. Flux-ratio calculation defined'
          print*,'but no solar file exists'
          stop
         endif
      endif


C     Open spectra file
      lspec=37
      CALL file(runname,runname,'spx')
      open(lspec,file=runname,status='old')
     
C     Open output file
      lout=38
      loutp=88
      lraw=36
      CALL file(runname,runname,'mre')
      open(lout,file=runname,status='unknown')
      CALL file(runname,runname,'raw')
      open(lraw,file=runname,status='unknown')
      CALL file(runname,plotname,'vre')
      open(loutp,file=plotname,status='unknown')
      write(lout,*)nspec,' ! Total number of retrievals'
      write(lraw,*)nspec,' ! Total number of retrievals'

      if(lin.gt.0)then
C      if previous retrieval to be considered, 
C      open previous raw retrieval file (copied to .pre)
       lpre=39
       CALL file(runname,runname,'pre')
       open(lpre,file=runname,status='old')
       read(lpre,*)nspecx
       if(nspec+ioff-1.gt.nspecx)then
        print*,'.pre file does not contain enough'
        print*,'retrievals'
        stop
       endif

       if(lin.eq.4)then
        lprx=40
        CALL file(runname,runname,'prx')
        open(lprx,file=runname,status='old')
        read(lprx,*)nspecx
        if(nspec+ioff-1.gt.nspecx)then
         print*,'.prx file does not contain enough'
         print*,'retrievals'
         stop
        endif
       endif

      endif


C     skip first ioff-1 spectra
      do ispec=1,ioff-1

       call readnextspavX(lspec,iform,woff,xlat,xlon,ngeom,nav,ny,y,se,
     & fwhm,nconv,vconv,angles,wgeom,flat,flon)

C      Look to see if previously retrieved information is to be used
C      and if so, skipped
       if(lin.gt.0)then
      
        call readraw(lpre,xlatx,xlonx,nprox,nvarx,varidentx,
     1    varparamx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,jloggx,
     2    jfracx,nxx,xnx,stx)
      
       endif

      enddo

      do 2999 ispec=ioff,ioff-1+nspec

C      Reading if there are previous retrieved variables
       if(lin.eq.4)then
        call readraw(lprx,xlatx1,xlonx1,nprox1,nvarx1,varidentx1,
     1    varparamx1,jsurfx1,jalbx1,jxscx1,jtanx1,jprex1,jradx1,jloggx1,
     2    jfracx1,nxx1,xnx1,stx1)
       endif

C     Read in measurement vector, obs. geometry and covariances
      call readnextspavX(lspec,iform,woff,xlat,xlon,ngeom,nav,ny,y,se,
     1  fwhm,nconv,vconv,angles,wgeom,flat,flon)

      if(ionpeel.eq.1.and.ngeom.gt.1)then
         if(idiag.gt.0)then
          print*,'Error in NemesisL - Onion-peeling method flag is set'
          print*,'but there is more than one geometry in .spx file'
         endif
      endif


C     Read in forward modelling errors
      call forwarderr(ename,ngeom,nconv,vconv,woff,rerr)

C     Add forward errors to measurement covariances
      k=0
      DO i=1,ngeom
       do j=1,nconv(i)
        k = k+1
        xerr=rerr(i,j)
        if(percbool.eqv..true.)xerr=rerr(i,j)*(y(j)/100)
        se(k)=se(k)+xerr**2
       enddo
      ENDDO

      if(ilbl.eq.0)then
C      Calculate the tabulated wavelengths of c-k look up tables
       do igeom=1,ngeom
        nconv1 = nconv(igeom)
        if(igeom.eq.1.and.idiag.gt.0)print*,nconv1
        do j=1,nconv1
         vconv1(j)=vconv(igeom,j)
         if(igeom.eq.1.and.idiag.gt.0)print*,vconv1(j)
        enddo
        CALL wavesetb(runname,vkstart,vkend,vkstep,nconv1,vconv1,fwhm,
     1   nwave1,vwave1)
        if(igeom.eq.1.and.idiag.gt.0)print*,nwave1
        do j=1,nwave1
         vwave(igeom,j)=vwave1(j)
         if(igeom.eq.1.and.idiag.gt.0)print*,vwave1(j)
        enddo
        nwave(igeom)=nwave1
       enddo
      endif

      if(ilbl.eq.2)then
       call file(runname,sfile,'sha')
        open(13,file=sfile,status='old')
        READ(13,*)ISHAPE
       close(13)

C      Calculate the tabulated wavelengths of lbl look up tables
       do igeom=1,ngeom
        nconv1=nconv(igeom)
        do j=1,nconv1
         vconv1(j)=vconv(igeom,j)
        enddo
        CALL wavesetc(runname,vkstart,vkend,vkstep,nconv1,vconv1,
     1  fwhm,ishape,nwave1,vwave1)
        do j=1,nwave1
         vwave(igeom,j)=vwave1(j)
        enddo
        nwave(igeom)=nwave1
       enddo
      endif

C     set up a priori of x and its covariance
      CALL readapriori(runname,lin,lpre,xlat,npro,nvar,varident,
     1  varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,jfrac,nx,xa,
     2  sa,lx)

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

      ilbl1=ilbl

      call coreretL(runname,ispace,iscat,ilbl,ica,kiter,phlimit,
     1  inum,fwhm,xlat,xlon,ngeom,nav,nwave,vwave,nconv,vconv,angles,
     2  npro,gasgiant,lin,lpre,nvar,varident,varparam,jsurf,jalb,jxsc,
     3  jtan,jpre,jrad,jlogg,jfrac,occult,ionpeel,wgeom,flat,flon,nx,
     4  lx,xa,sa,ny,y,se,xn,sm,sn,st,yn,kk,aa,dd)

C     Calculate retrieval errors.
C     Simple errors, set to sqrt of diagonal of ST
      do i=1,nx
       err1(i)=sqrt(abs(st(i,i)))
      enddo

C     write output

      if(occult.eq.2)iform=5
      if(occult.eq.3)iform=5
      if(occult.eq.4)iform=6
      CALL writeout(iform,runname,ispace,lout,ispec,xlat,xlon,npro,
     1 nvar,varident,varparam,nx,ny,y,yn,se,xa,sa,xn,err1,ngeom,
     2 nconv,vconv,gasgiant,jpre,jrad,jlogg,jfrac,iscat,lin)

      CALL writeoutp(iform,runname,ispace,loutp,ispec,xlat,xlon,npro,
     1 nvar,varident,varparam,nx,ny,y,yn,se,xa,sa,xn,err1,ngeom,
     2 nconv,vconv,gasgiant,jpre,jrad,jlogg,jfrac,iscat,lin)


      if(lin.eq.4)then
       CALL writerawx(lraw,ispec,xlat,xlon,npro,nvar,varident,
     1   varparam,nx,xn,st,nvarx1,varidentx1,varparamx1,nxx1,xnx1,stx1)
      else
       CALL writeraw(lraw,ispec,xlat,xlon,npro,nvar,varident,
     1   varparam,nx,xn,st)
      endif

      if(ica.eq.1)then
C       Write out all the error matrices if only one case retrieved
        call write_covariance(runname,npro,nvar,varident,varparam,
     1    nx,ny,sa,sm,sn,st,se,aa,dd,kk)
      endif

2999  continue

      close(lspec)
      close(lout)
      close(lraw)
      close(lpre)
      if(lin.eq.4)then
       close(lprx)
      endif

C     New compiler time
      call system_clock(time2)
      tot_time=(time2-time1)/rate

      if(idiag.gt.0)write(6,*)'Model run OK'
      if(idiag.gt.0)write(6,244)tot_time
244   FORMAT(' Elapsed time (s) = ',F8.1)


      END




