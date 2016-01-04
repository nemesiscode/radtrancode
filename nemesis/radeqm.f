      PROGRAM radeqm
C     $Id:
C     ******************************************************************
C
C     Pat Irwin	        Modified from NIMS retrieval code 21/3/00
C			Updated	4/4/01
C			Updated for continuous vmr profiles 7/10/03
C			Updated for FOV-averaging 9/2/04
C			Modifed from Nemesis 30/10/14
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

      character*100 buffer
      integer i,j,k,ioff,j1,j2
      integer iform
      integer npro,nvmr,ispace
      character*100 runname
      integer nwave,ngas,ncont
      real vwave(mwave),xmin,xmax
      real vkstart,vkend,vkstep,tstep
      integer idump,kiter
C     ********** Scattering variables **********************
      real xwave(maxsec),xf(maxcon,maxsec),xg1(maxcon,maxsec)
      real xg2(maxcon,maxsec)
      real tnco,twave,frac,tico
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


C     New compiler time
      CALL system_clock(count_rate=cr)
      CALL system_clock(count_max=cm)
      rate = REAL(cr)
      call system_clock(time1)

      CALL prompt('Enter run name : ')
      READ(5,1)buffer
1     FORMAT(a)
      runname = buffer(1:36)

      print*,'checking files'
C     Make sure input files are OK
      CALL checkfiles(runname)

      CALL readrefhead(runname,npro,nvmr,gasgiant)
      if(npro.gt.maxpro)then
       print*,'Error in Radeqm. npro > maxpro : ',npro,maxpro
      endif

      print*,'Files OK',gasgiant

      CALL file(runname,runname,'inp')
      OPEN(32,file=runname,status='old')

C     Read in whether to calculate with wavenumbers(0) or wavelength(1)
      READ(32,*)ispace

      CALL readkkhead(runname,vkstart,vkend,vkstep)

C     Read in time step and number of iterations
      READ(32,*)tstep,kiter 
      

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

      print*,'Enter wavelength/wavenumber range : '
      read(32,*)xmin,xmax

      close(32)

C     Calculate the tabulated wavelengths of c-k look up tables
      j1 = (xmin-vkstart)/vkstep
      j2 = (xmax-vkstart)/vkstep
      nwave=1+j2-j1
      print*,'nwave = ',nwave
      do j=1,nwave
       vwave(j)=vkstart+(j1+j-1)*vkstep
       print*,j,vwave(j)
      enddo

      idump=0	! flag for diagnostic print dumps


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
	
      call subradeqm(runname,ispace,nwave,vwave,gasgiant,kiter,
     1 tstep)


C     New compiler time
      call system_clock(time2)
      tot_time=(time2-time1)/rate

      write(6,*)'Model run OK'
      WRITE(6,244)tot_time
244   FORMAT(' Elapsed time (s) = ',F8.1)


      END




