************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C_TITLE:		PROGRAM CIRSDRV_WAVE
C
C
C_ARGS: Input Variables
C	DELV:REAL	The resolution of the k-table to be used.
C       DIST:REAL	Distance from Sun (as prompted by CIRSDRV_WAVE) in
C			units of AU.
C	FWHM:REAL	Full-Width-at-Half-Max
C	INTMOD:INT	Wavenumber interpolation mode
C	I:INT		Incrementor/Loop Counter for
C       INORMAL:INT     flag for ortho:para ratio (0=equilibrium =>1:1) 
C			(1=normal =>3:1)
C       ITYPE:INT	Value designating the chosen scattering routine
C	J:INT		Incrementor/Loop counter for
C	K:INT		Incrementor/Loop counter for
C	NCONV:INT
C	NCONVMIN:INT	Minimum value of returned
C	NCONVMAX:INT	Maximum value of returned
C       NWAVE:INT	The number of calculation wavenumbers.
C	NWAVEMIN:INT	Minimum value of returned calculated wavenumbers.
C	NWAVEMAX:INT	Maximum value of returned calculated wavenumbers.
C       NPATH:INT	Number of individual paths through the layers.  
C	NPOINTS:INT	The number of points for calculation puroses as
C			determined by ...
C			NPoints= ((WNumTop-WNumBot) / Resolution) + 1
C	OPFILE:CHARA*100	Operation filename
C	OUTFILE:CHARA*100	
C	PLANET:INT	The planet
C	RESOLUTION:REAL	
C	SPEC:REAL
C	USPEC:REAL
C	VCONV:REAL
C	VCONVMIN:INT	Minimum value of returned
C	VCONVMAX:INT	Maximum value of returned
C       VWAVE:REAL	Bin centres in wavenumber space.
C	VWAVEMIN:INT	Minimum value of returned
C	VWAVEMAX:INT	Maximum value of returned
C	WNUMBOT:REAL	Minimum wavenumber in selected range
C	WNUMTOP:REAL	Maximum wavenumber in selected range
C
C       Output variables
C       OUTPUT          Output values at each wavenumber for each output
C                       type for each path
C
C_DESC: 
C
C_HIST: 
C-----------------------------------------------------------------------

      PROGRAM CIRSdrv_wavePY

      IMPLICIT NONE

C Defines the maximum values for a series of variables (layers, bins,
C paths, etc.)
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../../nemesis/arraylen.f'

      integer :: intmod, nconv, nwave, I, npath, J, K
      integer :: itype, ispace, idump, ilbl
      integer :: NPoints, Planet, INormal, ichannel,nem,nphi
      integer ::IRAY,IPTF,IMIE,IMIE1,ishape,lcdr,isol,layint
      integer :: ngeom,ny,nmu,lowbc,nf,nlayer,laytyp,layht

      real :: xlat,xlon,fwhm,galb,woff,vkstart,vkend,vkstep
      integer :: nconv1(mgeom),nav(mgeom)
      real :: y(my),se(my),vconv1(mgeom,mconv),angles(mgeom,mav,3)
      real :: wgeom(mgeom,mav),flat(mgeom,mav),flon(mgeom,mav)
      real :: vconv(mconv)
      double precision mu(maxmu),wtmu(maxmu)
      logical :: gasgiant

      real ::vwave(maxbin), spec(maxout3),uspec(maxout3)
      real :: VREF,DELVK,VMIN,VMAX
      real :: Dist, WNumBot, WNumTop, Resolution, DelV
      real :: VWaveMin, VWaveMax, USpecMin, USpecMax
      real :: VConvMin, VConvMax, SpecMin, SpecMax

      REAL		tsurf,vem(MAXSEC),emissivity(MAXSEC)
      real tot_time
      real rate
      integer c1,c2,cr,time1,time2,cm
      CHARACTER*80	TEXT
      character (len=100) :: opfile, patfile, outfile, idlfile, ciafil
      character (len=100) :: runname,runname1,confil,sfile
      CHARACTER*1 ANS

      REAL xwave(MAXSEC),xf(MAXCON,MAXSEC),xg1(MAXCON,MAXSEC)
      REAL xg2(MAXCON,MAXSEC)
      REAL tnco,twave,frac,tico
      COMMON /hgphas/xwave,xf,xg1,xg2,tnco,twave,frac,tico
      COMMON /SCATDUMP/ IDUMP
      COMMON /IMIESCAT/ IMIE1
      COMMON /LBLTABLE/ ILBL

      INCLUDE '../includes/ciacom.f'
      INCLUDE '../includes/gascom.f'
      INCLUDE '../includes/planrad.f'

      CHARACTER*100 ANAME
      REAL DNU,DV
      INTEGER IPARA     

C     Solar spectrum variables and flags      
      integer iform,iread,solnpt
      real solwave(maxbin),solrad(maxbin),solradius
      character*100 solfile,solname
      logical solexist
      common/solardat/iread, iform, solradius, solwave, solrad,  solnpt



C-----------------------------------------------------------------------
C
C	Write a header
C
C-----------------------------------------------------------------------
      jradf=-1
      jloggf=-1


      CALL system_clock(count_rate=cr)
      CALL system_clock(count_max=cm)
      rate = REAL(cr)
      call system_clock(time1)
      xwave(1)=-1                       ! Reset to force read of hgphase*
C                                         files for scloud11-13

      WRITE(*,*)'           '
      WRITE(*,*)'           WELCOME TO CIRSDRV_WAVE'
      WRITE(*,*)' '


C     Read in reference gas information data
      CALL RESERVEGAS

      IDUMP=0
C-----------------------------------------------------------------------
C
C	Prompt the user for prescribed variables.
C
C-----------------------------------------------------------------------

      WRITE(*,*)'Give operation filename: '
      READ(*,1000)OpFile
1     FORMAT(A)

C-----------------------------------------------------------------------
C
C       Reading all required parameters by CIRSrtf_wave from .cdr file
C
C-----------------------------------------------------------------------

      lcdr=37
      CALL file(opfile,runname1,'cdr')
      open(lcdr,file=runname1,status='old')
      print*,runname1
      read(lcdr,*)dist
      read(lcdr,*)fwhm
      read(lcdr,*)ispace
      read(lcdr,*)ilbl
c      read(lcdr,*)nwave
c      do i=1,nwave
c        read(lcdr,*)vwave(i)
c      enddo
      read(lcdr,*)npath
      read(lcdr,*)nconv
      do i=1,nconv
        read(lcdr,*)vconv(i)
      enddo
      read(lcdr,*)nem
      do i=1,nem
        read(lcdr,*)vem(i),emissivity(i)
      enddo
      read(lcdr,*)tsurf
      close(lcdr)

      itype = 11    !scloud11wave
      
      print*,'nconv = ',nconv



C-----------------------------------------------------------------------
C
C       Calculating calculation wavenumbers/wavelengths
C
C-----------------------------------------------------------------------

      if(ilbl.eq.2)then  
        CALL readkklblhead(opfile,vkstart,vkend,vkstep)
        if(fwhm.gt.0.0)then  !if it is lower then we read .fil file
         call file(opfile,sfile,'sha')
         open(13,file=sfile,status='old')
         READ(13,*)ISHAPE
         close(13)
        endif
        print*,'fwhm = ',fwhm
        print*,'ishape = ',ishape
        call wavesetc(opfile,vkstart,vkend,vkstep,nconv,vconv,
     1   fwhm,ishape,nwave,vwave)
        print*,'nwave = ',nwave
c        do i=1,nwave
c         print*,vwave(i)
c        enddo 
c        pause
      endif


      if(ilbl.eq.0)then  
        CALL readkkhead(opfile,vkstart,vkend,vkstep)
        print*,'fwhm = ',fwhm
        call wavesetb(opfile,vkstart,vkend,vkstep,nconv,vconv,
     1   fwhm,nwave,vwave)             
        print*,'nwave = ',nwave
      endif




C-----------------------------------------------------------------------
C
C       Get scattering type and Ortho/Para information.
C
C-----------------------------------------------------------------------

      CALL FILE(OPFILE,CIAFIL,'cia')

      OPEN(12,FILE=CIAFIL,STATUS='OLD')
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

      CALL READFLAGS(OPFILE,INORMAL,IRAY,IH2O,ICH4,IO3,INH3,
     1 IPTF,IMIE,IUVSCAT)
      IMIE1=IMIE

C     See if there is a solar or stellar reference spectrum and read in
C     if present.  
      call file(opfile,solfile,'sol')
      inquire(file=solfile,exist=solexist)    
      if(solexist)then
         call opensol(solfile,solname)
         CALL init_solar_wave(ispace,solname)
      endif
      iform=0

C-----------------------------------------------------------------------
C
C	Call CIRSrtf_wave
C
C-----------------------------------------------------------------------
	
      WRITE(*,*)' CALLING CIRSrtf_wave'
      CALL cirsrtf_wave (opfile, Dist, INormal, Iray,FWHM, 
     1 ispace,vwave, nwave,
     2 npath, uspec, vconv, nconv,itype,nem,vem,
     3 emissivity,tsurf,spec)

      WRITE(*,*)' CIRSrtf_wave COMPLETE'
      WRITE(*,*)' '

C-----------------------------------------------------------------------
C
C	Find the minimum and maximum limits of the returned values/arrays.
C
C	Create and write to output files (both an .out and an .idl file)
C	The .out file is for comparison between the non-convolved and
C	  convolved spectrums.
C	The .idl file only contains the convolved output and should be
C	  used in comparison with Radtrans.
C
C-----------------------------------------------------------------------

      VWaveMin= WNumBot
      VWaveMax= WNumTop
      USpecMax= USpec(1)
      USpecMin= USpecMax

      VConvMin= WNumBot
      VConvMax= WNumTop
      SpecMax= Spec(1)
      SpecMin= SpecMax

      DO 100 i=1,NPoints
        USpecMin= MIN(USpecMin,USpec(i))
        USpecMax= MAX(USpecMax,USpec(i))

        SpecMin= MIN(SpecMin,Spec(i))
        SpecMax= MAX(SpecMax,Spec(i))
100   CONTINUE


      WRITE(*,*)' ********* Writing ASCII output *********'
      WRITE(*,*)' '

      CALL file (opfile, outfile, 'out')
      OPEN (UNIT= 2, FILE= outfile, STATUS= 'unknown')
      WRITE(2,*)npath,'     !! npath'
      CALL file (opfile, idlfile, 'idl')
      OPEN (UNIT= 3, FILE= idlfile, STATUS= 'unknown')
      WRITE(3,*)npath,'     !! npath'

      DO I= 1, npath
        WRITE(2,*)nwave,'     !! nwave'
        WRITE(2,*)'XLIM, YLIM :'
        WRITE(2,*)VWaveMin, VWaveMax, USpecMin, USpecMax
        WRITE(2,*)'XLABEL, YLABEL :'
        WRITE(2,*)'Wavenumbers [cm-1]'
        WRITE(2,*)'Radiance []'
        DO J=  1, nwave
          K= I + (J-1)*npath
          WRITE(2,*)vwave(J), uspec(K)
        ENDDO

        WRITE(2,*)nconv,'     !! nconv'
        WRITE(3,*)nconv,'     !! nconv'
        WRITE(2,*)' XLIM, YLIM :'
        WRITE(3,*)' XLIM, YLIM :'
        WRITE(2,*)VConvMin, VConvMax, SpecMin, SpecMax
        WRITE(3,*)VConvMin, VConvMax, SpecMin, SpecMax
        WRITE(2,*)' XLABEL, YLABEL :'
        WRITE(3,*)' XLABEL, YLABEL :'
        WRITE(2,*)'Wavenumbers [cm-1]'
        WRITE(3,*)'Wavenumbers [cm-1]'
        IF(ispace.eq.0)then
         WRITE(2,*)'Radiance []'
         WRITE(3,*)'Radiance []'
        ELSE
         WRITE(2,*)'Radiance []'
         WRITE(3,*)'Radiance []'
        ENDIF
        DO J= 1, nconv
          K= I + (J-1)*npath
          WRITE(2,*)vconv(J), spec(K)
          WRITE(3,*)vconv(J), spec(K)
        ENDDO
      ENDDO

      WRITE(*,*)'%CIRSDRV_WAVE.f :: calculation complete'
      call system_clock(time2)
      tot_time=(time2-time1)/rate
      WRITE(*,200)tot_time

C-----------------------------------------------------------------------
C
C	Format statements. Close file(s) [2= .out, 3= .idl files]
C	and End.
C
C-----------------------------------------------------------------------

        CLOSE (2)
        CLOSE (3)

200     FORMAT(' Elapsed time (sec)= ',F8.1)
1000    FORMAT (A60)

        END
 
************************************************************************
************************************************************************
