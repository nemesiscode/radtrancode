************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C_TITLE:		PROGRAM LBLDRV_WAVE
C
C
C_ARGS: Input Variables
C	DELV:REAL	The resolution of the k-table to be used.
C       DIST:REAL	Distance from Sun (as prompted by lbldrv_WAVE) in
C			units of AU.
C	FWHM:REAL	Full-Width-at-Half-Max
C	INTMOD:INT	Wavenumber interpolation mode
C	I:INT		Incrementor/Loop Counter for
C       INORMAL:INT     flag for ortho:para ratio (0=equilibrium =>1:1) 
C			(1=normal =>3:1)
C       ITYPE:INT	Value designating the chosen scattering routine
C			(currently only scloud8 through scloud11).
C	J:INT		Incrementor/Loop counter for
C	K:INT		Incrementor/Loop counter for
C	NCONV:INT
C	NCONVMIN:INT	Minimum value of returned
C	NCONVMAX:INT	Maximum value of returned
C       NWAVE:INT	The number of calculation wavenumbers.
C	NWAVEMIN:INT	Minimum value of returned calculated wavenumbers.
C	NWAVEMAX:INT	Maximum value of returned calculated wavenumbers.
C       NPATH:INT	Number of individual paths through the layers.  
C			required by scattering routine scloud8.
C	NPOINTS:INT	The number of points for calculation puroses as
C			determined by ...
C			NPoints= ((WNumTop-WNumBot) / Resolution) + 1
C	OPFILE:CHARA*100	Operation filename
C	OUTFILE:CHARA*100	
C	PLANET:INT	The planet
C	RESOLUTION:REAL	
C	SPEC:REAL
C	VCONV:REAL
C	VCONVMIN:INT	Minimum value of returned
C	VCONVMAX:INT	Maximum value of returned
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

      PROGRAM lbldrv_wave

      IMPLICIT NONE

C Defines the maximum values for a series of variables (layers, bins,
C paths, etc.)
      INCLUDE '../includes/arrdef.f'

      INTEGER		intmod, nconv, nwave, I, npath, J, K,
     1			itype, ispace, idump
      INTEGER		NPoints, Planet, INormal, Iray, ichannel,nem
      INTEGER		IPTF,IMIE,IMIE1
      REAL		vconv(maxbin), spec(maxout3)
      REAL		VREF,DELVK,VMIN,VMAX
      REAL		Dist, WNumBot, WNumTop, Resolution, DelV, FWHM,
     2			VConvMin, VConvMax, SpecMin, SpecMax

      REAL	x0,x1,wing,vrel,maxdv
      REAL		tsurf,vem(maxsec),emissivity(maxsec)
      real tot_time
      real rate
      integer c1,c2,cr,time1,time2,cm
      CHARACTER*80	TEXT
      CHARACTER*100	opfile, patfile, outfile, idlfile, ciafil
      CHARACTER*100     confil
      CHARACTER*1 ANS

      REAL xwave(maxsec),xf(maxcon,maxsec),xg1(maxcon,maxsec),
     1 xg2(maxcon,maxsec)
      REAL tnco,twave,frac,tico
      COMMON /hgphas/xwave,xf,xg1,xg2,tnco,twave,frac,tico
      COMMON /SCATDUMP/ IDUMP
      COMMON /IMIESCAT/IMIE1

      INCLUDE '../includes/ciacom.f'
      INCLUDE '../includes/gascom.f'
      INCLUDE '../includes/planrad.f'

      CHARACTER*100 ANAME
      REAL DNU
      INTEGER IPARA

C     Solar spectrum variables and flags      
      integer iform,iread,solnpt
      real solwave(maxbin),solrad(maxbin),solradius
      character*100 solfile,solname
      logical solexist
      common/solardat/iread, iform, solradius, solwave, solrad,  solnpt
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet


C-----------------------------------------------------------------------
C
C	Write a header
C
C-----------------------------------------------------------------------
      jradf=-1
      jloggf=-1

      idiag=1

      CALL system_clock(count_rate=cr)
      CALL system_clock(count_max=cm)
      rate = REAL(cr)
      call system_clock(time1)
      xwave(1)=-1                       ! Reset to force read of hgphase*
C                                         files for scloud10-11
      WRITE(*,*)'           '
      WRITE(*,*)'           WELCOME TO lbldrv_WAVE'
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
      CALL FILE(OpFile,patfile,'pat')
      OPEN(2,FILE=PATFILE,STATUS='OLD')
2     READ(2,1,END=9)TEXT
      CALL REMSP(TEXT)
      CALL UPCASE(TEXT)
      IF(TEXT(1:1).EQ.' ')THEN
C Skipping blank lines
        GOTO 2
      ELSE IF(TEXT(1:8).EQ.'INTERVAL')THEN
        READ(2,*)WnumBot,WnumTop,DELV,FWHM
      ENDIF
      GOTO 2

9     CONTINUE
      CLOSE(2)


C Calculate number of points in output array 
      NPoints= NINT((WNumTop-WNumBot)/DELV) + 1

      CALL PROMPT('Enter ISPACE (0=wavenumber,1=wavelength): ')
      READ*,ispace

      NConv= NPoints
      PRINT*,'Convolution wavelengths : ',NConv
      DO i= 1, NConv
        VConv(i)= WNumBot + ((i-1) * DELV)
        PRINT*,i,VConv(i)
      ENDDO
      print*,'OK (Y/N)?'
      print*,'N.B. if you answer no here you'
      print*,'will need to input the required'
      print*,'channel wavelengths/numbers'

      read(5,1)ans
      ichannel=0
      if(ans.ne.'y'.and.ans.ne.'Y')then
       call prompt('Enter nconv : ')
       read*,nconv
       print*,'Enter wavenumbers:'
       read*,(vconv(i),i=1,nconv)
       ichannel=1
       fwhm=0.0
      endif


10    WRITE(*,*)'Select a planet (1-9) to look up solar dist: '
      READ (*,*) Planet
      IF (Planet.EQ.1) Dist= 0.39	! Mercury
      IF (Planet.EQ.2) Dist= 0.72	! Venus
      IF (Planet.EQ.3) Dist= 1.00	! Earth
      IF (Planet.EQ.4) Dist= 1.50	! Mars
      IF (Planet.EQ.5) Dist= 5.20	! Jupiter
      IF (Planet.EQ.6) Dist= 9.50	! Saturn
      IF (Planet.EQ.7) Dist= 19.2	! Uranus
      IF (Planet.EQ.8) Dist= 30.1	! Neptune
      IF (Planet.EQ.9) Dist= 39.5	! Pluto
      IF (Planet.EQ.10) THEN
       print*,'Enter solar distance : '
       read*,Dist
      ENDIF
      IF ((Planet.GT.10).OR.(Planet.EQ.0)) GOTO 10



C-----------------------------------------------------------------------
C
C       Get scattering type and Ortho/Para information.
C
C-----------------------------------------------------------------------

40    WRITE(*,*)'Allowed scattering routines: '
      WRITE(*,*)'       11   scloud11wave '
      WRITE(*,*)'       12   scloud12wave '
      WRITE(*,*)'       13   scloud11flux '
      WRITE(*,*)' Enter option : '
      READ(*,*)IType
      IF((IType.LT.11).or.(IType.GT.13))THEN
        WRITE(*,*)' Parameter must lie in allowed range (11-13)'
        GOTO 40
      ENDIF

      CALL READFLAGS(OPFILE,INORMAL,IRAY,IH2O,ICH4,IO3,INH3,
     1  IPTF,IMIE, iuvscat)
      IMIE1=IMIE
      WRITE(*,*)'Enter surface temperature : '
      READ(*,*)tsurf

      nem=2
      vem(1)=-100.0
      vem(2)=1e7
      emissivity(1)=1.0
      emissivity(2)=1.0


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



C     See if there is a solar or stellar reference spectrum and read in
C     if present.  
      call file(opfile,solfile,'sol')
      inquire(file=solfile,exist=solexist)    
      if(solexist)then
         call opensol(solfile,solname)
         CALL init_solar_wave(ispace,solname)
         print*,'Solar file = ',solname
      endif
      iform=0

C-----------------------------------------------------------------------
C
C	Call lblrtf_wave
C
C-----------------------------------------------------------------------
	

      call file(opfile,opfile,'lbl')
      open(11,file=opfile,status='old')
       print*,'LBL data read in'
       read(11,*)x0,x1,delv
       print*,'Total calculation range and step (cm-1) : ',
     1   x0,x1
       read(11,*)wing,vrel,maxdv
       print*,'wing, vrel, maxdv : ',wing, vrel, maxdv
       if(maxdv.ne.vrel)then
        print*,'*Warning from lbldrv_wave.f: V_cutoff <> vrel.'
        print*,'*These should usually be set equal'
        print*,'V_cut_off = ',maxdv
        print*,'VREL = ',vrel
       endif
      close(11)

      WRITE(*,*)' CALLING lblrtf_wave'

      CALL lblrtf_wave (x0, x1, wing, vrel, maxdv,opfile, Dist, 
     1 INormal, Iray, DELV, FWHM, ispace, npath, vconv, nconv, 
     2 itype,nem,vem,emissivity,tsurf,spec,IPTF)

      WRITE(*,*)' lblrtf_wave COMPLETE'
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

      VConvMin= WNumBot
      VConvMax= WNumTop
      SpecMax= Spec(1)
      SpecMin= SpecMax

      DO 100 i=1,NPoints
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
         WRITE(2,*)'Radiance [W cm-2 sr-2 (cm-1)-1]'
         WRITE(3,*)'Radiance [W cm-2 sr-2 (cm-1)-1]'
        ELSE
         WRITE(2,*)'Radiance [W cm-2 sr-2 micron-1]'
         WRITE(3,*)'Radiance [W cm-2 sr-2 micron-1]'
        ENDIF
        DO J= 1, nconv

          K=nconv*(I-1)+J
          WRITE(2,*)vconv(J), spec(K)
          WRITE(3,*)vconv(J), spec(K)
 	ENDDO
      ENDDO

      WRITE(*,*)'%lbldrv_WAVE.f :: calculation complete'
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

200	FORMAT(' Elapsed time (sec)= ',F8.1)
1000	FORMAT (A60)

	END
		
************************************************************************
************************************************************************
