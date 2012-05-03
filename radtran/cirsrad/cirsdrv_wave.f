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

      PROGRAM CIRSdrv_wave

      IMPLICIT NONE

C Defines the maximum values for a series of variables (layers, bins,
C paths, etc.)
      INCLUDE '../includes/arrdef.f'

      INTEGER		intmod, nconv, nwave, I, npath, J, K,
     1			itype, ispace, idump
      INTEGER		NPoints, Planet, INormal, ichannel,nem
      INTEGER		IRAY,IPTF
      REAL		vconv(maxbin), vwave(maxbin), spec(maxout3),
     1                  uspec(maxout3)
      REAL		VREF,DELVK,VMIN,VMAX
      REAL		Dist, WNumBot, WNumTop, Resolution, DelV, FWHM,
     1			VWaveMin, VWaveMax, USpecMin, USpecMax,
     2			VConvMin, VConvMax, SpecMin, SpecMax

      REAL		TOT_TIME,tsurf,vem(MAXSEC),emissivity(MAXSEC)
      DOUBLE PRECISION	TIME,TIME1,TIME2

      CHARACTER*80	TEXT
      CHARACTER*100	opfile, patfile, outfile, idlfile, ciafil
      CHARACTER*100     confil
      CHARACTER*1 ANS

      REAL xwave(MAXSEC),xf(MAXCON,MAXSEC),xg1(MAXCON,MAXSEC)
      REAL xg2(MAXCON,MAXSEC)
      REAL tnco,twave,frac,tico
      COMMON /hgphas/xwave,xf,xg1,xg2,tnco,twave,frac,tico
      COMMON /SCATDUMP/ IDUMP


      INCLUDE '../includes/ciacom.f'
      INCLUDE '../includes/gascom.f'

      CHARACTER*100 ANAME
      REAL DNU
      INTEGER IPARA     


C-----------------------------------------------------------------------
C
C	Write a header
C
C-----------------------------------------------------------------------

      CALL GETTIME(TIME)
      TIME1= TIME
      xwave(1)=-1                       ! Reset to force read of hgphase*
C                                         files for scloud11-13

      WRITE(*,*)'           '
      WRITE(*,*)'           WELCOME TO CIRSDRV_WAVE'
      WRITE(*,*)' '

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
      print*,'N.B. If you want to use channel-integrated k-tables'
      print*,'you need to answer no here and input the required'
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
C       Two wavenumber arrays are considered; one for final output
C       (VConv), and one for actual use in calculating outputs (VWave).
C
C	NCONV	the number of data points for the prompted wavenumber
C		range and resolution.
C	VCONV	an array containing each wavenumber value for use when
C		convolving.
C
C	NWAVE	the number of points in the k-table that cover the
C		prompted wavenumber range.
C	VWAVE	an array of the points within this wavenumber range.
C
C	There are two ways of producing NWAVE & VWAVE: fast and slow.
C	In the slow mode, NWAVE & VWAVE are the same as the convolution
C	wavenumbers. In the fast mode, they are the same wavenumbers as
C	the K tables are calculated at.
C
C	Accuracy of methods is approximately the same.
C
C-----------------------------------------------------------------------

      if(ichannel.eq.1)then
       intmod=1
      else
30     WRITE(*,*)'Enter Wavenumber interpolation mode:'
       WRITE(*,*)'    (1) Fast post-calculation'
       WRITE(*,*)'    (2) Slow K-coeff regridding'
 
       WRITE(*,*)'Option: '
       READ(*,*)IntMod
       IF ((IntMod.LT.1).OR.(IntMod.GT.2)) GOTO 30
       WRITE(*,*)' '
       IF (IntMod.EQ.1) WRITE(*,*)'-> Fast post-calculation selected'
       IF (IntMod.EQ.2) WRITE(*,*)'-> Slow K-coeff regridding selected'
       WRITE(*,*)' '
      endif

      IF (IntMod.EQ.1) THEN
       if(ichannel.eq.0)then
        WRITE(*,*)'Enter reference wavenumber and step in ktables'
        READ (*,*) VREF,DelVK

        VMIN = VREF+DELVK*INT((WNumBot-0.5*FWHM-VREF)/DELVK)
        VMAX = VREF+DELVK*(NINT((WNumTop+0.5*FWHM-VREF)/DELVK)+1)


C Calculate number of points in tabulated wavenumber array plus one to
C account for the end point if FWHM <> 0.0

        NWave = 1 + NINT((VMAX-VMIN) / DelVK)

        IF(FWHM.EQ.0.0)Nwave=Nwave-1

        DO j= 1, NWave
          VWave(j)= VMIN + ((j-1) * DelVK)
        ENDDO
       else
        nwave=nconv
        DO j= 1, NWave
          VWave(j)= Vconv(j)
        ENDDO
        FWHM = 0.0
       endif        
      ELSE
        NWave= NConv + 2
        VWave(1)= VConv(1) - DELV
        VWave(NWave)= VConv(NConv) + DELV
        DO i= 1, NWave
          VWave(i+1)= VConv(i)
        ENDDO
      ENDIF

      PRINT*,'Calculation wavelengths: ',NWave
      DO I=1,NWave
        PRINT*,I,VWave(I)
      ENDDO

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
      IF((IType.lt.11).or.(IType.GT.13))THEN
        WRITE(*,*)' Parameter must lie in allowed range (11-13)'
        GOTO 40
      ENDIF

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

      CALL READFLAGS(OPFILE,INORMAL,IRAY,IH2O,ICH4,IPTF)



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
        WRITE(2,*)'Radiance [W cm-2 sr-2 (cm-1)-1]'
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
         WRITE(2,*)'Radiance [W cm-2 sr-2 (cm-1)-1]'
         WRITE(3,*)'Radiance [W cm-2 sr-2 (cm-1)-1]'
        ELSE
         WRITE(2,*)'Radiance [W cm-2 sr-2 micron-1]'
         WRITE(3,*)'Radiance [W cm-2 sr-2 micron-1]'
        ENDIF
        DO J= 1, nconv
          K= I + (J-1)*npath
          WRITE(2,*)vconv(J), spec(K)
          WRITE(3,*)vconv(J), spec(K)
 	ENDDO
      ENDDO

      WRITE(*,*)'%CIRSDRV_WAVE.f :: calculation complete'
      CALL GETTIME(TIME)
      TIME2= TIME
      TOT_TIME=TIME2-TIME1
      WRITE(*,200)TOT_TIME

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
