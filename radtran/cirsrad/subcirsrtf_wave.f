************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C_TITLE:		SUBROUTINE SUBCIRSRTF_WAVE
C
C_ARGS: Input Variables
C       
C       DELV:REAL     The resolution of the k-table to be used.
C       ispace           integer Indicates if wavelengths in vconv and 
C                               vwave are in wavenumbers(0) or
C                               wavelengths (1)
C       DIST:REAL     Distance from Sun (as prompted by CIRSDRV) in
C                       units of AU.
C       FWHM:REAL     Full-Width-at-Half-Max
C       INTMOD:INT      Wavenumber interpolation mode
C       I:INT           Incrementor/Loop Counter for
C       INORMAL:INT     flag for ortho:para ratio (0=equilibrium =>1:1)
C                       (1=normal =>3:1)
C       ITYPE:INT       Value designating the chosen scattering routine
C                       (currently only scloud8 through scloud11).
C       J:INT           Incrementor/Loop counter for
C       K:INT           Incrementor/Loop counter for
C       NCONV:INT
C       NWAVE:INT       The number of calculation wavenumbers.
C       NPATH:INT       Number of individual paths through the layers.
C       NPOINTS:INT     The number of points for calculation puroses as
C                       determined by ...
C                       NPoints= ((WNumTop-WNumBot) / Resolution) + 1 
C       OPFILE:CHARA*100        Operation filename
C       RESOLUTION:REAL
C       SPEC:REAL
C       USPEC:REAL
C       VCONV:REAL
C       VWAVE:REAL    Bin centres in wavenumber space.
C       WNUMBOT:REAL  Minimum wavenumber in selected range
C       WNUMTOP:REAL  Maximum wavenumber in selected range
C       
C
C       Output variables
C                       
C       OUTPUT          Output values at each wavenumber for each output
C                       type for each path
C
C_DESC: 
C
C_HIST: 
C-----------------------------------------------------------------------

	SUBROUTINE subcirsrtf_wave(opfile1, Dist, INormal, Iray, 
     1      ispace, vwave,nwave,npath1,itype1,nem,vem,emissivity,
     2      tsurf, output)

	IMPLICIT NONE

C	The variables defined here are those normally used when
C	calculating atmospheric paths. ! nptah. itype declared within
	INCLUDE '../includes/arrdef.f'
	INCLUDE '../includes/pathcom.f'
C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)

        INTEGER         lun, ulog, iphi, itype1
	PARAMETER       (lun=2, ulog=17)

	INTEGER		nwave, I, J, npath1, nout, ispace

	INTEGER		INormal,Iray,nalb,nem
	REAL		Dist,valb(maxsec),alb(maxsec),tsurf
        REAL		vem(maxsec),emissivity(maxsec)
	REAL		vwave(nwave), output(maxout3)
	LOGICAL		scatter, dust, solexist
        double precision mu1(maxmu), wt1(maxmu), galb

	CHARACTER*100	klist
	CHARACTER*100	opfile1,solfile,solname

	CHARACTER*100	logfil, drvfil, radfile, xscfil,albfile

        common/scatd/mu1, wt1, galb
        common/alb/nalb,valb,alb
C-----------------------------------------------------------------------
C
C	Begin Program.
C
C-----------------------------------------------------------------------

	opfile= opfile1		! Renamed and assigned here because
	iphi= itype1		! there was a conflict with similiar
				! declarations in pathcom.f which is
				! needed for other bits of this code.

C-----------------------------------------------------------------------
C
C	Run path file to create driver file.
C
C-----------------------------------------------------------------------

	CALL subpath(opfile)

	npath1= npath

C-----------------------------------------------------------------------
C
C	Open drive file and read in data.
C
C       rdlbld: Reads in arrays from LUN= 1 (i.e. DRVFIL) such as VMin,
C       DelV, DelH, NPoint, FWHM, Wing, VRel, and KEYFIL.
C
C-----------------------------------------------------------------------

	CALL file (opfile, drvfil, 'drv')
	OPEN (UNIT= 1, FILE= DRVFIL, STATUS= 'OLD')
	CALL RDLBLD
	CLOSE(1)

C-----------------------------------------------------------------------
C
C	Check sizes
C
C       NWAVE is the number of calculation wavelengths. This routine
C       used to use NPOINT, which is actually the number of points as
C       defined in the .pat file, which is irrelevant for NIMSRAD.
C
C       It would be nice to check that NCONV (the number of convolved
C       output wavelengths) is also within range since NCONV is usually 
C       larger than NWAVE. However, NCONV is not passed to this routine.
C
C-----------------------------------------------------------------------

     	NOUT=NPATH*NWAVE
        WRITE(*,*)'NPATH,NWAVE,NOUT = ',NPATH,NWAVE,NOUT

      	IF (NOUT.GT.MAXOUT3) THEN
        	WRITE(*,1080)
        	STOP
      	ENDIF

C-----------------------------------------------------------------------
C
C	Read the ktables.
C
C-----------------------------------------------------------------------

	CALL file (opfile, klist, 'kls')
	WRITE(*,1050)klist

	WRITE(*,*)'     CALLING read_klist'
	CALL read_klist (klist, ngas, idgas, isogas, nwave, vwave)
	WRITE(*,*)'     read_klist COMPLETE'
	WRITE(*,*)' '

C-----------------------------------------------------------------------
C
C	Now read the scattering files if required.
C
C-----------------------------------------------------------------------

	scatter= .false.
	DO I= 1, npath
		IF (imod(I).EQ.15.OR.imod(I).EQ.16) scatter= .true.
                IF (imod(I).EQ.21.OR.imod(I).eq.22) scatter = .true.
	ENDDO

	IF (scatter) THEN
		CALL file(opfile, radfile, 'sca')
		WRITE(*,*)'     CALLING get_scatter'
		CALL get_scatter(radfile,ncont)
		WRITE(*,*)'     get_scatter COMPLETE'
		WRITE(*,*)' '
                if(galb.lt.0.0)then
                 call file(opfile,albfile,'alb')
                 open(9,file=albfile,status='old')
                 read(9,*)nalb
                 print*,'Reading albedo file'
                 if(nalb.gt.maxsec)then
                  print*,'Error in subcirsrtf_wave nalb > maxsec'
                  print*,nalb,maxsec
                  stop
                 endif
                 do i=1,nalb
                  read(9,*)valb(i),alb(i)
                 enddo
                 close(9)
                endif
	ENDIF

C-----------------------------------------------------------------------
C
C	And the xsc files likewise
C
C-----------------------------------------------------------------------

	dust= .false.
	IF (ncont.gt.0) dust= .true.
	IF (dust) THEN
		CALL file(opfile, xscfil, 'xsc')
		WRITE(*,*)'Subcirsrtf_wave: CALLING get_xsec. ncont = ', ncont
		CALL get_xsec(xscfil, ncont)
		WRITE(*,*)'Subcirsrtf_wave: get_xsec COMPLETE'
c		WRITE(*,*)' '
		
	ENDIF
 
        call file(opfile1,solfile,'sol')
     
        inquire(file=solfile,exist=solexist)
        print*,'solexist = ',solexist

        if(solexist)then
         print*,'Initiating solar reference file'
         call opensol(solfile,solname)
         CALL init_solar_wave(ispace,solname)
        endif

C-----------------------------------------------------------------------
C
C	Call CIRSrad_wave.	
C
C-----------------------------------------------------------------------

        print*,'Calling cirsrad_wave'
	CALL cirsrad_wave (Dist, INormal, Iray, ispace, DelH, nlayer, 
     1    npath,ngas, maxlay, maxcon, press, temp, pp, amount, nwave,
     2    vwave, nlayin, maxinc, layinc, cont, scale, imod, idgas,
     3    isogas,emtemp,iphi,nem,vem,emissivity,tsurf,
     4    flagh2p,hfp,flagc, hfc, ifc, basep, baseh, output)
	WRITE(*,*)'Subcirsrtf_wave: cirsrad_wave COMPLETE'
c	WRITE(*,*)' '

C-----------------------------------------------------------------------
C
C	Wrap up.
C
C-----------------------------------------------------------------------

	WRITE(*,*)'     CALLING close_scat'
	CALL close_scat(ncont)
	WRITE(*,*)'     close_scat COMPLETE'
	WRITE(*,*)' '

C-----------------------------------------------------------------------
C
C	Formats, return, and end.
C
C-----------------------------------------------------------------------

1050	FORMAT (/, 6x, 'Klist filename: ', A)
1080    FORMAT (' subCIRSrtf_wave: output arrays too small - recompile')

	RETURN

	END

************************************************************************
************************************************************************
