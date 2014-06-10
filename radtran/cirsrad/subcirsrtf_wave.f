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
        INCLUDE '../includes/laycom.f'
C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)

        INTEGER         lun, ulog, iphi, itype1
	PARAMETER       (lun=2, ulog=17)

	INTEGER		nwave, I, J, npath1, nout, ispace

	INTEGER		INormal,Iray,nem
	REAL		Dist,tsurf,radius1,zheight(maxpro)
        REAL		vem(maxsec),emissivity(maxsec)
	REAL		vwave(nwave), output(maxout3),radextra
	LOGICAL		fscatter, fdust, solexist,fexist
        double precision mu1(maxmu), wt1(maxmu), galb

	CHARACTER*100	klist,buffer
	CHARACTER*100	opfile1,solfile,solname

	CHARACTER*100	logfil, drvfil, radfile, xscfil,albfile

        common/scatd/mu1, wt1, galb
C       Need simple way of passing planetary radius to nemesis
        INCLUDE '../includes/planrad.f'

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
                PRINT*,'MAXOUT3 = ',MAXOUT3
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

	fscatter= .false.
	DO I= 1, npath
		IF (imod(I).EQ.15.OR.imod(I).EQ.16) fscatter= .true.
                IF (imod(I).EQ.22.OR.imod(I).eq.23) fscatter = .true.
                IF (imod(I).EQ.24) fscatter = .true.
	ENDDO

	IF (fscatter) THEN
		CALL file(opfile, radfile, 'sca')
		WRITE(*,*)'     CALLING get_scatter'
		CALL get_scatter(radfile,ncont)
		WRITE(*,*)'     get_scatter COMPLETE'
	ENDIF

C-----------------------------------------------------------------------
C
C	And the xsc files likewise
C
C-----------------------------------------------------------------------

	fdust= .false.
	IF (ncont.gt.0) fdust= .true.
	IF (fdust) THEN
		CALL file(opfile, xscfil, 'xsc')
		WRITE(*,*)'Subcirsrtf_wave: CALLING get_xsec. ncont = ', ncont
		CALL get_xsec(xscfil, ncont)
		WRITE(*,*)'Subcirsrtf_wave: get_xsec COMPLETE'
c		WRITE(*,*)' '
		
	ENDIF
 
C-----------------------------------------------------------------------
C
C	Call CIRSrad_wave.	
C
C-----------------------------------------------------------------------

C       Pass radius of planet to cirsradg
C       radius2 is radius held in planrad common block. Pass this to
C       cirsrad_wave in case it's been updated.      
        radius1=radius2

C       Look to see if ipzen has been set to 2 and if so read in altitude of
C       top of atmosphere
        radextra=0.
        if(ipzen.eq.2)then
C          need the current height profile
           call readprfheight(opfile1,npro,zheight)
           radextra=zheight(npro)
        endif

        print*,'Calling cirsrad_wave'
	CALL cirsrad_wave (Dist, INormal, Iray, ispace, DelH, nlayer, 
     1    npath,ngas, maxlay, maxcon, totam, press, temp, pp, amount,
     2    nwave, vwave, nlayin, maxinc, layinc, cont, scale, imod, 
     3    idgas, isogas,emtemp,iphi,nem,vem,emissivity,tsurf,
     4    flagh2p,hfp,flagc, hfc, ifc, basep, baseh, RADIUS1,
     5    radextra,output)
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
