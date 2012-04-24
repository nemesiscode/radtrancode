************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C_TITLE:		SUBROUTINE lblrtf_wave
C
C_ARGS: Input Variables
C
C	CONVOUT
C       DIST:REAL       Distance from Sun (as prompted by CIRSDRV) in
C                       units of AU.
C       ispace           integer Indicates if wavelengths in vconv and
C                               vwave are in wavenumbers(0) or
C                               wavelengths (1) 
C       I:INT           Incrementor/Loop Counter for
C       INORMAL:INT     flag for ortho:para ratio (0=equilibrium =>1:1)
C                       (1=normal =>3:1)
C       ITYPE:INT       Value designating the chosen scattering routine
C                       (currently only scloud8 through scloud11).
C       J:INT           Incrementor/Loop counter for
C       K:INT           Incrementor/Loop counter for
C       NCONV:INT
C       NPATH:INT       Number of individual paths through the layers.
C       OPFILE:CHARA*100	Operation filename
C	VCONV:REAL
C	VWAVE:REAL      Bin centres in wavenumber space.
C	Y
C	YOUT
C
C
C       Output variables
C       
C       OUTPUT          Output values for each output type for each path
C
C_DESC:
C
C_HIST: 
C-----------------------------------------------------------------------

        SUBROUTINE lblrtf_wave(X0,X1,WING1,VREL1,MAXDV,opfile1, Dist, 
     1          INormal, IRay, DELV1, FWHM1, ispace,
     2          npath1, vconv, nconv, itype1, nem, vem, 
     3          emissivity,tsurf, convout)

	IMPLICIT NONE

C	The variables defined here are those normally used when
C	calculating atmospheric paths. ! nptah. itype declared within
	INCLUDE '../includes/arrdef.f'
	INCLUDE '../includes/pathcom.f'

C       Defines the maximum values for a series of variables (layers,
C         bins, paths, etc.)
	INCLUDE '../includes/contdef.f'

        INTEGER         nconv, npath1, itype1, I, J, K, mconv
 	INTEGER		MAXLIN1,ISHAPE
        PARAMETER       (mconv=2000)
	INTEGER		INormal,Iray, ispace,nem,IBS(2),IBD(2), IFLAG
	REAL		Dist, FWHM1,X0,X1,WING1,VREL1,MAXDV
        REAL		VMAX,VBOT,DELV1,VV,X
        REAL            vconv(nconv), convout(maxout3),
     1                  output(maxpat), yp(maxpat,2)
        REAL            yout(maxpat,mconv),tsurf,
     2			vem(maxsec),emissivity(maxsec), ynor(maxpat,mconv)
	REAL		AAMOUNT(maxlay,maxgas)

        INTEGER         lun, ulog, iphi,NLINR
        PARAMETER       (lun=2, ulog=17)

        INTEGER         nout, IB

        INTEGER         nalb,FSTREC
        REAL            valb(maxsec),alb(maxsec)
	REAL		XNEXT,XCOM
        LOGICAL         scatter, dust,solexist
        double precision mu1(maxmu), wt1(maxmu), galb

        CHARACTER*100   klist, opfile1, solname, solfile

        CHARACTER*100    logfil, drvfil, albfile, sfile

        common/scatd/mu1, wt1, galb
        common/alb/nalb,valb,alb

   
        INCLUDE '../includes/dbcom.f'
        INCLUDE '../includes/lincomc.f'

    
C-----------------------------------------------------------------------
C
C       Begin Program.
C
C-----------------------------------------------------------------------


        opfile= opfile1         ! Renamed and assigned here because
        iphi= itype1            ! there was a conflict with similiar
                                ! declarations in pathcom.f which is
                                ! needed for other bits of this code.

C-----------------------------------------------------------------------
C
C       Run path file to create driver file.
C
C-----------------------------------------------------------------------

        CALL subpath(opfile, npath)

        npath1= npath

        print*,'lblrtf_wave called'
        print*,X0,X1,WING1,VREL1,MAXDV
        print*,opfile1
        print*,Dist,INormal, IRay,DELV1, FWHM1
        print*,ispace,npath1,nconv
        do i=1,nconv
         print*,vconv(i)
        enddo
        print*,itype1,nem
        do i=1,nem
         print*,vem(i),emissivity(i)
        enddo
        print*,tsurf


C-----------------------------------------------------------------------
C
C       Open drive file and read in data.
C
C       rdlbld: Reads in arrays from LUN= 1 (i.e. DRVFIL) such as VMin,
C       DelV, DelH, NPoint, FWHM, Wing, VRel, and KEYFIL.
C
C-----------------------------------------------------------------------

        CALL file (opfile, drvfil, 'drv')
        OPEN (UNIT= 1, FILE= DRVFIL, STATUS= 'OLD')
        CALL RDLBLD
        CLOSE(1)

        CALL RDKEY(LUN)
        CALL RDGAS
        CALL RDISO


        OPEN(UNIT=DBLUN,FILE=DBFILE,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='DIRECT',RECL=DBRECL)

C-----------------------------------------------------------------------
C
C       Check sizes
C
C-----------------------------------------------------------------------

        DELV=DELV1
        FWHM=FWHM1
        WING=WING1
        VREL=VREL1


        DO I=1,NLAYER
         DO J=1,NGAS
          AAMOUNT(I,J)=AMOUNT(I,J)
         ENDDO
        ENDDO

C      Precompute temperature coeffients of layers
        CALL CALCTEMPPARAM(NLAYER,NGAS,PRESS,TEMP,AAMOUNT,
     1 IDGAS,ISOGAS)


        NOUT=NPATH*NCONV
        WRITE(*,*)'NPATH,NCONV,NOUT = ',NPATH,NCONV,NOUT

        IF(NCONV.GT.MCONV)THEN
         PRINT*,'Error in lblrtfwave, NCONV > MCONV'
         PRINT*,NCONV,MCONV
         STOP
        ENDIF

        IF (NOUT.GT.MAXOUT3) THEN
                PRINT*,'NOUT > MAXOUT3',NOUT,MAXOUT3
                STOP
        ENDIF

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
                  print*,'Error in lblrtf_wave nalb > maxsec'
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
C       And the xsc files likewise
C
C-----------------------------------------------------------------------

        dust= .false.
        IF (ncont.gt.0) dust= .true.
        IF (dust) THEN
             CALL file(opfile, xscfil, 'xsc')
             WRITE(*,*)'lblrtf_wave: CALLING get_xsec. ncont = ', ncont
             CALL get_xsec(xscfil, ncont)
             WRITE(*,*)'lblrtf_wave: get_xsec COMPLETE'

        ENDIF


C
C-----------------------------------------------------------------------
C
C      Initialise continuum
C
       VMIN=X0
       VMAX=X1

       PRINT*,'Calling init_cont',VMIN-VREL,VMAX+VREL,WING

       CALL INIT_CONT(VMIN-VREL,VMAX+VREL,WING)

       VBOT=VBIN(1)
       PRINT*,'VBOT = ',VBOT
       IBS(1)=1
       IBS(2)=2

C      Read in 2 arrays of lines
       IB=1
       CALL FNDWAV(VMIN-VREL)
       FSTREC = DBREC

       PRINT*,'FSTREC = ',FSTREC

       MAXLIN1=MAXLIN
       print*,'lblrtf_wave : maxlin = ',MAXLIN1
       CALL LOADBUFFER(VMIN-VREL,VMAX+VREL,FSTREC,MAXLIN1,MAXBIN,IB,
     1 NGAS,IDGAS,ISOGAS,VBOT,WING,NLINR,VLIN,SLIN,ALIN,ELIN,IDLIN,
     2 SBLIN,PSHIFT,DOUBV,TDW,TDWS,LLQ,NXTREC,FSTLIN,LSTLIN)

       NLINE(IB)=NLINR
       IBD(IB)=-1

       IB=2
       FSTREC=NXTREC


       CALL LOADBUFFER(VMIN-VREL,VMAX+VREL,FSTREC,MAXLIN,MAXBIN,IB,
     1 NGAS,IDGAS,ISOGAS,VBOT,WING,NLINR,VLIN,SLIN,ALIN,ELIN,IDLIN,
     2 SBLIN,PSHIFT,DOUBV,TDW,TDWS,LLQ,NXTREC,FSTLIN,LSTLIN)

       
       NLINE(IB)=NLINR
       IBD(IB)=-1

C-----------------------------------------------------------------------
C
C       Call subroutine lblrad_wave:
C
C-----------------------------------------------------------------------

        VV = VMIN-DELV
        IFLAG=0
        DO J=1,NCONV
         DO I=1,NPATH
          YOUT(I,J)=0.0
          YNOR(I,J)=0.0
         ENDDO
        ENDDO

        call file(opfile1,solfile,'sol')

        inquire(file=solfile,exist=solexist)

        if(solexist)then
         call opensol(solfile,solname)
         CALL init_solar_wave(ispace,solname)
        endif

        XNEXT=0.0

        call file(opfile,sfile,'sha')
        open(13,file=sfile,status='old')
        READ(13,*)ISHAPE
        close(13)

        print*,'ISHAPE = ',ISHAPE
134     VV=VV+DELV
        X=VV
        IF(ispace.eq.1)X=1e4/VV


        CALL lblrad_wave (X, WING, VMIN, VMAX, VREL, MAXDV, IBS, 
     1    IBD, Dist, INormal, Iray, ispace, DelH, nlayer, npath,
     1    ngas, maxlay, maxcon, press, temp, pp, amount,
     2    nlayin, maxinc, layinc, cont, scale, imod, idgas,
     3    isogas,iproc,emtemp,iphi,nem,vem,emissivity,tsurf,
     4    flagh2p,hfp,flagc, hfc, ifc, basep, baseh, output)




        DO I= 1, npath
         yp(I,1)=yp(I,2)
         yp(I,2)= output(I)
        ENDDO

        IF(IFLAG.EQ.1)THEN

           CALL lblconv(opfile,fwhm,ishape,npath,ispace,vv,delv,yp,
     1      nconv,vconv,yout,ynor)

        ENDIF

        IF(IFLAG.EQ.0)IFLAG=1

        XCOM = 100.0*(VV-VMIN)/(VMAX-VMIN)
        IF(XCOM.GE.XNEXT)THEN
            WRITE(*,1020)XCOM
            print*,(yout(1,j),J=1,nconv)
            XNEXT = XNEXT+10.0
        ENDIF
1020    FORMAT(' lblrtf_wave.f :: % complete : ',F5.1)


        IF(VV.LT.VMAX)GOTO 134

        DO I=1,NPATH
         DO J=1,NCONV
          YOUT(I,J)=YOUT(I,J)/YNOR(I,J)
          K=NCONV*(I-1)+J
          CONVOUT(K)=YOUT(I,J)
         ENDDO
        ENDDO


        CLOSE(DBLUN)

	RETURN

	END

************************************************************************
************************************************************************
