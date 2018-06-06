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
     3          emissivity,tsurf, convout, IPTF)

	IMPLICIT NONE

C	The variables defined here are those normally used when
C	calculating atmospheric paths. ! nptah. itype declared within
	INCLUDE '../includes/arrdef.f'
	INCLUDE '../includes/pathcom.f'
        INCLUDE '../includes/laycom.f'
        INCLUDE '../includes/lcocom.f'

C       Defines the maximum values for a series of variables (layers,
C         bins, paths, etc.)
	INCLUDE '../includes/contdef.f'

        INTEGER         nconv, npath1, itype1, I, J, K, mconv
 	INTEGER		MAXLIN1,ISHAPE, IPTF
        PARAMETER       (mconv=6000)
	INTEGER		INormal,Iray, ispace,nem,IBS(2),IBD(2), IFLAG
	REAL		Dist, FWHM1,X0,X1,WING1,VREL1,MAXDV
        REAL		VBOT,DELV1,RADIUS1,radextra
        REAL            vconv(nconv), convout(maxout3),zheight(maxpro)
        REAL            output(maxpat), yp(maxpat,2)                  
        REAL            yout(maxpat,mconv),tsurf,ynor(maxpat,mconv)
        REAL		vem(maxsec),emissivity(maxsec) 
	REAL		AAMOUNT(maxlay,maxgas),XX0
        DOUBLE PRECISION VV,XX,DX0,DX1
        REAL		VVR
        INTEGER         lun, ulog, iphi,NLINR,ipzen1
        PARAMETER       (lun=2, ulog=17)

        INTEGER         nout, IB

        INTEGER         FSTREC
	REAL		XNEXT,XCOM
        LOGICAL         fscatter, fdust,solexist
        double precision mu1(maxmu), wt1(maxmu), galb

        CHARACTER*100   klist, opfile1, solname, solfile, buffer

        CHARACTER*100   logfil, drvfil, sfile, FWHMFILE,lcofil
        INTEGER         NFWHM,MFWHM
        PARAMETER(MFWHM=1000)
        LOGICAL         FWHMEXIST,FEXIST
        REAL            VFWHM(MFWHM),XFWHM(MFWHM)


        common/scatd/mu1, wt1, galb

   
        INCLUDE '../includes/dbcom.f'
        INCLUDE '../includes/lincomc.f'

        common/defang/ipzen1
C       Need simple way of passing planetary radius to nemesis
        INCLUDE '../includes/planrad.f'

    
C-----------------------------------------------------------------------
C
C       Begin Program.
C
C-----------------------------------------------------------------------

        IJLCO=0
        opfile= opfile1         ! Renamed and assigned here because
        iphi= itype1            ! there was a conflict with similiar
                                ! declarations in pathcom.f which is
                                ! needed for other bits of this code.

C-----------------------------------------------------------------------
C
C       Run path file to create driver file.
C
C-----------------------------------------------------------------------

        CALL subpath(opfile)

        npath1= npath
        PRINT*,'LBLRTF_WAVE called.'

C       See if file is present forcing FWHM to vary with wavelength/wavenumber
        CALL FILE(OPFILE,FWHMFILE,'fwh')
        INQUIRE(FILE=FWHMFILE,EXIST=FWHMEXIST)
C       If such a file exists then read in the data
        IF(FWHMEXIST)THEN
         print*,'Reading FWHM information from : ',FWHMFILE
         OPEN(13,FILE=FWHMFILE,status='old')
          READ(13,*)NFWHM
          DO I=1,NFWHM
           READ(13,*)VFWHM(I),XFWHM(I)
          ENDDO
         CLOSE(13)
        ENDIF

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

        CALL FILE(opfile1,LCOFIL,'lco')
        print*,'LCOFIL = ',LCOFIL
        INQUIRE(FILE=LCOFIL,EXIST=FEXIST)
        IF(FEXIST)THEN
          print*,'Calling INIT_LCO'
          CALL INIT_LCO(LCOFIL)
        ENDIF

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
     1 IDGAS,ISOGAS,IPTF)


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

        fscatter= .false.
        DO I= 1, npath
                IF (imod(I).EQ.15.OR.imod(I).EQ.16) fscatter= .true.
                IF (imod(I).EQ.21.OR.imod(I).eq.22) fscatter = .true.
        ENDDO

        IF (fscatter) THEN
                CALL file(opfile, radfile, 'sca')
                WRITE(*,*)'     CALLING get_scatter'
                CALL get_scatter(radfile,ncont)
                WRITE(*,*)'     get_scatter COMPLETE'
                WRITE(*,*)' '
        ENDIF

C-----------------------------------------------------------------------
C
C       And the xsc files likewise
C
C-----------------------------------------------------------------------

        fdust= .false.
        IF (ncont.gt.0) fdust= .true.
        IF (fdust) THEN
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

       XX0 = VMIN-VREL
       IF(XX0.LT.0.)THEN
        XX0 = 0.0
       ENDIF
       PRINT*,'Calling init_cont',XX0,VMAX+VREL,WING

       CALL INIT_CONT(XX0,VMAX+VREL,WING)
       print*,'NBIN = ',NBIN
       do i=1,nbin
        print*,i,vbin(i)
       enddo
       VBOT=VBIN(1)
       PRINT*,'VBOT = ',VBOT
       IBS(1)=1
       IBS(2)=2

C      Read in 2 arrays of lines
       IB=1
       CALL FNDWAV(DBLE(XX0))
       FSTREC = DBREC

       PRINT*,'FSTREC = ',FSTREC

       MAXLIN1=MAXLIN
       print*,'lblrtf_wave : maxlin = ',MAXLIN1
       DX0 = DBLE(XX0)
       DX1 = DBLE(VMAX+VREL)
       print*,'lblrtf_wave : DX0,DX1 ',DX0,DX1
       CALL LOADBUFFER(DX0,DX1,FSTREC,MAXLIN1,MAXBIN,IB,
     1   NGAS,IDGAS,ISOGAS,VBOT,WING,NLINR,VLIN,SLIN,ALIN,ELIN,IDLIN,
     2 SBLIN,PSHIFT,DOUBV,TDW,TDWS,LLQ,NXTREC,FSTLIN,LSTLIN,LASTBIN)

       print*,'FSTREC,NXTREC,LASTBIN',fstrec,nxtrec,lastbin(ib)
       do i=1,nbin
        if(fstlin(ib,i).gt.0)then
         print*,ib,i,vbin(i),fstlin(ib,i),lstlin(ib,i),
     1    vlin(ib,fstlin(ib,i)),vlin(ib,lstlin(ib,i))
        else
         print*,ib,i,vbin(i),fstlin(ib,i),lstlin(ib,i)
        endif
       enddo
       NLINE(IB)=NLINR
       IBD(IB)=-1
       print*,'ib, nline(ib) = ',ib,nline(ib)

       IB=2
       FSTREC=NXTREC


       CALL LOADBUFFER(DX0,DX1,FSTREC,MAXLIN,MAXBIN,IB,
     1 NGAS,IDGAS,ISOGAS,VBOT,WING,NLINR,VLIN,SLIN,ALIN,ELIN,IDLIN,
     2 SBLIN,PSHIFT,DOUBV,TDW,TDWS,LLQ,NXTREC,FSTLIN,LSTLIN,LASTBIN)

       print*,'FSTREC,NXTREC,LASTBIN',fstrec,nxtrec,lastbin(ib)
       do i=1,nbin
        if(fstlin(ib,i).gt.0)then
         print*,ib,i,vbin(i),fstlin(ib,i),lstlin(ib,i),
     1   vlin(ib,fstlin(ib,i)),vlin(ib,lstlin(ib,i))
        else
         print*,ib,i,vbin(i),fstlin(ib,i),lstlin(ib,i)
        endif
       enddo

       NLINE(IB)=NLINR
       IBD(IB)=-1
       print*,'ib, nline(ib) = ',ib,nline(ib)

C-----------------------------------------------------------------------
C
C       Call subroutine lblrad_wave:
C
C-----------------------------------------------------------------------

C       Pass radius of planet to cirsradg
C       radius2 is radius held in planrad common block. Pass this to
C       cirsrad_wave in case it's been updated.
        radius1=radius2

C       Look to see if ipzen has been set to 2 and if so read in altitude of
C       top of atmosphere
        radextra=0.
        ipzen1=ipzen
        if(ipzen.eq.2)then
C          need the current height profile
           call readprfheight(opfile1,npro,zheight)
           radextra=zheight(npro)
        endif


        VV = DBLE(VMIN-DELV)
        IFLAG=0
        DO J=1,NCONV
         DO I=1,NPATH
          YOUT(I,J)=0.0
          YNOR(I,J)=0.0
         ENDDO
        ENDDO

        XNEXT=0.0

        call file(opfile,sfile,'sha')
        open(13,file=sfile,status='old')
        READ(13,*)ISHAPE
        close(13)

        print*,'ISHAPE = ',ISHAPE
C        open(37,file='raw.dat',status='unknown')
134     VV=VV+DBLE(DELV)
        XX=VV
        IF(ispace.eq.1)XX=1e4/VV
        VVR=SNGL(VV)

        CALL lblrad_wave (XX, WING, VMIN, VMAX, VREL, MAXDV, IBS, 
     1    IBD, Dist, INormal, Iray, ispace, DelH, nlayer, npath,
     1    ngas, maxlay, maxcon, totam, press, temp, pp, amount,
     2    nlayin, maxinc, layinc, cont, scale, imod, idgas,
     3    isogas,iproc,emtemp,iphi,nem,vem,emissivity,tsurf,
     4    flagh2p,hfp,flagc, hfc, ifc, basep, baseh, RADIUS1,
     5    radextra,output)




        DO I= 1, npath
         yp(I,1)=yp(I,2)
         yp(I,2)= output(I)
        ENDDO
C        write(37,*),vv,output(1)

        IF(IFLAG.EQ.1)THEN

           CALL lblconv(opfile,fwhm,ishape,npath,ispace,vv,delv,yp,
     1      nconv,dble(vconv),yout,ynor,FWHMEXIST,NFWHM,VFWHM,XFWHM)

        ENDIF

        IF(IFLAG.EQ.0)IFLAG=1

        XCOM = 100.0*(VVR-VMIN)/(VMAX-VMIN)
        IF(XCOM.GE.XNEXT)THEN
            WRITE(*,1020)XCOM
            XNEXT = XNEXT+10.0
        ENDIF
1020    FORMAT(' lblrtf_wave.f :: % complete : ',F5.1)
    

        IF(VV.LT.DBLE(VMAX))GOTO 134
C        close(37)
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
