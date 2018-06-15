************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C_TITLE:		SUBROUTINE CIRSRAD_WAVE
C
C_DESCR:
C
C_ARGS:	Input Variables
C
C	AvgCONTMP:REAL	The average of the three points returned by CIACON
C	INORMAL1:INT     flag for ortho:para ratio (0=equilibrium =1:1)
C			(1=normal   =3:1)
C       DELH(NLAYER):REAL Physical thickness of each layer (km).
C       ispace           integer Indicates if wavelengths in vconv and
C                               vwave are in wavenumbers(0) or
C                               wavelengths (1)
C       DIST:REAL Distance from Sun (as prompted by CIRSDRV) in
C			units of AU.
C	NLAYER		Number of layers to consider.
C	NPATH		Number of individual paths through the layers.
C       NGAS		Number of gases to consider.
C	LIMLAY
C	LIMCONT
C       PRESS		total pressure for each layer (atm).
C       TEMP		temperature (Kelvin) of each layer.
C	PP		partial pressure of each gas (atm).
C       AMOUNT		number of molecules/cm2 for each layer
C                       and gas (normal incidence).
C	NWAVE		The number of calculation wavenumbers.
C	VWAVE		Bin centres in wavenumber space.
C       NLAYIN		Each layer defined above can be included in the 
C			calculation for each path several times. This 
C			is the number of layers included in the 
C			calculation for each path (if you include it 
C			twice count it twice).
C	INCDIM		First dimension of SCALE and LAYINC. 
C                       INCDIM >= all NLAYIN
C       LAYINC		The layer numbers in order of inclusion for each 
C			path. The order also defines the processing of 
C			each layer depending upon the model.
C	CONT
C       SCALE		Scaling factor (eg cosine of incidence) for each 
C			included layer in each of the paths. Note: when 
C			a layer is included more than once different 
C			scaling factors can be defined.
C       IMOD		Type of model to use to compute OUTPUT from 
C                       transmissions for each path
C       IDGAS		The LOCAL gas identifier for each gas. The local 
C			gas identifier agrees with HITRAN id's as far as 
C			possible (Jan 1988) with extra id's for gases not
C			in the HITRAN compilation. eg.those in the GEISA 
C			compilation.
C       ISOGAS		The local isotopic identifier, if zero all 
C			isotopes of the gas are included. Isotope id's 
C			also agree as far as possible with HITRAN id's. 
C			Similar (1,2,3...n) id's have been defined for 
C			additional gases.
C       EMTEMP		Temperatures to use in emission calculations.
C       ITYPE		Value designating the chosen scattering routine
C       NEM             integer Number of points in surface
C                                       emissivity spectrum    
C       VEM(NEM)        Wavelengths of emissivity spectrum
C       EMISSIVITY(NEM) Tabulated emissivities
C       TSURF           Surface temperature (K)     
C       RADIUS1         Planetary radius (km) at reference pressure (usually 
C								      1 bar)      
C	RADEXTRA	Any correction to RADIUS1 needed for disc-integration
C       FLAGH2P         FLAGH2P=1 if a para-H2 fraction
C                                 profile has been read.
C	HFP(NLAYER)	Para-H2 fraction in each layer (if FLAGH2P=1)
C       FLAGC		FLAGC=1 if variable fractional cloud cover profile
C			        is defined
C       HFC(NLAYER)     Fractional cloud cover in each layer (if FLAGC=1)
C	IFC(NLAYER)	integer to state if cloud is uniform or broken.
C	BASEP(NLAYER)	Base pressure of each layer (atm)
C	BASEH(NLAYER)	Base height of each layer (km)
C
C
C
C	Output variables
C
C       OUTPUT		Output values at each wavenumber for each output 
C			type for each path
C
C_HIST:	
C-----------------------------------------------------------------------

      SUBROUTINE cirsrad_wave (Dist, INormal, Iray, ispace,DelH,nlayer,
     1    npath, ngas, limlay, limcont, totam, press, temp, pp, amount, 
     2    nwave, vwave, nlayin, incdim, layinc, cont, scale, imod,
     3    idgas, isogas, emtemp, itype, nem, vem, emissivity, tsurf, 
     4    flagh2p,hfp, flagc, hfc, ifc, basep, baseh, RADIUS1, 
     5    radextra,output)

      IMPLICIT NONE


C		Internal dimensions

C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)
       INCLUDE '../includes/arrdef.f'

C		Passed variables

	INTEGER	Inormal,ik,iray,flagh2p,flagc,ispace,jf
	INTEGER	nlayer, npath, ngas, limlay, limcont, nwave, 
     1		nlayin(npath), incdim, layinc(incdim, npath), 
     2          idgas(ngas), isogas(ngas), imod(npath),igas,
     3 		infr,ish,irayx

	REAL	press(nlayer), temp(nlayer), pp(limlay,ngas), 
     1		amount(limlay,ngas), vwave(nwave), cont(limcont,nlayer), 
     2		scale(incdim,npath), emtemp(incdim,npath), hfp(nlayer), 
     3		output(npath,nwave),bb(maxlay,maxpat),xf,tsurf,dv,
     4          hfc(nlayer),vwavef(maxbin),basep(nlayer),baseh(nlayer),
     5          esurf,radground,totam(nlayer),RADIUS1,xdist,xfac,
     6		eoutput(npath,nwave)
        REAL    vem(MAXSEC),emissivity(MAXSEC),interpem
        REAL	xmu,dtr,xt1,radextra
        INTEGER ifc(limcont,nlayer),nem,j0,jnh3,jch4,jh2,jhe,ILBL

	REAL	XCOM,XNEXT,FPARA,vv,XRAY,RAYLEIGHJ,RAYLEIGHA,RAYLEIGHV
        REAL    RAYLEIGHLS,fheh2,fch4h2,fh2,fhe,fch4,fnh3
C		Dust variables

	INTEGER	nsec, ncont
	REAL	vsec(maxsec), xsec(2,maxcon,maxsec)

C		K table variables

	INTEGER	lun(maxbin,maxgas), ireck(maxbin,maxgas), lun0, irec0, 
     1		npk, ntk, ng, intmod(maxbin,maxgas), ntab,
     2		lunlbl(maxgas),npklbl,ntklbl
	REAL	xmink(maxbin,maxgas), delk(maxbin,maxgas), xmin, delx, 
     1		k_g(maxg), k_g1(maxg), k_g2(maxg), q1, q2, 
     2		pk(maxk), tk(maxk), g_ord(maxg), t2k(maxk,maxk), 
     3		delg(maxg), kl_g(maxg,maxlay),xminklbl,
     4		kout(maxlay,maxgas,maxg),p1,p2,frack(maxbin,maxgas)
        REAL    dkoutdt(maxlay,maxgas,maxg)
        REAL    dkoutdtlbl(maxlay,maxgas,maxg)   
	REAL    pklbl(maxk),tklbl(maxk),t2klbl(maxk,maxk)
	REAL	delvklbl,koutlbl(maxlay,maxgas),delklbl
        REAL    basehf(maxlay),basepf(maxlay),basehS(maxlay)
        REAL    scaleS(maxlay),Jsource(maxlay),delhs(maxlay)
C		Scattering variables

	INTEGER	liscat(maxcon), isol, nmu, imu, lowbc, npointk,
     1          lncons(maxcon), ibaseH(maxlay), ic, ilays,
     2		kf,nf, nfX, lnorm(maxcon), itype
	REAL	lfrac(maxcon,maxlay), solar, pi, theta0, dist, 
     1		bnu(maxscatlay), aphi, radg(maxmu), solar1, 
     2		lcons(maxcon,maxscatpar), rad1, omega(maxscatlay),
     3		tsun, v, frcscat(maxcon,maxlay), sol_ang, emiss_ang
        REAL 	solara
        REAL    bnuS(maxlay),fcover(maxscatlay),xsolar
        REAL    umif(maxmu,maxscatlay,maxf),sum
        REAL    uplf(maxmu,maxscatlay,maxf)	
        REAL    gmif(maxmu,maxscatlay,maxf)
        REAL    gplf(maxmu,maxscatlay,maxf)	
        REAL    gmifa(maxmu,maxscatlay)
        REAL    gplfa(maxmu,maxscatlay)	
        REAL    gmifc(maxmu,maxscatlay)
        REAL    gplfc(maxmu,maxscatlay)	
        REAL    gmifb(maxmu,maxscatlay,maxf)
        REAL    gplfb(maxmu,maxscatlay,maxf)	
        REAL    umift(maxmu,maxmu,maxscatlay,maxf)
        REAL    uplft(maxmu,maxmu,maxscatlay,maxf)	
        DOUBLE PRECISION mu1(maxmu), wt1(maxmu), galb, galb1, 
     1		pplsto(0:40,maxmu,maxmu,maxscatlay),
     2		pmisto(0:40,maxmu,maxmu,maxscatlay),
     3          eps(maxscatlay),epsS(maxlay),omegas(maxscatlay),
     4		epsb(maxscatlay),omegasb(maxscatlay)
	LOGICAL	scatter, single, sphsingle

C		Internal variables

	INTEGER	I, J, K, L, Ipath, Ig, nlays, lstcel,Jpath
	INTEGER irec0lbl
	REAL	utotl(maxlay), qh(maxlay), qhe(maxlay),
     1		frac(maxlay,maxgas), qh_he(maxlay), dist1,
     2		x, taucon(maxlay)
        REAL    taugas(maxlay), tauscat(maxlay), p, t
        REAL    taugasc(maxlay),xp
        REAL    tau, tau2, asec(maxsec), bsec(maxsec),
     1          tausc(maxcon), taus(maxlay),
     5		tautmp(maxlay), corkout(maxpat,maxg),
     7		ecorkout(maxpat,maxg),ctrans(mpoint,maxg),
     6		planck_wave, tmp, error, asec2(maxsec), 
     3          bsec2(maxsec), xsec2(2,maxcon,maxsec), muemiss,
     3          muinc, f(maxlay), fint(maxlay), tmp1, 
     4          intscat, pastint, phase(maxcon), fup(maxlay,maxg),
     5          tmp2, taug(maxlay),ssfac
        REAL    taucloud(maxcon,maxlay),tauclscat(maxcon,maxlay)
        REAL    taucl(maxcon,maxlay),taucs(maxcon,maxlay)
        integer icloud(maxcon,maxlay)
        REAL	fdown(maxlay,maxg),fwhmk,
     1          delvk,tauray(maxlay),f1(maxlay),g11(maxlay),
     2          g21(maxlay),taur(maxlay) 
	REAL ftop(maxg)
	DOUBLE PRECISION	tr, trold, taud, dtmp1, dtmp2, dpi, 
     1          dphase, calpha, dmuinc, dmuemiss, draddeg, dinc, demiss,
     2          tausun

C		Misc variables or continuum variables

	INTEGER	ii, id1,J1,J2,L1,L2,IABSORB(5)
	REAL	DelH(nlayer), AAmount(maxgas), PPP(maxgas),
     1		XLen, AvgCONTMP, CONTMP(IORDP1), DABSORB(7)

        REAL   rad(maxlay+1),refrac(MAXLAY),r0(3),vec(3),r1(3)
        REAL   ethick,scaleH,svec(3),vec1(3),ff,sproduct,rtmp
        REAL   trace(maxpat,3),disttrace(maxpat),scaletrace(maxpat)
        REAL   scaletraceLOS(maxpat),traceLOS(maxpat,3),tauLOS
        REAL   tausol,vectrace(maxpat,3),vectraceLOS(maxpat,3)
        integer jtrace(maxlay),startlay,npathLOS,ipzen,icont
        integer jtraceLOS(maxlay),npathsol,JLAY,ILAY,KLAY
        common/defang/ipzen
        common/lbltable/ILBL

        INTEGER LUNIS,IRECL,IOFF,NLAYERF,NMUF,NWAVEF,NGF,IFLUX,NFF
        integer fintrad,first
        character*100 fintname
C		Common blocks and parameters

C       Solar reference spectrum common block
        real swave(maxbin),srad(maxbin),solrad
        integer iread,nspt,iform
        common /solardat/iread,iform,solrad,swave,srad,nspt


	common/dust/vsec,xsec,nsec,ncont
	common/interpk/lun, ireck, xmink, delk, frack, pk, npk, tk, 
     1      t2k, ntk, ng, delvk, fwhmk, g_ord, delg, kout, dkoutdt
	common/interpklbl/lunlbl,irec0lbl,xminklbl,delklbl,npointk, 
     1    pklbl,npklbl,tklbl,t2klbl,ntklbl,koutlbl,dkoutdtlbl
	common/scatd/mu1, wt1, galb
	common/scatter1/nmu, isol, dist1, lowbc, liscat, lnorm,
     1		lncons, lcons, sol_ang, emiss_ang, aphi, nf
	common/phasesto/pplsto,pmisto
        common/intrad/fintrad,fintname
        common/interp_phase_initial/first

	PARAMETER	(tsun=5900.,theta0=4.65241e-3, pi=3.141593,
     1			error=5.e-5,LUNIS=61,infr=62) 

C-----------------------------------------------------------------------
C
C	Check input parameters for possible problems.
C
C-----------------------------------------------------------------------



	if (nwave.gt.maxbin) then
		write (*,*) ' CIRSRAD_WAVE: Too many bins'
		write (*,*) ' NWave = ',nwave,' Maxbin = ',maxbin
		stop
	endif

	if (nlayer.gt.maxlay) then
		write (*,*) ' CIRSRAD_WAVE: Too many layers'
		write (*,*) ' Nlayer = ',nlayer,' maxlay = ',maxlay
		stop
	endif

	if (npath.gt.maxpat) then
		write (*,*) ' CIRSRAD_WAVE: Too many paths'
		write (*,*) ' Npath = ',npath,' Maxpat = ',maxpat
		stop
	endif

	if (ngas.gt.maxgas) then
		write (*,*) ' CIRSRAD_WAVE: Too many gases'
		write (*,*) ' Ngas = ',ngas,' Maxgas = ',maxgas
		stop
	endif

	if (ncont.gt.maxcon) then
		write (*,*) ' CIRSRAD_WAVE: Too many dust continua'
		write (*,*) ' Ncont = ',ncont,' maxcon = ',maxcon
		stop
	endif

	if (nsec.gt.maxsec) then
		write (*,*) ' CIRSRAD_WAVE: Too many dust continua pts'
		write (*,*) ' Nsec = ',nsec,' Maxsec = ',maxsec
		stop
	endif

        IF(ILBL.EQ.2)THEN
         NG=1
         DELG(1)=1.
        ELSE
 	  if (ng.gt.maxg) then
		write (*,*) ' CIRSRAD_WAVE: Too many g ordinates'
		write (*,*) ' Ng = ',ng,' Maxg = ',maxg
		stop
	  endif
        ENDIF

	if ((npk.gt.maxk).or.(abs(ntk).gt.maxk)) then
		write (*,*) ' CIRSRAD_WAVE: Too many P/T points in K',
     1			' tables' 
		write (*,*) ' Npk = ',npk,' Maxk = ',maxk
		write (*,*) ' Ntk = ',ntk,' Maxk = ',maxk
		stop
	endif

	if ((npklbl.gt.maxk).or.(abs(ntklbl).gt.maxk)) then
		write (*,*) ' CIRSRAD_WAVE: Too many P/T points in K',
     1			' tables' 
		write (*,*) ' Npklbl = ',npklbl,' Maxk = ',maxk
		write (*,*) ' Ntklbl = ',ntklbl,' Maxk = ',maxk
		stop
	endif

	do Ipath = 1, npath
		if (nlayin(Ipath).gt.incdim) then
		   write (*,*) ' CIRSRAD_WAVE: Nlayin exceeds Incdim'
			write (*,*) ' Ipath = ', Ipath,' Nlayin = ',
     1				nlayin(Ipath), ' Incdim = ', incdim
			stop
		endif
		do I = 1, nlayin(Ipath)
			if (layinc(I,Ipath).gt.nlayer) then
			write (*,*) ' CIRSRAD_WAVE: Layinc exceeds',
     1					' nlayer'
				write (*,*) ' Ipath = ', Ipath,
     1					' Layer = ', I, ' Layinc = ',
     2					layinc(I,Ipath), ' Nlayer = ',
     3					nlayer
				stop
			endif
		enddo
	enddo


C-----------------------------------------------------------------------
C
C     	Initialise some variables.
C
C		imod = 15		Multiple scattering
C		imod = 16		Single scattering approx
C		imod = 22		Limb scattering
C		imod = 23		Limb scattering (precomputed flux)
C
C-----------------------------------------------------------------------

C       Set flag for code to read in modified PHASE*.DAT files if required
        first=0

	dpi = 4.0d0 * datan(1.0d0)
	draddeg = 2.0d0 * dpi/360.0d0

	scatter=.FALSE.
	single=.FALSE.
	do Ipath = 1, npath
		if (imod(Ipath).eq.15.or.
     1	  (imod(Ipath).ge.22.and.imod(Ipath).le.27)) then
        		scatter=.true.
		end if
		if (imod(Ipath).eq.16) then
			single = .true.
			dinc = dble(sol_ang)
			demiss = dble(emiss_ang)
			dmuinc = cos(dinc * draddeg)
			dmuemiss = cos(demiss * draddeg) 
			muemiss = sngl(dmuemiss)
			muinc = sngl(dmuinc)
                        ssfac = muinc/(muinc+muemiss)
                        do I = 1, nlayer
                                scale(I,1) = sngl(1./dmuemiss +
     1                                  1./dmuinc)
                        enddo

			
			calpha = sin(sol_ang*pi/180.0)*
     1				 sin(emiss_ang*pi/180.0)*
     2				 cos((aphi-pi)*pi/180.0) - 
     3				 muemiss*muinc


		endif
                if(imod(ipath).eq.28)then
                  sphsingle=.TRUE.
                else
                  sphsingle=.FALSE.
                endif
	end do

c	Isec = min(4,nsec)

        if(ng.gt.1)then
 	 ntab = npk * abs(ntk) * ng
        else
         ntab = npklbl*abs(ntklbl)
        endif
C-----------------------------------------------------------------------
C
C	Precompute volume fractions. The factor of 1e-20 applied to
C	utotl and totam arises for the following reason. In order to avoid 
C	underflow, RADTRAN (and CIRSRAD) routines multiply line strengths
C	by a factor of 1.e47. Thus, when final calculation routines get 
C	hold of the strengths, they must correct by 1.e-47. This is broken
C	into two factors, 1.e-27 applied directly in calculating the K
C	tables, and 1.e-20 applied to the amounts where they are used in
C	conjunction with the K tables.
C
C-----------------------------------------------------------------------

        JH2 = -1
        JNH3 = -1
        JHE = -1
        JCH4 = -1
        DO IGAS=1,NGAS
          IF(IDGAS(IGAS).EQ.39.AND.
     1     (ISOGAS(IGAS).EQ.0.OR.ISOGAS(IGAS).EQ.1))THEN
            JH2 = IGAS
          ENDIF
          IF(IDGAS(IGAS).EQ.40.AND.
     1     (ISOGAS(IGAS).EQ.0.OR.ISOGAS(IGAS).EQ.1))THEN
            JHE = IGAS
          ENDIF
          IF(IDGAS(IGAS).EQ.6.AND.
     1     (ISOGAS(IGAS).EQ.0.OR.ISOGAS(IGAS).EQ.1))THEN
            JCH4 = IGAS
          ENDIF
          IF(IDGAS(IGAS).EQ.11.AND.
     1     (ISOGAS(IGAS).EQ.0.OR.ISOGAS(IGAS).EQ.1))THEN
            JNH3 = IGAS
          ENDIF
        ENDDO




	DO I= 1, nlayer
		utotl(I)= 0.
		DO J= 1, ngas
			frac(I,J) = pp(I,J) / press(I)
			utotl(I) = utotl(I) + amount(I,J)
		ENDDO
		utotl(I) = utotl(I) * 1.e-20
		totam(I) = totam(I) * 1.e-20

	ENDDO

C-----------------------------------------------------------------------
C
C	Begin loop over wavenumber bins. At each wavenumber, first 
C	initialise arrays over the levels, then step over the g 
C	ordinates, calculating the relevant output.
C
C	xsec(1....) is the extinction coefficient, xsec(2....) the 
C	single scattering albedo.
C
C-----------------------------------------------------------------------

	IF (nsec.gt.1) THEN
		DO K = 1, ncont
		 	DO L = 1, nsec
				asec(L) = xsec(1,K,L)
				bsec(L) = xsec(2,K,L)
		 	ENDDO
                     
		 	CALL cspline (vsec,asec,nsec,1.e30,1.e30,asec2)
		 	CALL cspline (vsec,bsec,nsec,1.e30,1.e30,bsec2)

		 	DO L = 1, nsec
				xsec2(1,K,L) = asec2(L)
				xsec2(2,K,L) = bsec2(L)
		 	ENDDO
		ENDDO
	ENDIF

        XNEXT=0.0
        IRECL = 4

        IFLUX=0
        do I=1,npath
         if(imod(i).eq.15.and.itype.eq.13)then
          iflux=1
         endif
         if(imod(i).eq.23)then 
          iflux=2
         endif
        ENDDO

        print*,'IFLUX = ',IFLUX
        IOFF = 0

C       Look to see if there is a scattering net-flux calculation. If so
C       open output file for internal radiation field output. 
        ish=0
        do ipath=1,npath
         if(imod(ipath).eq.24)ish=1
        enddo
        print*,'ish = ',ish
        if(ish.eq.1)then
         print*,'opening ish.dat'
         open(infr,file='ish.dat',status='unknown')
         write(infr,*)nlayer
         write(infr,*)nmu
         write(infr,*)nwave
         write(infr,*)nf
         write(infr,*)sol_ang
         write(infr,*)(vwave(i),i=1,nwave)
         write(infr,*)(basep(i),i=1,nlayer)
         write(infr,*)(baseh(i),i=1,nlayer)
        endif
  

        IF(IFLUX.EQ.1)THEN

         OPEN(LUNIS,FILE=fintname,STATUS='UNKNOWN',
     1    ACCESS='DIRECT',RECL=IRECL)

         WRITE(LUNIS,REC=1)nlayer
         WRITE(LUNIS,REC=2)nmu
         WRITE(LUNIS,REC=3)nwave
         WRITE(LUNIS,REC=4)ng
         WRITE(LUNIS,REC=5)nf

         DO I=1,nwave
          WRITE(LUNIS,REC=5+I)vwave(i)
         ENDDO
         DO I=1,nlayer
          WRITE(LUNIS,REC=5+NWAVE+I)basep(i)
         ENDDO
         DO I=1,nlayer
          WRITE(LUNIS,REC=5+NWAVE+NLAYER+I)baseh(i)
         ENDDO

         IOFF = 5+NWAVE+2*NLAYER+1

        ENDIF

        IF(IFLUX.EQ.2)THEN
         IF(ISOL.EQ.0)THEN
          FINTNAME='internal.fld'
         ELSE
          FINTNAME='internalsol.fld'
         ENDIF

         OPEN(LUNIS,FILE=FINTNAME,STATUS='OLD',
     1    ACCESS='DIRECT',RECL=IRECL)

         READ(LUNIS,REC=1)nlayerF
         READ(LUNIS,REC=2)nmuF
         READ(LUNIS,REC=3)nwaveF
         READ(LUNIS,REC=4)ngF
         READ(LUNIS,REC=5)nfF
         print*,'LUNIS',LUNIS,nlayerF,nmuF,nwaveF,ngF,nfF

         if(ngF.ne.ng)then
          print*,'Error in cirsrad_wave: ngf <> ng'
          print*,ngf,ng
          stop
         endif

         print*,nwaveF
         DO I=1,nwaveF
          READ(LUNIS,REC=5+I)vwaveF(i)
          print*,vwaveF(i)
         ENDDO

         DO I=1,nlayerF
          READ(LUNIS,REC=5+NWAVEF+I)basepF(i)
         ENDDO

         DO I=1,nlayerF
          READ(LUNIS,REC=5+NWAVEF+NLAYERF+I)basehF(i)
          print*,basehF(i)
         ENDDO

         IOFF = 5+NWAVEF+2*NLAYERF+1

        ENDIF


        DO I = 1, nwave
		x = vwave(I)
                esurf = interpem(nem,vem,emissivity,x)
C               Set vv to the current WAVENUMBER
                vv = x
                if(ispace.eq.1)then
                  vv=1e4/x 
                endif
		IF ((scatter).OR.(single).OR.(sphsingle)) THEN
			IF (isol.EQ.1) THEN
			  CALL get_solar_wave(x,dist,solar)
			ELSE
			  solar = 0.
			ENDIF
		ENDIF

		IF (single) THEN
			DO J = 1, ncont
				CALL get_hg(x,calpha,ncont,J,dphase)
			  	phase(J) = sngl(dphase)
			ENDDO	
C                       Calculate Rayleigh scattering too
                        phase(ncont+1)=sngl(0.75*(1+calpha*calpha))
		ELSE
			DO J = 1, ncont
				phase(J) = 0.
			ENDDO
		ENDIF
C-----------------------------------------------------------------------
C
C	Begin loop over layers to initialise various arrays. Set values 
C	of continuum and scattering optical depths, as well as getting 
C	K coefficients for each layer.
C
C	Taucon is the continuum optical depth due to both gas and
C	particle contributions. The particle component arises from the 
C	extinction coefficient.
C	Tauscat is the scattering optical depth derived from the the
C	particle distribution and the single scattering albedo. For 
C	multiple scattering, this is passed to scattering routines.
C
C-----------------------------------------------------------------------

		intscat = 0.
		pastint = 0.

C     read in k-coefffients for each gas and each layer
                IF(ILBL.EQ.0)THEN
                  CALL GET_K(NLAYER,PRESS,TEMP,NGAS,I,X)
                ELSE
                  CALL GET_KLBL(NLAYER,PRESS,TEMP,NGAS,X)
                ENDIF

 
                FNH3=0.
                FH2=0.
                FHe=0.
                FCH4=0.

		DO J = 1, nlayer
		 taucon(J) = 0.
		 taugasc(J) = 0.
		 tauscat(J) = 0.
		 f(J) = 0.
		 p = press(J)
		 t = temp(J)

                 IF(JNH3.GT.0.)FNH3 = FRAC(J,JNH3)
                 IF(JCH4.GT.0.)FCH4 = FRAC(J,JCH4)
                 IF(JH2.GT.0.)FH2 = FRAC(J,JH2)
                 IF(JHE.GT.0.)FHE = FRAC(J,JHE)

C-----------------------------------------------------------------------
C
C	Add contributions from dust absorption and scattering into the 
C	continuum and scattering optical depths. Note that dust (or 
C	equivalent) is defined for each layer by an array of number of 
C	dust particles/cm2 (CONT(M,LAYER)) for each M particle type. The 
C	dust cross sections are held in XSEC(2,M,NSEC) where the 
C	corresponding wavenumbers are VSEC(NSEC). The first array
C	of XSEC (i.e. XSEC(1,.....)) contains extinction X-sections,
C	the second (1.e. XSEC(2,....)) is the single scattering albedo
C	multiplied by the extinction X-section (performed in READ_XSEC).  
C
C	Interpolations to the correct wavelength are done using 
C	cubic splines. The required derivatives are calculated earlier
C	outside of the loop over levels.
C
C	Calculate fractional contribution of each dust type to overall
C	scattering optical depth.	
C
C	f is used in the single scattering approximation and is the
C	product of single scattering albedo, number density, extinction
C	X-section, phase function, and solar radiance divided by 4pi.
C
C-----------------------------------------------------------------------

			DO K = 1, ncont
			 tausc(K) = 0.
			 IF (cont(K,J).GT.0.) THEN
			   IF (nsec.EQ.1) THEN
				tau = cont(K,J) * xsec(1,K,1)
				tau2 = cont(K,J) * xsec(2,K,1)
			   ELSE
				DO L = 1, nsec
				 asec(L) = xsec(1,K,L)
				 bsec(L) = xsec(2,K,L)
				 asec2(L) = xsec2(1,K,L)
				 bsec2(L) = xsec2(2,K,L)

				ENDDO

				CALL csplint (vsec, asec, asec2, 
     1					nsec, x, tau)
				CALL csplint (vsec, bsec, bsec2,
     1					nsec, x, tau2)

				if(tau.lt.0)then
C             			 print*,'tau lt 0: Particle type ',K
C                                 print*,nsec,vsec(1),vsec(nsec)
C                                 print*,nsec,asec(1),asec(nsec)
C                                 print*,x,tau
C                                 print*,'Do linear interpolation'
                                 jf=-1
                                 do l=1,nsec-1
                         if(x.ge.vsec(l).and.x.lt.vsec(l+1))then
                          xf = (x-vsec(l))/(vsec(l+1)-vsec(l))
                          tau = (1.0-xf)*asec(l)+xf*asec(l+1)
                          jf=l
                         endif
                                 enddo
                                 if(x.le.vsec(1))tau=asec(1)
                                 if(x.ge.vsec(nsec))tau=asec(nsec)
C                                 print*,'new tau ',tau
C				 if(jf.gt.0)then
C                                  print*,jf,vsec(jf),vsec(jf+1),xf
C                                  print*,asec(jf),asec(jf+1)
C                                 endif
                                endif


				if(tau2.lt.0)then
C              			  print*,'tau2 lt 0: Particle type ',K
C                                 print*,nsec,vsec(1),vsec(nsec)
C                                 print*,nsec,bsec(1),bsec(nsec)
C                                 print*,x,tau2
C                                 print*,'Do linear interpolation'
                                 jf=-1
                                 do l=1,nsec-1
                         if(x.ge.vsec(l).and.x.lt.vsec(l+1))then
                          xf = (x-vsec(l))/(vsec(l+1)-vsec(l))
                          tau2 = (1.0-xf)*bsec(l)+xf*bsec(l+1)
                          jf=l
                         endif
                                 enddo
                                 if(x.le.vsec(1))tau2=bsec(1)
                                 if(x.ge.vsec(nsec))tau2=bsec(nsec)
C                                 print*,'new tau2 ',tau2
C				 if(jf.gt.0)then
C				  print*,jf,vsec(jf),vsec(jf+1),xf
C				  print*,bsec(jf),bsec(jf+1)
C				 endif
                                endif

                                if(cont(K,J).lt.0)then
                                 print*,'CONT,K,J',cont(k,j),k,j
                                endif

				tau = tau * cont(K,J)
				tau2 = tau2 * cont(K,J)

			   ENDIF	
			   taucon(J) = taucon(J) + tau
			   tauscat(J) = tauscat(J) + tau2
                           taucloud(K,J) = tau
                           tauclscat(K,J) = tau2   
			   tausc(K) = tau2
C                          ##### f(J) is average phase function ######
  			   f(J) = f(J)+phase(K)*tau2

			 ENDIF
			ENDDO
			intscat = intscat + tauscat(J)
			pastint = intscat

			DO K = 1, ncont
				IF (tauscat(J).eq.0.) THEN
					frcscat(K,J) = 0.0
				ELSE
					frcscat(K,J) = tausc(K)/
     1						tauscat(J)
				ENDIF
			ENDDO

C-----------------------------------------------------------------------
C
C	Add continuum contributions from hydrogen and other gases into
C	the layer optical depths. 
C
C	NGASCON is a RADTRANS routine which
C	calculates continuum contributions from various gases
C
C-----------------------------------------------------------------------

C		        WRITE (*,*) '       CALLING NGascon for each gas'
			DO K = 1, ngas

C	Computes a polynomial approximation to any known continuum spectra 
C       for a particular gas over a defined wavenumber region.

			 CALL ngascon(vv,idgas(K),isogas(K),
     1			  amount(J,K),pp(J,K),p,t,AvgCONTMP)
                          if(AvgCONTMP.ge.0)then
C			    Can occasionally get -ve AvCONTMPs at
C			    edge of Karkoschka CH4 data.
 			    taucon(J) = taucon(J) + AvgCONTMP
 			    taugasc(J) = taugasc(J) + AvgCONTMP
                          endif

			ENDDO
C-----------------------------------------------------------------------
C
C	CIACON: to compute gaseous continuum spectra from a variety of gas
C	pairs due to Collisionally Induced-Absorption.
C
C	Here compute gaseous continua not covered by wings in LBL or by
C	band data in other codes and add to the continuum polynomials.
C
C	Computed in the same way as filter profile interpolation but in
C	a different subroutine to allow easy modification. 
C
C-----------------------------------------------------------------------

C	Reduce the 2d arrays for gas amounts (no./unit vol.) and partial
C	pressure to 1d for passing into subroutine
			DO 320 ii=1, NGAS
				AAmount(ii) = Amount(J,ii)
				PPP(ii) = PP(J,ii)
320			CONTINUE

C	ID1= 0 means no diagnostic print statements while in CIACON.
			id1 = 0
			XLEN = DELH(J)

C       This goes off to do the collision induced-absorption
C                        print*,'FLAGH2P = ',FLAGH2P
	             	IF(FLAGH2P.EQ.1)THEN
                	   FPARA = HFP(J)
 			   CALL NPARACON(vv,P,T,
     1			     NGas,idgas,isogas,AAmount,PPP,
     2			     FPARA,XLen,AvgCONTMP,IABSORB,DABSORB,
     3                       id1)
             		ELSE
			   CALL NCIACON(vv,P,T,INormal,
     1			     NGas,idgas,isogas,AAmount,PPP,
     2			     XLen,AvgCONTMP,IABSORB,DABSORB,id1)
             		ENDIF
C                        print*,'FFF',J,taucon(J)
			taucon(J)= taucon(J) + AvgCONTMP
			taugasc(J)= taugasc(J) + AvgCONTMP

C	The code below is if Rayleigh scattering is considered.
C        print*,'IRAY, ITYPE = ',IRAY,ITYPE
	IF(IRAY.GT.0.AND.ITYPE.NE.0)THEN
C               Set vv to the current WAVENUMBER
                vv = x
                if(ispace.eq.1)then
                  vv=1e4/x 
                endif
                if(IRAY.EQ.1)THEN
                 xray = RAYLEIGHJ(vv,p,t)*1E20
                elseif(IRAY.EQ.2)then
                 xray = RAYLEIGHV(vv,p,t)*1E20
                elseif(IRAY.EQ.3)then
                 xray = RAYLEIGHA(vv,p,t)*1E20
                else
                 if(FH2.GT.0.0)THEN
                  fheh2=FHE/FH2
                  fch4h2=FCH4/FH2
                 else
                  fheh2=0.
                  fch4h2=0.
                 endif
                 xray = RAYLEIGHLS(vv,fheh2,fch4h2,fnh3)*1E20
		endif

                avgcontmp = totam(j)*xray
                tauray(J)=avgcontmp

                if(AvgCONTMP.ge.0.0)then
 		 taucon(J)= taucon(J) + AvgCONTMP
                 taugasc(J)= taugasc(J) + AvgCONTMP
                 tauscat(J)=tauscat(J) + AvgCONTMP
	        endif

C               Calculate single-scattering contribution
                f(J) = f(J) + phase(ncont+1)*tauray(J)
                           
        ELSE
               	tauray(J)=0.0
	ENDIF


C       ### renormalise f(J) to be tau-averaged phase function
        if(tauray(J)+tauscat(J).gt.0)then
           f(J) = f(J)/tauscat(J)
        else
           f(J) = 0.0
        endif


C-----------------------------------------------------------------------
C
C	For each gas, look up the K coefficients. Stepping through the 
C	gases, combine each set of coeffs with the combination of all
C	previous gases using OVERLAP. Finally, store the layer dependent
C	summed K coefficients. End the loop over each layer.
C
C-----------------------------------------------------------------------


        	IF(ILBL.EQ.0)THEN
	  	  DO K = 1, ngas
C                         print*,'I,K',I,K,lun(I,K)
			 IF (lun(I,K).LT.0) THEN
				DO L = 1, ng
					k_g2(L) = 0.
				ENDDO
			 ELSE
				DO L = 1,NG
                          		k_g2(L)=KOUT(J,K,L)
				ENDDO
			 ENDIF
C                         print*,'k_g2',k_g2(5)

		 	 IF (K.EQ.1) THEN
				q1 = frac(J,K)
				DO L = 1, ng
					k_g(L) = k_g2(L)
				ENDDO
C                                print*,'AA',q1,k_g(5)
			 ELSE
				q2 = frac(J,K)
				DO L = 1, ng
					k_g1(L) = k_g(L)
				ENDDO
C                                print*,'BB',q2,k_g1(5)
C	This subroutine combines the absorption coefficient distributions
C	of two overlapping gases. The overlap is implicitly assumed to be 
C	random and the k-distributions are assumed to have NG-1 mean
C	values and NG-1 weights. Correspondingly there are NG ordinates in
C	total.
                                if(q2.gt.0.and.q1.gt.0)then
				 CALL overlap(delg,ng,k_g1,q1,
     1					k_g2,q2,k_g)
                                else
                                 if(q2.gt.0)then
                                  do L=1,ng
                                   k_g(l)=k_g2(l)
                                  enddo                                  
                                 endif
C                                if(q1.gt.0)then leave k_g unchanged.
                                endif

				q1 = q1 + q2
C                                print*,'CC',q1,k_g(5)
			 ENDIF

C                        print*,'Final'
			DO L = 1, ng
				kl_g(L,J) = k_g(L)
C                                print*,L,j,kl_g(L,j)
			ENDDO
		  ENDDO

	        ELSE

                        KL_G(1,J)=0.

			DO K = 1, ngas

                         KL_G(1,J)=KL_G(1,J)+KOUTLBL(J,K)*FRAC(J,K)

			ENDDO

	        ENDIF

        ENDDO
C-----------------------------------------------------------------------
C
C	For each bin, we must step over each g-ordinate.
C
C	At each g ordinate, the optical depth of the layer is given by
C	the continuum value + k u where u is the integral over the 
C	absorber number along the path (.i.e. the sums of the 
C	constituent absorber numbers per unit area). Units of K are
C	(molecules number)**-1 cm**2 x 1.e20 (the later is taken care
C	of earlier in the code).
C
C-----------------------------------------------------------------------

C                print*,'NG = ',NG
		DO Ig = 1, ng
			DO J = 1, nlayer
				tautmp(J) = taucon(J) + kl_g(Ig,J) 
     1					* utotl(J)
				taugas(J) = taugasc(J) + kl_g(Ig,J) 
     1					* utotl(J)
			ENDDO
C-----------------------------------------------------------------------
C
C	Step through the number of paths, calculating the required 
C	output properties. These are
C
C	      Imod
C		0	(Atm) Pure transmission
C		1	(Atm) Absorption (useful for small transmissions)
C		2	(Atm) Emission. Planck function evaluated at each
C				wavenumber. NOT SUPPORTED HERE.
C		3	(Atm) Emission. Planck function evaluated at bin 
C				center.
C		8	(Combined Cell,Atm) The product of two
C				previous output paths.
C		11	(Atm) Contribution function.
C		13	(Atm) SCR Sideband
C		14	(Atm) SCR Wideband
C		15	(Atm) Multiple scattering (multiple models)
C		16	(Atm) Single scattering approximation.
C		21	(Atm) Net flux calculation (thermal)
C		22	(Atm) Limb scattering calculation
C		23	(Atm) Limb scattering calculation using precomputed
C			      internal radiation field.
C		24	(Atm) Net flux calculation (scattering)
C		25	(Atm) Upwards flux (internal) calculation (scattering)  
C		26	(Atm) Upwards flux (top) calculation (scattering)  
C		27	(Atm) Downwards flux (bottom) calculation (scattering)  
C		28	(Atm) Single scattering approximation (spherical).
C	then end the loop over the g ordinate.
C
C-----------------------------------------------------------------------

			DO Ipath = 1, npath
				nlays = nlayin(Ipath)
				corkout(Ipath,Ig) = 0.
				DO J = 1, nlays

		taus(J) = tautmp(layinc(J,ipath)) * scale(J,ipath)
		taur(J) = tauray(layinc(J,ipath)) * scale(J,ipath)


C               New arrays for scloud12wave. taucl holds the extinction
C		optical depth of each cloud type. taucs holds the 
C  		scattering optical depth of each cloud type. taug holds
C		the total gas extinction optical depth (excluding 
C		Rayleigh scattering).
                do k=1,ncont
                 taucl(k,j)=taucloud(k,layinc(J,ipath))*scale(J,ipath)
                 taucs(k,j)=tauclscat(k,layinc(J,ipath))*scale(J,ipath)
                enddo
 
		taug(J) = taugas(layinc(J,ipath)) * scale(J,ipath)
C               Subtract Rayleigh scattering
                taug(J)=taug(J)-taur(J)

                if(single)then
                  fint(j)=f(layinc(J,ipath))
		endif

		DO K = 1, ncont
			lfrac(K,J) = frcscat(K,layinc(J,ipath)) 
		ENDDO
				ENDDO

C        print*,'!! imod,ipath = ',imod,ipath
	
	IF (imod(ipath).EQ.0) THEN

C		WRITE(*,*) 'CIRSRAD_WAVE: Imod= 0 =Pure Transmission,',
C     1			' creating output'

		DO J = 1, nlays
			corkout(Ipath,Ig) = corkout(Ipath,Ig) + taus(J)
C                        print*,'J,tau',J,taus(J)
		ENDDO
		corkout(Ipath,Ig) = exp(-corkout(Ipath,Ig))

C               If a solar file exists and iform=4 we should multiply the
C               transmission by the solar flux
                if(iread.eq.999.and.iform.eq.4)then
                 call get_solar_wave(x,dist,xsolar)
                 corkout(Ipath,Ig)=corkout(Ipath,Ig)*xsolar
                endif


	ELSEIF (imod(ipath).EQ.1) THEN	

cc         WRITE(*,*) 'CIRSRAD_WAVE: Imod= 1 =Absorption, creating',
cc     1                  ' output'
 
		DO J = 1, nlays
			corkout(Ipath,Ig) = corkout(Ipath,Ig) + taus(J)
		enddo
		corkout(Ipath,Ig) = 1.- exp(-corkout(Ipath,Ig))

	ELSEIF (imod(ipath).EQ.3) THEN	

C		WRITE(*,*) 'CIRSRAD_WAVE: Imod= 3 =Emission, creating',
C     1                  ' output'

C               The code below calculates the radiance
C               spectrum in units of W cm-2 sr-1 (cm-1)-1 or W cm-2 sr-1 um-1.
C
                xfac=1.
C               If iform = 1 or iform = 3 we need to calculate the total
C               spectral flux from the planet.
                if(iform.eq.1.or.iform.eq.3)then
                  xfac=xfac*pi*4.*pi*((RADIUS1+radextra)*1e5)**2
                endif
C               If a solar file exists and iform=1 we should divide the planet
C               flux by the solar flux
                if(iread.eq.999.and.iform.eq.1)then
C                Set dist to -1 to get total power spectrum of star
                 xdist=-1.0
                 call get_solar_wave(x,xdist,xsolar)
                 xfac=xfac/xsolar
                endif
C               If doing integrated flux from planet need a factor to stop the
C               matrix inversion crashing
                if(iform.eq.3)xfac=xfac*1e-18

		taud = 0.
		trold = 1.
		do J = 1, nlays
			taud = taud + taus(J)
			tr = dexp(-taud)
                        if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace,x,
     1				emtemp(J,Ipath))
                        endif

		        corkout(Ipath,Ig) = corkout(Ipath,Ig) + 
     1                  xfac*sngl((trold-tr)) * bb(J,Ipath)
 			trold = tr
		enddo

C               Check to see if this is a limb path
                j1=int(0.5*nlays)

                p1=press(layinc(j1,ipath))
                p2=press(layinc(nlays,ipath))


C               If not a limb path, add on surface radiance
                if(p2.gt.p1)then


                 if(tsurf.le.0.0)then
                  radground = bb(nlays,Ipath)
                 else
                  radground = esurf*planck_wave(ispace,x,tsurf)
                 endif
                 corkout(Ipath,Ig) = corkout(Ipath,Ig) +
     1                  xfac*sngl(trold)*radground
                endif

        ELSEIF (imod(ipath).eq.8) THEN
C             model 8, product of two path outputs

 		corkout(ipath,Ig)=corkout(layinc(1,Ipath),Ig)*
     1              corkout(layinc(2,Ipath),Ig)

	ELSEIF (imod(ipath).EQ.11) THEN	

cc                WRITE(*,*) 'CIRSRAD_WAVE: Imod= 11 =Contribution Function,',
cc     1                  ' creating output'

		trold = 1.
		do J = 1, nlays
			taud = taus(J)
			tr = dexp(-taud)
                        if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace, 
     1			    x,emtemp(J,Ipath))
                        endif
			corkout(Ipath,Ig) = corkout(Ipath,Ig) + 
     1                  sngl((trold-tr)) * bb(J,Ipath)
 			trold = tr
                enddo

        ELSEIF (imod(ipath).eq.13) THEN

C          model 13, SCR sideband transmission (1-cell transmission)
			taud = taus(J)

           corkout(Ipath,Ig)=1.0-exp(-taus(1))

           LSTCEL=IPATH

        ELSEIF (imod(ipath).eq.14) THEN

C          model 14 SCR wideband transmission

           corkout(Ipath,Ig)=0.5*(1.0+exp(-taus(1)))

           LSTCEL=IPATH

	ELSEIF (imod(ipath).EQ.15) THEN

C           WRITE(*,*) 'CIRSRAD_WAVE: Imod= 15 =Multiple Scattering, ',
C     1                  ' creating output'

		do J = 1, nlays
			if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace,x,
     1			   emtemp(J,Ipath))
                        endif
			bnu(J) = bb(J,Ipath)
			IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
				dtmp1 = dble(tauscat(layinc(J,Ipath)))
				dtmp2 = dble(tautmp(layinc(J,Ipath)))
                                omegas(J)=dtmp1/dtmp2
				eps(J) = 1.0d00 - omegas(J)
		        ELSE
  				EPS(J)=1.0
         			omegaS(J)=0.
         		END IF


			fcover(J) = HFC(layinc(J,Ipath))
                        do k=1,ncont
                         icloud(k,j)=IFC(k,layinc(J,Ipath))
                        enddo
		enddo


C                         bb(J,Ipath)=planck_wave(ispace,x,
C     1                          emtemp(J,Ipath))

               
C               If tsurf > 0, then code assumes bottom layer is at the
C               surface and so uses tsurf to compute radiation upwelling.
C               Otherwise, code assumes optical depth is large and sets
C               upwelling radiation field to local temperature.

                if(tsurf.le.0.0)then
                 radground = bb(nlays,Ipath)
                else
                 radground = esurf*planck_wave(ispace,x,tsurf)
                endif

                galb1=galb

                if(galb1.lt.0)then
                         galb1 = dble(1.-esurf)
                endif

                do J=1,nmu
                 radg(J)=radground
                enddo

C               The codes below calculates the radiance
C               spectrum in units of W cm-2 sr-1 (cm-1)-1 or W cm-2 sr-1 um-1.
C
                xfac=1.
C               If iform = 1 or iform = 3 we need to calculate the total
C               spectral flux from the planet.
                if(iform.eq.1.or.iform.eq.3)then
                 xfac=xfac*pi*4.*pi*((RADIUS1+radextra)*1e5)**2
                endif
C               If a solar file exists and iform=1 we should divide the planet
C               flux by the solar flux
                if(iread.eq.999.and.iform.eq.1)then
C                Set dist to -1 to get total power spectrum of star
                 xdist=-1.0
                 call get_solar_wave(x,xdist,xsolar)
                 xfac=xfac/xsolar
                endif
C               If doing integrated flux from planet need a factor to stop the
C               matrix inversion crashing
                if(iform.eq.3)xfac=xfac*1e-18


                IF (itype.EQ.11) THEN

C                        WRITE (*,*) '     CALLING scloud11wave'
C                        WRITE (*,*)'IRAY,INORMAL = ',IRAY,INORMAL
C                        print*,sol_ang,solar
C                        print*,emiss_ang,aphi
C                        print*,radg
C                        print*,lowbc,galb1
C                        print*,nmu
C                        print*,(mu1(j),j=1,nmu)
C                        print*,(wt1(j),j=1,nmu)
C                        if(taur(1).gt.10000.0)print*,'taurA',
C     1  (taur(j),j=1,nlays)

C                        do j=1,nlays
C                         print*,j,eps(j),bnu(j),taus(j),taur(j),
C     1                           taus(j)-taur(j)
C                        enddo
C                        print*,'iray = ',iray

C                        print*,'galb1,solar,emiss_ang',galb1,solar,
C     1  emiss_ang
      		  	call scloud11wave(rad1, sol_ang, emiss_ang,
     1                          aphi, radg, solar, lowbc, galb1, iray,
     2				mu1, wt1, nmu,   
     3				nf, Ig, x, vv, eps, omegas,bnu, taus, 
     4				taur,nlays, ncont,lfrac)
C                        if(taur(1).gt.10000.0)print*,'taurB',
C     1 (taur(j),j=1,nlays)



 		  	corkout(Ipath,Ig) = xfac*rad1


                ELSEIF (itype.EQ.12) THEN

C                        WRITE (*,*) '     CALLING scloud12wave'
C                        WRITE (*,*)'IRAY,INORMAL = ',IRAY,INORMAL
C                        print*,'sol_ang,solar',sol_ang,solar
C                        print*,'emiss_ang,aphi',emiss_ang,aphi
C                        print*,'radg',radg
C                        print*,'lowbc,galb1',lowbc,galb1

			if(FLAGC.NE.1)THEN
                         print*,'Need to define a fractional cloud'
                         stop
                        endif

		  	call scloud12wave(rad1, sol_ang, emiss_ang, aphi, 
     1				radg, solar, lowbc, galb1, iray,
     2				mu1, wt1, nmu,   
     3				nf, Ig, x, vv, bnu, taucl, taucs,
     4                          icloud, fcover, taug, taur, nlays, 
     5                          ncont)

 		  	corkout(Ipath,Ig) = xfac*rad1

                ELSEIF (itype.EQ.13) THEN


   	  	   call scloud11flux(radg, solar, sol_ang, 
     1               lowbc, galb1, iray, mu1, wt1, nmu, nf, Ig, x,
     2               vv, eps, omegas,bnu, taus, taur, 
     3               nlays, ncont, lfrac, umif, uplf)
                    
                     open(48,file='test.dat',status='unknown',
     1                 form='unformatted')
                      write(48)umif
                      write(48)uplf
                     close(48)

                     tmp=1e10
                     j0=0
                     dtr = 3.1415927/180.0
                     xmu = cos(sol_ang*dtr)
                     do j=1,nmu
                      xt1 = abs(xmu-sngl(mu1(j)))
                      if(xt1.lt.tmp)then
                       j0=j
                       tmp=xt1
                      endif
                     enddo

                     solar1 = solar/(2.0*3.1415927*sngl(wt1(j0)))
                     
                     call dumpflux(LUNIS,ioff,nlays,nmu,nf,radg,umif,
     1                uplf,I,x,nwave,Ig,ng,j0,xmu,solar1,taus)

 		     corkout(Ipath,Ig) = umif(1,1,1)

                ELSE

			WRITE (*,*) ' Scattering ITYPE =',itype,
     1				' not defined for CIRSRAD_WAVE.'
			STOP
                ENDIF

	ELSEIF (imod(ipath).EQ.16) THEN

C		WRITE(*,*) 'CIRSRAD_WAVE: Imod= 16 =Single Scat Approx,',
C     1                  ' creating output'

		IF (Ig.EQ.1) WRITE(*,*)' SINGLE SCATTERING IN USE '

                galb1 = galb
                        
                if(galb1.lt.0)then
                   galb1 = dble(1.-esurf)
                  if(Ig.eq.1)print*,'x,galb1',x,galb1
                endif

                taud = 0.
                tausun=0.0
                trold = 1.


C               This code below calculates the radiance
C               spectrum in units of W cm-2 sr-1 (cm-1)-1 or W cm-2 sr-1 um-1.
C
                xfac=1.
C               If iform = 1 or iform = 3 we need to calculate the total
C               spectral flux from the planet.
                if(iform.eq.1.or.iform.eq.3)then
                 xfac=xfac*pi*4.*pi*((RADIUS1+radextra)*1e5)**2
                endif
C               If a solar file exists and iform=1 we should divide the planet
C               flux by the solar flux
                if(iread.eq.999.and.iform.eq.1)then
C                Set dist to -1 to get total power spectrum of star
                 xdist=-1.0
                 call get_solar_wave(vwave,xdist,xsolar)
                 xfac=xfac/xsolar
                endif
C               If doing integrated flux from planet need a factor to stop the
C               matrix inversion crashing
                if(iform.eq.3)xfac=xfac*1e-18


		DO J= 1, nlays

                 IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
	          tmp1 = tauscat(layinc(J,Ipath))
		  tmp2 = tautmp(layinc(J,Ipath))
		  omega(J) = tmp1/tmp2
	         ELSE
  		  omega(J)=0.0
                 END IF

                 taud = taud + taus(J)
                 tr = dexp(-taud)


		 corkout(Ipath,Ig) = corkout(Ipath,Ig) + 
     1			xfac*sngl((trold-tr))*ssfac*fint(J)*
     2			omega(J)*solar/(4*pi)


                 if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace,x,
     1				emtemp(J,Ipath))
                 endif

	         corkout(Ipath,Ig) = corkout(Ipath,Ig) + 
     1                  xfac*sngl((trold-tr)) * bb(J,Ipath)


                 trold = tr
		ENDDO

C               Add in reflectance from the ground
                corkout(Ipath,Ig)=corkout(Ipath,Ig) + 
     1                 xfac*sngl(trold*solar*muinc*galb1)/pi

C                print*,'galb,solar,muinc,xfac,trold',galb1,
C     1  solar,muinc,xfac,trold
                
C		WRITE (*,*) ' Calculated: ', Ig, corkout(Ipath,Ig)

	ELSEIF (imod(ipath).EQ.21)THEN

cc		WRITE(*,*) 'CIRSRAD_WAVE: Imod= 21 =Net flux calculation,',
cc     1        'thermal emission only. Creating output'

C		Computes up and down flux at BOTTOM of each layer
C		Compute once for Ipath=Npath
C		Increase pathlength by 5/3 to approximate for
C               Hemispherical integration (See Houghton, 1986, p13)

		if(Ipath.eq.npath)then
                  do j=1,nlayin(npath)
                   fup(Ig,j)=0.0
                   fdown(Ig,j)=0.0
                   if(scale(J,npath).gt.1.01)then
                    print*,'Error in cirsrad_wave, layer properties'
		    print*,'should calculated at zero emission'
   		    print*,'angle.'
                    stop
                   endif
                   taus(layinc(J,npath))=tautmp(layinc(J,npath)) *
     &			scale(J,npath)*5./3.
		   if(Ig.eq.1)then
		     bb(layinc(j,npath),1)=planck_wave(ispace,
     1                  x,emtemp(j,npath))
		   endif
		  enddo

C                 If tsurf > 0, then code assumes bottom layer is at the
C                 surface and so uses tsurf to compute radiation upwelling.
C                 Otherwise, code assumes optical depth is large and sets
C                 upwelling radiation field to local temperature.
                  if(tsurf.le.0.0)then
                   radground = bb(1,1)
                  else
                   radground = esurf*planck_wave(ispace,x,tsurf)
                  endif

	   	  do j=1,nlayin(npath)
           	   if(j.eq.1)then
		    fup(j,Ig)=pi*radground
		   else
                    taud=taus(j-1)
		    tr = dexp(-taud)
                    fup(j,Ig)=fup(j-1,Ig)*sngl(tr)+sngl((1.0-tr))*
     1			pi*bb(j,1)
		   endif
		  enddo 

	   	  do j=nlayin(npath),1,-1
                   taud=taus(j)
                   tr=dexp(-taud)
            	   if(j.eq.nlays)then
		    fdown(j,Ig)=sngl((1.0-tr))*pi*bb(j,1)
		   else
                    fdown(j,Ig)=sngl((1.0-tr))*pi*bb(j,1)+
     1			sngl(tr)*fdown(j+1,Ig)
		   endif
		  enddo 

C                  open(12,file='netflux_therm.dat',status='unknown')
C                   write(12,*)Ig
C                   do j=1,nlayin(npath)
C                    write(12,*)j,taus(j),fup(J,Ig),fdown(J,Ig)
C                   enddo
C                 close(12)

C                 if(Ig.eq.5)stop

                 do Jpath=1,Npath
 	          corkout(Jpath,Ig) = fup(Jpath,Ig) - fdown(Jpath,Ig)
                 enddo

                endif 

            
 	
	ELSEIF (imod(ipath).EQ.22) THEN

cc                WRITE(*,*) 'CIRSRAD_WAVE: Imod= 22 =Limb Scattering, ',
cc     1                  ' creating output'
C                print*,ispace,x
C       ************************************************************
C       First we need to compute the internal field at ALL altitudes
C       ************************************************************
		do J1 = 1, nlayer
                        J = nlayer+1-J1
			if(Ig.eq.1)then
                         bnu(J) = TEMP(J1)
                         bnu(J) = planck_wave(ispace,x,TEMP(J1))
                        endif
			IF(TAUSCAT(J1).GT.0.0) THEN
				dtmp1 = dble(tauscat(J1))
				dtmp2 = dble(tautmp(J1))
                                omegas(J)=dtmp1/dtmp2
				eps(J) = 1.0d00 - omegas(J)
		        ELSE
  				EPS(J)=1.0
				omegas(J)=0.
         		END IF
		enddo

C               If tsurf > 0, then code assumes bottom layer is at the
C               surface and so uses tsurf to compute radiation upwelling.
C               Otherwise, code assumes optical depth is large and sets
C               upwelling radiation field to local temperature.
                if(tsurf.le.0.0)then
                 radground = bnu(nlayer)
                else
                 radground = esurf*planck_wave(ispace,x,tsurf)
                endif

                galb1=galb

                if(galb1.lt.0)then
                         galb1 = dble(1.-esurf)
                endif

                do J=1,nmu
                 radg(J)=radground
                enddo


C                WRITE (*,*) '     CALLING scloud11flux, nf = ',nf
		call scloud11flux(radg, solar, sol_ang, 
     1             lowbc, galb1, iray, mu1, wt1, nmu, nf, Ig, x,
     2             vv, eps, omegas, bnu, taus, taur, 
     3             nlayer, ncont, lfrac, umif, uplf)

C                do k=1,nlays
C                 print*,k,press(k),(umif(j,k,1),j=1,nmu)
C                 print*,k,press(k),(uplf(j,k,1),j=1,nmu)
C                enddo

                do J = 1, nlays
                        if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace,x,
     1                     emtemp(J,Ipath))
                        endif
                        bnu(J) = bb(J,Ipath)
                        IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
                                dtmp1 = dble(tauscat(layinc(J,Ipath)))
                                dtmp2 = dble(tautmp(layinc(J,Ipath)))
				omegas(J)=dtmp1/dtmp2
                                eps(J) = 1.0d00 - omegas(J)
                        ELSE
                                EPS(J)=1.0
				omegas(J)=0.
                        END IF
                enddo

                call hgaverage(nlays,ncont,taus,eps,lfrac,
     1           taur,x,f1,g11,g21)

                call limbscatter(nlays,nmu,wt1,mu1,eps,bnu,
     1                taus,f1,g11,g21,umif,uplf,nf,rad1)

 		corkout(Ipath,Ig) = rad1

	ELSEIF (imod(ipath).EQ.23) THEN

                WRITE(*,*) 'CIRSRAD_WAVE: Imod= 23 = modified source, ',
     1                  ' creating output'
C                print*,ispace,x
C       ************************************************************
C       First need to read in the internal field at ALL altitudes
C       ************************************************************

                if(isol.eq.0) then
                 print*,'Importing azimuth-integrated flux'
                 call impflux(LUNIS,ioff,nlayerf,nmuf,nff,umif,uplf,
     1                x,nwavef,vwavef,Ig,ng)
                else
                 print*,'Importing total internal flux field'
                 call impfluxsol(LUNIS,ioff,nlayerf,nmuf,nff,umift,
     1                uplft,x,nwavef,vwavef,Ig,ng)


                endif

                do J = 1, nlays
                        if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace,x,
     1                     emtemp(J,Ipath))
                        endif
                        bnuS(J) = bb(J,Ipath)
                        IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
                                dtmp1 = dble(tauscat(layinc(J,Ipath)))
                                dtmp2 = dble(tautmp(layinc(J,Ipath)))
                                epsS(J) = 1.0d00 - dtmp1/dtmp2
                        ELSE
                                epsS(J)=1.0
                        END IF

                        basehS(J) = baseh(layinc(J,Ipath))
                        ibaseH(J) = layinc(J,Ipath) 
                        scaleS(J) = scale(J,Ipath)
			delhs(J) = delh(layinc(J,Ipath))

                enddo

                call hgaverage(nlays,ncont,taus,epsS,lfrac,
     1           taur,x,f1,g11,g21)

                if(isol.eq.0)then
                 call scatsource(nlays,nmuf,wt1,mu1,epsS,bnuS,
     1            taus,basehS,scaleS,f1,g11,g21,nlayerf,basehf,
     2            umif,uplf,nff,Jsource,Ig,I)
                else
C                 print*,'Calling scatsourcesol'
                 call scatsourcesol(nlays,nmuf,wt1,mu1,epsS,bnuS,
     1            basehS,ibaseH,tautmp,emiss_ang,sol_ang,aphi,solar,
     2            f1,g11,g21,nlayerf,basehf,umift,uplft,nff,Jsource,
     3            Ig,I)
                endif

                taud = 0.
                trold = 1.
                do J = 1, nlays
                        taud = taud + taus(J)
                        tr = dexp(-taud)
                        corkout(Ipath,Ig) = corkout(Ipath,Ig) +
     1                  sngl((trold-tr)) * Jsource(J)
                        trold = tr
                enddo
C               If surface temperature defined and bottom layer
C               at zero altitude, assume we have to add on planck
C               emission from ground.
                if(tsurf.gt.0.0.and.basehS(nlays).eq.0.0)then
                 radground = esurf*planck_wave(ispace,x,tsurf)
                 corkout(Ipath,Ig) = corkout(Ipath,Ig) +
     1            sngl(trold) * radground
                endif

	ELSEIF (imod(ipath).EQ.24) THEN

cc		WRITE(*,*) 'CIRSRAD_WAVE: Imod= 24 =Net flux calculation,',
cc     1        ' multiple scattering'. Creating output'

                if(Ipath.eq.Npath)then

                 nlays = nlayin(Ipath)

		 do J = 1, nlays
			if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace,x,
     1			   emtemp(J,Ipath))
                        endif
			bnu(J) = bb(J,Ipath)
			IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
				dtmp1 = dble(tauscat(layinc(J,Ipath)))
				dtmp2 = dble(tautmp(layinc(J,Ipath)))
				omegas(J)=dtmp1/dtmp2
				eps(J) = 1.0d00 - omegas(J)
		        ELSE
  				EPS(J)=1.0
				omegas(J)=0.
         		END IF
                        EPSB(J)=1.0
                        omegasb(J)=0.

			fcover(J) = HFC(layinc(J,Ipath))
                        do k=1,ncont
                         icloud(k,j)=IFC(k,layinc(J,Ipath))
                        enddo
		  enddo


C                 If tsurf > 0, then code assumes bottom layer is at the
C                 surface and so uses tsurf to compute radiation upwelling.
C                 Otherwise, code assumes optical depth is large and sets
C                 upwelling radiation field to local temperature.

                  if(tsurf.le.0.0)then
                    radground = bb(nlays,Ipath)
                  else
                    radground = esurf*planck_wave(ispace,x,tsurf)
                  endif
                  print*,'tsurf,radgound',tsurf,radground

                  galb1=galb

                  if(galb1.lt.0)then
                         galb1 = dble(1.-esurf)
                  endif

                  do J=1,nmu
                    radg(J)=radground
                  enddo

C		  Calculate radiation field for thermal emission only, 
C                  but scattering clouds
	          solarA=0.
                  nfx=0

   	    	  call scloud11flux(radg, solarA, sol_ang, 
     1               lowbc, galb1, iray, mu1, wt1, nmu, nfx, Ig, x,
     2               vv, eps, omegas, bnu, taus, taur, 
     3               nlays, ncont, lfrac, umif, uplf)

                  do imu=1,nmu
                   do ilays=1,nlays
                     if(ig.eq.1)then
                      gplfa(imu,ilays)=delg(ig)*uplf(imu,ilays,1)
                      gmifa(imu,ilays)=delg(ig)*umif(imu,ilays,1)
                     else
                      gplfa(imu,ilays)=gplfa(imu,ilays)+
     1                  delg(ig)*uplf(imu,ilays,1)
                      gmifa(imu,ilays)=gmifa(imu,ilays)+
     1                  delg(ig)*umif(imu,ilays,1)
                     endif
                   enddo
                  enddo

C		  Calculate radiation field for thermal emission only, 
C                  but non-scattering clouds
	          solarA=0.
                  nfx=0
                  irayx=0

   	    	  call scloud11flux(radg, solarA, sol_ang, 
     1               lowbc, galb1, irayx, mu1, wt1, nmu, nfx, Ig, x,
     2               vv, epsb, omegasb, bnu, taus, taur, 
     3               nlays, ncont, lfrac, umif, uplf)

                  do imu=1,nmu
                   do ilays=1,nlays
                     if(ig.eq.1)then
                      gplfc(imu,ilays)=delg(ig)*uplf(imu,ilays,1)
                      gmifc(imu,ilays)=delg(ig)*umif(imu,ilays,1)
                     else
                      gplfc(imu,ilays)=gplfc(imu,ilays)+
     1                  delg(ig)*uplf(imu,ilays,1)
                      gmifc(imu,ilays)=gmifc(imu,ilays)+
     1                  delg(ig)*umif(imu,ilays,1)
                     endif
                   enddo
                  enddo


C		  Calculate radiation field with sunlight on but for
C                 for non-scattering conditions to get direct solar component.

   	    	  call scloud11flux(radg, solar, sol_ang, 
     1               lowbc, galb1, irayx, mu1, wt1, nmu, nf, Ig, x,
     2               vv, epsb, omegasb, bnu, taus, taur, 
     3               nlays, ncont, lfrac, umif, uplf)

                  do imu=1,nmu
                   do ilays=1,nlays
                    do ic=1,nf+1
                     if(ig.eq.1)then
                      gplfb(imu,ilays,ic)=delg(ig)*uplf(imu,ilays,ic)
                      gmifb(imu,ilays,ic)=delg(ig)*umif(imu,ilays,ic)
                     else
                      gplfb(imu,ilays,ic)=gplfb(imu,ilays,ic)+
     1                  delg(ig)*uplf(imu,ilays,ic)
                      gmifb(imu,ilays,ic)=gmifb(imu,ilays,ic)+
     1                  delg(ig)*umif(imu,ilays,ic)
                     endif
                    enddo
                   enddo
                  enddo


C	          Finally, calculate radiation field for nominal 
C		   scattering case


   	    	  call scloud11flux(radg, solar, sol_ang, 
     1               lowbc, galb1, iray, mu1, wt1, nmu, nf, Ig, x,
     2               vv, eps, omegas, bnu, taus, taur, 
     3               nlays, ncont, lfrac, umif, uplf)

                  do imu=1,nmu
                   do ilays=1,nlays
                    do ic=1,nf+1
                     if(ig.eq.1)then
                      gplf(imu,ilays,ic)=delg(ig)*uplf(imu,ilays,ic)
                      gmif(imu,ilays,ic)=delg(ig)*umif(imu,ilays,ic)
                     else
                      gplf(imu,ilays,ic)=gplf(imu,ilays,ic)+
     1                  delg(ig)*uplf(imu,ilays,ic)
                      gmif(imu,ilays,ic)=gmif(imu,ilays,ic)+
     1                  delg(ig)*umif(imu,ilays,ic)
                     endif
                    enddo
                   enddo
                  enddo


                  call streamflux(nlays,nmu,mu1,wt1,radg,galb1,
     1                ig,umif,uplf,fup,fdown,ftop)                    



                  do Jpath=1,Ipath
  	           corkout(Jpath,Ig) = fup(Jpath,Ig) - fdown(Jpath,Ig)
                  enddo

                endif

	ELSEIF (imod(ipath).EQ.25) THEN

cc		WRITE(*,*) 'CIRSRAD_WAVE: Imod= 25 = Upwards flux calculation,',
cc     1        ' multiple scattering'. Creating output'

                if(Ipath.eq.Npath)then

                 nlays = nlayin(Ipath)

		 do J = 1, nlays
			if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace,x,
     1			   emtemp(J,Ipath))
                        endif
			bnu(J) = bb(J,Ipath)
			IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
				dtmp1 = dble(tauscat(layinc(J,Ipath)))
				dtmp2 = dble(tautmp(layinc(J,Ipath)))
				omegas(J)=dtmp1/dtmp2
				eps(J) = 1.0d00 - omegas(J)
		        ELSE
  				EPS(J)=1.0
				omegas(J)=0.
         		END IF


			fcover(J) = HFC(layinc(J,Ipath))
                        do k=1,ncont
                         icloud(k,j)=IFC(k,layinc(J,Ipath))
                        enddo
		  enddo


C                 If tsurf > 0, then code assumes bottom layer is at the
C                 surface and so uses tsurf to compute radiation upwelling.
C                 Otherwise, code assumes optical depth is large and sets
C                 upwelling radiation field to local temperature.

                  if(tsurf.le.0.0)then
                    radground = bb(nlays,Ipath)
                  else
                    radground = esurf*planck_wave(ispace,x,tsurf)
                  endif

                  galb1=galb

                  if(galb1.lt.0)then
                         galb1 = dble(1.-esurf)
                  endif

                  do J=1,nmu
                    radg(J)=radground
                  enddo

   	    	  call scloud11flux(radg, solar, sol_ang, 
     1               lowbc, galb1, iray, mu1, wt1, nmu, nf, Ig, x,
     2               vv, eps, omegas, bnu, taus, taur, 
     3               nlays, ncont, lfrac, umif, uplf)

                  call streamflux(nlays,nmu,mu1,wt1,radg,galb1,
     1                ig,umif,uplf,fup,fdown,ftop)                    


                  do Jpath=1,Ipath
  	           corkout(Jpath,Ig) = fup(Jpath,Ig)
                  enddo

                endif

	ELSEIF (imod(ipath).EQ.26) THEN
cc		WRITE(*,*) 'CIRSRAD_WAVE: Imod= 26 = Upwards flux at top,',
cc     1        ' multiple scattering'. Creating output'


		do J = 1, nlays
			if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace,x,
     1			   emtemp(J,Ipath))
                        endif
			bnu(J) = bb(J,Ipath)
			IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
				dtmp1 = dble(tauscat(layinc(J,Ipath)))
				dtmp2 = dble(tautmp(layinc(J,Ipath)))
                                omegas(J)=dtmp1/dtmp2
				eps(J) = 1.0d00 - omegas(J)
		        ELSE
  				EPS(J)=1.0
         			omegaS(J)=0.
         		END IF


			fcover(J) = HFC(layinc(J,Ipath))
                        do k=1,ncont
                         icloud(k,j)=IFC(k,layinc(J,Ipath))
                        enddo
		enddo



               
C               If tsurf > 0, then code assumes bottom layer is at the
C               surface and so uses tsurf to compute radiation upwelling.
C               Otherwise, code assumes optical depth is large and sets
C               upwelling radiation field to local temperature.

                if(tsurf.le.0.0)then
                 radground = bb(nlays,Ipath)
                else
                 radground = esurf*planck_wave(ispace,x,tsurf)
                endif

                galb1=galb

                if(galb1.lt.0)then
                         galb1 = dble(1.-esurf)
                endif

                do J=1,nmu
                 radg(J)=radground
                enddo

C               The codes below calculates the upwards flux
C               spectrum in units of W cm-2 (cm-1)-1 or W cm-2 um-1.
C
                xfac=1.
C               If iform = 1 or iform = 3 we need to calculate the total
C               spectral flux from the planet.
                if(iform.eq.1.or.iform.eq.3)then
                 xfac=xfac*4.*pi*((RADIUS1+radextra)*1e5)**2
                endif
C               If a solar file exists and iform=1 we should divide the planet
C               flux by the solar flux
                if(iread.eq.999.and.iform.eq.1)then
C                Set dist to -1 to get total power spectrum of star
                 xdist=-1.0
                 call get_solar_wave(x,xdist,xsolar)
                 xfac=xfac/xsolar
                endif
C               If doing integrated flux from planet need a factor to stop the
C               matrix inversion crashing
                if(iform.eq.3)xfac=xfac*1e-18


   	   	call scloud11flux(radg, solar, sol_ang, 
     1               lowbc, galb1, iray, mu1, wt1, nmu, nf, Ig, x,
     2               vv, eps, omegas, bnu, taus, taur, 
     3               nlays, ncont, lfrac, umif, uplf)

                call streamflux(nlays,nmu,mu1,wt1,radg,galb1,
     1                ig,umif,uplf,fup,fdown,ftop)                    

 		corkout(Ipath,Ig) = xfac*ftop(ig)


	ELSEIF (imod(ipath).EQ.27) THEN
cc		WRITE(*,*) 'CIRSRAD_WAVE: Imod= 27 = Downwards flux at bottom,',
cc     1        ' multiple scattering'. Creating output'


		do J = 1, nlays
			if(Ig.eq.1)then
                         bb(J,Ipath)=planck_wave(ispace,x,
     1			   emtemp(J,Ipath))
                        endif
			bnu(J) = bb(J,Ipath)
			IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
				dtmp1 = dble(tauscat(layinc(J,Ipath)))
				dtmp2 = dble(tautmp(layinc(J,Ipath)))
                                omegas(J)=dtmp1/dtmp2
				eps(J) = 1.0d00 - omegas(J)
		        ELSE
  				EPS(J)=1.0
         			omegaS(J)=0.
         		END IF


			fcover(J) = HFC(layinc(J,Ipath))
                        do k=1,ncont
                         icloud(k,j)=IFC(k,layinc(J,Ipath))
                        enddo
		enddo



               
C               If tsurf > 0, then code assumes bottom layer is at the
C               surface and so uses tsurf to compute radiation upwelling.
C               Otherwise, code assumes optical depth is large and sets
C               upwelling radiation field to local temperature.

                if(tsurf.le.0.0)then
                 radground = bb(nlays,Ipath)
                else
                 radground = esurf*planck_wave(ispace,x,tsurf)
                endif

                galb1=galb

                if(galb1.lt.0)then
                         galb1 = dble(1.-esurf)
                endif

                do J=1,nmu
                 radg(J)=radground
                enddo

C               The codes below calculates the downwards flux
C               spectrum in units of W cm-2 (cm-1)-1 or W cm-2 um-1.
C
                xfac=1.
C               If iform = 1 or iform = 3 we need to calculate the total
C               spectral flux from the planet.
                if(iform.eq.1.or.iform.eq.3)then
                 xfac=xfac*4.*pi*((RADIUS1+radextra)*1e5)**2
                endif
C               If a solar file exists and iform=1 we should divide the planet
C               flux by the solar flux
                if(iread.eq.999.and.iform.eq.1)then
C                Set dist to -1 to get total power spectrum of star
                 xdist=-1.0
                 call get_solar_wave(x,xdist,xsolar)
                 xfac=xfac/xsolar
                endif
C               If doing integrated flux from planet need a factor to stop the
C               matrix inversion crashing
                if(iform.eq.3)xfac=xfac*1e-18

   	   	call scloud11flux(radg, solar, sol_ang, 
     1               lowbc, galb1, iray, mu1, wt1, nmu, nf, Ig, x,
     2               vv, eps, omegas, bnu, taus, taur, 
     3               nlays, ncont, lfrac, umif, uplf)

                call streamflux(nlays,nmu,mu1,wt1,radg,galb1,
     1                ig,umif,uplf,fup,fdown,ftop)                    

 		corkout(Ipath,Ig) = xfac*fdown(1,ig)

	ELSEIF (imod(ipath).EQ.28) THEN

		WRITE(*,*) 'CIRSRAD_WAVE: Imod= 28 =Single Scat Approx,',
     1                  ' (Spherical atmosphere) creating output'

C		IF (Ig.EQ.1) WRITE(*,*)' SINGLE SCATTERING IN USE '

                galb1 = galb
                        
                if(galb1.lt.0)then
                   galb1 = dble(1.-esurf)
                   if(Ig.eq.1)print*,'x,galb1',x,galb1
                endif

                taud = 0.
                tausun=0.0
                trold = 1.


C               This code below calculates the radiance
C               spectrum in units of W cm-2 sr-1 (cm-1)-1 or W cm-2 sr-1 um-1.
C
                xfac=1.
C               If iform = 1 or iform = 3 we need to calculate the total
C               spectral flux from the planet.
                if(iform.eq.1.or.iform.eq.3)then
                 xfac=xfac*pi*4.*pi*((RADIUS1+radextra)*1e5)**2
                endif
C               If a solar file exists and iform=1 we should divide the planet
C               flux by the solar flux
                if(iread.eq.999.and.iform.eq.1)then
C                Set dist to -1 to get total power spectrum of star
                 xdist=-1.0
                 call get_solar_wave(vwave,xdist,xsolar)
                 xfac=xfac/xsolar
                endif
C               If doing integrated flux from planet need a factor to stop the
C               matrix inversion crashing
                if(iform.eq.3)xfac=xfac*1e-18

		DO J= 1, nlays
                 IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
	          tmp1 = tauscat(layinc(J,Ipath))
		  tmp2 = tautmp(layinc(J,Ipath))
		  omega(J) = tmp1/tmp2
	         ELSE
  		  omega(J)=0.0
                 END IF
                 refrac(J)=1.0
                 rad(J)=RADIUS1+baseh(J)
                ENDDO
                rad(nlays+1)=RADIUS1+baseh(nlays)+delH(nlays)
                ethick=(rad(2)-rad(1))/10.

C               Compute path properties along line of sight
                startlay=nlays+1
                j=nlays/2
                scaleH=-delH(j)/log(basep(j+1)/basep(j))

C               Set up initial ray position and vector. Also set up 
C               direction to Sun
                call raystart(emiss_ang,sol_ang,aphi,ipzen,
     1           nlays,rad,r0,vec,svec)
              
                call traceray(nlays,rad,refrac,startlay,r0,vec,ethick,
     1             scaleH,npathLOS,traceLOS,jtraceLOS,disttrace,
     2             scaletraceLOS,vectraceLOS)


C               Compute transmission along line of sight
                tauLOS=0.
                trold=1.
                DO J=1,npathLOS-1
                 J1=jtraceLOS(J)
                 J2=jtraceLOS(J+1)
C                If ray going downwards need to take opacity of layer 
C                with ordinate equal to boundary at end of this subsection.
C                If ray going upwardswards need to take opacity of layer 
C                with ordinate equal to boundary at start of this subsection.
 
                 IF(J1.GE.J2)THEN
                  JLAY=J2
                 ELSE
                  JLAY=J1
                 ENDIF                 
                 ILAY=1+NLAYS-JLAY

               
                 tauLOS=tauLOS+tautmp(JLAY)*
     1			   scaletraceLOS(J+1)

C                Compute solar transmission to this point.

                 do L=1,3
                   r0(L)=traceLOS(J+1,L)
                   vec1(L)=vectraceLOS(J+1,L)
                   vec(L)=svec(L)
                 enddo
C                Find local scattering angle
                 calpha=sproduct(vec1,svec)

                 ff=0.
                 sum=0
                 do icont=1,ncont
                       CALL get_hg(x,calpha,ncont,icont,dphase)
                       phase(icont) = sngl(dphase)
                       ff=ff+phase(icont)*tauclscat(icont,JLAY)
                       sum=sum+tauclscat(icont,JLAY)
                 enddo
                 IF(IRAY.GT.0)THEN
C                  Calculate Rayleigh scattering too
                   phase(ncont+1)=sngl(0.75*(1+calpha*calpha)) 
                   ff=ff+phase(ncont+1)*tauray(JLAY)
                   sum=sum+tauray(JLAY)
                 ENDIF
                 if(sum.gt.0)then
                  ff=ff/sum
                 else
                  ff=0.
                 endif
                 startlay=jtraceLOS(J+1)

                 call traceray(nlays,rad,refrac,startlay,r0,vec,ethick,
     1             scaleH,npathsol,trace,jtrace,disttrace,scaletrace,
     2		   vectrace)
                 tausol=0.
                 DO L=1,NPATHSOL-1
                   L1=JTRACE(L)
                   L2=JTRACE(L+1)
                   IF(L1.GE.L2)THEN
                    KLAY=L2                 
                   ELSE
                    KLAY=L1
                   ENDIF
                   tausol=tausol+tautmp(KLAY)*
     1				scaletrace(L+1)
                 ENDDO
                 tr=exp(-(tauLOS+tausol)) 
                 DO L=1,3
                  r1(L)=traceLOS(J+1,L)
                 ENDDO
                 rtmp=sqrt(sproduct(r1,r1))                 
                 print*,'rtmp-RADIUS1 = ',rtmp-RADIUS1
                 
                 muinc = abs(sproduct(svec,r1)/rtmp)
                 muemiss = abs(sproduct(vec1,r1)/rtmp)
                 print*,muinc,muemiss
                 dtr=3.1415927/180.0
                 print*,J,acos(muinc)/dtr,acos(muemiss)/dtr
                 ssfac = muinc/(muinc+muemiss)


 		 corkout(Ipath,Ig) = corkout(Ipath,Ig) + 
     1			xfac*sngl((trold-tr))*ff*ssfac*
     2			omega(ILAY)*solar/(4*pi)

                 if(Ig.eq.1)then
                          bb(J,Ipath)=planck_wave(ispace,x,
     1				emtemp(ILAY,Ipath))
C                          print*,JLAY,ILAY,emtemp(ILAY,Ipath)
                 endif

 	         corkout(Ipath,Ig) = corkout(Ipath,Ig) + 
     1                  xfac*sngl((trold-tr)) * bb(J,Ipath)


                 trold = tr
                 
                ENDDO                 

C               Add in reflectance from the ground

                IF(npathLOS.EQ.nlays+1)THEN

C                print*,'npathLOS',npathLOS
C                print*,'galb,solar,muinc,xfac,trold',galb1,
C     1  solar,muinc,xfac,trold

                 corkout(Ipath,Ig)=corkout(Ipath,Ig) + 
     1                 xfac*sngl(trold*solar*muinc*galb1)/pi
                ENDIF
C		WRITE (*,*) ' Calculated: ', Ig, corkout(Ipath,Ig)


	ELSE
		WRITE (*,*) ' Imod = ', imod(ipath), ' is not a valid',
     1			' CIRSRAD_WAVE option'
		STOP

	ENDIF
			ENDDO
		ENDDO

C-----------------------------------------------------------------------
C
C	Now integrate over g ordinates and then end loop over bins.
C
C-----------------------------------------------------------------------


		DO Ipath = 1, npath
			output(Ipath,I) = 0.
			DO Ig = 1, ng
			 output(Ipath,I) = output(Ipath,I) + 
     1					corkout(Ipath,Ig) * delg(Ig)
			ENDDO
C		 	print*,'Ipath, I, output',Ipath,I,
C     1				output(Ipath,I)
		ENDDO

                if(ish.eq.1)then
                 write(infr,*)solar
                 write(infr,*)radground
                 write(infr,*)galb1
                 do imu=1,nmu
                  do ilays=1,nlays
                   do ic=1,nf+1
                     write(infr,*)gplf(imu,ilays,ic)
                   enddo
                   do ic=1,nf+1
                     write(infr,*)gmif(imu,ilays,ic)
                   enddo
                   do ic=1,nf+1
                     write(infr,*)gplfb(imu,ilays,ic)
                   enddo
                   do ic=1,nf+1
                     write(infr,*)gmifb(imu,ilays,ic)
                   enddo
                   write(infr,*)gplfa(imu,ilays)
                   write(infr,*)gmifa(imu,ilays)
                   write(infr,*)gplfc(imu,ilays)
                   write(infr,*)gmifc(imu,ilays)
                  enddo
                 enddo
                endif

		XCOM = 100.0*FLOAT(I)/FLOAT(NWAVE)
		IF(XCOM.GE.XNEXT)THEN
C			WRITE(*,1000)
C			WRITE(*,1010) I, nwave, x
C			WRITE(*,1020)XCOM

			XNEXT = XNEXT+10.0
		ENDIF

	ENDDO

C-----------------------------------------------------------------------
C
C	Return and end.
C
C-----------------------------------------------------------------------

1000	FORMAT (/,'***************************************************')
1010	FORMAT (' Ending wavenumber loop',I6,' of',I6,' : ',
     1   F8.3,' cm^-1')
1020	FORMAT(' CIRSRAD_WAVE.f :: % complete : ',F5.1)

        IF(IFLUX.GT.0)THEN
         CLOSE(LUNIS)
        ENDIF

        if(ish.eq.1)then
         close(infr)
        endif

	RETURN

	END

************************************************************************
************************************************************************
