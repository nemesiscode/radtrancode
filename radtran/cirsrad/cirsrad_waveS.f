************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C_TITLE:		SUBROUTINE CIRSRAD_WAVES
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
C       DIST:REAL*4 Distance from Sun (as prompted by CIRSDRV) in
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
C			(currently only scloud8 through scloud11).
C
C
C	Dust variables
C
C	NCONT		Number of dust types included
C       CONT		The number of dust particles/cm2 through each 
C			layer (vertical path).
C       NSEC		Number of wavelengths for which the dust cross-
C			sections are defined.
C	XSEC		Dust cross sections for each dust type
C	VSEC		Corresponding wavelengths for dust x-sections
C
C
C	Output variables
C
C       OUTPUT		Output values at each wavenumber for each output 
C			type for each path
C
C_HIST:	Fall1999 PDP	Modified original version from NIMSrad.
C	20.2.2000 PDP	Replace the h2con program with ciacon and the
C			associated changes that allowed ciacon to work
C			properly.
C	26nov02	PDP	Defined "single=.false." so that code didn't go
C                       off and look for any Henyey-Greenstein phase
C                       files automatically. Required if using the Intel
C			FORTRAN compiler.
C-----------------------------------------------------------------------

      SUBROUTINE cirsrad_waveS (Dist, INormal,iray,ispace,DelH,nlayer,
     1    npath, ngas, limlay, limcont, press, temp, pp, amount, nwave,
     2    vwave, nlayin, incdim, layinc, cont, scale, imod,
     3    idgas, isogas, emtemp, itype, nem, vem, emissivity, tsurf, 
     4    flagh2p, hfp, flagc, hfc, ifc, basep, baseh, output)

      IMPLICIT NONE


C		Internal dimensions

C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)
       INCLUDE '../includes/arrdef.f'

C		Passed variables

	INTEGER	Inormal,ik,iray,flagh2p,flagc,ispace,jf
	INTEGER	nlayer, npath, ngas, limlay, limcont, nwave, 
     1		nlayin(npath), incdim, layinc(incdim, npath), 
     2          idgas(ngas), isogas(ngas), imod(npath)

	REAL	press(nlayer), temp(nlayer), pp(limlay,ngas), 
     1		amount(limlay,ngas), vwave(nwave), cont(limcont,nlayer), 
     2		scale(incdim,npath), emtemp(incdim,npath), hfp(nlayer), 
     3		output(npath,nwave),bb(maxlay,maxpat),xf,tsurf,dv,
     4          hfc(nlayer),vwavef(maxbin),basep(nlayer),baseh(nlayer),
     5          jb(maxlay,maxpat),esurf
        REAL    vem(MAXSEC),emissivity(MAXSEC),interpem
        INTEGER ifc(limcont,nlayer),nem
        LOGICAL modsource

	REAL	XCOM,XNEXT,FPARA,vv,XRAY,RAYLEIGHJ,RAYLEIGHA,RAYLEIGHV
        double precision get_albedo

C		Dust variables

	INTEGER	nsec, ncont
	REAL	vsec(maxsec), xsec(2,maxcon,maxsec)

C		K table variables

	INTEGER	lun(maxbin,maxgas), ireck(maxbin,maxgas), lun0, irec0, 
     1		npk, ntk, ng, intmod(maxbin,maxgas), ntab
	REAL	xmink(maxbin,maxgas), delk(maxbin,maxgas), xmin, delx, 
     1		k_g(maxg), k_g1(maxg), k_g2(maxg), q1, q2, 
     2		pk(maxk), tk(maxk), g_ord(maxg), 
     3		delg(maxg), kl_g(maxg,maxlay),
     4		kout(maxlay,maxgas,maxg)

        REAL    basehf(maxlay),basepf(maxlay),basehS(maxlay)
        REAL    scaleS(maxlay),jsource(maxlay)
C		Scattering variables

	INTEGER	liscat(maxcon), isol, nmu, lowbc, 
     1          lncons(maxcon), 
     2		nf, lnorm(maxcon), itype
	REAL	lfrac(maxcon,maxlay), solar, pi, theta0, dist, 
     1		bnu(maxscatlay), aphi, radg(maxmu),radground, 
     2		lcons(maxcon,maxscatpar), rad1, omega(maxscatlay),
     3		tsun, v, frcscat(maxcon,maxlay), sol_ang, emiss_ang
        REAl    bnuS(maxlay)
        REAL    press1(maxscatlay),temp1(maxscatlay),utot1(maxscatlay)
        REAL	fcover(maxscatlay)
        REAL    umif(maxmu,maxscatlay,maxf)
        REAL    uplf(maxmu,maxscatlay,maxf)	
        DOUBLE PRECISION mu1(maxmu), wt1(maxmu), galb, galb1, 
     1		pplsto(0:40,maxmu,maxmu,maxscatlay),
     2		pmisto(0:40,maxmu,maxmu,maxscatlay),
     3          epsS(maxlay),eps(maxscatlay)
	LOGICAL	scatter, single

C		Internal variables

	INTEGER	I, J, K, L, Ipath, Ig, nlays
	REAL	utotl(maxlay), qh(maxlay), qhe(maxlay),
     1		frac(maxlay,maxgas), qh_he(maxlay), dist1,
     2		totamh(maxlay), x, taucon(maxlay)
        REAL    taugas(maxlay), tauscat(maxlay), p, t
        REAL    taugasc(maxlay),xp
        REAL    tau, tau2, asec(maxsec), bsec(maxsec),
     1          tausc(maxcon), taus(maxlay),
     5		tautmp(maxlay), corkout(maxpat,maxg), 
     6		planck_wave, tmp, error, asec2(maxsec), 
     3          bsec2(maxsec), xsec2(2,maxcon,maxsec), muemiss,
     3          muinc, ssfac, f(maxlay), fint(maxlay), tmp1, 
     4          intscat, pastint, phase(maxcon), fup(maxlay,maxg),
     5          tmp2, taug(maxlay)
        REAL    taucloud(maxcon,maxlay),tauclscat(maxcon,maxlay)
        REAL    taucl(maxcon,maxlay),taucs(maxcon,maxlay)
        integer icloud(maxcon,maxlay)
        REAL	fdown(maxlay,maxg),fwhmk,
     1          delvk,tauray(maxlay),f1(maxlay),g11(maxlay),
     2          g21(maxlay),taur(maxlay) 
	DOUBLE PRECISION	tr, trold, taud, dtmp1, dtmp2, dpi, 
     1          dphase, calpha, dmuinc, dmuemiss, draddeg, dinc, demiss

C		Misc variables or continuum variables

	INTEGER	ii, id1,j1, IABSORB(5)
	REAL	DelH(nlayer), AAmount(maxgas), PPP(maxgas),
     1		XLen, AvgCONTMP, CONTMP(IORDP1), DABSORB(7)

        INTEGER LUNIS,IRECL,IOFF,NLAYERF,NMUF,NWAVEF,NGF,IFLUX,NFF

C		Common blocks and parameters

	common/dust/vsec,xsec,nsec,ncont
	common/interpk/lun, ireck, xmink, delk, pk, npk, tk, ntk, ng, 
     1		delvk, fwhmk, g_ord, delg, kout
	common/scatd/mu1, wt1, galb
	common/scatter1/nmu, isol, dist1, lowbc, liscat, lnorm,
     1		lncons, lcons, sol_ang, emiss_ang, aphi, nf
	common/phasesto/pplsto,pmisto

	PARAMETER	(tsun=5900.,theta0=4.65241e-3, pi=3.141593,
     1			error=5.e-5,LUNIS=61) 

C-----------------------------------------------------------------------
C
C	Check input parameters for possible problems.
C
C-----------------------------------------------------------------------


C	PRINT*,'CIRSRAD_WAVES calling input parameters'
C	PRINT*,'Dist = ',Dist
C	PRINT*,'INormal1 = ',Inormal
C	PRINT*,'nlayer,npath,ngas',nlayer,npath,ngas
C	PRINT*,'limlay,limcont,ncont',limlay,limcont,ncont
C	PRINT*,'press,temp,delh,cont'
C        DO i=1,nlayer 
C         print*,i,press(i),temp(i),delh(i),(cont(j,i),j=1,ncont)
C        ENDDO 
C	PRINT*,'pp, amount'
C        DO i=1,nlayer
C         print*,(pp(i,j),j=1,ngas)
C         print*,(amount(i,j),j=1,ngas)
C        ENDDO
C	PRINT*,'nwave',nwave,(vwave(i),i=1,nwave)
C	PRINT*,'nlayin',(nlayin(i),i=1,npath)
C	PRINT*,'imod',(imod(i),i=1,npath)
C	PRINT*,'incdim',incdim
C	PRINT*,'layinc,scale,emtemp'
C        do i=1,npath
C         print*,'ipath = ',i
C         do j=1,nlayin(i)
C          print*,layinc(j,i),scale(j,i),emtemp(j,i)
C         enddo
C        enddo
C	PRINT*,'idgas',(idgas(i),i=1,ngas)
C	PRINT*,'isogas',(isogas(i),i=1,ngas)
C	PRINT*,'itype',itype
C	PRINT*,'tsurf',tsurf
C        print*,'flagh2p',flagh2p
C        print*,'hfp : ',(hfp(i),i=1,nlayer)
C        print*,'flagc : ',flagc
C        print*,'hfc : ',(hfc(i),i=1,nlayer)
C        print*,'ifc, basep, baseh'
C        do i=1,nlayer
C         print*,(ifc(j,i),j=1,ncont),basep(i),baseh(i)
C        enddo

	if (nwave.gt.maxbin) then
		write (*,*) ' CIRSRAD_WAVES: Too many bins'
		write (*,*) ' NWave = ',nwave,' Maxbin = ',maxbin
		stop
	endif

	if (nlayer.gt.maxlay) then
		write (*,*) ' CIRSRAD_WAVES: Too many layers'
		write (*,*) ' Nlayer = ',nlayer,' maxlay = ',maxlay
		stop
	endif

	if (npath.gt.maxpat) then
		write (*,*) ' CIRSRAD_WAVES: Too many paths'
		write (*,*) ' Npath = ',npath,' maxpat = ',maxpat
		stop
	endif

	if (ngas.gt.maxgas) then
		write (*,*) ' CIRSRAD_WAVES: Too many gases'
		write (*,*) ' Ngas = ',ngas,' Maxgas = ',maxgas
		stop
	endif

	if (ncont.gt.maxcon) then
		write (*,*) ' CIRSRAD_WAVES: Too many dust continua'
		write (*,*) ' Ncont = ',ncont,' maxcon = ',maxcon
		stop
	endif

	if (nsec.gt.maxsec) then
		write (*,*) ' CIRSRAD_WAVES: Too many dust continua pts'
		write (*,*) ' Nsec = ',nsec,' Maxsec = ',maxsec
		stop
	endif

	if (ng.gt.maxg) then
		write (*,*) ' CIRSRAD_WAVES: Too many g ordinates'
		write (*,*) ' Ng = ',ng,' Maxg = ',maxg
		stop
	endif

	if (nmu.gt.maxmu) then
		write (*,*) ' CIRSRAD_WAVES: Too many zenith angles'
		write (*,*) ' nmu = ',nmu,' maxnmu = ',maxmu
		stop
	endif

	if ((npk.gt.maxk).or.(ntk.gt.maxk)) then
		write (*,*) ' CIRSRAD_WAVES: Too many P/T points in K',
     1			' tables' 
		write (*,*) ' Npk = ',npk,' Maxk = ',maxk
		write (*,*) ' Ntk = ',ntk,' Maxk = ',maxk
		stop
	endif

	do Ipath = 1, npath
		if (nlayin(Ipath).gt.incdim) then
		   write (*,*) ' CIRSRAD_WAVES: Nlayin exceeds Incdim'
			write (*,*) ' Ipath = ', Ipath,' Nlayin = ',
     1				nlayin(Ipath), ' Incdim = ', incdim
			stop
		endif
		do I = 1, nlayin(Ipath)
			if (layinc(I,Ipath).gt.nlayer) then
			write (*,*) ' CIRSRAD_WAVES: Layinc exceeds',
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
C		imod = 21		Limb scattering
C
C-----------------------------------------------------------------------
cc	print*,'TEST: cirsrad_waveS: nwave,nmu=',nwave,nmu


	dpi = 4.0d0 * datan(1.0d0)
	draddeg = 2.0d0 * dpi/360.0d0

	scatter=.FALSE.
	single=.FALSE.
	do Ipath = 1, npath
		if (imod(Ipath).eq.15.or.imod(Ipath).eq.21) then
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
                        print*,'ssfac',muinc,muemiss,ssfac
			do I = 1, nlayer
				scale(I,1) = sngl(1./dmuemiss + 
     1					1./dmuinc)
			enddo
			
			calpha = sin(sol_ang*pi/180.0)*
     1				 sin(emiss_ang*pi/180.0)*
     2				 cos((aphi-pi)*pi/180.0) - 
     3				 muemiss*muinc

			print*,'calpha',calpha

		endif
	end do

c	Isec = min(4,nsec)

	ntab = npk * ntk * ng


C-----------------------------------------------------------------------
C
C	Precompute volume fractions. The factor of 1e-20 applied to
C	utotl arises for the following reason. In order to avoid 
C	underflow, RADTRAN (and CIRSRAD) routines multiply line strengths
C	by a factor of 1.e47. Thus, when final calculation routines get 
C	hold of the strengths, they must correct by 1.e-47. This is broken
C	into two factors, 1.e-27 applied directly in calculating the K
C	tables, and 1.e-20 applied to the amounts where they are used in
C	conjunction with the K tables.
C
C	Totamh is the mass (number?) of molecules per cm^2 for H + He
C	and is used later in calculating the CIA opacity of H2. Because
C	it is not used with K tables, no 1.e-20 factor is needed.
C
C-----------------------------------------------------------------------

	DO I= 1, nlayer
		utotl(I)= 0.
		DO J= 1, ngas
			frac(I,J) = pp(I,J) / press(I)
			utotl(I) = utotl(I) + amount(I,J)
C                        print*,i,j,amount(i,j),utotl(i)
		ENDDO
		utotl(I) = utotl(I) * 1.e-20
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
                     
		 	CALL spline (vsec,asec,nsec,1.e30,1.e30,asec2)
		 	CALL spline (vsec,bsec,nsec,1.e30,1.e30,bsec2)

		 	DO L = 1, nsec
				xsec2(1,K,L) = asec2(L)
				xsec2(2,K,L) = bsec2(L)
		 	ENDDO
		ENDDO
	ENDIF

        XNEXT=0.0
        IRECL = 4

        IFLUX=2

        print*,'IFLUX = ',IFLUX

        IF(IFLUX.EQ.2)THEN
         OPEN(LUNIS,FILE='internal.fld',STATUS='OLD',
     1    ACCESS='DIRECT',RECL=IRECL)

         READ(LUNIS,REC=1)nlayerF
         READ(LUNIS,REC=2)nmuF
         READ(LUNIS,REC=3)nwaveF
         READ(LUNIS,REC=4)ngF
         READ(LUNIS,REC=5)nfF

         if(ngF.ne.ng)then
          print*,'Error in cirsrad_waveS: ngf <> ng'
          print*,ngf,ng
          stop
         endif

         DO I=1,nwaveF
          READ(LUNIS,REC=5+I)vwaveF(i)
         ENDDO

         DO I=1,nlayerF
          READ(LUNIS,REC=5+NWAVEF+I)basepF(i)
         ENDDO

         DO I=1,nlayerF
          READ(LUNIS,REC=5+NWAVEF+NLAYERF+I)basehF(i)
         ENDDO

         IOFF = 5+NWAVEF+2*NLAYERF+1

        ENDIF

        print*,'IOFF = ',IOFF

        DO I = 1, nwave
		x = vwave(I)
                esurf = interpem(nem,vem,emissivity,x)
C               Set vv to the current WAVENUMBER
                vv = x
                if(ispace.eq.1)then
                  vv=1e4/x 
                endif
		IF ((scatter).OR.(single)) THEN
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
                        phase(ncont+1)=0.75*sngl((1.+calpha*calpha))
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
                CALL GET_K(NLAYER,PRESS,TEMP,NGAS,I,X)
 
		DO J = 1, nlayer
			taucon(J) = 0.
			taugasc(J) = 0.
			tauscat(J) = 0.
			f(J) = 0.
			p = press(J)
			t = temp(J)

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

				CALL splint (vsec, asec, asec2, 
     1					nsec, x, tau)
				CALL splint (vsec, bsec, bsec2,
     1					nsec, x, tau2)

				if(tau.lt.0)then
             			 print*,'tau lt 0: Particle type ',K
                                 print*,nsec,vsec(1),vsec(nsec)
                                 print*,nsec,asec(1),asec(nsec)
                                 print*,x,tau
                                 print*,'Do linear interpolation'
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
                                 print*,'new tau ',tau
				 if(jf.gt.0)then
                                  print*,jf,vsec(jf),vsec(jf+1),xf
                                  print*,asec(jf),asec(jf+1)
                                 endif
                                endif


				if(tau2.lt.0)then
             			 print*,'tau2 lt 0: Particle type ',K
                                 print*,nsec,vsec(1),vsec(nsec)
                                 print*,nsec,bsec(1),bsec(nsec)
                                 print*,x,tau2
                                 print*,'Do linear interpolation'
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
                                 print*,'new tau2 ',tau2
				 if(jf.gt.0)then
				  print*,jf,vsec(jf),vsec(jf+1),xf
				  print*,bsec(jf),bsec(jf+1)
				 endif
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
C                           print*,'J,phase(K),F(J)',J,phase(K),f(J)

			 ENDIF
			ENDDO
			intscat = intscat + tauscat(J)
c			if (intscat.gt.pastint) then
c				write (*,*) ' Layer ', J,':',
c     1					' Integrated scattering: ',
c     2					intscat
c			endif
			pastint = intscat

			DO K = 1, ncont
				IF (tauscat(J).eq.0.) THEN
					frcscat(K,J) = 1./ncont
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
C	GASCON is a RADTRANS routine which
C	calculates continuum contributions from various gases
C
C-----------------------------------------------------------------------

cc		        WRITE (*,*) '       CALLING NGascon for each gas'
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
C	a different subroutine to allow easy modification. GASCON returns
C	a polynomial approximation to the continuum of a particular gas
C	over a particular bin.
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
	             	IF(FLAGH2P.EQ.1)THEN
                	   FPARA = HFP(J)
 			   CALL NPARACON(vv,P,T,
     1			     NGas,idgas,isogas,AAmount,PPP,
     2			     FPARA,XLen,AvgCONTMP,IABSORB,DABSORB,
     3			     id1)
             		ELSE
			   CALL NCIACON(vv,P,T,INormal,
     1			     NGas,idgas,isogas,AAmount,PPP,
     2			     XLen,AvgCONTMP,IABSORB,DABSORB,id1)
             		ENDIF

			taucon(J)= taucon(J) + AvgCONTMP
			taugasc(J)= taugasc(J) + AvgCONTMP

C	The code below is if Rayleigh scattering is considered.

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
                else
                 xray = RAYLEIGHA(vv,p,t)*1E20
                endif
                avgcontmp = utotl(j)*xray
                tauray(J)=avgcontmp

                if(AvgCONTMP.ge.0.0)then
 		 taucon(J)= taucon(J) + AvgCONTMP
                 taugasc(J)= taugasc(J) + AvgCONTMP
                 tauscat(J)=tauscat(J) + AvgCONTMP
	        endif

C               Calculate single-scattering contribution
                f(J) = f(J) + phase(ncont+1)*tauray(J)
C                print*,'single',J,f(J),phase(ncont+1)
                           
        ELSE
               	tauray(J)=0.0
	ENDIF

C       ### renormalise f(J) to be tau-averaged phase function
        if(tauray(J)+tauscat(J).gt.0)then
           f(J) = f(J)/(tauray(J)+tauscat(J))
        else
           f(J) = 0.0
        endif
C        print*,'Final f(j) ',f(j)

C-----------------------------------------------------------------------
C
C	For each gas, look up the K coefficients. Stepping through the 
C	gases, combine each set of coeffs with the combination of all
C	previous gases using OVERLAP. Finally, store the layer dependent
C	summed K coefficients. End the loop over each layer.
C
C-----------------------------------------------------------------------

			DO K = 1, ngas
			 IF ((idgas(K).EQ.39).OR.(idgas(K).EQ.40)
     1				.OR.
     2				(lun(I,K).LT.0)) THEN
				DO L = 1, ng
					k_g2(L) = 0.
				ENDDO
			 ELSE
				DO L = 1,NG
                          		k_g2(L)=KOUT(J,K,L)
				ENDDO
			 ENDIF

		 	 IF (K.EQ.1) THEN
				q1 = frac(J,K)
				DO L = 1, ng
					k_g(L) = k_g2(L)
				ENDDO
			 ELSE
				q2 = frac(J,K)
				DO L = 1, ng
					k_g1(L) = k_g(L)
				ENDDO

C	This subroutine combines the absorption coefficient distributions
C	of two overlapping gases. The overlap is implicitly assumed to be 
C	random and the k-distributions are assumed to have NG-1 mean
C	values and NG-1 weights. Correspondingly there are NG ordinates in
C	total.
				CALL overlap(delg,ng,k_g1,q1,
     1					k_g2,q2,k_g)
				q1 = q1 + q2
			 ENDIF
			ENDDO

			DO K = 1, ng
				kl_g(K,J) = k_g(K)
			ENDDO
C			print*,'Layer : ',J
C			print*,'K = ',(kl_g(K,J),K=1,NG)
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

		DO Ig = 1, ng
			DO J = 1, nlayer
				tautmp(J) = taucon(J) + kl_g(Ig,J) 
     1					* utotl(J)
				taugas(J) = taugasc(J) + kl_g(Ig,J) 
     1					* utotl(J)

C                                print*,J,Ig,kl_g(Ig,J)*utotl(J),
C     1                             taucon(J),utotl(J)


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
C				previous output paths. NOT SUPPORTED HERE.
C		11	(Atm) Contribution function.
C		15	(Atm) Multiple scattering (multiple models)
C		16	(Atm) Single scattering approximation.
C		20	(Atm) Net flux calculation
C	then end the loop over the g ordinate.
C
C-----------------------------------------------------------------------

                modsource=.true.


                if(modsource)then
                   call impflux(LUNIS,ioff,nlayerf,nmuf,nff,umif,uplf,
     1                x,nwavef,vwavef,Ig,ng)
C                   print*,'Fluxes imported'
                endif


			DO Ipath = 1, npath
				nlays = nlayin(Ipath)
				corkout(Ipath,Ig) = 0.
				DO J = 1, nlays

		taus(J) = tautmp(layinc(J,ipath)) * scale(J,ipath)
		taur(J) = tauray(layinc(J,ipath)) * scale(J,ipath)
                if(Ig.eq.1)then
                    bb(J,Ipath)=planck_wave(ispace,x,emtemp(J,Ipath))
                    jb(J,Ipath)=bb(J,Ipath)
                    bnuS(J) = bb(J,Ipath)
                endif

                if(modsource)then

C                        print*,J,taus(J),tautmp(layinc(J,Ipath)),
C     1                    tauscat(layinc(J,Ipath))
                        IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
                                dtmp1 = dble(tauscat(layinc(J,Ipath)))
                                dtmp2 = dble(tautmp(layinc(J,Ipath)))
C                                print*,J,dtmp1,dtmp2
                                 epsS(J) = 1.0d00 - dtmp1/dtmp2
C                                print*,layinc(J,IPath)
C                                print*,dtmp1/dtmp2
C                                print*,epsS(J)
                        ELSE
                                epsS(J)=1.0
                        END IF

                        basehS(J) = baseh(layinc(J,Ipath))
                        scaleS(J) = scale(J,Ipath)

                endif



C                print*,'taus',J,taus(J)

C               New arrays for scloud12wave. taucl holds the extinction
C		optical depth of each cloud type. taucs holds the 
C  		scattering optical depth of each cloud type. taug holds
C		the total gas extinction optical depth (excluding 
C		Rayleigh scattering).
                do k=1,ncont
                 taucl(k,j)=taucloud(k,layinc(J,ipath))*scale(J,ipath)
                 taucs(k,j)=tauclscat(k,layinc(J,ipath))*scale(J,ipath)
                enddo
 
C                print*,'taucl',(taucl(k,j),k=1,ncont)
C                print*,'taucs',(taucs(k,j),k=1,ncont)

		taug(J) = taugas(layinc(J,ipath)) * scale(J,ipath)
C               Subtract Rayleigh opdepth from taug
                taug(J) = taug(J)-taur(J)

C                print*,'taug',taug(J)

                if(single)then
C                  print*,'single'
                  fint(j)=f(layinc(J,ipath))
C                  print*,j,ipath,layinc(j,ipath),fint(j),scale(j,ipath)
C                  print*,taus(J),press(layinc(J,ipath))
		endif

		DO K = 1, ncont
			lfrac(K,J) = frcscat(K,layinc(J,ipath)) 
		ENDDO
				ENDDO


C        stop


        if(modsource)then


          call hgaverage(nlays,ncont,taus,epsS,lfrac,tauray,x,
     1     f1,g11,g21)

          call scatsource(nlays,nmuf,wt1,mu1,epsS,bnuS,taus,basehS,
     1      scaleS,f1,g11,g21,nlayerf,basehf,umif,uplf,nff,Jsource)

          do j=1,nlays
           jb(J,Ipath)=Jsource(J)
           print*,'scat',ipath,j,bb(j,ipath),jb(j,ipath)
          enddo

        endif


C       ***************************************************************	
	IF (imod(ipath).EQ.0) THEN

cc		WRITE(*,*) 'CIRSRAD_WAVES: Imod= 0 =Pure Transmission,',
cc     1			' creating output'

		DO J = 1, nlays
			corkout(Ipath,Ig) = corkout(Ipath,Ig) + taus(J)
		ENDDO
		corkout(Ipath,Ig) = exp(-corkout(Ipath,Ig))

	ELSEIF (imod(ipath).EQ.1) THEN	

cc         WRITE(*,*) 'CIRSRAD_WAVES: Imod= 1 =Absorption, creating',
cc     1                  ' output'
 
		DO J = 1, nlays
			corkout(Ipath,Ig) = corkout(Ipath,Ig) + taus(J)
		enddo
		corkout(Ipath,Ig) = 1.- exp(-corkout(Ipath,Ig))

	ELSEIF (imod(ipath).EQ.3) THEN	

C		WRITE(*,*) 'CIRSRAD_WAVES: Imod= 3 =Emission, creating',
C     1                  ' output'

		taud = 0.
		trold = 1.
		do J = 1, nlays
			taud = taud + taus(J)
			tr = dexp(-taud)
		        corkout(Ipath,Ig) = corkout(Ipath,Ig) + 
     1                  sngl((trold-tr)) * jb(J,Ipath)
 			trold = tr
		enddo

	ELSEIF (imod(ipath).EQ.11) THEN	

cc                WRITE(*,*) 'CIRSRAD_WAVES: Imod= 11 =Contribution Function,',
cc     1                  ' creating output'

		trold = 1.
		do J = 1, nlays
			taud = taus(J)
			tr = dexp(-taud)
			corkout(Ipath,Ig) = corkout(Ipath,Ig) + 
     1                  sngl((trold-tr)) * jb(J,Ipath)
 			trold = tr
                enddo

	ELSEIF (imod(ipath).EQ.15) THEN

cc                WRITE(*,*) 'CIRSRAD_WAVES: Imod= 15 =Multiple Scattering, ',
cc     1                  ' creating output'

		do J = 1, nlays
			bnu(J) = jb(J,Ipath)
			IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
				dtmp1 = dble(tauscat(layinc(J,Ipath)))
				dtmp2 = dble(tautmp(layinc(J,Ipath)))
				eps(J) = 1.0d00 - dtmp1/dtmp2
		        ELSE
  				EPS(J)=1.0
         		END IF
                        temp1(J) = TEMP(layinc(J,Ipath))
                        press1(J) = PRESS(layinc(J,Ipath))
			utot1(J) = UTOTL(layinc(J,Ipath))
			fcover(J) = HFC(layinc(J,Ipath))
                        do k=1,ncont
                         icloud(k,j)=IFC(k,layinc(J,Ipath))
                        enddo
C                        print*,J,layinc(J,Ipath),bnu(J),eps(j),
C     1     taus(j),temp1(j),press1(j),utot1(j),lfrac(1,j),liscat(1)
		enddo


               
C               If tsurf > 0, then code assumes bottom layer is at the
C               surface and so uses tsurf to compute radiation upwelling.
C               Otherwise, code assumes optical depth is large and sets
C               upwelling radiation field to local temperature.

                if(tsurf.le.0.0)then
                 radground = jb(nlays,Ipath)
                else
                 radground = esurf*planck_wave(ispace,x,tsurf)
                endif

                galb1=galb

                if(galb1.lt.0)then
                         galb1 = get_albedo(nem,vem,emissivity,x)
                endif

                do J=1,nmu
                 radg(J)=radground*sngl(1.0-galb1)
                enddo

                IF (itype.EQ.11) THEN

C                        WRITE (*,*) '     CALLING scloud11wave'
C                        WRITE (*,*)'IRAY,INORMAL = ',IRAY,INORMAL
C                        print*,sol_ang,solar
C                        print*,emiss_ang,aphi
C                        print*,radg
C                        print*,lowbc,galb1

		  	call scloud11wave(rad1, sol_ang, emiss_ang,
     1                          aphi, radg, solar, lowbc, galb1, iray,
     2				mu1, wt1, nmu,   
     3				nf, Ig, x, vv, eps, bnu, taus, taur,
     4                          nlays,ncont,lfrac)

 		  	corkout(Ipath,Ig) = rad1

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

 		  	corkout(Ipath,Ig) = rad1

                ELSEIF (itype.EQ.13) THEN

                   if(nf.ne.0)then
                    print*,'Error in cirsrad_waveS, nf <> 0',nf
                   endif

   	  	   call scloud11flux(radg, solar, sol_ang, 
     1               lowbc, galb1, iray, mu1, wt1, nmu, nf, Ig, x,
     2               vv, eps, bnu, taus, taur, 
     3               nlays, ncont, lfrac, umif, uplf)
                    
                     open(48,file='test.dat',status='unknown')
                     write(48,*)nlays
                     do k=1,nlays
                      write(48,*)k,press1(k),(umif(j,k,1),j=1,nmu)
                     enddo
                     do k=1,nlays
                      write(48,*)k,press1(k),(uplf(j,k,1),j=1,nmu)
                     enddo
                     close(48)

                     call dumpflux(LUNIS,ioff,nlays,nmu,nf,radg,umif,
     1                uplf,I,x,nwave,Ig,ng)

 		     corkout(Ipath,Ig) = umif(1,1,1)

                ELSE

			WRITE (*,*) ' Scattering ITYPE =',itype,
     1				' not defined for CIRSRAD_WAVES.'
			STOP
                ENDIF

	ELSEIF (imod(ipath).EQ.16) THEN

cc		WRITE(*,*) 'CIRSRAD_WAVES: Imod= 16 =Single Scat Approx,',
cc     1                  ' creating output'

C		IF (Ig.EQ.1) WRITE(*,*)' SINGLE SCATTERING IN USE '

                galb1 = galb
                        
                if(galb1.lt.0)then
                   galb1 = get_albedo(nem,vem,emissivity,x)
C                  if(Ig.eq.1)print*,x,galb1
                endif

                taud = 0.
                trold = 1.
		DO J= 1, nlays

  		        IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
			 tmp1 = tauscat(layinc(J,Ipath))+
     1				  tauray(layinc(J,Ipath))

 
			 tmp2 = tautmp(layinc(J,Ipath))
			 omega(J) = tmp1/tmp2
		        ELSE
  			 omega(J)=0.0
         	        END IF

                        taud = taud + taus(J)
                        tr = dexp(-taud)
			corkout(Ipath,Ig) = corkout(Ipath,Ig) + 
     1				sngl(trold-tr)*ssfac*fint(J)*
     2				omega(J)*solar/(4*pi)

                        trold = tr

		ENDDO
C               Add in reflectance from the ground
                corkout(Ipath,Ig)=corkout(Ipath,Ig) + 
     1                          sngl(trold*solar*galb1)/pi

C		WRITE (*,*) ' Calculated: ', Ig, corkout(Ipath,Ig)

	ELSEIF (imod(ipath).EQ.20)THEN

cc		WRITE(*,*) 'CIRSRAD_WAVES: Imod= 20 =Net flux calculation,',
cc     1                  ' creating output'

C		Computes up and down flux at BOTTOM of each layer
C		Compute once for Ipath=1

		if(Ipath.eq.1)then
                  do j=1,nlayin(npath)
                   fup(Ig,j)=0.0
                   fdown(Ig,j)=0.0
                   taus(layinc(J,npath))=tautmp(layinc(J,npath)) *
     &			scale(J,npath)
		  enddo

	   	  do j=1,nlayin(npath)
           	   if(j.eq.1)then
                    j1 = layinc(j,npath)
		    fup(j,Ig)=bb(j1,1)
		   else
		    tr = exp(-taus(j-1))
               fup(j,Ig)=fup(j-1,Ig)*sngl(tr) + sngl(1.0-tr)*bb(j1,1)
C		    print*,j,bb(j1,1),tr,fup(j,Ig)	
		   endif
		  enddo 

	   	  do j=nlayin(npath),1,-1
                   tr=exp(-taus(j))
            	   if(j.eq.nlays)then
		    fdown(j,Ig)=sngl(1.0-tr)*bb(j1,1)
		   else
             fdown(j,Ig)=sngl(1.0-tr)*bb(j1,1)+sngl(tr)*fdown(j+1,Ig)
C		    print*,j,bb(j1,1),tr,fdown(j,Ig)	
		   endif
		  enddo 

                  if(Ig.eq.1)then
                   do j=1,nlayin(npath)
C		    print*,j,fup(j,Ig),fdown(j,Ig),
C     &			fup(j,Ig)-fdown(j,Ig)			
		   enddo
                  endif
                endif 
            
	        corkout(Ipath,Ig) = fup(Ipath,Ig) - fdown(Ipath,Ig) 	

	ELSE
		WRITE (*,*) ' Imod = ', imod(ipath), ' is not a valid',
     1			' CIRSRAD_WAVES option'
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

		XCOM = 100.0*FLOAT(I)/FLOAT(NWAVE)
		IF(XCOM.GE.XNEXT)THEN
			WRITE(*,1000)
			WRITE(*,1010) I, nwave, x
			WRITE(*,1020)XCOM

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
1020	FORMAT(' CIRSRAD_WAVES.f :: % complete : ',F5.1)

        IF(IFLUX.GT.0)THEN
         CLOSE(LUNIS)
        ENDIF

C        print*,'cirsrad_waveS'
C        do ipath=1,npath
C         ioff = (ipath-1)*nwave
C         print*,(output(ipath,i),i=1,nwave)
C        enddo

	RETURN

	END

************************************************************************
************************************************************************
