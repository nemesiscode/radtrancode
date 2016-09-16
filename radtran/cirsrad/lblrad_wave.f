************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C_TITLE:		SUBROUTINE LBLrad_wave
C
C_DESCR:
C
C_ARGS:	Input Variables
C
C	AvgCONTMP:REAL	The average of the three points returned by CIACON
C	INORMAL:INT     flag for ortho:para ratio (0=equilibrium =1:1)
C			(1=normal   =3:1)
C       DELH(NLAYER):REAL Physical thickness of each layer (km).
C       ispace           integer Indicates if wavelengths in vconv
C                               are in wavenumbers(0) or
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
C	IORDER		The order of the continuum polynomial= 3
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
C_HIST:	
C-----------------------------------------------------------------------

      SUBROUTINE LBLrad_wave (X, WING, VMIN, VMAX, VREL, MAXDV, 
     1    IBS, IBD, Dist, INormal, IRay, ispace, DelH, nlayer,
     2    npath, ngas, limlay, limcont, totam, press, temp, pp, 
     3    amount, nlayin, incdim, layinc, cont, scale, imod,
     4    idgas, isogas, iproc, emtemp, itype, nem, vem, emissivity, 
     5    tsurf, 
     6    flagh2p,hfp, flagc, hfc, ifc, basep, baseh, RADIUS1, radextra, 
     7    output)

      IMPLICIT NONE


C		Internal dimensions


C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/contdef.f'

      INCLUDE '../includes/lincomc.f'
C		Passed variables

      INTEGER	Inormal,ik,iray,flagh2p,flagc,ispace,jf
      INTEGER	nlayer, npath, ngas, limlay, limcont, 
     1		nlayin(npath), incdim, layinc(incdim, npath), 
     2          idgas(ngas), isogas(ngas), imod(npath),
     3		iproc(ngas)
      INTEGER	MAXLIN1

      REAL	press(nlayer), temp(nlayer), pp(limlay,ngas), 
     1		amount(limlay,ngas), cont(limcont,nlayer), 
     2		scale(incdim,npath), emtemp(incdim,npath), hfp(nlayer), 
     3		bb(maxlay,maxpat),xf,tsurf,dv,bsurf,
     4          hfc(nlayer),basep(nlayer),baseh(nlayer),
     5          esurf,radground, totam(nlayer),p1,p2,RADIUS1,radextra
      REAL    vem(maxsec),emissivity(maxsec),interpem
      REAL	xmu,dtr,xt1,xfac
      INTEGER ifc(limcont,nlayer),nem,j0,IBS(2),IBD(2)

	REAL	XCOM,XNEXT,FPARA,vv,XRAY,RAYLEIGHJ,RAYLEIGHA,RAYLEIGHV
        REAL    WING,MAXDV,VMIN,VMAX

        double precision get_albedo

C		Dust variables

	INTEGER	nsec, ncont, IFCONT, IB, FSTREC
	REAL	vsec(maxsec), xsec(2,maxcon,maxsec), xk

C		K table variables

        REAL    basehf(maxlay),basepf(maxlay),basehS(maxlay)
        REAL    scaleS(maxlay),Jsource(maxlay),delhs(maxlay)
C		Scattering variables

	INTEGER	liscat(maxcon), isol, nmu, lowbc, 
     1          lncons(maxcon), ibaseH(maxlay),
     2		kf,nf, lnorm(maxcon), itype
	REAL	lfrac(maxcon,maxlay), solar, pi, theta0, dist, 
     1		bnu(maxscatlay), aphi, radg(maxmu), solar1, 
     2		lcons(maxcon,maxscatpar), rad1, omega(maxscatlay),
     3		tsun, v, frcscat(maxcon,maxlay), sol_ang, emiss_ang
        REAL    bnuS(maxlay),fcover(maxscatlay)
        REAL    umif(maxmu,maxscatlay,maxf)
        REAL    uplf(maxmu,maxscatlay,maxf)	
        REAL    umift(maxmu,maxmu,maxscatlay,maxf)
        REAL    uplft(maxmu,maxmu,maxscatlay,maxf)	
        DOUBLE PRECISION mu1(maxmu), wt1(maxmu), galb, galb1, 
     1		pplsto(0:maxscatlay,maxmu,maxmu,maxscatlay),
     2		pmisto(0:maxscatlay,maxmu,maxmu,maxscatlay),
     3          eps(maxscatlay),epsS(maxlay),
     4          omegas(maxscatlay),omegasS(maxlay)
	LOGICAL	scatter, single

C		Internal variables

	INTEGER	I, J, K, L, Ipath, Ig, nlays
	REAL	qh(maxlay), qhe(maxlay),
     1		frac(maxlay,maxgas), qh_he(maxlay), dist1,
     2		x, taucon(maxlay)
        REAL    taugas(maxlay), tauscat(maxlay), p, t
        REAL    taugasc(maxlay),xp,VREL
        REAL    tau, tau2, asec(maxsec), bsec(maxsec),
     1          tausc(maxcon), taus(maxlay),
     5		tautmp(maxlay), output(maxpat), 
     6		planck_wave, tmp, error, asec2(maxsec), 
     3          bsec2(maxsec), xsec2(2,maxcon,maxsec), muemiss,
     3          muinc, f(maxlay), fint(maxlay), tmp1, 
     4          intscat, pastint, phase(maxcon), fup(maxlay),
     5          tmp2, taug(maxlay),ssfac,xsolar,xdist
        REAL    taucloud(maxcon,maxlay),tauclscat(maxcon,maxlay)
        REAL    taucl(maxcon,maxlay),taucs(maxcon,maxlay)
        integer icloud(maxcon,maxlay)
        REAL	fdown(maxlay),fwhmk,
     1          delvk,tauray(maxlay),f1(maxlay),g11(maxlay),
     2          g21(maxlay),taur(maxlay) 
	DOUBLE PRECISION	tr, trold, taud, dtmp1, dtmp2, dpi, 
     1          dphase, calpha, dmuinc, dmuemiss, draddeg, dinc, demiss,
     2          tausun

C		Misc variables or continuum variables

	INTEGER	ii, id1,j1,LAYER,NLINR,IABSORB(5)
	REAL	DelH(nlayer), AAmount(maxgas), PPP(maxgas),
     1		XLen, AvgCONTMP, VBOT,DABSORB(7)

        INTEGER LUNIS,IRECL,IOFF,NLAYERF,NMUF,NWAVEF,NGF,NFF
        integer fintrad
        character*100 fintname
C		Common blocks and parameters


C       Solar reference spectrum common block
        real swave(maxbin),srad(maxbin),solrad
        integer iread,nspt,iform
        common /solardat/iread,iform,solrad,swave,srad,nspt

	common/dust/vsec,xsec,nsec,ncont
	common/scatd/mu1, wt1, galb
	common/scatter1/nmu, isol, dist1, lowbc, liscat, lnorm,
     1		lncons, lcons, sol_ang, emiss_ang, aphi, nf
	common/phasesto/pplsto,pmisto
        common/intrad/fintrad,fintname


	PARAMETER	(tsun=5900.,theta0=4.65241e-3, pi=3.141593,
     1			error=5.e-5,LUNIS=61) 


C       IBS(2) specifies which of the two buffers is buffer 1 and buffer 2
C              at start IBS = [1,2], but when a new set of lines is read in
C	       it becomes [2,1] and then keeps switching as lines are read.
C	IBD(2) specifies if the lines in a buffer need to be added to the
C	       continuum. -1 indicates that the lines are new and need adding
C	       and +1 indicates that the lines have already been added to the
C	       continuum.
C-----------------------------------------------------------------------
C
C	Check input parameters for possible problems.
C
C-----------------------------------------------------------------------

C        print*,X,WING,VREL,MAXDV,IBS(1),IBS(2),IBD(1),IBD(2)
C        print*,Dist, INormal1, ispace,nlayer

C        print*,'npath, ngas, ncont'
C        print*,npath, ngas, ncont
C        print*,'incdim,itype',incdim,itype
C        print*,'tsurf,flagh2p,flagc'
C        print*,tsurf,flagh2p,flagc
C        print*,'imod',(imod(i),i=1,npath)
C        do i=1,ngas
C         print*,i,idgas(i),isogas(i),iproc(i)
C        enddo


C        do i=1,nlayer
C         print*,DelH(i),press(i),temp(i),basep(i),baseh(i)
C         print*,(pp(i,j),j=1,ngas)
C         print*,(amount(i,j),j=1,ngas)
C         print*,(cont(j,i),j=1,ncont)
C         print*,(ifc(j,i),j=1,ngas)
C         print*,hfp(i),hfc(i)
C        enddo

C        print*,'path'
C        do i=1,npath     
C         print*,nlayin(i)
C         print*,imod(i)
C         print*,(layinc(j,i),j=1,nlayin(i))
C         print*,(scale(j,i),j=1,nlayin(i))
C         print*,(emtemp(j,i),j=1,nlayin(i))
C        enddo
        
C        print*,'emissivity',nem
C        do i=1,nem
C         print*,vem(i),emissivity(i)
C        enddo



        
	if (nlayer.gt.maxlay) then
		write (*,*) ' LBLrad_wave: Too many layers'
		write (*,*) ' Nlayer = ',nlayer,' maxlay = ',maxlay
		stop
	endif

	if (npath.gt.maxpat) then
		write (*,*) ' LBLrad_wave: Too many paths'
		write (*,*) ' Npath = ',npath,' maxpat = ',maxpat
		stop
	endif

	if (ngas.gt.maxgas) then
		write (*,*) ' LBLrad_wave: Too many gases'
		write (*,*) ' Ngas = ',ngas,' Maxgas = ',maxgas
		stop
	endif

	if (ncont.gt.maxcon) then
		write (*,*) ' LBLrad_wave: Too many dust continua'
		write (*,*) ' Ncont = ',ncont,' maxcon = ',maxcon
		stop
	endif

	if (nsec.gt.maxsec) then
		write (*,*) ' LBLrad_wave: Too many dust continua pts'
		write (*,*) ' Nsec = ',nsec,' Maxsec = ',maxsec
		stop
	endif


	do Ipath = 1, npath
		if (nlayin(Ipath).gt.incdim) then
		   write (*,*) ' LBLrad_wave: Nlayin exceeds Incdim'
			write (*,*) ' Ipath = ', Ipath,' Nlayin = ',
     1				nlayin(Ipath), ' Incdim = ', incdim
			stop
		endif
		do I = 1, nlayin(Ipath)
			if (layinc(I,Ipath).gt.nlayer) then
			write (*,*) ' LBLrad_wave: Layinc exceeds',
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


	dpi = 4.0d0 * datan(1.0d0)
	draddeg = 2.0d0 * dpi/360.0d0

	scatter=.FALSE.
	single=.FALSE.
	do Ipath = 1, npath
		if (imod(Ipath).eq.15.or.imod(Ipath).eq.21.
     1            or.imod(Ipath).eq.22) then
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
C                        print*,'ssfac',muinc,muemiss,ssfac
                        do I = 1, nlayer
                                scale(I,1) = sngl(1./dmuemiss +
     1                                  1./dmuinc)
                        enddo

			
			calpha = sin(sol_ang*pi/180.0)*
     1				 sin(emiss_ang*pi/180.0)*
     2				 cos((aphi-pi)*pi/180.0) - 
     3				 muemiss*muinc

C			print*,'calpha',calpha

		endif
	end do

c	Isec = min(4,nsec)


C-----------------------------------------------------------------------
C
C	Precompute volume fractions. 
C
C-----------------------------------------------------------------------

	DO I= 1, nlayer
		DO J= 1, ngas
			frac(I,J) = pp(I,J) / press(I)
		ENDDO
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

        esurf = interpem(nem,vem,emissivity,x)

C       Set vv to the current WAVENUMBER
        vv = x
        if(ispace.eq.1)then
          vv=1e4/x 
        endif

C        print*,'lblrad_wave ispace,vv=',ispace,vv

        solar=0.
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
C               Calculate Rayleigh scattering too
                phase(ncont+1)=0.75*sngl(1.+calpha*calpha)
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

        VBOT=VBIN(1)
 
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

			CALL csplint (vsec, asec, asec2, 
     1				nsec, x, tau)
			CALL csplint (vsec, bsec, bsec2,
     1				nsec, x, tau2)

			if(tau.lt.0)then
C             		 print*,'tau lt 0: Particle type ',K
C                         print*,nsec,vsec(1),vsec(nsec)
C                         print*,nsec,asec(1),asec(nsec)
C                         print*,x,tau
C                         print*,'Do linear interpolation'
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
C                         print*,'new tau ',tau
			  if(jf.gt.0)then
C                            print*,jf,vsec(jf),vsec(jf+1),xf
C                            print*,asec(jf),asec(jf+1)
                          endif
                        endif


			if(tau2.lt.0)then
C             		 print*,'tau2 lt 0: Particle type ',K
C                         print*,nsec,vsec(1),vsec(nsec)
C                         print*,nsec,bsec(1),bsec(nsec)
C                         print*,x,tau2
C                         print*,'Do linear interpolation'
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
C                         print*,'new tau2 ',tau2
		          if(jf.gt.0)then
C			  print*,jf,vsec(jf),vsec(jf+1),xf
C			  print*,bsec(jf),bsec(jf+1)
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
C                  ##### f(J) is average phase function ######
  		   f(J) = f(J)+phase(K)*tau2
C                  print*,'J,phase(K),F(J)',J,phase(K),f(J)
C                   print*,'X1',taucon(j)
	 	 ENDIF
		ENDDO
		intscat = intscat + tauscat(J)
c		if (intscat.gt.pastint) then
c			write (*,*) ' Layer ', J,':',
c     1				' Integrated scattering: ',
c     2				intscat
c		endif
		pastint = intscat

		DO K = 1, ncont
			IF (tauscat(J).eq.0.) THEN
				frcscat(K,J) = 0.0
			ELSE
				frcscat(K,J) = tausc(K)/
     1					tauscat(J)
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

cc	        WRITE (*,*) '       CALLING NGascon for each gas'
		DO K = 1, ngas

C	Computes a polynomial approximation to any known continuum spectra 
C       for a particular gas over a defined wavenumber region.

		AvgCONTMP=0.0
		 CALL ngascon(vv,idgas(K),isogas(K),
     1		  amount(J,K),pp(J,K),p,t,AvgCONTMP)

                 if(AvgCONTMP.ge.0)then
C		    Can occasionally get -ve AvCONTMPs at
C		    edge of Karkoschka CH4 data.
 		    taucon(J) = taucon(J) + AvgCONTMP
 		    taugasc(J) = taugasc(J) + AvgCONTMP
C                    print*,'X2',taucon(j)
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
320		CONTINUE

C	ID1= 0 means no diagnostic print statements while in CIACON.
		id1 = 0
		XLEN = DELH(J)

C       This goes off to do the collision induced-absorption

		AvgCONTMP=0.0
	        IF(FLAGH2P.EQ.1)THEN
                   FPARA = HFP(J)
 		   CALL NPARACON(vv,P,T,
     1		     NGas,idgas,isogas,AAmount,PPP,
     2		     FPARA,XLen,AvgCONTMP,IABSORB,DABSORB,id1)
             	ELSE
		   CALL NCIACON(vv,P,T,INormal,
     1		     NGas,idgas,isogas,AAmount,PPP,
     2		     XLen,AvgCONTMP,IABSORB,DABSORB,id1)
             	ENDIF

                if(AvgCONTMP.gt.0.0) then
 		 taucon(J)= taucon(J) + AvgCONTMP
		 taugasc(J)= taugasc(J) + AvgCONTMP
C		 print*,'X3',taucon(j)
		endif

C	The code below is if Rayleigh scattering is considered.

		IF(IRAY.GT.0.AND.ITYPE.NE.0)THEN
                  if(IRAY.EQ.1)THEN
                   xray = RAYLEIGHJ(vv,p,t)
 		  elseif(IRAY.EQ.2)then
                   xray = RAYLEIGHV(vv,p,t)
		  else
                   xray = RAYLEIGHA(vv,p,t)
		  endif
                  avgcontmp = totam(J)*xray
                  tauray(J)=avgcontmp
C                  print*,J,totam(j),xray,tauray(J)

                  if(AvgCONTMP.ge.0.0)then
 	   	   taucon(J)= taucon(J) + AvgCONTMP
                   taugasc(J)= taugasc(J) + AvgCONTMP
                   tauscat(J)=tauscat(J) + AvgCONTMP
C		   print*,'X4',taucon(j)
	          endif

C                 Calculate single-scattering contribution
                  f(J) = f(J) + phase(ncont+1)*tauray(J)
C                 print*,'single',J,f(J),phase(ncont+1)
                           
        	ELSE
               	  tauray(J)=0.0
		ENDIF       

C       ### renormalise f(J) to be tau-averaged phase function
        	if(tauray(J)+tauscat(J).gt.0)then
C           	 f(J) = f(J)/(tauray(J)+tauscat(J))
            	 f(J) = f(J)/tauscat(J)
        	else
          	 f(J) = 0.0
        	endif

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

C check bins in layer 1        
        J=1
        CALL CALC_FINE(VV,WING,MAXDV,J,NLAYER,NGAS,PRESS,TEMP,
     1       FRAC,IDGAS,ISOGAS,IPROC,IBS,IFCONT,XK)

C          print*,'calc_fine OK'

        IF(IFCONT.EQ.1.AND.J.EQ.1)THEN
C          If we are in the 1st layer and CALC_FINE indicates that we have run
C          out of lines then we need to relabel <buffer 2> as <buffer 1> and
C     read in a new <buffer 2>.
           IB=IBS(1)
           FSTREC=NXTREC
           print*,'Loading more lines into what was <buffer 1>, '
           print*,'but will now be <buffer 2>'
           print*,'FSTREC = ',FSTREC
           MAXLIN1 = MAXLIN
           print*,'lblrad_wave: MAXLIN = ',MAXLIN1
           print*,IBS(1),IBS(2)
           print*,IBD(1),IBD(2)
           print*,VMIN,VMAX,VREL
           
           CALL LOADBUFFER(VMIN-VREL,VMAX+VREL,FSTREC,MAXLIN1,MAXBIN,
     1          IB,NGAS,IDGAS,ISOGAS,VBOT,WING,NLINR,VLIN,SLIN,ALIN,
     3          ELIN,IDLIN,SBLIN,PSHIFT,DOUBV,TDW,TDWS,LLQ,NXTREC,
     3          FSTLIN,LSTLIN,LASTBIN)
           
           NLINE(IB)=NLINR
           IBD(IB)=-1
           
C     Swap round the IBS and IBD arrays
           J1=IBS(1)
           IBS(1)=IBS(2)
           IBS(2)=J1
           
           J1=IBD(1)
           IBD(1)=IBD(2)
           IBD(2)=J1
           
           
           print*,IBS(1),IBS(2)
           print*,IBD(1),IBD(2)
           

C           print*,'IB,NLINE(IB),IBD',IB,NLINE(IB),IBD(IB)
C           DO I=1,NBIN
C     PRINT*,I,FSTLIN(IBS(1),I),LSTLIN(IBS(1),I),
C     1        FSTLIN(IBS(2),I),LSTLIN(IBS(2),I)
C     ENDDO
           


        ENDIF

C now loop over layers
      	DO J = 1, nlayer
C	  First compute continuum contributions as necessary in
C         all continuum bins and all layers.
          DO I=1,2
           IF(IBD(I).LT.0)THEN
            print*,'Precalculating continuum from array : ',I
            CALL CALC_CONT(WING,MAXDV,NLAYER,NGAS,PRESS,TEMP,
     1       FRAC,IDGAS,ISOGAS,IPROC,IBS(I))
            IBD(I)=1
           ENDIF
          ENDDO
C	  Now calculate fine structure for this layer and wavenumber
C          print*,'Calculating fine structure. VV,Layer = ',VV,J

          CALL CALC_FINE(VV,WING,MAXDV,J,NLAYER,NGAS,PRESS,TEMP,
     1 FRAC,IDGAS,ISOGAS,IPROC,IBS,IFCONT,XK)

	  tautmp(J) = taucon(J) + xk
	  taugas(J) = taugasc(J) + xk
C          if(j.eq.1)print*,'X5',tautmp(j)

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
C		21	(Atm) Limb scattering calculation
C	then end the loop over the g ordinate.
C
C-----------------------------------------------------------------------

      DO Ipath = 1, npath
	nlays = nlayin(Ipath)
	output(Ipath) = 0.
	DO J = 1, nlays

		taus(J) = tautmp(layinc(J,ipath)) * scale(J,ipath)
		taur(J) = tauray(layinc(J,ipath)) * scale(J,ipath)

C                print*,'J,taus,taur',J,taus(J),taur(J)
C                print*,layinc(J,ipath),tautmp(layinc(J,ipath)),
C     1 tauray(layinc(J,ipath)),scale(J,ipath)


C                if(j.eq.41)print*,'X6',ipath,j,tautmp(layinc(J,ipath)),
C     &			taus(j),taur(j)


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
C                print*,'X7',taug(j)
C               Subtract Rayleigh scattering
                taug(J)=taug(J)-taur(J)

C                print*,'X8 taug',taug(J)

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

        Ig=1
	
	IF (imod(ipath).EQ.0) THEN

cc		WRITE(*,*) 'LBLrad_wave: Imod= 0 =Pure Transmission,',
cc     1			' creating output'

		DO J = 1, nlays
			output(Ipath) = output(Ipath) + taus(J)
		ENDDO
      		output(Ipath) = exp(-output(Ipath))
C                print*,'Ipath, output',Ipath,output(Ipath)

C               If a solar file exists and iform=4 we should multiply the
C               transmission by the solar flux
                if(iread.eq.999.and.iform.eq.4)then
                 call get_solar_wave(x,dist,xsolar)
                 output(Ipath)=output(Ipath)*xsolar
                endif



	ELSEIF (imod(ipath).EQ.1) THEN	

cc         WRITE(*,*) 'LBLrad_wave: Imod= 1 =Absorption, creating',
cc     1                  ' output'
 
		DO J = 1, nlays
			output(Ipath) = output(Ipath) + taus(J)
		enddo
		output(Ipath) = 1.- exp(-output(Ipath))

	ELSEIF (imod(ipath).EQ.3) THEN	

C		WRITE(*,*) 'LBLrad_wave: Imod= 3 =Emission, creating',
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
 
                         bb(J,Ipath)=planck_wave(ispace,x,
     1				emtemp(J,Ipath))
		        output(Ipath) = output(Ipath) + 
     1                  xfac*sngl(trold-tr) * bb(J,Ipath)
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
                 output(Ipath) = output(Ipath) +
     1                  xfac*sngl(trold)*radground

                endif


        ELSEIF (imod(ipath).eq.8) THEN
C             model 8, product of two path outputs
                output(Ipath)=output(layinc(1,Ipath))*
     1              output(layinc(2,Ipath))


	ELSEIF (imod(ipath).EQ.11) THEN	

cc                WRITE(*,*) 'LBLrad_wave: Imod= 11 =Contribution Function,',
cc     1                  ' creating output'

		trold = 1.
		do J = 1, nlays
			taud = taus(J)
			tr = dexp(-taud)

                         bb(J,Ipath)=planck_wave(ispace, 
     1			    x,emtemp(J,Ipath))

			output(Ipath) = output(Ipath) + 
     1                  sngl(trold-tr) * bb(J,Ipath)
 			trold = tr
                enddo


        ELSEIF (imod(ipath).eq.13) THEN

C          model 13, SCR sideband transmission (1-cell transmission)
                        taud = taus(J)

           output(Ipath)=1.0-exp(-taus(1))

        ELSEIF (imod(ipath).eq.14) THEN

C          model 14 SCR wideband transmission

           output(Ipath)=0.5*(1.0+exp(-taus(1)))


	ELSEIF (imod(ipath).EQ.15) THEN

C           WRITE(*,*) 'LBLrad_wave: Imod= 15 =Multiple Scattering, ',
C     1                  ' creating output'


		do J = 1, nlays

                         bb(J,Ipath)=planck_wave(ispace,x,
     1			   emtemp(J,Ipath))

			bnu(J) = bb(J,Ipath)
C                        print*,J,taus(J),taur(j)
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


C                         bb(J,Ipath)=planck_wave(ispace,x,
C     1                          emtemp(J,Ipath))

               
C               If tsurf > 0, then code assumes bottom layer is at the
C               surface and so uses tsurf to compute radiation upwelling.
C               Otherwise, code assumes optical depth is large and sets
C               upwelling radiation field to local temperature.

                if(tsurf.le.0.0)then
                 radg(1) = bb(nlays,Ipath)
                else
                 radg(1) = esurf*planck_wave(ispace,x,tsurf)
                endif

                galb1=galb

                if(galb1.lt.0)then
                         galb1 = get_albedo(nem,vem,emissivity,x)
C                         print*,x,galb1
                endif

                do J=1,nmu
                 radg(J)=radg(1)*sngl(1.0-galb1)
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
C                        print*,(mu1(i),i=1,nmu)
C                        print*,(wt1(i),i=1,nmu)
 
		  	call scloud11wave(rad1, sol_ang, emiss_ang,
     1                          aphi, radg, solar, lowbc, galb1, iray,
     2				mu1, wt1, nmu,   
     3				nf, Ig, x, vv, eps, omegas,bnu, taus, 
     4				taur,nlays, ncont,lfrac)

 		  	output(Ipath) = xfac*rad1


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

 		  	output(Ipath) = xfac*rad1

                ELSEIF (itype.EQ.13) THEN

C                   if(nf.ne.0)then
C                    print*,'Error in LBLrad_wave, nf <> 0',nf
C                   endif

   	  	   call scloud11flux(radg, solar, sol_ang, 
     1               lowbc, galb1, iray, mu1, wt1, nmu, nf, Ig, x,
     2               vv, eps, omegas, bnu, taus, taur, 
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
C                     print*,'sol_ang, xmu: ',sol_ang,xmu
                     do j=1,nmu
                      xt1 = abs(xmu-sngl(mu1(j)))
                      if(xt1.lt.tmp)then
                       j0=j
                       tmp=xt1
                      endif
C                      print*,j,j0,mu1(j)
                     enddo

                     solar1 = solar/(2.0*3.1415927*sngl(wt1(j0)))
C                     print*,'sol',j0,wt1(j0),solar,solar1
                     
C                     call dumpflux(LUNIS,ioff,nlays,nmu,nf,radg,umif,
C     1                uplf,I,x,nwave,Ig,Ig,j0,xmu,solar1,taus)

 		     output(Ipath) = xfac*umif(1,1,1)

                ELSE

			WRITE (*,*) ' Scattering ITYPE =',itype,
     1				' not defined for LBLrad_wave.'
			STOP
                ENDIF

	ELSEIF (imod(ipath).EQ.16) THEN

		WRITE(*,*) 'LBLrad_wave: Imod= 16 =Single Scat Approx,',
     1                  ' creating output'

		WRITE(*,*)' SINGLE SCATTERING IN USE '

                galb1 = galb
                        
                if(galb1.lt.0)then
                   galb1 = get_albedo(nem,vem,emissivity,x)
C                   print*,x,galb1
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
                 call get_solar_wave(x,xdist,xsolar)
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

C                 print*,J,taus(J),omega(J),fint(J),emtemp(J,ipath)
C                 print*,muemiss,muinc
C                 taud = taud + taus(J)*(1.0-omega(J))

                 taud = taud + taus(J)
                 tr = dexp(-taud)
C                 print*,taud,tausun,tr

		 output(Ipath) = output(Ipath) + 
     1			sngl(trold-tr)*ssfac*fint(J)*
     2			omega(J)*solar/(4*pi)



                 bb(J,Ipath)=planck_wave(ispace,x,
     1				emtemp(J,Ipath))


	         output(Ipath) = output(Ipath) + 
     1                  xfac*sngl(trold-tr) * bb(J,Ipath)


                 trold = tr
		ENDDO

C               Add in reflectance and emission from the ground
		bsurf=planck_wave(ispace,x,Tsurf)
                output(Ipath)=output(Ipath) + 
     1                          xfac*sngl(trold*solar*galb1)/pi +
     2				xfac*sngl(trold*esurf*bsurf)


	ELSEIF (imod(ipath).EQ.20)THEN

cc		WRITE(*,*) 'LBLrad_wave: Imod= 20 =Net flux calculation,',
cc     1                  ' creating output'

C		Computes up and down flux at BOTTOM of each layer
C		Compute once for Ipath=1

		if(Ipath.eq.1)then
                  do j=1,nlayin(npath)
                   fup(j)=0.0
                   fdown(j)=0.0
                   taus(layinc(J,npath))=tautmp(layinc(J,npath)) *
     &			scale(J,npath)

		    bb(layinc(j,npath),1)=planck_wave(ispace,
     1                  x,emtemp(j,npath))
		  enddo

	   	  do j=1,nlayin(npath)
           	   if(j.eq.1)then
		    fup(j)=bb(j,1)
		   else
		    tr = exp(-taus(j-1))
               fup(j)=fup(j-1)*sngl(tr) + sngl(1.0-tr)*bb(j,1)
		   endif
		  enddo 

	   	  do j=nlayin(npath),1,-1
                   tr=exp(-taus(j))
            	   if(j.eq.nlays)then
		    fdown(j)=sngl(1.0-tr)*bb(j,1)
		   else
               fdown(j)=sngl(1.0-tr)*bb(j,1)+sngl(tr)*fdown(j+1)
C		    print*,j,bb(j,1),tr,fdown(j)	
		   endif
		  enddo 

                endif 
            
	        output(Ipath) = fup(Ipath) - fdown(Ipath)
 	
	ELSEIF (imod(ipath).EQ.21) THEN

cc                WRITE(*,*) 'LBLrad_wave: Imod= 21 =Limb Scattering, ',
cc     1                  ' creating output'
C                print*,ispace,x
C       ************************************************************
C       First we need to compute the internal field at ALL altitudes
C       ************************************************************

		do J1 = 1, nlayer
                        J = nlayer+1-J1
                        bnu(J) = TEMP(J1)
                        bnu(J) = planck_wave(ispace,x,TEMP(J1))

			IF(TAUSCAT(J1).GT.0.0) THEN
				dtmp1 = dble(tauscat(J1))
				dtmp2 = dble(tautmp(J1))
				omegas(J) = dtmp1/dtmp2
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
                 radg(1) = bnu(nlayer)
                else
                 radg(1) = esurf*planck_wave(ispace,x,tsurf)
                endif

                galb1=galb

                if(galb1.lt.0)then
                         galb1 = get_albedo(nem,vem,emissivity,x)
                endif

                do J=1,nmu
                 radg(J)=radg(1)*sngl(1.0-galb1)
                enddo


C                WRITE (*,*) '     CALLING scloud11flux, nf = ',nf
		call scloud11flux(radg, solar, sol_ang, 
     1             lowbc, galb1, iray, mu1, wt1, nmu, nf, Ig, x,
     2             vv, eps,omegas, bnu, taus, taur, 
     3             nlayer, ncont, lfrac, umif, uplf)

C                do k=1,nlays
C                 print*,k,press(k),(umif(j,k,1),j=1,nmu)
C                 print*,k,press(k),(uplf(j,k,1),j=1,nmu)
C                enddo

                do J = 1, nlays
 
                        bb(J,Ipath)=planck_wave(ispace,x,
     1                     emtemp(J,Ipath))

                        bnu(J) = bb(J,Ipath)
                        IF(TAUSCAT(layinc(J,Ipath)).GT.0.0) THEN
                                dtmp1 = dble(tauscat(layinc(J,Ipath)))
                                dtmp2 = dble(tautmp(layinc(J,Ipath)))
                                omegas(J) = dtmp1/dtmp2
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

 		output(Ipath) = rad1

	ELSE
		WRITE (*,*) ' Imod = ', imod(ipath), ' is not a valid',
     1			' LBLrad_wave option'
		STOP

	ENDIF
      ENDDO

C-----------------------------------------------------------------------
C
C	Return and end.
C
C-----------------------------------------------------------------------

	RETURN

	END

************************************************************************
************************************************************************
