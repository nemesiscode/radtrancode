      SUBROUTINE cirsradg_wave (dist,inormal,iray,delh,nlayer,npath,
     1 ngas,press,temp,pp,amount,iwave,ispace,AMFORM,vwave,nlayin,
     2 layinc,cont,scale,imod,idgas,isogas,emtemp,itype,
     3 nem, vem, emissivity, tsurf,gtsurf,RADIUS1,flagh2p,hfp,
     4 radextra,nsw,isw,output,doutputdq)
C***********************************************************************
C_TITL:	CIRSRADG
C
C_DESC:
C
C_ARGS:	Input variables:
C	dist			REAL	Distance from Sun [AU].
C	inormal			INTEGER	Flag for ortho:para ratio
C					(0=equilibrium =1:1) (1=normal =>
C					3:1) Only used if FLAGH2P=0.
C	iray			INTEGER	Rayleigh optical depth flag
C       delh(nlayer)		REAL	Physical thickness [KM] of each
C					layer.
C	nlayer			INTEGER	Number of layers to consider.
C	npath			INTEGER	Number of individual paths
C					through the layers.
C       ngas			INTEGER	Number of gases to consider.
C       press(nlayer)		REAL	Total pressure [atm] for each
C					layer.
C       temp(nlayer)		REAL	Temperature [Kelvin] of each
C					layer.
C	pp(nlayer,ngas)		REAL	Partial pressure [atm] of each
C					gas.
C       amount(nlayer,ngas)	REAL	Number of molecules/cm2 for each
C					layer and gas (normal incidence).
C	iwave			INTEGER	Which wavenumber (need to
C					interpolate k-table correctly).

C       ispace          	INTEGER Indicates if wavelengths in vconv and
C                               	vwave are in wavenumbers(0) or
C                               	wavelengths (1)
C	AMFORM			INTEGER	Indicates if profile is one where
C					the sum of vmrs=1 (1). Equals 0 otherwise
C	vwave 			REAL 	Required wavenumber.
C	nlayin(npath)		INTEGER	Each layer defined above can be
C					included in the calculation for
C					each path several times. This is
C					the number of layers included in
C					the calculation for each path (if
C					you include it twice count it
C					twice).
C	layinc(nlayer,npath)	INTEGER	The layer numbers in order of
C					inclusion for each path. The order
C					also defines the processing of
C					each layer depending upon the
C					model.
C	cont(ncont,nlayer)	REAL	The number of dust particles/cm2
C					through each layer (vertical
C					path).
C	scale(nlayer,npath)	REAL	Scaling factor (eg one over
C					cosine of incidence) for each
C					included layer in each of the
C					paths. NOTE: when a layer is 
C					included more than once different
C					scaling factors can be defined.
C	imod(npath)		INTEGER	Type of model to use to compute
C					OUTPUT from transmissions for
C					each path.
C	idgas(ngas)		INTEGER	The LOCAL gas identifier for each
C					gas. The local gas identifier
C					agrees with HITRAN id's as far as
C					possible (Jan 1988) with extra
C					IDs for gases not in the HITRAN
C					compilation. eg.those in the GEISA
C					compilation.
C	isogas(ngas)		INTEGER	The local isotopic identifier, if
C					zero all isotopes of the gas are
C					included. Isotope IDs also agree
C					as far as possible with HITRAN
C					IDs. Similar (1,2,3...n) id's
C					have been defined for additional
C					gases.
C	emtemp(nlayer,npath)	REAL	Temperatures to use in emission
C					calculations.
C	itype			INTEGER	Value designating the chosen
C					scattering routine (currently not
C					implemented).
C       nem			integer	Number of points in surface
C					emissivity spectrum
C	vem(nem)		real	Wavelengths of emissivity spectrum
C	emissivity(nem)		real	Tabulated emissivities
C	tsurf			real	Surface temperature (K)
C	RADIUS1			real	Planetary radius (km) at reference
C					pressure (usually 1 bar)
C	RADEXTRA		real	Any correction to RADIUS1 needed for
C					disc-integration.
C	flagh2p			INTEGER	FLAGH2P=1 if a para-H2 fraction
C					profile has been read.
C	hfp(nlayer)		REAL	Para-H2 fractions in each layer if
C					FLAGH2P=1.
C	nsw			INTEGER	Number of variable elements for
C					which gradients are required.
C	isw(nsw)		INTEGER	ID of required variable
C					gradients. Code is:
C					1 to NGAS: Gas amount in each
C					  layer as listed in the .drv
C					  file.
C					NGAS+1: Temperature in each layer.
C					NGAS+2 to NGAS+2+NCONT-1: Aerosol
C					  OD in each layer.
C					NGAS+2+NCONT: para-H2 fraction in
C					  each layer.
C
C	Output variables:
C       output(npath)		REAL	Calculated Output values for each
C					path.
C	doutputdq(npath,nlayer,ngas+1+ncont+flagf2p)	...
C				REAL	Calculated rate of change of 
C					output with gas vmr, T, dust and
C					para-H2 fraction in each level.
C	gtsurf(npath)		real	Rate of change of output with surface
C					temperature for each path
C
C_CALL:	check_limits	Checks that array sizes are not exceeded.
C	spline		Returns an array that contains the second
C			derivatives of the interpolating function.
C	get_solar_wave	Calculates solar flux using a look up table.
C	get_hg		I presume, if scattering is set, gets the
C			Heney-Greenstein phase function .dat files (1
C			through 8) -- not documented very well.
C	get_kg		Interpolates/Returns the k-value at a selected
C			pressure and temperature.
C	splint		Returns a cubic spline interpolated value of Y
C			using  the result of spline.f.
C	ngascon		Computes a continuum absorbtion to any known
C			continuum spectra for a particular gas over a
C			defined wavenumber region.
C	nparacon	Computes gaseous continuum spectra from H2-H2 and
C			H2-He only but with variable para-H2 fraction.
C	nciacon		Computes gaseous continuum spectra from a variety
C			of gas pairs due to Collisionally
C			Induced-Absorption.
C	rayleighj	Evaluates the Rayleigh scattering cross section
C			for Jovian air using Allen (1976) astrophysical
C			quantities.
C	noverlapg	Combines the absorption coefficient distributions
C			of two overlapping gases.
C	planck_wave	Calculates the Planck function.
C	planckg_wave	Calculate the dB/dT of the Planck function.
C
C_HIST:	
C	???????	PGJI	Conversion to gradient version.
C	7aug03	NT	corrected bug on line511 cont(j,k) -> cont(k,j).
C	14aug03	NT	variable 'solar' initialised properly
C	21jun05	NT	spline interpolation replaced by linear interp
C			for the aerosol xsc interp (splines were unstable for
C			small numbers of points)
C***************************** VARIABLES *******************************

      IMPLICIT NONE


C The include file ...
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'
C Defines the maximum values for a series of variables (layers, bins,
C paths, etc.)


C Definition of input and output variables ...
      INTEGER inormal,iwave,ispace,itype,flagh2p
      INTEGER nlayer,npath,ngas,nsw,AMFORM,nlay1
      INTEGER idgas(maxgas),isogas(maxgas),isw(maxgas+2+maxcon)
      INTEGER layinc(maxlay,maxpat),nlayin(maxpat),imod(maxpat)
      REAL xdist,dist,vwave,delh(nlayer),esurf,radground
      REAL cont(maxcon,maxlay),tsurf,bsurf,gtsurf(maxpat)
      REAL tempgtsurf(maxpat,maxg),gradground
      REAL press(maxlay),temp(maxlay),hfp(maxlay)
      REAL scale(maxlay,maxpat),emtemp(maxlay,maxpat)
      REAL pp(maxlay,maxgas),amount(maxlay,maxgas)
      REAL output(maxpat),doutputdq(maxpat,maxlay,maxgas+2+maxcon)
      REAL rad1,xfac,RADIUS1,p1,p2


C Definition of general variables ...
      INTEGER i,j,k,l,ipath,ig,iray,ipath1,lstcel,ipath2
      INTEGER k1,nlays,nparam,ILBL
      INTEGER jnh3,jch4,jh2,jhe,igas
      REAL ppp(maxgas),aamount(maxgas),vv
      REAL dbdt(maxlay,maxpat),bb(maxlay,maxpat)
      REAL fpara,xray
      REAL rayleighj,rayleigha,rayleighv,planck_wave,planckg_wave
      REAL rayleighls,fheh2,fch4h2,fh2,fhe,fch4,fnh3
      REAL p,t,tau,tau2,dpexp,dist1
      REAL muemiss,intscat,pastint
      REAL utotl(maxlay),taucon(maxlay),tauscat(maxlay)
      REAL taus(maxlay),tautmp(maxlay),f(maxlay),fint(maxlay)
      REAL asec(maxsec),asec2(maxsec),bsec(maxsec),bsec2(maxsec)
      REAL xsec2(2,maxcon,maxsec),tauray(maxlay),x
      REAL tausc(maxcon),phase(maxcon),taur(maxlay)
C PHASE: Resulting phase function from get_hg?
      REAL frac(maxlay,maxgas)
      REAL tmp,corkout(maxpat,maxg)
C TMP: ??? (relates to corkout)

      REAL dtolddq(maxlay,maxgas+2+maxcon),
     1 dtrdq(maxlay,maxgas+2+maxcon),
     2 dtaudq(maxlay,maxgas+2+maxcon),
     3 de3olddq(maxlay,maxgas+2+maxcon),
     4 de3dq(maxlay,maxgas+2+maxcon)
      DOUBLE PRECISION e3old,e3new,e3grad
      REAL dcoutdq(maxpat,maxg,maxlay,maxgas+2+maxcon)
      DOUBLE PRECISION tr,trold,taud,dpi,dphase,calpha,tlayer

C DPHASE: ??? Outputted from get_hg.f.
C CALPHA: ??? Input into get_hg.f.
      DOUBLE PRECISION dmuinc,dmuemiss,draddeg,dinc,demiss

C Dust variables ...
      INTEGER nsec,ncont,npro
      real zheight(maxpro)
C NSEC: Number of wavelengths for which dust cross-sections are defined.
C NCONT: Number of dust types included.

      REAL vsec(maxsec),xsec(2,maxcon,maxsec)
C VSEC: Dust cross sections for each dust type.
C XSEC: Corresponding wavelengths for dust x-sections.


C K-table variables ...
      INTEGER npk,ntk,ng
C NPK: Number of pressures in the k-table.
C NTK: Number of temperatures in the k-table.
C NG: Number of Gauss-Legendre ordinates in k-distribution.
      INTEGER lun(maxbin,maxgas),ireck(maxbin,maxgas) 
      INTEGER lunlbl(maxgas),irec0lbl,npointk,npklbl,ntklbl
      REAL fwhmk,delvk,xminklbl,delklbl
      REAL pklbl(maxk),tklbl(maxk),t2klbl(maxk,maxk)
      REAL koutlbl(maxlay,maxgas)
      REAL dkoutdtlbl(maxlay,maxgas)
C FWHMK: Full-Width Half Maximum of the k-table.
      REAL pk(maxk),tk(maxk),t2k(maxk,maxk)
C PK: K-table pressures [atm].
C TK: K-table temperatures [Kelvin].
      REAL xmink(maxbin,maxgas),delk(maxbin,maxgas) 
C XMINK: K-table wavenumber minimum [cm-1].
C DELK: K-table point-spacing [cm-1].
      REAL amo(maxgas)
      REAL k_gn(maxg,maxgas),dkgndt(maxg,maxgas)
      REAL k_g(maxg),delg(maxg),g_ord(maxg)
      REAL frack(maxbin,maxgas)
C K_G: Calculated k-distribution.
C DEL_G: Gauss-Legendre weights for integration.
C G_ORD: Gauss-Legendre ordinates for calculating the k-distribution.
      REAL kl_g(maxg,maxlay)
      REAL kout(maxlay,maxgas,maxg),dkoutdt(maxlay,maxgas,maxg)
      REAL dkdq(maxg,maxgas+1)
      REAL dkdql(maxg,maxgas+1,maxlay)
      REAL dtaucondq(maxlay,maxgas+2+maxcon),
     1 dtautmpdq(maxlay,maxgas+2+maxcon),
     2 dtauscadq(maxlay,maxgas+2+maxcon)


C Scattering variables ...
      INTEGER isol,nmu,lowbc,nf,first
C ISOL: 0=Sun off, 1=Sun on.
C NMU: Number of zenith ordinates.
C LOWBC: Lower boundary condition: thermal=0, Lambert=1.
C NF: Number of Fourier azimuth coeffs.
      INTEGER liscat(maxcon),lncons(maxcon),lnorm(maxcon)

      REAL solar,sol_ang,emiss_ang,aphi
      DOUBLE PRECISION GALB,galb1,get_albedo
C GALB: Ground albedo.
C SOL_ANG: Solar zenith angle.
C EMISS_ANG: Viewing zenith angle.
C APHI: Azimuth angle.
      REAL bnu(maxlay),bnug(maxlay),radg(maxmu)
      REAL tauscatl(maxcon,maxlay),tauscat1(maxcon,maxlay)
      REAL lcons(maxcon,maxscatpar),xsolar
      INTEGER NALB
      real alb(maxsec),valb(maxsec)
      INTEGER nem,nbase,ioff
      REAL vem(maxsec),emissivity(maxsec),interpem,radextra


      DOUBLE PRECISION mu1(maxmu),wt1(maxmu)
      DOUBLE PRECISION pplsto(0:maxscatlay,maxmu,maxmu,maxscatlay),
     1 pmisto(0:maxscatlay,maxmu,maxmu,maxscatlay)

      LOGICAL scatter,single


C Misc variables or continuum variables ...
      INTEGER ii,idump,j1,iabsorb(5)
C IDUMP: If =1 then print diagnostic print statements.
C IABSORB: Flag set to the gas numbers corresponding to the active gas
C gradients in DABSROB.

      REAL xlen,avgcontmp,dabsorb(7)
C XLEN: ? (relates to delh)
C AVGCONTMP: Calculated optical depth.
C DABSORB: Rate of change of optical depth with: 1: H2, 2: He, 3: N2, 
C 4: CH4 and 5: CO2 amounts, 6: temperature, 7: para-H2 fraction.

C     Solar reference spectrum common block
      real swave(maxbin),srad(maxbin),solrad
      integer iread,nspt,iform
      common /solardat/iread,iform,solrad,swave,srad,nspt

C Common blocks ...
      COMMON /dust/ vsec,xsec,nsec,ncont
      COMMON /interpk/ lun,ireck,xmink,delk,frack,pk,npk,tk,t2k,ntk,
     1 ng, delvk,fwhmk,g_ord,delg,kout,dkoutdt
      common/interpklbl/lunlbl,irec0lbl,xminklbl,delklbl,npointk,
     1    pklbl,npklbl,tklbl,t2klbl,ntklbl,koutlbl,dkoutdtlbl
      common /lbltable/ILBL

      COMMON /scatd/ mu1,wt1,galb 
      common/alb/nalb,valb,alb
      COMMON /scatter1/ nmu,isol,dist1,lowbc,liscat,lnorm,
     1 lncons,lcons,sol_ang,emiss_ang,aphi,nf
      COMMON /phasesto/ pplsto,pmisto
      COMMON /initial/ first


C********************************* CODE ********************************

      CALL check_limits(nlayer,npath,ngas,ncont,nsec,
     1 ng,npk,ntk,nlayin,layinc)

      IF(ILBL.EQ.2)THEN
       NG=1
       DELG(1)=1.0
      ENDIF


C=======================================================================
C
C     	Initialise some variables.
C
C		imod = 15		Multiple scattering
C		imod = 16		Single scattering approx
C
C=======================================================================
      dpi = 4.0d0*datan(1.0d0)
      draddeg = 2.0d0*dpi/360.0d0

      scatter = .FALSE.
      single = .FALSE.
      DO ipath = 1, npath
        IF(imod(ipath).EQ.15)THEN
          scatter = .TRUE.
          first = 0                       ! Interp_phase will initialise
        ENDIF
        IF(imod(ipath).EQ.16)THEN
          single = .TRUE.
          dinc = DBLE(sol_ang)
          demiss = DBLE(emiss_ang)
          dmuinc = COS(dinc * draddeg)
          dmuemiss = COS(demiss * draddeg) 
          muemiss = SNGL(dmuemiss)
          DO i=1,nlayer
            scale(i,1) = SNGL(1./dmuemiss + 1./dmuinc)
          ENDDO
          calpha = 1.0*dmuinc*dmuemiss - SQRT(1.0d0 - dmuinc**2)*
     1    SQRT(1.0d0 - dmuemiss**2)*COS(aphi*draddeg)
        ENDIF
      ENDDO

      nparam = ngas + 1 + ncont
      IF(flagh2p.EQ.1)nparam = nparam + 1

C     Set vv to the current WAVENUMBER
      vv = vwave
      x = vv
      esurf = interpem(nem,vem,emissivity,x)
      if(ispace.eq.1)vv=1e4/vwave

C=======================================================================
C
C	Precompute volume fractions. The factor of 1e-20 applied to
C	utotl arises because in order to avoid 
C	underflow, the CIRSRAD routines multiply k-coefficients in
C	the k-tables by 1e20
C
C=======================================================================
 
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



      DO i=1,nlayer
        utotl(i) = 0.0
        DO j=1,ngas
          frac(i,j) = pp(i,j) / press(i)
          utotl(i) = utotl(i) + amount(i,j)
C          print*,i,j,amount(i,j),utotl(i)
        ENDDO
        utotl(i) = utotl(i)*1.e-20
      ENDDO

C=======================================================================
C
C	Begin loop over wavenumber bins. At each wavenumber, first 
C	initialise arrays over the levels, then step over the g 
C	ordinates, calculating the relevant output.
C
C	xsec(1....) is the extinction coefficient, xsec(2....) the 
C	single scattering albedo.
C
C=======================================================================

      IF(nsec.GT.1)THEN
        DO k=1,ncont
          DO l=1,nsec
            asec(k) = xsec(1,k,l)
            bsec(l) = xsec(2,k,l)
          ENDDO

          if(vwave.lt.vsec(1).or.
     &          vwave.gt.vsec(nsec))then
           print*,'Error in cirsradg_wave'
           print*,'wavelengths in xsc file not consistent with'
           print*,'calculation wavelength : ',vwave
           print*,'vsec : ',vsec(1),vsec(nsec)
           stop
          endif
cc Nick T; sack this bit off and use linear interpolation instead as
cc more stable for widely spaced points
cc          CALL spline (vsec,asec,nsec,1.e30,1.e30,asec2)
cc          CALL spline (vsec,bsec,nsec,1.e30,1.e30,bsec2)
cc
cc          DO l=1,nsec
cc            xsec2(1,k,l) = asec2(l)
cc            xsec2(2,k,l) = bsec2(l)
cc          ENDDO
        ENDDO
      ENDIF

      IF( ((scatter).OR.(single)) .and. (isol.eq.1)) THEN
          CALL get_solar_wave(vwave,dist,solar)
      ELSE
          solar = 0.0
      ENDIF

      IF(single)THEN
        DO j=1,ncont
          CALL get_hg(vwave,calpha,ncont,j,dphase)
          phase(j) = SNGL(dphase)
        ENDDO
      ELSE
        DO j=1,ncont
          phase(j) = 0.0
        ENDDO
      ENDIF

C=======================================================================
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
C=======================================================================

      intscat = 0.0
      pastint = 0.0

C Read in k-coefffients for each gas and each layer
      IF(ILBL.EQ.0)THEN
        CALL get_kg(nlayer,press,temp,ngas,iwave,vwave)
      ELSE
        CALL get_klblg(nlayer,press,temp,ngas,vwave)
      ENDIF
 
C      print*,'OK here'  
 
      DO j=1,nlayer
        taucon(j) = 0.0
        tauscat(j) = 0.0
        DO k=1,(ngas+ncont+2)
          dtaucondq(j,k) = 0.0
        ENDDO
        f(j) = 0.0
        p = press(j)
        t = temp(j)

        IF(JNH3.GT.0.)FNH3 = FRAC(J,JNH3)
        IF(JCH4.GT.0.)FCH4 = FRAC(J,JCH4)
        IF(JH2.GT.0.)FH2 = FRAC(J,JH2)
        IF(JHE.GT.0.)FHE = FRAC(J,JHE)

C=======================================================================
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
C=======================================================================

        DO k=1,ncont
          tausc(k) = 0.0
          IF (nsec.EQ.1) THEN
            tau = xsec(1,k,1)
            tau2 = xsec(2,k,1)
          ELSE
            DO l=1,nsec
              asec(l) = xsec(1,k,l)
              bsec(l) = xsec(2,k,l)
ccNick T: use linear interpolation instead
cc              asec2(l) = xsec2(1,k,l)
cc              bsec2(l) = xsec2(2,k,l)
            ENDDO
cc            CALL splint (vsec, asec, asec2, nsec, vwave, tau)
cc            CALL splint (vsec, bsec, bsec2, nsec, vwave, tau2)
            call interp(vsec, asec, nsec, tau, vwave)
            call interp(vsec, bsec, nsec, tau2, vwave)
          ENDIF
          taucon(j) = taucon(j) + tau*cont(k,j)
          tauscat(j) = tauscat(j) + tau2*cont(k,j)
          tausc(k) = tau2*cont(k,j)
	  f(j) = f(j) + tau2*cont(k,j)*phase(k)*solar/(4.*pi)
          dtaucondq(j,ngas+1+k) = tau
          dtauscadq(j,ngas+1+k) = tau2
          tauscat1(k,j)=tausc(k)
        ENDDO
        intscat = intscat + tauscat(j)
        pastint = intscat


C=======================================================================
C
C	Add continuum contributions from hydrogen and other gases into
C	the layer optical depths. 
C
C	GASCON is a RADTRANS routine which calculates continuum
C	contributions from various gases.
C
C=======================================================================

        DO k=1,ngas
C Computes a polynomial approximation to any known continuum spectra 
C for a particular gas over a defined wavenumber region.
          CALL ngascon(vv,idgas(k),isogas(k),amount(j,k),
     1    pp(j,k),press(j),temp(j),avgcontmp)

          taucon(j) = taucon(j) + avgcontmp

C Compute rate of change of radiance with amount. Neglect temperature
C dependences for now (K=ngas+1)
          IF(avgcontmp.GT.0.AND.amount(j,k).GT.0.0)THEN
            dtaucondq(j,k) = avgcontmp/amount(j,k)
          ENDIF
        ENDDO

C=======================================================================
C
C	NCIACON: to compute gaseous continuum spectra from a variety of
C	gas pairs due to Collisionally Induced-Absorption.
C
C	NPARACON: to compute gaseous continuum spectra from H2-H2 and
C	H2-He only but with variable para-H2 fraction.
C
C=======================================================================

C Reduce the 2d arrays for gas amounts (no./unit vol.) and partial
C pressure to 1d for passing to subroutine
        DO 320 ii=1,ngas
          aamount(ii) = amount(j,ii)
          ppp(ii) = pp(j,ii)
320     CONTINUE

C idump= 0 means no diagnostic print statements while in NCIACON.
        idump = 0
        xlen = delh(j)

        IF(flagh2p.EQ.1)THEN
          fpara = hfp(j)
          CALL NPARACON_ALL(vv,p,t,ngas,idgas,isogas,aamount,
     1    ppp,fpara,xlen,avgcontmp,iabsorb,dabsorb,idump)
        ELSE
          CALL NCIACON(vv,p,t,inormal,ngas,idgas,isogas,aamount,
     1    ppp,xlen,avgcontmp,iabsorb,dabsorb,idump)
        ENDIF
c	idump=0
c	print*,'IABSORB:',iabsorb
	
        taucon(j) = taucon(j) + avgcontmp

        DO ii=1,5
          k = iabsorb(ii)
          IF(k.GT.0)THEN
            dtaucondq(j,k) = dtaucondq(j,k) +
     1      dabsorb(ii)/(utotl(j)*1.0e20)
          ENDIF
        ENDDO


C Temperature:
        dtaucondq(j,ngas+1) = dtaucondq(j,ngas+1) + dabsorb(6)
C Para-H2 fraction:
        IF(flagh2p.EQ.1)THEN
          dtaucondq(j,nparam) = dtaucondq(j,nparam) + dabsorb(7)
        ENDIF

        IF(iray.GT.0.AND.itype.NE.0)THEN
          IF(iray.EQ.1)then
           xray = RAYLEIGHJ(vv,press(j),temp(j))*1E20
          ELSEIF(iray.eq.2)then
           xray = RAYLEIGHV(vv,press(j),temp(j))*1E20
          ELSEIF(iray.eq.3)then
           xray = RAYLEIGHA(vv,press(j),temp(j))*1E20
          ELSE
           if(FH2.GT.0.0)THEN
             fheh2=FHE/FH2
             fch4h2=FCH4/FH2
           else
             fheh2=0.
             fch4h2=0.
           endif
           xray = RAYLEIGHLS(vv,fheh2,fch4h2,fnh3)*1E20
          ENDIF

          avgcontmp = utotl(j)*xray
          tauray(j) = avgcontmp

          taucon(j) = taucon(j) + tauray(j)

C Compute rate of change of radiance with gas amounts. Neglect temperature
C gradient for now (i.e. k = ngas + 1)
          DO k=1,ngas
            dtaucondq(j,k) = dtaucondq(j,k) + xray*1.0e-20
          ENDDO
        ENDIF


C=======================================================================
C
C	For each gas, look up the k-coefficients. Stepping through the 
C	gases, combine each set of coeffs with the combination of all
C	previous gases using OVERLAP. Finally, store the layer dependent
C	summed K coefficients. End the loop over each layer.
C
C=======================================================================

        DO k=1,ngas
          amo(k) = amount(j,k)*1.0e-20
          if(ILBL.EQ.0)THEN
           IF(lun(iwave,k).LT.0) THEN
            DO l=1,ng
              k_gn(l,k) = 0.0
              dkgndt(l,k) = 0.0
            ENDDO
           ELSE
            DO l=1,ng
              k_gn(l,k) = kout(j,k,l)
              dkgndt(l,k) = dkoutdt(j,k,l)
            ENDDO
           ENDIF
          else
           if(lunlbl(k).lt.0)then
                k_gn(1,k) = 0.0
                dkgndt(1,k) = 0.0
           else
                k_gn(1,k) = koutlbl(j,k)
                dkgndt(1,k) = dkoutdtlbl(j,k)
           endif
          endif

C          if(ng.eq.1)then
C             print*,K,k_gn(1,k)
C          else
C             print*,K,k_gn(ng/2,k)
C          endif

        ENDDO

        IF(ILBL.EQ.0)THEN
         CALL noverlapg(idump,delg,ng,ngas,amo,k_gn,dkgndt,k_g,
     1    dkdq)
        ELSE
         CALL noverlapg1(idump,delg,ng,ngas,amo,k_gn,dkgndt,k_g,
     1    dkdq)
        ENDIF


        DO l=1,ng
          kl_g(l,j) = k_g(l)
C Need to correct for the 1e20 term again to get the right gradient with
C absorber amount.
          DO k=1,ngas
            dkdql(l,k,j) = dkdq(l,k)*1.0e-20
          ENDDO
C Temperature
          dkdql(l,ngas+1,j) = dkdq(l,ngas+1)
        ENDDO
C        if(ng.eq.1)then
C         print*,'K',J,press(J),kl_g(1,j)
C        else
C         print*,'K',J,press(J),kl_g(ng/2,j)
C        endif
C        stop
      ENDDO

C=======================================================================
C
C	Now we must step over each g-ordinate.
C
C	At each g ordinate, the optical depth of the layer is given by
C	the continuum value + k u where u is the integral over the 
C	absorber number along the path (.i.e. the sums of the 
C	constituent absorber numbers per unit area). Units of K are
C	(molecules number)**-1 cm**2 x 1.e20 (the later is taken care
C	of earlier in the code).
C
C=======================================================================


      DO ig=1,ng
        DO j=1,nlayer
          tautmp(j) = taucon(j) + kl_g(ig,j) 
C          print*,j,press(j),taucon(j),kl_g(ig,j),tautmp(j)
C Continuum component to gradients:
          DO k=1,nparam
            dtautmpdq(j,k) = dtaucondq(j,k)
          ENDDO

C K-data contribution to gradients:
          DO k=1,ngas+1 
            dtautmpdq(j,k) = dtautmpdq(j,k) + dkdql(ig,k,j)
          ENDDO

        ENDDO


C      ***********************************************************
C      If(AMFORM.EQ.1)then we need to adjust the gradients slightly, since
C      the sum of vmrs at all levels must add up to 1.

C       IF(AMFORM.EQ.1)THEN
C        CALL ADJUSTGRAD(dtautmpdq,amount,ngas,nlayer)
C       ENDIF


C      ***********************************************************

C=======================================================================
C
C	Step through the number of paths, calculating the required 
C	output properties. These are
C
C	      Imod
C		0	(Atm) Pure transmission
C		1	(Atm) Absorption (useful for small transmissions)
C		3	(Atm) Emission. Planck function evaluated at bin 
C               8       (Combined Cell,Atm) The product of two
C                               previous output paths.                    
C               13      (Atm) SCR Sideband
C               14      (Atm) SCR Wideband
C		15	(Atm) Multiple scattering.
C		18	(Atm) Emission into hemisphere (for disc-averaging)
C
C	then end the loop over the g-ordinate.
C
C=======================================================================

        DO ipath=1,npath
          nlays = nlayin(ipath)
          corkout(ipath,ig) = 0.0
          tempgtsurf(ipath,ig) = 0.0
          DO j=1,nlays
            DO k=1,nparam
              dcoutdq(ipath,ig,j,k) = 0.0
            ENDDO
          ENDDO
          DO j=1,nlays
            taus(j) = tautmp(layinc(j,ipath)) * scale(j,ipath)
            taur(j) = tauray(layinc(j,ipath)) * scale(j,ipath)
C            print*,'tau',j,ipath,layinc(j,ipath),scale(j,ipath),
C     1		tautmp(layinc(j,ipath)),taus(j)

            do k=1,ncont
             tauscatl(k,j)=tauscat1(k,layinc(j,ipath)) * scale(j,ipath)
            enddo
C            print*,J,taus(J),tautmp(layinc(J,Ipath)),
C     1                    tauscat(layinc(J,Ipath)),utotl(J)
            DO k=1,nparam
              dtaudq(j,k) = dtautmpdq(layinc(j,ipath),k)*scale(j,ipath)
            ENDDO
            IF(single)THEN
              fint(j) = f(layinc(j,ipath))/muemiss
            ENDIF

          ENDDO


          IF(imod(ipath).EQ.0)THEN
cc	      WRITE(*,*)' CIRSRADG.f :: Imod= 0 ==> pure transmission'
cc	      WRITE(*,*)' creating output ...'
            DO j=1,nlays
              corkout(ipath,ig) = corkout(ipath,ig) + taus(j)
            ENDDO
            tmp = corkout(ipath,ig)
            corkout(ipath,ig) = dpexp(-tmp)

C           If a solar file exists and iform=4 we should multiply the
C           transmission by the solar flux
            xfac=1.
            if(iread.eq.999.and.iform.eq.4)then
              call get_solar_wave(vwave,dist,xsolar)
              xfac=xsolar
            endif

            DO j=1,nlays
              DO k1=1,nsw
                k=isw(k1)
                dcoutdq(ipath,ig,j,k)=-tmp*corkout(ipath,iG)*
     1          dtaudq(j,k)*xfac
              ENDDO
            ENDDO

            corkout(ipath,ig)=corkout(ipath,ig)*xfac

          ELSE IF(imod(ipath).EQ.1)THEN
cc	      WRITE(*,*)' CIRSRADG.f :: Imod= 1 ==> absorption calc.'
cc	      WRITE(*,*)' creating output ...'
            DO j=1,nlays
              corkout(ipath,ig) = corkout(ipath,ig) + taus(j)
            ENDDO
            tmp = corkout(ipath,ig)
            corkout(ipath,ig) = 1. - dpexp(-tmp)
            DO j=1,nlays
              DO k1=1,nsw
		k=isw(k1)
                dcoutdq(ipath,ig,j,k) = tmp*
     1          (1 - corkout(ipath,ig))*dtaudq(j,k)
              ENDDO
            ENDDO
          ELSE IF(imod(ipath).EQ.3)THEN
cc	      WRITE(*,*)' CIRSRADG.f :: Imod= 3 ==> emission calc.'
cc	      WRITE(*,*)' creating output ...'
            trold = 1.
            taud = 0.0

            DO j=1,nlays
              DO k=1,nparam
                dTolddq(j,k) = 0.0
                dtrdq(j,k) = 0.0
              ENDDO
            ENDDO


C           The code below calculates the radiance
C           spectrum in units of W cm-2 sr-1 (cm-1)-1 or W cm-2 sr-1 um-1.
C
            xfac=1.
C           If iform = 1 or iform = 3 we need to calculate the total 
C           spectral flux from the planet.
            if(iform.eq.1.or.iform.eq.3)then
              xfac=xfac*pi*4.*pi*((RADIUS1+radextra)*1e5)**2
            endif
C           If a solar file exists and iform=1 we should divide the planet
C           flux by the solar flux
            if(iread.eq.999.and.iform.eq.1)then
C            Set dist to -1 to get total power spectrum of star
             xdist=-1.0
             call get_solar_wave(vwave,xdist,xsolar)
             xfac=xfac/xsolar
            endif
C           If doing integrated flux from planet need a factor to stop the
C           matrix inversion crashing
            if(iform.eq.3)xfac=xfac*1e-18

            DO J=1,nlays
              IF(ig.EQ.1)THEN
                 bb(j,ipath) = planck_wave(ispace,vwave,
     1					emtemp(j,ipath))
                 dbdt(j,ipath) = planckg_wave(ispace,vwave,
     1					emtemp(j,ipath))
              ENDIF

              tlayer = dexp(-dble(taus(j)))
              taud = taud + taus(j)
              tr = trold*tlayer

              corkout(ipath,ig) = corkout(ipath,ig) + 
     1        xfac*sngl(trold - tr)*bb(J,Ipath)


              DO k1=1,nsw
                k = isw(k1)
                DO j1=1,j
                  IF(j1.LT.j)THEN
                    dtrdq(j1,k) = dtolddq(j1,k)*sngl(tlayer)
                  ELSE
                    dtrdq(j1,k) = -dtaudq(j1,k)*sngl(tlayer*trold)
                  ENDIF
                  dcoutdq(ipath,ig,j1,k) =
     1                 dcoutdq(ipath,ig,j1,k) + xfac*bb(j,ipath)*
     2                     (dtolddq(j1,k) - dtrdq(j1,k))
                ENDDO
                IF(k.EQ.ngas+1)THEN
                  dcoutdq(ipath,ig,j,k) =
     1       dcoutdq(ipath,ig,j,k)+xfac*sngl(trold-tr)*dbdt(j,ipath)
               ENDIF
              ENDDO

              trold = tr
              DO j1=1,j
                DO k1=1,nsw
                  k = isw(k1)
                  dtolddq(j1,k) = dtrdq(j1,k)
                ENDDO
              ENDDO

            ENDDO

C            if(ig.eq.5)then
C            DO j=1,nlays
C               DO k1=1,nsw
C                 k = isw(k1)
C                 print*,'rad',ipath,ig,j,k1,k,dcoutdq(ipath,ig,j,k)
C               ENDDO
C            ENDDO
C            endif
 
C           Check to see if this is a limb path
            j1=int(0.5*nlays)

            p1=press(layinc(j1,ipath))
            p2=press(layinc(nlays,ipath))


C           If not a limb path, add on surface radiance
            if(p2.gt.p1)then

             if(tsurf.le.0.0)then
               radground = bb(nlays,Ipath)
               gradground = dbdt(nlays,Ipath)
             else
               radground = esurf*planck_wave(ispace,vwave,tsurf)
               gradground = esurf*planckg_wave(ispace,vwave,tsurf)
             endif

             corkout(Ipath,Ig) = corkout(Ipath,Ig) +
     1               xfac*sngl(trold)*radground

C            gradient wrt surface temperature
             tempgtsurf(ipath,ig) = xfac*sngl(trold)*gradground

C            gradient with respect to atmospheric properties
             DO j=1,nlays
                DO k1=1,nsw
                  k = isw(k1)
                   dcoutdq(ipath,ig,j,k) =  
     1        dcoutdq(ipath,ig,j,k)+xfac*radground*dtolddq(j,k)
C                 print*,'radg: j,k1,k,ig,ipath,',j,k1,k,ig,ipath,
C     1			dcoutdq(ipath,ig,j,k)
                ENDDO
             ENDDO
            endif

C            if(ig.eq.5)then
C            DO j=1,nlays
C               DO k1=1,nsw
C                 k = isw(k1)
C                 print*,'radX',ipath,ig,j,k1,k,dcoutdq(ipath,ig,j,k)
C               ENDDO
C            ENDDO
C            endif


          ELSEIF(imod(ipath).eq.8)then

C           model 8, product of two path outputs

            ipath1=layinc(1,Ipath)
            ipath2=layinc(2,Ipath)

C           Assume one of these must be an atmospheric calculation
C	    and one must be a cell calculation. Assume that first
C           one is a cell transmission and the second one is an
C           atmosphere

C            print*,'ipath1, imod',ipath1,imod(ipath1)
C            print*,'ipath2, imod',ipath2,imod(ipath2)            

            if(imod(ipath1).lt.13.or.imod(ipath1).gt.14)then
             print*,'Error in cirsradg_wave: cell and atm'
             print*,'paths not in expected order. Stopping'
             stop
            endif
           
C           Multiply two outputs together
            corkout(ipath,Ig)=corkout(ipath1,Ig)*
     1              corkout(ipath2,Ig)

C            print*,layinc(1,ipath),layinc(2,ipath)
C            print*,corkout(ipath1,ig),
C     1		corkout(ipath2,Ig)


            nlays=nlayin(ipath2)

            DO j=1,nlays
             do k1=1,nsw
                k = isw(k1)
                dcoutdq(ipath,ig,j,k) = corkout(ipath1,Ig)*
     1		  dcoutdq(ipath2,ig,j,k)
C                  if(ig.eq.5)then
C                  print*,'mod 8',j,k1,k,ig,ipath,ipath2,
C     1		  corkout(ipath1,Ig),
C     1		  dcoutdq(ipath2,ig,j,k),dcoutdq(ipath,ig,j,k)
C                  endif
             enddo
            ENDDO


C           gradient wrt surface temperature
            tempgtsurf(ipath,ig) = tempgtsurf(ipath2,ig)*
     1		corkout(ipath1,Ig)

          ELSEIF (imod(ipath).eq.13) THEN

C          model 13, SCR sideband transmission (1-cell transmission)
                        taud = taus(J)

           corkout(Ipath,Ig)=1.0-exp(-taus(1))
C           print*,'SIDEBAND: ',taus(1),1.0-exp(-taus(1))
           LSTCEL=IPATH

          ELSEIF (imod(ipath).eq.14) THEN

C          model 14 SCR wideband transmission

           corkout(Ipath,Ig)=0.5*(1.0+exp(-taus(1)))
C           print*,'WIDEBAND: ',taus(1),0.5*(1.0+exp(-taus(1)))

           LSTCEL=IPATH


          ELSEIF(imod(ipath).EQ.15)then

C           WRITE(*,*) 'CIRSRADG_WAVE: Imod= 15 =Multiple Scattering, ',
C     1                  ' creating output'

            print*,'CIRSRADG_WAVE. Multiple scattering not implemented.'

            STOP


          ELSE IF(imod(ipath).EQ.18)THEN
C         Disc-integrated emission
cc	      WRITE(*,*)' CIRSRADG.f :: Imod= 18 ==> hemisphere emission'
cc	      WRITE(*,*)' creating output ...'

            trold = 1.0
            e3old = 0.5

            DO j=1,nlays
              DO k=1,nparam
                dTolddq(j,k) = 0.0
                dtrdq(j,k) = 0.0
                de3olddq(j,k) = 0.0
                de3dq(j,k) = 0.0
              ENDDO
            ENDDO

C           The code below calculates the hemispherically-integrated emission 
C           spectrum in units of W cm-2 (cm-1)-1 or W cm-2 um-1.
C
            xfac=1.
C           If iform = 1 or iform = 3 we need to calculate the total 
C           spectral flux from the planet.
            if(iform.eq.1.or.iform.eq.3)then
             xfac=xfac*4.*pi*((RADIUS1+radextra)*1e5)**2
            endif
C           If a solar file exists and iform=1 we should divide the planet
C           flux by the solar flux
            if(iread.eq.999.and.iform.eq.1)then
C            Set dist to -1 to get total power spectrum of star
             xdist=-1.0
             call get_solar_wave(vwave,xdist,xsolar)
             xfac=xfac/xsolar
            endif
C           If doing integrated flux from planet need a factor to stop the
C           matrix inversion crashing
            if(iform.eq.3)xfac=xfac*1e-18

            DO J=1,nlays
              IF(ig.EQ.1)THEN
                 bb(j,ipath) = planck_wave(ispace,vwave,
     1					emtemp(j,ipath))
                 dbdt(j,ipath) = planckg_wave(ispace,vwave,
     1					emtemp(j,ipath))
                 bsurf=planck_wave(ispace,vwave,Tsurf)
              ENDIF

              tlayer = dexp(-dble(taus(j)))
              tr = trold*tlayer

              call e3interp(tr,e3new,e3grad)

              corkout(ipath,ig) = corkout(ipath,ig) + 
     1         xfac*2.0*pi*sngl((e3old - e3new))*bb(J,Ipath)

              DO k1=1,nsw
                k = isw(k1)
                DO j1=1,j
                  IF(j1.LT.j)THEN
                    dtrdq(j1,k) = dtolddq(j1,k)*sngl(tlayer)
                  ELSE
                    dtrdq(j1,k) = -dtaudq(j1,k)*sngl(tr)
                  ENDIF
                  de3dq(j1,k)=sngl(e3grad)*dtrdq(j1,k)

                  dcoutdq(ipath,ig,j1,k) =
     1                 dcoutdq(ipath,ig,j1,k) + bb(j,ipath)*
     2                     xfac*2.0*pi*(de3olddq(j1,k) - de3dq(j1,k))
                ENDDO
                IF(k.EQ.ngas+1)THEN
                 dcoutdq(ipath,ig,j,k) =
     1      dcoutdq(ipath,ig,j,k) + xfac*2.0*pi*sngl((e3old-e3new))*
     2           dbdt(j,ipath)
               ENDIF
              ENDDO

              e3old = e3new
              trold = tr

              DO j1=1,j
                DO k1=1,nsw
                  k = isw(k1)
                  dtolddq(j1,k) = dtrdq(j1,k)
                  de3olddq(j1,k) = de3dq(j1,k)
                ENDDO
              ENDDO

            ENDDO

C           We need to add on any radiation from the surface
C *** Need to add in code to compute the gradient also and report back.

            corkout(ipath,ig) = corkout(ipath,ig) +
     1        xfac*2.0*pi*sngl(e3old)*(esurf*bsurf)

C           gradient wrt surface temperature
            tempgtsurf(ipath,ig) = xfac*2.0*pi*sngl(e3old)*esurf*
     1        planckg_wave(ispace,vwave,Tsurf)

C            print*,'XXX',vwave,Tsurf,
C     1       planckg_wave(ispace,vwave,Tsurf),gradtsurf(ipath)

          ELSE
            WRITE(*,*)' CIRSRADG.f :: Error: model not defined.'
            WRITE(*,*)' Stopping program.'
            WRITE(*,*)' '
            WRITE(*,*)' Imod = ',imod(ipath),' is not a currently' 
            WRITE(*,*)' valid CIRSRADG option.'
            STOP
          ENDIF
        ENDDO
      ENDDO

C=======================================================================
C
C	Now integrate over g-ordinates and then end loop over bins.
C
C=======================================================================

      DO ipath=1,npath
        output(ipath) = 0.0
        gtsurf(ipath) = 0.0
        DO j=1,maxlay
          DO k=1,nparam
            doutputdq(ipath,j,k) = 0.0
          ENDDO
        ENDDO

C       Find correct number of atmospheric layers for IMOD=8 
C           (i.e. SCR radiometer).
        if(imod(ipath).eq.8)then
           nlay1=nlayin(1)
        else
           nlay1=nlayin(ipath)
        endif

        DO ig=1,ng
          output(ipath) = output(ipath) + corkout(ipath,ig)*delg(ig)
          gtsurf(ipath) = gtsurf(ipath) + 
     1			tempgtsurf(ipath,ig)*delg(ig)

          DO j=1,nlay1
            DO k1=1,nsw
              k=isw(k1)
              doutputdq(ipath,j,k) = doutputdq(ipath,j,k) +
     1        dcoutdq(ipath,ig,j,k)*delg(ig)
            ENDDO
          ENDDO


C        if(ig.eq.5)then
C        print*,'CCALC',ipath,output(ipath),gtsurf(ipath)
C        do k1=1,nsw
C          k=isw(k1)
C          do j=1,nlay1
C           print*,'CGRAD',k,j,doutputdq(ipath,j,k)
C          enddo
C        enddo
C        print*,ipath,output(ipath)
C        endif

        
        ENDDO


      ENDDO

      RETURN

      END
C***********************************************************************
C***********************************************************************
