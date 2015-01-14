      SUBROUTINE GENRADS(NLAYER,NPATH,AMOUNT,PP,PRESS,TEMP,
     1 DELH,NLAYIN,LAYINC,SCALE,EMTEMP,NGAS,IDGAS,ISOGAS,INORMAL,
     2 IRAY,IPTF,IPROC,VMIN,DELV,NPOINT,FWHM,ICONV,WINGIP,VRELIP,MAXDV,
     3 IMOD,INLTE,NFILT,FILTER,VFILT,NCONT,CONT,NSEC,VSEC,
     4 XSEC,PRESSKTA,TEMPKTA,NPK,NTK,IRECK,GABSC,DELG,NG,
     5 DOP,ERRLIM,FLAGH2P,HFP,OUTPUT,ERROR)
C-------------------------------------------------------------------------
C     $Id: genrads.f,v 1.22 2011-09-06 15:33:16 irwin Exp $
C-------------------------------------------------------------------------
C_TITLE: GENRAD: to compute multiple path spectra at "infinite" or finite
C        resolution for several gases.
C
C_ARGS:  NLAYER:INTEGER  number of layers (including cell layers) to consider
C        NPATH:INTEGER    number of individual paths through these layers
C	 DELH(MAXLAY):REAL Height of each layer (km)
C        AMOUNT(MAXLAY,NGAS):REAL   number of molecules/cm2 for each layer
C                         and gas (normal incidence)
C        PP(MAXLAY,NGAS):REAL   partial pressure of each gas (atm)
C        PRESS(NLAYER):REAL     total pressure for each layer (atm)
C        TEMP(NLAYER):REAL      temperature (Kelvin) of each layer
C        NLAYIN(NPATH):INTEGER    each layer defined above can be included
C                         in the calculation for each path several times. This
C                         is the number of layers included in the calculation
C                         for each path (if you include it twice count
C                         it twice)
C        LAYINC(MAXINC,NPATH):INTEGER    the layer numbers in order of
C                         inclusion for each path. The order also defines the
C                         processing of each layer depending upon the model
C                         ( e.g. the 1st layer may be a gas cell)
C        SCALE(MAXINC,NPATH):REAL   scaling factor (eg cosine of incidence)
C                         for each included layer in each of the paths
C                         note: when a layer is included more than once
C                         different scaling factors can be defined.
C        EMTEMP(MAXINC,NPATH):REAL temperatures to use in emission calculations
C        NGAS:INTEGER     number of gases to consider
C        IDGAS(NGAS):INTEGER   the LOCAL gas identifier for each gas
C                         The local gas identifier agrees with HITRAN id's
C                         as far as possible (Jan 1988) with extra id's
C                         for gases not in the HITRAN compilation. eg.
C                         those in the GEISA compilation
C        ISOGAS(NGAS):INTEGER   the local isotopic identifier, if zero all
C                         isotopes of the gas are included.
C                         Isotope id's also agree as far as possible with
C                         HITRAN id's. Similar (1,2,3...n) id's have been
C                         defined for additional gases.
C                         If zero then line strengths are used as tabulated
C                         (i.e. corrected for normal terrestrial distribution)
C                         If a specific isotope is selected then its strength
C                         is scaled to the pure isotope using the HITRAN
C                         data base values
C        IPROC(NGAS):INTEGER   type of absorption coefficient calculation
C                         to use for each gas
C                          0 = default Voigt broadening
C        VMIN:REAL        minimum wavenumber for output
C        DELV:REAL        wavenumber spacing for output
C        FWHM:REAL	  resolution
C        NPOINT:INTEGER   the number of output points
C                         i.e. each wavenumber, v=vmin+(i-1)*delv for
C                         i=1 to npoint
C        ICONV:INTEGER    convolution required for output
C                         0 = none
C                         1 = square with width delv CENTERED on each
C                         wavenumber
C        WINGIP:REAL      lines more than wing wavenumbers from current
C                         integration wavenumber may be treated as continuum
C        VRELIP:REAL      lines more than VREL away are ignored
C                         note: both WING and VREL can be modified by the
C                         program but the dummy input parameters WINGIP and
C                         VRELIP are unchanged
C        IMOD(NPATH):INTEGER   type of model to use to compute OUTPUT from
C                         transmissions for each path
C                         0 = transmission
C                         1 = absorption (gives improved accuracy
C                             for low absorption)
C                         2 = emission
C                         3 = emission but uses source function at centre
C                             of current bin (note TROLD is set to total
C                             transmission in emission calculations)
C                         4 = the total transmission computed during the
C                             last emission calculation (TROLD)
C                         5 = Single cell transmission
C                         6 = PMR SB transmission ( transmission layer 1 -
C                             transmission layer 2) only 2 layers
C                         7 = the last 'cell' output path multiplied
C                             by the last 'atmospheric' output path
C                         8 = the output is the product of two previous
C                             output paths. The two path numbers are stored
C                             in the LAYINC array instead of layer numbers
C                         9 = the output is the product of a single layer
C                             transmission path and the output of the previous
C                             path (this is useful for speeding up calculation
C                             of transmission weighting functions)
C                        10 = as 2 but assumes transmissions calculated by
C                             single curtis-godson path
C                        11 = as 10 but uses source function at centre
C                             of bin
C			 12 = PMR WB transmission ( transmission layer 1 +
C			      transmission layer 2)/2
C			 13 = SCR SB
C			 14 = SCR WB
C			 15 = Scattering
C			 16 = Single scattering
C			 17 = Thermal emission into hemisphere
C			 18 = Thermal emission into hemisphere but using
C			     source function at centre of current bin
C			 19 = As 17 but using CG
C			 20 = As 18 but using CG
C			 21 = Limb scattering
C
C        FILTER(NFILT):REAL   filter transmissions
C        VFILT(NFILT):REAL   wavenumbers for each filter transmission
C        NFILT:INTEGER    number of filter points provided (if 1 ignores
C                         filter) note: must be non zero
C        NCONT:INTEGER    number of dust types included
C        CONT(MAXCON,NLAYER):REAL   the number of dust particles/cm2 through
C			  each layer (vertical path).
C        NSEC:INTEGER     number of wavenumbers for which the dust cross-
C			  sections are defined.
C	 XSEC(2,MAXCON,MAXSEC):REAL	dust cross sections for each dust type
C	 VSEC(MAXSEC):REAL	corresponding wavenumbers for dust x-sections
C        DOP(NLAYER)      doppler shift for each layer. DOP is effectively
C                         ADDED to each line position
C        ERRLIM(NPATH):REAL   fractional integration error acceptable
C                         for each path
C        OUTPUT(NPATH,NPOINT):REAL   output values at each wavenumber for
C                         each output type for each path
C        ERROR(NPATH,NPOINT):REAL     the estimated error in the calculated
C                         output (not used for infinite resolution spectra)
C
C        Additional values may be required by some models. e.g. a
C        source function for each layer for emmission calculations
C        These are either passed via common blocks or read in from some
C        input channel. The parameters are described in the relevant input
C        sections as they are added to the program
C
C_KEYS:  SUBR,ATMO,SPEC,VMS
C
C_DESCR: General purpose transmission calculation routine in three parts:
C
C	 1) LBL
C        The spectral region is divided into bins WING wavenumbers wide.
C        Line data is stored in blocks corresponding to the WING sized
C        bins. Each monocromatic calculation first checks that the appropriate
C        bins are in memory then swaps bins in and out if necessary.
C        For a monochromatic calculation lines from the bin at that wavenumber
C        and the two adjacent bins are treated explicitly.
C        Absorption by lines from other bins is calculated from
C        precomputed continuum polynomials based upon fits to the line
C        wings.
C        Note - the integration steps and the bin sizes are not related
C        so an integration may cover many bins or just one.
C        Values can be output as an infinite resolution spectrum
C        or averaged using an adaptive Simpson's rule integration.
C        In averaging over the instrument function each integration is
C        performed from the last upper integration limit to the next line
C        centre or over DVMAX wavenumbers whichever is smaller. This ensures
C        that complex regions have a large number of small integrations
C        while regions with few lines are treated seperately with larger
C        integration steps.
C        The program computes monochromatic transmission for each layer (for
C        several paths through the layers) then these transmissions are
C        used to compute the output value which may be just the total
C        transmission or the total emmission etc. The output value computed
C        is determined by the parameter IMOD
C        note: The layers input are rearranged into paths depending upon
C        the array LAYINC. The original order of layer input is not
C        significant.
C
C	 2) Band
C
C	 3) Correlated-K
C
C
C_FILES  LUN - line data base files
C        none unless used by a particular OUTPUT model. See code for
C        descriptions of model dependent parameters and input
C
C_CALLS:
C
C_BUGS:
C
C_HIST:  19jun86 SBC ORIGINAL VERSION
C         3feb88 SBC modified to include new line data base routines and
C                replacing NAG routines with numerical recipies routines.
C        12may89 SBC adding 'D lines' diagnostic code to check line binning
C        15mar90 SBC restructured, revised integration scheme and added
C                explicit error output
C        20sep91 SBC modified to use revised line data routines
C        18nov92 SBC added EMTEMP as a parameter, changed BB calulations to use
C                this
C        23nov92 SBC chaned function of models 10 and 11 to conform to new
C                paths computed by limb.f
C        15feb93 PGJI debugged dust continuum calculation
C        26feb93 PGJI add option to choose CO2 line shape
C        26feb93 PGJI add stimulated emission term
C	 24oct94 PGJI modified to combine RT models.
C        10feb97 CAN  modified to call new CIACON instead of HYDCON,
C                to incorporate the new CIA gas codes for Titan
C	18/3/04	NT added iproc=6 option, Rosenkrantz-Ben-Reuven lineshape
C			with voigt lineshape as default (NH3 only)
C------------------------------------------------------------------------------
C     Variables:
C------------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      
C     Input variables
      INTEGER NLAYER,NPATH,NLAYIN(MAXPAT)
      INTEGER LAYINC(MAXINC,MAXPAT),NGAS,IDGAS(MAXGAS),ISOGAS(MAXGAS)
      INTEGER IPROC(MAXGAS),NPOINT,ICONV,IMOD(MAXPAT),NFILT,NCONT
      INTEGER INLTE,IRAY,IPTF
      REAL AMOUNT(MAXLAY,MAXGAS),PP(MAXLAY,MAXGAS),PRESS(MAXLAY),NU
      REAL HFP(MAXLAY),LINECONTRIB,FH2,FNH3
      INTEGER JH2,JNH3
      REAL TEMP(MAXLAY),SCALE(MAXINC,MAXPAT),EMTEMP(MAXINC,MAXPAT)
      REAL VMIN,DELV,FWHM,WINGIP,VRELIP,FILTER(NFILT),VFILT(NFILT)
      REAL CONT(MAXCON,MAXLAY),DOP(MAXLAY),ERRLIM(MAXPAT),DELH(MAXLAY)
      REAL OUTPUT(NPATH,NPOINT),ERROR(NPATH,NPOINT)
      REAL MAXDV,DV1,E3,E3OLD,E3NEW
      REAL SOL_ANG,EMISS_ANG,APHI
      INTEGER NF
      CHARACTER*100 RADFILE1
C------------------------------------------------------------------------------

      INCLUDE '../includes/parcom.f'
C------------------------------------------------------------------------------
C     Adaptive integration variables
      REAL XNEW(MAXPAT),XOLD(MAXPAT),ERR(MAXPAT),EMAX,ELMAX
      REAL SUMEND(MAXPAT),SUMEVEN(MAXPAT),SUMODD(MAXPAT)
      INTEGER ISIMP,FLAGH2P,KK
      LOGICAL LCALCH,KMETH,SKIP,ADDBAND
C------------------------------------------------------------------------------
C     Include line record variables:
      INCLUDE '../includes/lincom.f'
C------------------------------------------------------------------------------
C     Include spectral bin variables
      INCLUDE '../includes/bincom.f'
C------------------------------------------------------------------------------
C     Continuum variables
      INTEGER ISUM
C     IORDER is the order of the continuum polynomial
      REAL CONTINK(IORDP1,MAXLAY,MAXBIN)
      REAL CONTIN(IORDP1,MAXLAY,MAXBIN),CONTMP(IORDP1)
      REAL CONVAL(NWAV),CONWAV(NWAV),FILCON(IORDP1,MAXBIN)
      REAL MATRIX(IORDP1,IORDP1),UNIT(IORDP1,IORDP1)
      DOUBLE PRECISION DMATRIX(IORDP1,IORDP1),DUNIT(IORDP1,IORDP1)
      REAL CONSCA(MAXCON,IORDP1,MAXLAY,MAXBIN)
C     CONTIN holds continuum polynomial for each bin, CONTMP is a
C      temporary store for these.
C     CONVAL holds the absorption coefficients at wavenumbers CONWAV prior to
C      fitting continuum.
C     FILCON holds the continuum polynomial for the filter profile.
C     MATRIX and UNIT are both used for matrix inversion when computing
C     polynomials where insufficient points for least squares.
      COMMON /CONCOM/CONTIN,CONTINK,CONTMP,CONVAL,CONWAV,FILCON,
     &MATRIX,UNIT
C------------------------------------------------------------------------------
C     Dust variables
      INTEGER NSEC
      REAL VSEC(MAXSEC),XSEC(2,MAXCON,MAXSEC),TSEC(MAXSEC)
C------------------------------------------------------------------------------
C     General variables
      INTEGER BBOFF,MAXBBBIN
      PARAMETER (MAXBBBIN = 750000)
      REAL BB(MAXINC,MAXPAT),BBBIN(MAXBBBIN)
      REAL TAUTMP(MAXLAY),FRAC(MAXLAY,MAXGAS),TAUSCA(MAXCON,MAXLAY)
      REAL TAUSCAT(MAXLAY)
      INTEGER LSTATM,LSTCEL

C     BB is the (black body) source function for each layer
C     BBBIN is the source function precomputed at each BIN centre
C     TAUTMP holding the normal incidence case optical depth of each layer
C     FRAC holds volume fraction of each gas in each layer
C     LSTATM holds the path number of the last "atmospheric" path calculated
C     LSTCEL the path number of the last cell path calculated
C------------------------------------------------------------------------------
C     Misc variables or variables defined in the code
      INTEGER I,ii,J,K,L,LINE
      INTEGER LAYER,IPATH
      INTEGER NSAMP,NCALC
      REAL DOPMAX,WING,VREL,VMAX,VTOP,VBOT,RANGE,XMASS,FILCUR
      REAL V1,VTMP,TAU,TR,TROLD,GETMASS
      DOUBLE PRECISION V,V1DP,V2DP,DVSAMP,VCENT
      REAL ABSCO,LNABSCO,AD,X,Y,MINPR,MXMASS,MXWID,XWID,FPARA
      REAL PARTF,AAMOUNT(MAXGAS),PPP(MAXGAS),XLEN
      INTEGER LABEL
      LOGICAL NOFILT
      REAL M, TRAT
C------------------------------------------------------------------------------
C     Include line data base variables
      INCLUDE '../includes/dbcom.f'
C------------------------------------------------------------------------------
C     Partition function variables
      REAL TCORS1(MAXLAY,MAXGAS),TCORS2(MAXLAY),TRATIO(MAXLAY)
      REAL TCORDW(MAXLAY,MAXGAS),NLTE(MAXLAY),CALC_NLTE,PLANCK
      REAL TSTIM,TS1,TS2,TPART(4)
C     TCORS1 is the partition function temperature dependence times absorber
C     amount
C     TCORS2 is the temperature dependence for the Boltzman distribution
C     TCORDW is the temperature dependence of the doppler width
C     TRATIO is just 296/T for use in calculating T dependence of line width
C     TSTIM is the correction for stimulated emission
C------------------------------------------------------------------------------
C     Spectral model identifiers
      LOGICAL LINEDATA,LBLCALC,BANDCALC,CORKCALC
C------------------------------------------------------------------------------
C     Variables for Overlapping lines and Random Band model
C     Rodgers numerical data
      DOUBLE PRECISION RODGERS(4,0:7)
C     Look up table of voigt equivalent widths
      REAL VWIDTH(117,161)
      COMMON /VOIGTDATA/VWIDTH
C     Codes
      INTEGER WID_MOD,CALC_MOD,COM_MOD,LOR_MOD,DOP_MOD
      REAL EQW_LORENTZ,EQW_DOPPLER,EQW,SUMEQW
      REAL MALKMUS_LOR,GOODY_LOR,GODSON_LOR
C------------------------------------------------------------------------------
C     Parband variables
      REAL ALCORR,TAU_GOODY_VOIGT1,TAU_MG_LOR,TAU_EKS,TGV3
      LOGICAL MG,VOIGT
      REAL KNU,KNU0,DELAD,Y0,YV,EL,E1,E2,AA,SFB,P,T,U,Q,QROT
C------------------------------------------------------------------------------
C     Correlated-K variables
      INTEGER NG,IGPATH,ID1
      PARAMETER (IGPATH=200)
      REAL K_G(MAXG),KL_G(MAXG,MAXLAY),G_ORD(MAXG),UTOT,DELG(MAXG)
      REAL PI,UTOTL(MAXLAY),K_G1(MAXG),K_G2(MAXG)
      PARAMETER (PI=3.1415927)
      REAL CORKOUT(IGPATH,MAXG),QACT(MAXLAY),GABSC(MAXG)
      DOUBLE PRECISION SL,BL,SS,TAUD,BB1,PAR(MAXGAS,3)
      REAL PRESSKTA(MAXK),TEMPKTA(MAXK),PIN,TIN,VIN
      INTEGER IRECK(10),NPK,NTK,ULOG
      PARAMETER (ULOG=17)
C------------------------------------------------------------------------------
C     GASCON switches
      INCLUDE '../includes/gascom.f'
C
C     Scattering variables
C
      DOUBLE PRECISION mu1(maxmu), wt1(maxmu),galb
      real rad, radg(maxmu), eps(maxscatlay), bnu(maxscatlay),
     1 taus(maxscatlay),flux, dff, solar, lfrac(maxcon,maxscatlay),
     2 theta0,dist,tsun, lcons(maxcon,10),lfrac1(maxcon,maxscatlay)
      parameter(theta0=4.65241E-03,tsun=5900.0)
      integer lncons(maxcon),lnorm(maxcon),nmu,liscat(maxcon),
     1 lowbc,nlays
      real press1(maxscatlay),temp1(maxscatlay),utot1(maxscatlay)
      character*30 scatfile(maxcon)
      logical scatter

      DOUBLE PRECISION PPLSTO(0:25,MAXMU,MAXMU,maxscatlay)
      DOUBLE PRECISION PMISTO(0:25,MAXMU,MAXMU,maxscatlay)
C

      COMMON/PHASESTO/PPLSTO,PMISTO

      REAL P0,T0
      PARAMETER(P0 = 1.0, T0 = 296.0)

      character*512 buffer
      real cthet(50),thet(50),phas(10,3,50),head(10,3,3)
      real dv
      integer irec1(10),nphas,ncons,maxrec(10)

      common /phase/thet,cthet,phas,head,irec1,nphas,maxrec

C added orth-para Hydrogen ratio: bimodal (but in future could be continuosly
C     variable)
C
C 0=equilibrium H2 (ortho:para = 2:1)
C 1='normal' H2 (orth:para = 3:1)

      integer INORMAL,isol,iunit,npoint1,irec,igas,igdist,ng1
      integer miter,itmp,ioff,ntot,kl,jl,iodd,jbin,i1,lk
      integer j1,lun0,irec0,il,nphi
      real v0,vw1,xpw,xpd,vbmin,vbmax,del,vv,width,delwav,sums
      real sumw,dpexp,al,calc_lor_width,calc_dop_width
      real trantmp,xcorr,c1,c2,x1,vmin1,wing1,delv1,qu,sst,bbt
      real vrel1,sumk,sumy,q1,q2,vc

      REAL BANDQ(42)
      DATA BANDQ/1.5,1.0,1.5,1.0,1.0,1.5,1.0,1.0,1.5,1.5,
     &           1.5,1.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.5,
     &		 1.5,1.0,1.0,1.5,1.5,1.0,1.5,1.5,1.0,1.5,
     &		 1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.0,1.0,
     &		 1.5,1.5/

C------------------------------------------------------------------------------
C
C
C
C------------------------------------------------------------------------------
C     CODE
C------------------------------------------------------------------------------
C
C
C
C------------------------------------------------------------------------------
C     checking input parameters
      IF(NLAYER.GT.MAXLAY)THEN
        WRITE(*,1001)
1001    FORMAT(' %too many layers for GENRAD')
        print*,'MAXLAY, NLAYER = ',MAXLAY,NLAYER
        STOP
      END IF
      IF(NPATH.GT.MAXPAT)THEN
        WRITE(*,1002)
1002    FORMAT(' %too many paths for GENRAD')
        print*,'MAXPAT,NPATH = ',MAXPAT,NPATH
        STOP
      END IF
      IF(NGAS.GT.MAXGAS)THEN
        WRITE(*,1003)
1003    FORMAT(' %too many gases for GENRAD')
        print*,'MAXGAS,NGAS = ',MAXGAS,NGAS
        STOP
      END IF

      DO 53 IPATH=1,NPATH
       IF(NLAYIN(IPATH).LE.MAXINC)GOTO 79
       WRITE(*,78)
78     FORMAT(' %number of layers defined is inconsistent with arrays')
       print*,'IPATH,NLAYIN(IPATH),MAXINC = ',IPATH,NLAYIN(IPATH),MAXINC
       STOP
79     CONTINUE
       IF(IMOD(IPATH).EQ.8)GOTO 53
       DO 55 LAYER=1,NLAYIN(IPATH)
        IF(LAYINC(LAYER,IPATH).LE.NLAYER)GOTO 55
        WRITE(*,73)
73     FORMAT(' %attempt to include an undefined layer in path')
       print*,'LAYINC(LAYER,IPATH),NLAYER = ',LAYINC(LAYER,IPATH),NLAYER
        STOP
55     CONTINUE
53    CONTINUE
C------------------------------------------------------------------------------
C     Defining some useful variables
      LINMAX=MAXLIN
C-----------------------------------------------------------------------------
C     Reading in scattering information
      scatter=.false.
      do ipath=1,npath
       if(imod(ipath).eq.15.or.imod(ipath).eq.21)then
        scatter=.true.
       end if
       if(imod(ipath).eq.16)scatter=.true.
      end do

      if(scatter)then

       radfile1 = radfile

       call read_scatter1(radfile,nmu,mu1,wt1,isol,dist,lowbc,
     1  galb,ncons,liscat,lnorm,lncons,lcons,scatfile,sol_ang,
     2  emiss_ang,aphi,nf)

       do i=1,ncont

        if(liscat(i).eq.4)then
         iunit=10+i
         open(unit=iunit,file=scatfile(i),status='old',access='direct',
     1    recl=512,form='formatted')
         read(iunit,1,rec=1)buffer
1        format(a512)
         if(buffer(2:2).eq.'w')then
          read(buffer(12:512),*)v0,v1,dv,npoint1,nphas
         else
          read(buffer,*)v0,v1,dv,npoint1,nphas
         end if

         maxrec(i)=npoint1+3
         irec1(i)=3

         if(i.eq.1)then
          irec=3
          read(iunit,1,rec=irec)buffer
          read(buffer,*)(thet(j),j=1,nphas)
          do j=1,nphas
           cthet(j)=cos(thet(j)*0.0174532)
          end do
         end if
        end if

       end do
      end if

C-----------------------------------------------------------------------------
C
C
C
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

C     Identify which gas is NH3 and H2 (if applicable) so that the
C     weird new NH3 lineshape can be used.
      JH2 = -1
      JNH3 = -1       
      DO IGAS=1,NGAS
       IF(IDGAS(IGAS).EQ.39.AND.
     1  (ISOGAS(IGAS).EQ.0.OR.ISOGAS(IGAS).EQ.1))THEN
        JH2 = IGAS
       ENDIF
       IF(IDGAS(IGAS).EQ.11.AND. 
     1  (ISOGAS(IGAS).EQ.0.OR.ISOGAS(IGAS).EQ.1))THEN
        JNH3 = IGAS
       ENDIF
      ENDDO

C
C     Setting up Spectral model identifiers
C
C-----------------------------------------------------------------------------

      WRITE(ULOG,*)' '

      IF(ICONV.LT.10)THEN
       LBLCALC=.TRUE.
       IGDIST = 1
       WRITE(ULOG,*)'Line by line model'
      ELSE
       LBLCALC=.FALSE.
      ENDIF

      IF(ICONV.GE.10.AND.ICONV.LT.20)THEN
       BANDCALC=.TRUE.
       WRITE(ULOG,*)'Band model'

       IF(ICONV.EQ.10.OR.ICONV.EQ.11)THEN
        CALL READ_RODGERS(RODGERS)
        CALL READ_VOIGT
        WID_MOD = INT(WINGIP)
        CALC_MOD = INT(VRELIP)

        COM_MOD=INT(WID_MOD/9)
        WID_MOD=WID_MOD - 9*COM_MOD
        LOR_MOD=INT(WID_MOD/3)
        DOP_MOD=INT(WID_MOD - 3*LOR_MOD)
       END IF

       IF(ICONV.EQ.10)THEN
        WRITE(ULOG,*)'Non-overlapping lines'
        WRITE(ULOG,*)'COM_MOD = ',COM_MOD
        WRITE(ULOG,*)'LOR_MOD = ',LOR_MOD
        WRITE(ULOG,*)'DOP_MOD = ',DOP_MOD

       ELSE IF(ICONV.EQ.11)THEN
        IF(CALC_MOD.LT.2)THEN
         WRITE(ULOG,*)'General random band model'
        ELSE IF(CALC_MOD.EQ.2)THEN
         WRITE(ULOG,*)'Malkmus-Lorentz model'
        ELSE IF(CALC_MOD.EQ.3)THEN
         WRITE(ULOG,*)'Goody-Lorentz model'
        ELSE IF(CALC_MOD.EQ.4)THEN
         WRITE(ULOG,*)'Godson-Lorentz model'
        END IF
        WRITE(ULOG,*)'COM_MOD = ',COM_MOD
        WRITE(ULOG,*)'LOR_MOD = ',LOR_MOD
        WRITE(ULOG,*)'DOP_MOD = ',DOP_MOD

       ELSE IF(ICONV.EQ.12)THEN
C       Malkmus Lorentz
        WRITE(ULOG,*)'Malkmus-Lorentz model'
        VOIGT=.FALSE.
        MG=.TRUE.
        ALCORR=0.25
       ELSE IF(ICONV.EQ.13)THEN
C       Goody-Lorentz
        WRITE(ULOG,*)'Goody-Lorentz model'
        VOIGT=.FALSE.
        MG=.FALSE.
        ALCORR=1.
       ELSE IF(ICONV.EQ.14)THEN
C       Godson-Lorentz
        WRITE(ULOG,*)'Godson-Lorentz model'
        PRINT*,'Godson-Lorentz model not ready'
        STOP
       ELSE IF(ICONV.EQ.15)THEN
C       Goody-Voigt
        WRITE(ULOG,*)'Goody-Voigt model'
        VOIGT=.TRUE.
       ELSE
        PRINT*,'Band model undefined'
        print*,'ICONV = ',ICONV
        STOP
       END IF
      ELSE
       BANDCALC=.FALSE.
      END IF


      IF(ICONV.GE.20)THEN
       WRITE(ULOG,*)'Correlated-K model'
       CORKCALC=.TRUE.
       IF(ICONV.NE.24)THEN
C       Read in G-intervals
        WRITE(ULOG,*)'Reading in g - intervals from',DGFILE
        OPEN(15,FILE=DGFILE,STATUS='OLD')
        READ(15,*)NG
        DO 234 I=1,NG
         READ(15,*)GABSC(I),DELG(I)
234     CONTINUE
        G_ORD(1)=0
        DO 324 I=1,NG
         G_ORD(I+1)=G_ORD(I) + DELG(I)
324     CONTINUE
        CLOSE(15)
       ELSE
        WRITE(ULOG,*)'K-table look up'
        WRITE(ULOG,*)'Number of pressure points = ',NPK
        DO I=1,NPK
         WRITE(ULOG,*)PRESSKTA(I)
        END DO
        WRITE(ULOG,*)'Number of temperature points = ',NTK
        DO I=1,NTK
         WRITE(ULOG,*)TEMPKTA(I)
        END DO
        G_ORD(1)=0.
        DO I=1,NG
         G_ORD(I+1)=G_ORD(I)+DELG(I)
        END DO
       END IF

       WRITE(ULOG,*)'NG = ',NG
       NG1=NG+1
       WRITE(ULOG,*)'G_ORD(I)    GABSC(I)       DELG(I)'
       DO I=1,NG
        WRITE(ULOG,*)G_ORD(I),GABSC(I),DELG(I)
       END DO
       WRITE(ULOG,*)G_ORD(NG+1)


       IF(ICONV.EQ.21)THEN
        CALC_MOD = INT(WINGIP)
        WRITE(ULOG,*)'Malkmus-Lorentz corr-k CALC_MOD = ',CALC_MOD
       END IF
       IF(ICONV.EQ.22)THEN
        CALC_MOD = INT(WINGIP)
        LCALCH=.FALSE.
        WRITE(ULOG,*)'Goody-Voigt corr-k'
        IF(CALC_MOD.EQ.1)LCALCH=.TRUE.
        WRITE(ULOG,*)'CALC_MOD, LCALCH= ',CALC_MOD,LCALCH
       END IF
       IF(ICONV.EQ.23)THEN
        WRITE(ULOG,*)'Exponential-sum'
        CALC_MOD = INT(WINGIP)
        WRITE(ULOG,*)'CALC_MOD = ',CALC_MOD
        MITER=100
       END IF
       IF(ICONV.EQ.24)THEN
        WRITE(ULOG,*)'Correlated-K Lookup table'
       END IF
      ELSE
       CORKCALC=.FALSE.
      END IF

      IF(ICONV.LT.12.OR.ICONV.EQ.20)THEN
       WRITE(ULOG,*)'Line data base keyfile : ',KEYFIL
       LINEDATA=.TRUE.
      ELSE
       IF(ICONV.NE.24)THEN
        WRITE(ULOG,*)'Band data base keyfile : ',KEYFIL
       ELSE
        WRITE(ULOG,*)'Corrk Look-up table keyfile : ',KEYFIL
       END IF
       LINEDATA=.FALSE.
      END IF


      IF(IRAY.GT.0) THEN
       PRINT*,'Setting Rayleigh scattering odepths on'
      ELSE
       PRINT*,'No Rayleigh scattering optical depths applied' 
      ENDIF

      PRINT*,'INORMAL,FLAGH2P,IRAY = ',INORMAL,FLAGH2P,IRAY

CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

   
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
C     precalculating layer properties
C
C------------------------------------------------------------------------------
C     computing volume fractions
      WRITE(ULOG,*)'Number of gases = ',NGAS
      WRITE(ULOG,*)'Number of layers = ',NLAYER
      WRITE(ULOG,*)'NLTE Code (INLTE) = ',INLTE
      DO 251 I=1,NLAYER
       QACT(I)=0.
       UTOTL(I)=0.
       DO 241 J=1,NGAS
        FRAC(I,J)=PP(I,J)/PRESS(I)
        QACT(I)=QACT(I)+FRAC(I,J)
        UTOTL(I)=UTOTL(I)+AMOUNT(I,J)*1.E-20
241    CONTINUE
       WRITE(ULOG,*)'Active Fraction for layer ',I,' = ',QACT(I)
251   CONTINUE

      MXMASS = -20

      DO 16 I=1,NGAS
C     Check isotopes included in model for each gas and set mass for doppler
C     width calculation

       XMASS=GETMASS(IDGAS(I),ISOGAS(I))

       WRITE(ULOG,420)IDGAS(I),ISOGAS(I),GASNAM(IDGAS(I)),XMASS
420    FORMAT(1X,'gas:',I2,'/',I2,2X,1A6,' mass=',F7.2)

       MXMASS = MAX(MXMASS,XMASS)

       DO 16 J=1,NLAYER
C      Calculate temperature tempendance parameters for the gas lines
C      (Note: TCORS1 includes factor of 1.E-47 for scaling of stored line
C       strengths. Scaling is applied in two stages to avoid numerical
C       overflow)

        TCORS1(J,I)=PARTF(IDGAS(I),ISOGAS(I),TEMP(J),IPTF)*
     1   (AMOUNT(J,I)*1.E-20)
        TCORS1(J,I)=1.E-27*TCORS1(J,I)
        TCORDW(J,I)=4.301E-7*SQRT(TEMP(J)/XMASS)
16    CONTINUE

      VW1 = VMIN + 0.5*(NPOINT-1)*DELV

      IF(MXMASS.LE.0)THEN
       PRINT*,'Warning. No Gases present. Setting MXMASS to 100'
       MXMASS=100
      END IF
      XPW = 0.1*PRESS(1)*SQRT(296.0/TEMP(1))
      XPD = 4.301E-7*VW1*SQRT(TEMP(1)/MXMASS)
      MXWID = SQRT(XPD**2 + XPW**2)


      DO 17 J=1,NLAYER
       TRATIO(J)=296./TEMP(J)
       TCORS2(J)=1.439*(TEMP(J)-296.)/(296.*TEMP(J))
C      Include any empirical function to account for nonLTE source function
       NLTE(J)=CALC_NLTE(PRESS(J),INLTE)
       WRITE(ULOG,*)'LAYER, NLTE Correction : ',J,NLTE(J)
       MINPR=MIN(MINPR,PRESS(J))
       XPW= 0.1*PRESS(J)*SQRT(TRATIO(J))
       XPD = 4.301E-7*VW1*SQRT(TEMP(J)/MXMASS)
       XWID = SQRT(XPD**2 + XPW**2)
       MXWID = MIN(MXWID,XWID)
17    CONTINUE


CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL


      IF(LBLCALC.OR.ICONV.EQ.20)THEN
       IF(WINGIP.EQ.0)THEN
         WINGIP=MXWID*3.0
         WINGIP=0.01*INT(100*WINGIP)
         WINGIP=MAX(WINGIP,1.0)
         PRINT*,'Default setting WING to ',WINGIP
       END IF
       IF(VRELIP.EQ.0)THEN
         VRELIP = 20*WINGIP
         VRELIP=0.01*INT(100*VRELIP)
         PRINT*,'Default setting VREL to ',VRELIP
       ENDIF
      END IF


CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
C
C     Initialising pointers and bins for calculations
C
C------------------------------------------------------------------------------
C
C--------------------------------------------------
      IF(LBLCALC)THEN
       WING=WINGIP
C      Increase WING by maximum doppler shift so that even for shifted
C      layers, enough lines are read in for accurate calculation
       DOPMAX=0.
       DO 171 LAYER=1,NLAYER
        DOPMAX=MAX(ABS(DOP(LAYER)),DOPMAX)
171    CONTINUE
       WING=WING+DOPMAX
       VREL=VRELIP

C      VMAX is the maximum wavenumber for output
       VMAX=VMIN+DELV*FLOAT(NPOINT-1)
C      VBOT-VTOP are range for calculation

       VBOT=VMIN-VREL

       IF(VBOT.LT.0.0)THEN
		VBOT=0.0
       ENDIF

       VTOP=VMAX+VREL
       IF(ICONV.EQ.1)THEN
C       allow for finite resolution
        VBOT=VBOT-DELV/2.
        VTOP=VTOP+DELV/2.
	IF (VBOT.LT.0.0) THEN
	 VBOT=VBOT+DELV
	ENDIF
       END IF

       RANGE=VTOP-VBOT
       NBIN=INT(RANGE/WING)+1
       DELWAV=WING


      ELSE
C--------------------------------------------------

C      Set up 'bins' for BAND or CK calculations
       VBOT=VMIN-0.5*FWHM
       VMAX=VMIN+(NPOINT-1)*DELV
       VTOP=VMAX + 0.5*FWHM

       IF(ICONV.EQ.20)THEN
        VBOT=VBOT-VRELIP
        VTOP=VTOP+VRELIP
       END IF
       DELWAV=DELV

       IF(LINEDATA)THEN
        IF(ICONV.EQ.10.OR.ICONV.EQ.11)THEN
          WING=FWHM
        ELSE
          WING=WINGIP
        END IF
       END IF


       IF(LINEDATA)THEN
        RANGE=VTOP-VBOT
        NBIN=INT(RANGE/WING)+1
       ELSE
        NBIN=0
       END IF


      END IF
C--------------------------------------------------


      IF(NBIN.GT.MAXBIN)THEN
        WRITE(*,*)'NBIN'
        WRITE(*,*)'MAXBIN = ',MAXBIN
        WRITE(*,11)
11      FORMAT(' %too many bins - must recompile GENRAD')
        STOP
      END IF

      IF(LINEDATA)THEN
       DO 127 I=1,NBIN
        VBIN(I)=VBOT+(I-1)*WING
127     CONTINUE
      END IF

      IF(LBLCALC)THEN
        DO 137 I=1,NBIN
          VBINB(I)=VBIN(I)
137      CONTINUE
        NCALC=NBIN
        DV=0.5*WING
       ELSE
        DO 147 I=1,NPOINT
          VBINB(I)=VBOT+(I-1)*DELV
          IF(ICONV.EQ.20)VBINB(I)=VBINB(I)+VRELIP
147      CONTINUE
        NCALC=NPOINT
        DV=0.5*FWHM
       ENDIF


       IF(NCALC.GT.MAXBIN)THEN
        PRINT*,'NCALC > MAXBIN',NCALC,MAXBIN
        STOP
       END IF

       DO 157 I=1,NCALC
        V1=VBINB(I)+DV
        DO 85 IPATH=1,NPATH
        ITMP = IMOD(IPATH)

        IF(ITMP.EQ.3.OR.ITMP.EQ.11.OR.ITMP.EQ.18.OR.ITMP.EQ.20)THEN
         DO 96 LAYER=1,NLAYIN(IPATH)
          IOFF = BBOFF(NPATH,NLAYIN,NCALC,IPATH,LAYER,I,MAXBBBIN)
          IF(EMTEMP(LAYER,IPATH).GT.0.0)THEN
           BBBIN(IOFF)=NLTE(LAYINC(LAYER,IPATH))*
     1       PLANCK(V1,EMTEMP(LAYER,IPATH))
          ELSE
           BBBIN(IOFF) = 0.0
          ENDIF
96       CONTINUE
        ENDIF
85      CONTINUE

157   CONTINUE

      NTOT=MAX(NCALC,NBIN)


C------------------------------------------------------------------------------
C     initialising continuum polynomial
      DO 14 I=1,NCALC
      DO 14 K=1,NLAYER
      DO 14 J=1,IORDP1
       CONTIN(J,K,I)=0.
      DO 14 KL=1,NCONT
       CONSCA(KL,J,K,I)=0
14    CONTINUE

CBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB


      IF(LINEDATA)THEN
        CALL LOADBINS(WING,NGAS,IDGAS,ISOGAS)
      END IF

CFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
C
C     Now computing filter profile as polynomial by fitting to input profile
C     for each bin
C         filter(v)=filcon(0)+filcon(1).dv+filcon(2).dv**2......
C     where dv=v-vbin and vbin is the lower wavenumber limit
C     Filter profile is defined by NFILT points which have no relation
C     to bin size.
C
C------------------------------------------------------------------------------
      NOFILT=.TRUE.
      IF(NFILT.GE.2)THEN
       NOFILT=.FALSE.
       IF(VFILT(1).GT.VBOT.OR.VFILT(NFILT).LT.VTOP)THEN
          WRITE(*,115)
115       FORMAT(' %filter profile doesn"t cover computation range')
       STOP
       END IF
       WRITE(*,111)
111    FORMAT(' %fitting to filter profile')

       DO 112 I=1,NCALC
C      calculate filter polynomial for each bin
       IF(LBLCALC)THEN
         VBMIN=VBIN(I)
         VBMAX=VBIN(I)+WING
       ELSE
         VBMIN=VBINB(I)
         VBMAX=VBINB(I)+FWHM
       END IF

       CALL CALC_PCOEFF(NFILT,FILTER,VFILT,VBMIN,VBMAX,CONTMP)

       DO 122 L=1,IORDP1
        FILCON(L,I)=CONTMP(L)
122    CONTINUE
112   CONTINUE
      END IF

CFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Now step through each bin computing continuum polynomial of line wings
C     from other bins if LBL calculation. Only need to compute for central
C     bins if monochromatic calc.
C
C     For other models which use direct line data, read in the lines for all
C     the bins.
C
C------------------------------------------------------------------------------


C------------------------------------------------------------------------------
      FSTBIN=1
      LSTBIN=NBIN

      IF(LBLCALC.OR.ICONV.EQ.20)THEN
      IF(ICONV.EQ.0)THEN
       FSTBIN=INT((VMIN-VBOT)/WING)
       LSTBIN=INT((VMAX-VBOT)/WING)+2
      ELSE IF(ICONV.EQ.1)THEN
       FSTBIN=INT((VMIN - 0.5*DELV - VBOT)/WING)
       LSTBIN=INT((VMAX + 0.5*DELV - VBOT)/WING)+2
      ELSE IF(ICONV.EQ.20)THEN
       FSTBIN=INT((VMIN-0.5*FWHM - VBOT)/WING)
       LSTBIN=INT((VMAX+0.5*FWHM - VBOT)/WING)+2
      END IF
      WRITE(*,71)NBIN
      WRITE(ULOG,71)NBIN
71    FORMAT(' %computing continuum polynomial from ',I4,' bins')

      DEL = VBIN(2)-VBIN(1)

      DO 444 K=1,IORDP1
C     Wavenumbers to compute continuum at
C     Note: polynomials are in wavenumbers from start of bin
       CONWAV(K)=FLOAT(K-1)*WING/FLOAT(IORDER)
444   CONTINUE
C     setting up matrix of wavenumbers for use in polynomial calculation
      DO 19 K=1,IORDP1
      MATRIX(1,K)=1.
      DMATRIX(1,K)=1.
      DO 19 J=2,IORDP1
      MATRIX(J,K)=MATRIX(J-1,K)*CONWAV(K)
      DMATRIX(J,K)=DMATRIX(J-1,K)*CONWAV(K)
19    CONTINUE
 
      L=IORDP1
      CALL DMATINV(DMATRIX,L,L,DUNIT)

      DO K=1,IORDP1
       DO J=1,IORDP1
        UNIT(J,K)=SNGL(DUNIT(J,K))
       ENDDO
      ENDDO

      DO 13 I=1,NBIN
C       computing continuum for all except adjacent bins
        PRINT*,'Continuum bin ',I,' of ',NBIN
        print*,'FSTLIN,LSTLIN',FSTLIN(I),LSTLIN(I)
        DO 15 J=FSTBIN,LSTBIN
          IF(ABS(I-J).LE.1)GOTO 15
C         for each layer
          DO 18 LAYER=1,NLAYER
C           for each line
            DO 22 K=1,IORDP1
              CONVAL(K)=0.
22          CONTINUE

            DO 20 LINE=FSTLIN(I),LSTLIN(I)
C            compute absorption coefficient for normal incidence

             DO 21 K=1,IORDP1
              VV=CONWAV(K)+VBIN(J)
              DV1=(VV-VLIN(LINE))

C             *********** CIRS Line processing Parameters ******************
	      IF(ABS(DV1).LE.MAXDV)THEN
C	        Ignore lines more than MAXDV widths away

               IF(IDGAS(IDLIN(LINE)).EQ.1.AND.IH2O.GT.0.AND.
     &           ABS(DV1).GT.25.0)THEN
C                 Don't calc continuum more than 25cm-1 from H2O lines
C                 if H2O continuum is turned on
C                 print*,'Ignoring water line: ',VLIN(LINE)
               ELSE

                  J1=IDLIN(LINE)
                 
                  FNH3=-1.
                  FH2=-1.
                  IF(JNH3.GT.0.) FNH3 = FRAC(LAYER,JNH3)
                  IF(JH2.GT.0.) FH2 = FRAC(LAYER,JH2)

                  CONVAL(K) = CONVAL(K)+ LINECONTRIB(IPROC(J1),
     1       IDGAS(J1),VV,TCORDW(LAYER,J1),TCORS1(LAYER,J1),
     2       TCORS2(LAYER),PRESS(LAYER),TEMP(LAYER),FRAC(LAYER,J1),
     3       VLIN(LINE),SLIN(LINE),ELIN(LINE),ALIN(LINE),SBLIN(LINE),
     4       TDW(LINE),TDWS(LINE),LLQ(LINE),DOUBV(LINE),FNH3,FH2)

               ENDIF

              ENDIF
C             *********************************************************

21           CONTINUE
20          CONTINUE
            DO 23 K=1,IORDP1
             DO 314 L=1,IORDP1
              CONTIN(K,LAYER,J)=CONTIN(K,LAYER,J)+UNIT(L,K)*CONVAL(L)
314          CONTINUE
23          CONTINUE
18        CONTINUE

15      CONTINUE



13    CONTINUE


      IF(ICONV.EQ.20)THEN
       DO 8881 ISUM=1,IORDP1
        DO 8882 LAYER=1,NLAYER
         DO 8883 J=1,NBIN
          CONTINK(ISUM,LAYER,J)=CONTIN(ISUM,LAYER,J)/UTOTL(LAYER)
          CONTIN(ISUM,LAYER,J)=0
8883     CONTINUE
8882    CONTINUE
8881   CONTINUE
      END IF
C -------------------------------------------
      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
C     Now compute the continuum effects due to dust and add to wing continua
C     as filter but added to contin for each layer
C     note that dust (or equivalent) is defined for each layer by an array
C     of number of dust particles/cm2 for each particle type CONT(M,LAYER)
C     The dust cross sections are held in XSEC(2,M,NSEC) where the
C     corresponding wavenumbers are VSEC(NSEC)
C
C-----------------------------------------------------------------------------
      WRITE(ULOG,*)'Dust extinct. and scat. X-sections (cm2) : '
      WRITE(ULOG,*)'NSEC = ',NSEC
      DO K=1,NSEC
       WRITE(ULOG,*)VSEC(K),(XSEC(1,J,K),J=1,NCONT)
       WRITE(ULOG,*)'              ',(XSEC(2,J,K),J=1,NCONT)
      END DO




      IF(NSEC.LT.2)GOTO 600
      IF(VSEC(1).GT.VBOT.OR.VSEC(NSEC).LT.VTOP)THEN
        WRITE(*,601)
601     FORMAT(' dust x-sections don"t cover computational range')
        PRINT*,'VBOT, VSEC(1)    = ',VBOT,VSEC(1)
        PRINT*,'VTOP, VSEC(NSEC) = ',VTOP,VSEC(NSEC)
        STOP
      END IF

      WRITE(*,602)
602   FORMAT(' %fitting to supplied dust continua')

      DO 603 I=1,NCALC
        IF(LBLCALC)THEN
         VBMIN=VBIN(I)
         VBMAX=VBIN(I)+WING
        ELSE
         VBMIN=VBINB(I)
         VBMAX=VBINB(I)+FWHM
        END IF

        DO 621 J=1,NCONT
         DO 620 LAYER=1,NLAYER
          DO 622 L=1,NSEC
           TSEC(L)=XSEC(1,J,L)*CONT(J,LAYER)
622       CONTINUE
          CALL CALC_PCOEFF(NSEC,TSEC,VSEC,VBMIN,VBMAX,CONTMP)
          DO 617 L=1,IORDP1
           CONTIN(L,LAYER,I)=CONTIN(L,LAYER,I)+CONTMP(L)
617       CONTINUE
620      CONTINUE
621     CONTINUE


        DO 671 J=1,NCONT
         DO 670 LAYER=1,NLAYER
          DO 672 L=1,NSEC
           TSEC(L)=XSEC(2,J,L)*CONT(J,LAYER)
672       CONTINUE
          CALL CALC_PCOEFF(NSEC,TSEC,VSEC,VBMIN,VBMAX,CONTMP)
          DO 677 L=1,IORDP1
           CONSCA(J,L,LAYER,I)=CONTMP(L)
677       CONTINUE
670      CONTINUE
671     CONTINUE

603   CONTINUE
600   CONTINUE

CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD



CGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
C
C     Here compute gaseous continua not covered by wings in LBL or by band
C     data in other codes and add to the continuum polynomials.
C
C     Computed in the same way as filter profile interpolation but in
C     a different subroutine to allow easy modification.
C     GASCON returns a polynomial approximation to the continuum of a
C     particular gas over a particular bin.
C
C-------------------------------------------------------------------------

      DO 312 K=1,NCALC
        IF(LBLCALC)THEN
          VBMIN=VBIN(K)
          WIDTH=WING
        ELSE
          VBMIN=VBINB(K)
          WIDTH=FWHM
        END IF
        DO 310 J=1,NGAS
          DO 310 I=1,NLAYER
            CALL GASCON(VBMIN,WIDTH,IDGAS(J),ISOGAS(J),
     1       AMOUNT(I,J),PP(I,J),PRESS(I),TEMP(I),CONTMP)
            DO 311 L=1,IORDP1
              CONTIN(L,I,K)=CONTIN(L,I,K)+CONTMP(L)
311         CONTINUE
310    CONTINUE


       DO 317 I=1,NLAYER
C reduce the 2d arrays for gas amounts (no./unit vol.)
C and partial pressure to 1d for passing to subroutine


             do 42 ii=1, NGAS
                AAMOUNT(ii) = AMOUNT(I,ii)
                PPP(ii) = PP(I,ii)
 42          continue

             XLEN = DELH(I)
             id1 = 0				! diagnostic print


             IF(FLAGH2P.EQ.1)THEN
  	        FPARA = HFP(I)
                CALL FPARACON(VBMIN,WIDTH,PRESS(I),TEMP(I),
     1            NGAS,IDGAS,ISOGAS,AAMOUNT,PPP,
     2		  FPARA,XLEN,CONTMP,ID1)
             ELSE
                 DO KK=1,IORDP1
                  CONTMP(KK)=0.0
                 ENDDO

                 CALL CIACON(VBMIN,WIDTH,PRESS(I),TEMP(I),
     1            INORMAL,NGAS,IDGAS,ISOGAS,AAMOUNT,PPP,
     2		  XLEN,CONTMP,ID1)

             ENDIF

             DO 313 L=1,IORDP1
                CONTIN(L,I,K)=CONTIN(L,I,K)+CONTMP(L)
 313         CONTINUE

             IF(IRAY.GT.0)THEN
              CALL CONRAY(IRAY,VBMIN,WIDTH,PRESS(I),TEMP(I),
     1		UTOTL(I),CONTMP)

              DO 315 L=1,IORDP1
                CONTIN(L,I,K)=CONTIN(L,I,K)+CONTMP(L)
315           CONTINUE

             ENDIF


 317   CONTINUE

 312  CONTINUE


CGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG



CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C     Now performing spectral calculation
C
C     The calculation may have to be performed in several parts in LBL
C     since all lines cannot be stored in memory simultaneously.
C     The number of steps is calculated assuming that the number of
C     bins rather than the number of lines limits the calculation size. This
C     feature does not apply to other models
C
C-----------------------------------------------------------------------------

      WRITE(*,72)
72    FORMAT(' %now beginning spectral calculation')



CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
C     Now stepping through different models (ICONV) and setting up calculation
C     accordingly
C
C-------------------------------------------------------------------------
C
C     Line by Line Models:


      IF(ICONV.EQ.0)THEN
C      infinite resolution
       DO 103 I=1,NPOINT
        V=DBLE(VMIN+FLOAT(I-1)*DELV)
        VV=SNGL(V)
        ASSIGN 2001 TO LABEL
        GOTO 2000
2001   CONTINUE
103    CONTINUE
      ELSEIF(ICONV.EQ.1) THEN
C-------------------------------------------------------------------------
C      finite resolution LBL
       ELMAX=0.0
       DO 1017 K=1,NPATH
        ELMAX=MAX(ELMAX,ERRLIM(K))	!Find maximum allowable percentage error
1017   CONTINUE
       ELMAX = ELMAX*100

       DO 104 I=1,NPOINT
        VCENT=DBLE(VMIN+(I-1)*DELV)
C       the integration limits for this interval are
        V1DP=VCENT-DBLE(DELV)/2.
        V2DP=VCENT+DBLE(DELV)/2.
        WRITE(*,1105)I,V1DP,V2DP
1105   FORMAT(' interval',I7,'  from',F10.3,' to',F10.3)

        XWID = MAX(0.001,MXWID)
C Multiply by three so as to have the center point, and one at both  VSTART
C (= VMIN - 0.5*FWHM) and VEND (= VSTART + FWHM).
        NSAMP=3*INT(0.5*SNGL(V2DP-V1DP)/XWID)

        IF(NSAMP.LE.4)NSAMP=3
C        NSAMP=20
        DVSAMP=(V2DP-V1DP)/DBLE(NSAMP)
        WRITE(*,1102)NSAMP,DVSAMP

1102    FORMAT(' sub dividing interval into ',I4,' intervals of ',
     1   F9.4,'cm^-1')

C       -----------------------------------------------------------
C       initialise integration variables

        DO 233 K=1,NPATH
          SUMEND(K)=0.
          SUMEVEN(K)=0.
          SUMODD(K)=0.
233     CONTINUE

        DO 105 J=0,NSAMP
         V=V1DP+DBLE(J)*DVSAMP
         VV=SNGL(V)
         ASSIGN 2002 TO LABEL
         GOTO 2000
2002     CONTINUE


         IODD = J - 2*INT(0.5*J)		!Flag for odd, even ordinates

         DO 211 K=1,NPATH
           IF(J.EQ.0.OR.J.EQ.NSAMP)THEN
            SUMEND(K)=SUMEND(K)+OUTPUT(K,I)	!Sum of end values
           ELSE IF(IODD.EQ.1)THEN
            SUMODD(K)=SUMODD(K)+OUTPUT(K,I)	!Sum of odd values
           ELSE
            SUMEVEN(K)=SUMEVEN(K)+OUTPUT(K,I)	!Sum of even values
           ENDIF
211      CONTINUE

105     CONTINUE

C       Calculate averaged spectrum for each path
        DO 222 K=1,NPATH
         XOLD(K) = SNGL((DVSAMP/3.0)*
     &   (SUMEND(K)+4*SUMODD(K)+2*SUMEVEN(K))/(V2DP-V1DP))

222     CONTINUE

C       -----------------------------------------------------------
C       The variables are now initialised for a SIMPSON integration

        ISIMP=1
C       -----------------------------------------------------------
C       iteration loop
107     CONTINUE
        ISIMP=ISIMP+1
        WRITE(*,1103)ISIMP,DVSAMP/2
1103    FORMAT(' Integration "index":',I5,'DVSAMP : ',F9.4)

C       Now halve the sample distance. Doing this means that the new sum
C       of the even points is just the sum of the even AND odd points from
C       the last step
        DO K=1,NPATH
         SUMEVEN(K) = SUMEVEN(K)+SUMODD(K)
         SUMODD(K)=0.
        ENDDO

C       Calculate the new odd points
        DO 108 J=0,NSAMP-1
         V=V1DP+DBLE(J+0.5)*DVSAMP
         VV=SNGL(V)
         ASSIGN 2004 TO LABEL
         GOTO 2000
2004     CONTINUE


         DO K=1,NPATH
          SUMODD(K)=SUMODD(K)+OUTPUT(K,I)
         END DO

108     CONTINUE

        DVSAMP=DVSAMP*0.5	!halve the sample distance

C       Calculate new averaged outputs

        DO K=1,NPATH
         XNEW(K) = SNGL((DVSAMP/3.0)*
     &    (SUMEND(K)+4*SUMODD(K)+2*SUMEVEN(K))/(V2DP-V1DP))
         IF(XOLD(K).NE.0)THEN
          ERR(K) = 100*ABS(XNEW(K)-XOLD(K))/XOLD(K)	!% difference
         ELSEIF(XNEW(K).EQ.0)THEN			!from last step
          ERR(K)= 0.0
         ELSE
          ERR(K) = 100.0
         ENDIF
         IF(K.EQ.1)THEN
          EMAX = ERR(1)
         ELSE
          EMAX = MAX(EMAX,ERR(K))
         ENDIF
        END DO

C       If not converged, double the number of points and do again.
        IF(EMAX.GT.ELMAX)THEN
         NSAMP=NSAMP*2
         IF(NSAMP.GT.MAXINT)THEN
          PRINT*,'NSAMP > MAXINT',NSAMP,MAXINT
	  PRINT*,'Stopping now.'
          GOTO 277
         ENDIF
         DO K=1,NPATH
          XOLD(K)=XNEW(K)
         END DO
         GOTO 107
        END IF

277     DO K=1,NPATH
         OUTPUT(K,I)=XNEW(K)
         ERROR(K,I)=XNEW(K)-XOLD(K)
        END DO

104     CONTINUE
      END IF

C-------------------------------------------------------------------------
C     End of ICONV=0,1 If Statement of GENLBL

C-------------------------------------------------------------------------

C     Band models

      IF(BANDCALC)THEN
C      BAND MODELS
       DO 333 I=1,NPOINT
        V=DBLE(VMIN+FLOAT(I-1)*DELV)
        VV=SNGL(V)
         PRINT*,'BANDCALC: V = ',SNGL(V)
         ASSIGN 2009 TO LABEL
         GOTO 2000
2009    CONTINUE
333    CONTINUE
      END IF

C-------------------------------------------------------------------------
C
C     Correlated-K models
C
C-------------------------------------------------------------------------

      IF(CORKCALC)THEN

       DO 431 I=1,NPOINT

       V=DBLE(VMIN+FLOAT(I-1)*DELV)
       VV=SNGL(V)

       PRINT*,'CORKCALC: V = ',SNGL(V)
       DO 433 IGDIST=1,NG
        ASSIGN 2010 TO LABEL
        GOTO 2000
2010    CONTINUE
        DO 534 IPATH=1,NPATH
         CORKOUT(IPATH,IGDIST)=OUTPUT(IPATH,I)
534     CONTINUE
433    CONTINUE
       DO 434 IPATH=1,NPATH
        OUTPUT(IPATH,I)=CORKOUT(IPATH,1)*DELG(1)
        DO 435 IGDIST=2,NG
         OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+CORKOUT(IPATH,IGDIST)
     1*DELG(IGDIST)
435     CONTINUE
434    CONTINUE

431    CONTINUE
      END IF

CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM



CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      RETURN
C=============================================================================
C
C     This section calculates the transmission etc at wavenumber V.
C     The output values are stored in the array OUTPUT(,I)
C     This would be better as a subroutine, but there are so many parameters
C     to pass that it would take forever so entry is via a GOTO and exit
C     via assigned GOTO
2000  CONTINUE

CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
C     First compute normal path optical depth through each layer (TAUTMP)
C
C----------------------------------------------------------------------------
      DO 50 LAYER=1,NLAYER
C      Compute absorption coefficient for LAYER
C      First calculate continuum contribution. CURBIN is the current BIN
C      for which calculations are being made.
       IF(LBLCALC)THEN
        CURBIN=INT((SNGL(V-DBLE(VBOT)))/WING)+1
       ELSE
        CURBIN=I
       END IF

C      Calculate the continuum absorption via the IORDP1 polynomial
C      coefficients held in CONTIN

       TAUTMP(LAYER)=CONTIN(1,LAYER,CURBIN)
       VTMP=SNGL(V-DBLE(VBINB(CURBIN)))-DOP(LAYER)
       DO 51 ISUM=2,IORDP1
        TAUTMP(LAYER)=TAUTMP(LAYER)+CONTIN(ISUM,LAYER,CURBIN)*VTMP
        VTMP=VTMP*VTMP
51     CONTINUE

       IF(TAUTMP(LAYER).LT.0.0)THEN
        PRINT*,'Error in continuum TAUTMP. TAU < 0.0!'
        PRINT*,'LAYER = ',LAYER
        PRINT*,'V = ',V
        PRINT*,'I,CURBIN,TAUTMP = ',I,CURBIN,TAUTMP(LAYER)
        PRINT*,'Resetting to zero'
        TAUTMP(LAYER)=0.0
       ENDIF

       TAUSCAT(LAYER)=0.
       DO 793 JL=1,NCONT
       TAUSCA(JL,LAYER)=CONSCA(JL,1,LAYER,CURBIN)
       VTMP=SNGL(V-DBLE(VBIN(CURBIN)))-DOP(LAYER)
       DO 571 ISUM=2,IORDP1
        TAUSCA(JL,LAYER)=TAUSCA(JL,LAYER)+
     1    CONSCA(JL,ISUM,LAYER,CURBIN)*VTMP
        VTMP=VTMP*VTMP
571     CONTINUE
        TAUSCAT(LAYER)=TAUSCAT(LAYER)+TAUSCA(JL,LAYER)
793    CONTINUE

       IF(SCATTER)THEN

        IF(NCONT.GT.maxcon)THEN
         PRINT*,'Error in GENRADS. NCONT > maxcon'
         PRINT*,NCONT,maxcon
         STOP
        ENDIF
        IF(LAYER.GT.maxscatlay)THEN
         PRINT*,'Error in GENRADS. LAYER > maxscatlay'
         PRINT*,LAYER,maxscatlay
         STOP
        ENDIF

        IF(TAUSCAT(LAYER).EQ.0.)THEN
         DO JL=1,NCONT
          LFRAC(JL,LAYER)=1.0/REAL(MAX(1,NCONT))
         END DO
        ELSE
         DO JL=1,NCONT
          LFRAC(JL,LAYER)=TAUSCA(JL,LAYER)/TAUSCAT(LAYER)
         END DO
        END IF

       ENDIF

CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
C      Step though and select the spectral model required
C
C****************************************************************************
       IF(LBLCALC)THEN

       IF(CURBIN.LT.2)CURBIN=2
       DO 507 JBIN=CURBIN-1, CURBIN+1


        DO 52 LINE=FSTLIN(JBIN),LSTLIN(JBIN)
C        compute absorption coefficient for normal incidence

         DV1=SNGL(DBLE(V-VLIN(LINE)))-DOP(LAYER)


C        *********** Line processing Parameters ******************
	 IF(ABS(DV1).LT.MAXDV)THEN
C	  Ignore lines more than MAXDV away

          IF(IDGAS(IDLIN(LINE)).EQ.1.AND.IH2O.GT.0.AND.
     &      ABS(DV1).GT.25.0)THEN
C          Don't calc continuum more than 25cm-1 from H2O lines if
C          IH2O is turned on.

          ELSE

            J1=IDLIN(LINE)

            FNH3=-1.
            FH2=-1.
            IF(JNH3.GT.0.) FNH3 = FRAC(LAYER,JNH3)
            IF(JH2.GT.0.) FH2 = FRAC(LAYER,JH2)

            TAUTMP(LAYER) = TAUTMP(LAYER)+ LINECONTRIB(IPROC(J1),
     1       IDGAS(J1),VV,TCORDW(LAYER,J1),TCORS1(LAYER,J1),
     2       TCORS2(LAYER),PRESS(LAYER),TEMP(LAYER),FRAC(LAYER,J1),
     3       VLIN(LINE),SLIN(LINE),ELIN(LINE),ALIN(LINE),SBLIN(LINE),
     4       TDW(LINE),TDWS(LINE),LLQ(LINE),DOUBV(LINE),FNH3,FH2)

          ENDIF
         ENDIF
C        *********************************************************

52      CONTINUE


        IF(TAUTMP(LAYER).LT.0.0)THEN
         PRINT*,'Error in LBL TAUTMP. TAU < 0.0!'
         PRINT*,'LAYER = ',LAYER
         PRINT*,'V = ',V
         PRINT*,'I,CURBIN = ',I,CURBIN
         PRINT*,'TAUTMP = ',TAUTMP(LAYER)
        ENDIF


507     CONTINUE

C****************************************************************************
       ELSE IF(ICONV.EQ.10.OR.ICONV.EQ.11)THEN
C       BAND: overlapping lines, general random band models


        SUMEQW=0.
        SUMS=0.
        SUMW=0.
        DELWAV=FWHM

        DO 502 LINE=FSTLIN(CURBIN),LSTLIN(CURBIN)
C         Compute stimulated emission term
          IF(VLIN(LINE).LT.10000)THEN
           TS1 = 1.0-DPEXP(-1.439*VLIN(LINE)/TEMP(LAYER))
           TS2 = 1.0-DPEXP(-1.439*VLIN(LINE)/296.0)
           TSTIM=1.0
           IF(TS2.NE.0) TSTIM = TS1/TS2
          ELSE
           TSTIM=1.
          END IF

          LNABSCO=LOG(SLIN(LINE))+LOG(TCORS1(LAYER,IDLIN(LINE)))+
     &        (TCORS2(LAYER)*ELIN(LINE))+LOG(TSTIM)
          ABSCO=EXP(LNABSCO)

          AD=TCORDW(LAYER,IDLIN(LINE))*VLIN(LINE)

          M = FRAC(LAYER, IDLIN(LINE))
          TRAT = TRATIO(LAYER)
          AL = (ALIN(LINE) * (1 - M) * TRAT ** TDW(LINE) +
     $         (ALIN(LINE) - SBLIN(LINE)) * M * TRAT ** TDWS(LINE)) *
     $         PRESS(LAYER)


          IF(CALC_MOD.GE.2)THEN

            SUMW=SUMW+ABSCO
            SUMS=SUMS+SQRT(ABSCO*AL)

          ELSE

            IF(COM_MOD.EQ.1.AND.ABSCO.GT.0.)THEN
C           **** VOIGT WIDTH BY APPROX. USING AL AND AD *******

             EQW_LORENTZ = CALC_LOR_WIDTH(ABSCO,AL,LOR_MOD,RODGERS)
             EQW_DOPPLER = CALC_DOP_WIDTH(ABSCO,AD,DOP_MOD,RODGERS)
             EQW = EQW_LORENTZ**2 + EQW_DOPPLER**2
             EQW = EQW - (EQW_LORENTZ*EQW_DOPPLER/ABSCO)**2
             EQW = SQRT(MAX(EQW,0.0))

            ELSEIF(COM_MOD.EQ.2.AND.ABSCO.GT.0)THEN
C           **** VOIGT WIDTH USING LOOK-UP TABLE *************

             CALL VOIGT_INTERP(ABSCO,AL,AD,EQW,RODGERS)
            ELSE
C           **** LORENTZ WIDTH ONLY **************************
             EQW = CALC_LOR_WIDTH(ABSCO,AL,LOR_MOD,RODGERS)

            END IF

            SUMEQW=SUMEQW + EQW

          END IF

502     CONTINUE

        IF(CALC_MOD.GE.2)THEN
          SUMW=SUMW/DELWAV
          SUMS=SUMS*2./DELWAV
          IF(CALC_MOD.EQ.2)THEN
           SUMEQW = MALKMUS_LOR(SUMS,SUMW,DELWAV)
          ELSEIF(CALC_MOD.EQ.3)THEN
           SUMEQW = GOODY_LOR(SUMS,SUMW,DELWAV)
          ELSEIF(CALC_MOD.EQ.4)THEN
           SUMEQW = GODSON_LOR(SUMS,SUMW,DELWAV)
          END IF
        END IF


        IF(ICONV.EQ.10)THEN
C       ************** NON OVERLAPPING LINES MODEL **************************
          TRANTMP = 1. - SUMEQW/DELWAV
          IF(TRANTMP.LT.0.)THEN
           PRINT*,'Approximation error, TRANTMP < 0. at layer',LAYER
           PRINT*,TRANTMP
           PRINT*,'Change .drv to a band model driver'
           STOP
          END IF

          TAUTMP(LAYER) = TAUTMP(LAYER)-LOG(TRANTMP)
C         Remember, TAUTMP(LAYER) is the OPTICAL DEPTH of the layer.

        ELSE
C       ****************** RANDOM BAND MODELS  ******************************
          TAUTMP(LAYER)=TAUTMP(LAYER)+SUMEQW/DELWAV
        ENDIF







C****************************************************************************
       ELSE IF(ICONV.GE.12.AND.ICONV.LT.20)THEN
C      Band models using Band data

        DO 56 J=1,NGAS
         P = PRESS(LAYER)
         T = TEMP(LAYER)
         U = AMOUNT(LAYER,J)*1E-20
         Q = FRAC(LAYER,J)
         QROT = BANDQ(IDGAS(J))
         XCORR = 0.0
         IF(BANDTYP(J).LT.2)THEN
          KNU0 = BANDPAR(I,J,1)
          DELAD = BANDPAR(I,J,2)
          Y0 = BANDPAR(I,J,3)
          EL = BANDPAR(I,J,4)
          SFB = BANDPAR(I,J,5)
          IF(BANDTYP(J).EQ.1)THEN
           C1 = ABS(BANDPAR(I,J,6)*1E-4)
           C2 = BANDPAR(I,J,7)             
           XCORR = C1*U*(P/P0)*(T/T0)**C2
          ENDIF
         ELSE
          KNU0 = BANDPAR(I,J,1)
          DELAD = BANDPAR(I,J,2)
          Y0 = BANDPAR(I,J,3)
          E1 = BANDPAR(I,J,4)
          SFB = BANDPAR(I,J,5)
          E2 = BANDPAR(I,J,6)
          AA = BANDPAR(I,J,7)
         ENDIF
C	 turn off methane band models above 9550cm-1 (to avoid
C	 overlap with Karkoschka continuum data
         KMETH = .FALSE.
         IF(IDGAS(J).EQ.6.AND.(ISOGAS(J).EQ.0.OR.ISOGAS(J).EQ.1))THEN
          KMETH=.TRUE.
         ENDIF
C        IF gas is methane then KMETH=.TRUE.
C        If the gas continuum is included in gascon, then ICH4=1
         ADDBAND=.TRUE.
         IF(V.GE.9550.0.AND.KMETH.AND.ICH4.EQ.1)THEN
            ADDBAND=.FALSE.
         ENDIF

         IF(ADDBAND)THEN
          IF(VOIGT)THEN
           IF(BANDTYP(J).LT.2)THEN
            X1 = TAU_GOODY_VOIGT1(KNU0,DELAD,Y0,EL,SFB,QROT,P,T,U,Q)
           ELSE
            IF(BANDTYP(J).EQ.2)THEN
             I1 = 1
            ELSE
             I1 = 2
             IF(ISOGAS(J).EQ.0)THEN
              K=0
             ELSE
              K=ISOGAS(J)
             END IF
             TPART(1)=DBQTA(K,IDGAS(J))
             TPART(2)=DBQTB(K,IDGAS(J))
             TPART(3)=DBQTC(K,IDGAS(J))
             TPART(4)=DBQTD(K,IDGAS(J))
            ENDIF

            X1 = TGV3(I1,QROT,TPART,KNU0,DELAD,Y0,AA,E1,E2,SFB,P,T,
     1                  U,Q)
           ENDIF
            
           TAUTMP(LAYER) = TAUTMP(LAYER) + X1 + XCORR
             
          ELSE

           IF(DELAD.NE.0.)THEN
            YV = ALCORR*Y0/DELAD
            TAUTMP(LAYER) = TAUTMP(LAYER) + TAU_MG_LOR(MG,KNU0,YV,EL,
     1      SFB,QROT,P,T,U,Q)
           END IF

          END IF

         ENDIF

56      CONTINUE







C****************************************************************************
       ELSE IF(ICONV.EQ.20)THEN
C       CORRK: K-calculated from linedata
        IF(IGDIST.EQ.1)THEN
         VMIN1=VBINB(I)
         WING1=WINGIP
         VREL1=VRELIP

         P=PRESS(LAYER)
         T=TEMP(LAYER)
         XPW= 0.1*P
         XPD = 4.301E-7*VMIN1*SQRT(T/MXMASS)
         XWID = SQRT(XPD**2 + XPW**2)
         DELV1 = XWID/5.
         IF(DELV1.LT.0.001)DELV1 = 0.001
         NPOINT1 = 1 + INT(FWHM/DELV1)
         IF(NPOINT1.GT.1)THEN
          DELV1 = FWHM/(1.0*(NPOINT1-1))
         ENDIF

         DO LK=1,NGAS
          PPP(LK)=PP(LAYER,LK)
         END DO


         CALL LBL_KDISTS(PPP,P,T,NGAS,IDGAS,ISOGAS,IPROC,VMIN1,DELV1,
     1 NPOINT1,GABSC,DELG,K_G,NG,LAYER,IPTF)


         DO LK=1,NG
          KL_G(LK,LAYER)=K_G(LK)
         END DO

        END IF

        TAUTMP(LAYER)=TAUTMP(LAYER)+UTOTL(LAYER)*KL_G(IGDIST,LAYER)

C****************************************************************************
       ELSE IF(ICONV.EQ.21)THEN
C      Correlated-K using  - Malkmus-Lorentz Approximation

        IF(IGDIST.EQ.1)THEN
C       only need to calculate k-distribution once!

        SS=0.
        BB1=0.
        QU=QACT(LAYER)

        DO 516 J=1,NGAS
         QROT=BANDQ(IDGAS(J))
         KNU0 = BANDPAR(I,J,1)
         DELAD = BANDPAR(I,J,2)
         Y0 = BANDPAR(I,J,3)
         EL = BANDPAR(I,J,4)
         SFB = BANDPAR(I,J,5)
         P = PRESS(LAYER)
         T = TEMP(LAYER)
         U = AMOUNT(LAYER,J)*1.E-20
         Q = FRAC(LAYER,J)

C        Extract Lorentz line width parameter

         IF(DELAD.NE.0)THEN
          YV = 0.25*Y0/DELAD

          IF (KNU0.GT.0.AND.Q.GT.0.)THEN
           CALL ML_LACIS(KNU0,YV,EL,SFB,QROT,P,T,Q,SL,BL)
           PAR(J,1)=SL
           PAR(J,2)=BL
           PAR(J,3)=AMOUNT(LAYER,J)*1.E-20
          ELSE
           PAR(J,1)=0.
           PAR(J,2)=0.
           PAR(J,3)=0.
          END IF
          SS=SS + SL*Q
          BB1=BB1 + SQRT(SL*BL*Q)
         END IF
516     CONTINUE

        SST=SNGL(SS)
        BBT=SNGL(BB1)

        IF(CALC_MOD.LT.3)THEN

         IF(QU.GT.0.AND.SS.GT.0.0)THEN
          SS=SS/QU
          BB1=BB1/SQRT(QU)
          BB1=BB1*BB1/SS
         END IF

        ELSE

         TAUD=0.
         SS=0.

         DO 653 J=1,NGAS
           SL=PAR(J,1)
           BL=PAR(J,2)
           U =SNGL(PAR(J,3))
           IF(SL.GT.0.AND.BL.GT.0)THEN
            SS=SS+SL*U
            TAUD=TAUD+0.5*PI*BL*(DSQRT(1. + 4.*SL*U/(PI*BL)) - 1.)
           END IF
653      CONTINUE

         IF(TAUD.EQ.0.)THEN
          WRITE(*,*)'Warning, TAUD=0. setting S,B to be as PAR1=1'

          IF(QU.GT.0.AND.SST.GT.0.0)THEN
           SS=SST/QU
           BB1=BBT/SQRT(QU)
           BB1=BB1*BB1/SS
          END IF

         ELSE
         IF(SS.EQ.0.OR.(SS*UTOTL(LAYER)).LE.TAUD)THEN
          WRITE(*,*)'Warning, SS=0. or SS*UTOTL(LAYER)<TAUD'
          WRITE(*,*)'Setting S,B to be as PAR1=1'
          IF(QU.GT.0.AND.SST.GT.0.0)THEN
           SS=SST/QU
           BB1=BBT/SQRT(QU)
           BB1=BB1*BB1/SS
          ELSE
           BB1=0.
           SS=0.
          END IF
         ELSE
           SS=SS/UTOTL(LAYER)
           BB1=(TAUD*TAUD/PI)/(SS*UTOTL(LAYER) - TAUD)
         END IF
        END IF

        END IF


C       SS and BB1 now calculated for the layer.
C       now calculate k(g)
        IF(CALC_MOD.EQ.1)THEN
         CALL CALC_MLK1_K(SS,BB1,G_ORD,K_G,NG1)
        ELSE
         UTOT=UTOTL(LAYER)
         CALL CALC_MLK2_K(SS,BB1,UTOT,G_ORD,K_G,NG1)
        END IF

        DO J=1,NG
         KL_G(J,LAYER)=K_G(J)
        END DO

        END IF

        TAUTMP(LAYER)=TAUTMP(LAYER)+KL_G(IGDIST,LAYER)*UTOTL(LAYER)


C****************************************************************************
       ELSE IF(ICONV.EQ.22)THEN
C	Correlated-K using Band data: Goody-Voigt Model

        IF(IGDIST.EQ.1)THEN

        SUMK=0.
        SUMY=0.
        QU=QACT(LAYER)

        DO 416 J=1,NGAS
         QROT=BANDQ(IDGAS(J))
         KNU0 = BANDPAR(I,J,1)
         DELAD = BANDPAR(I,J,2)
         Y0 = BANDPAR(I,J,3)
         EL = BANDPAR(I,J,4)
         SFB = BANDPAR(I,J,5)
         P = PRESS(LAYER)
         T = TEMP(LAYER)
         U = AMOUNT(LAYER,J)*1.E-20
         Q = FRAC(LAYER,J)


C        calculate y=alpha_L/alpha_D at path conditions
         Y=Y0*P*(17.20465053/T)*(Q + (1-Q)/SFB)

C        Calculate absorption coefficient at path temperature
         KNU=(KNU0*(296./T)**QROT)*dpexp(1.439*EL*(1/296. - 1/T))

         SUMK=SUMK + KNU*Q
         SUMY=SUMY + SQRT(KNU*Y*Q)

C         QU=QU+Q

416     CONTINUE


        IF(QU.GT.0.AND.SUMK.GT.0)THEN
          KNU=SUMK/QU
          SUMY=SUMY/SQRT(QU)
          Y=SUMY*SUMY/KNU
        END IF

        UTOT=UTOTL(LAYER)

        CALL CALC_GVK_K(KNU,DELAD,Y,T,UTOT,LCALCH,G_ORD,K_G,NG1)

        DO J=1,NG
         KL_G(J,LAYER)=K_G(J)
        END DO

        END IF

        TAUTMP(LAYER)=TAUTMP(LAYER)+KL_G(IGDIST,LAYER)*UTOTL(LAYER)

C****************************************************************************
       ELSE IF(ICONV.EQ.23)THEN
C      Exponential sum

        IF(IGDIST.EQ.1)THEN

        SUMK=0.
        SUMY=0.
        QU=QACT(LAYER)
        UTOT=UTOTL(LAYER)

        DO 496 J=1,NGAS
         QROT=BANDQ(IDGAS(J))
         KNU0 = BANDPAR(I,J,1)
         DELAD = BANDPAR(I,J,2)
         Y0 = BANDPAR(I,J,3)
         EL = BANDPAR(I,J,4)
         SFB = BANDPAR(I,J,5)
         P = PRESS(LAYER)
         T = TEMP(LAYER)
         U = AMOUNT(LAYER,J)*1.E-20
         Q = FRAC(LAYER,J)


C        calculate y=alpha_L/alpha_D at path conditions
         Y=Y0*P*SQRT(296.0/T)*(Q + (1-Q)/SFB)

C        Calculate absorption coefficient at path temperature
         KNU=(KNU0*(296./T)**QROT)*dpexp(1.439*EL*(1/296. - 1/T))

         SUMK=SUMK + KNU*Q
         SUMY=SUMY + SQRT(KNU*Y*Q)
C         QU=QU+Q

496     CONTINUE


        IF(QU.GT.0.AND.SUMK.GT.0)THEN
          KNU=SUMK/QU
          SUMY=SUMY/SQRT(QU)
          Y=SUMY*SUMY/KNU
        END IF

        CALL CALC_ESUM2(KNU,DELAD,Y,T,CALC_MOD,G_ORD,DELG,K_G,NG,
     1   MITER)

        END IF

        TAUTMP(LAYER)=TAUTMP(LAYER)+K_G(IGDIST)*UTOTL(LAYER)

C****************************************************************************
       ELSE IF(ICONV.EQ.24) THEN
C      Correlated-K Look up Table

        IF (IGDIST.EQ.1)THEN


         PIN = PRESS(LAYER)
         TIN = TEMP(LAYER)
         VIN = SNGL(V)
         DO 243 J=1,NGAS
          LUN0 = 40 + J -1
          IREC0 = IRECK(J)
          CALL INTERP_KTABLE(PRESSKTA,TEMPKTA,NPK,NTK,NG,LUN0,IREC0,
     1     PIN,TIN,VIN,VMIN,DELV,K_G2)

          IF(J.EQ.1)THEN
           DO K=1,NG
            K_G(K) = K_G2(K)
           END DO
           Q1=FRAC(LAYER,J)
          ELSE
           Q2=FRAC(LAYER,J)
           DO K=1,NG
            K_G1(K)=K_G(K)
           END DO
           CALL OVERLAP(G_ORD,NG1,K_G1,Q1,K_G2,Q2,K_G)
           Q1 = Q1 + Q2
          END IF
243      CONTINUE

         DO J=1,NG
          KL_G(J,LAYER)=K_G(J)
         END DO

        END IF

        TAUTMP(LAYER)=TAUTMP(LAYER) + KL_G(IGDIST,LAYER)*UTOTL(LAYER)

C****************************************************************************
       ELSE
        PRINT*,'ICONV NOT DEFINED'
        STOP
       END IF

50    CONTINUE

CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
CLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL


CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
C     Layer optical depths determined.
C     Now calculate the path transmissions and relevant outputs
C
C----------------------------------------------------------------------------
      DO 61 IPATH=1,NPATH
      IF(IMOD(IPATH).EQ.0)THEN
C       model 0,  transmissions only output
        OUTPUT(IPATH,I)=0.
        DO 62 LAYER=1,NLAYIN(IPATH)
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+TAUTMP(LAYINC(LAYER,IPATH))*
     1  SCALE(LAYER,IPATH)
62      CONTINUE
        OUTPUT(IPATH,I)=DPEXP(-OUTPUT(IPATH,I))
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.1)THEN
C       model 1, absorption
        OUTPUT(IPATH,I)=0.
        DO 76 LAYER=1,NLAYIN(IPATH)
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+TAUTMP(LAYINC(LAYER,IPATH))*
     1  SCALE(LAYER,IPATH)
76      CONTINUE
C        IF(OUTPUT(IPATH,I).GT.0.01)OUTPUT(IPATH,I)=
C     1  1.-DPEXP(-OUTPUT(IPATH,I))
        OUTPUT(IPATH,I)=1.-DPEXP(-OUTPUT(IPATH,I))
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.2)THEN
C       model 2, emission. NOTE planck function recomputed at each wavenumber
C       which is very inefficient
        TAU=0.
        TROLD=1.
        OUTPUT(IPATH,I)=0.
        DO 83 LAYER=1,NLAYIN(IPATH)
        BB(LAYER,IPATH)=NLTE(LAYINC(LAYER,IPATH))*
     1   PLANCK(SNGL(V),EMTEMP(LAYER,IPATH))
        TAU=TAU+TAUTMP(LAYINC(LAYER,IPATH))*SCALE(LAYER,IPATH)
        TR=DPEXP(-TAU)
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+
     1  BB(LAYER,IPATH)*(TROLD-TR)
        TROLD=TR
83      CONTINUE
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.3)THEN
C       model 3, emission as 2 but uses planck function at bin centre
        TAU=0.
        TROLD=1.
        OUTPUT(IPATH,I)=0.
        DO 87 LAYER=1,NLAYIN(IPATH)
        TAU=TAU+TAUTMP(LAYINC(LAYER,IPATH))*SCALE(LAYER,IPATH)
        TR=DPEXP(-TAU)
        IOFF = BBOFF(NPATH,NLAYIN,NCALC,IPATH,LAYER,CURBIN,MAXBBBIN) 
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+
     1  BBBIN(IOFF)*(TROLD-TR)
        TROLD=TR
87      CONTINUE
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.4)THEN
C       model 4, outputs last transmission computed by an emission calculation
C       (to allow surface emission to be added etc)
        OUTPUT(IPATH,I)=TROLD
       ELSE IF(IMOD(IPATH).EQ.5)THEN
C       model 5, mean cell transmission (cell transmission)
C       note: there is no check that only one layer is included but only
C       one is used
        OUTPUT(IPATH,I)=
     1  DPEXP(-TAUTMP(LAYINC(1,IPATH))*SCALE(1,IPATH))
        LSTCEL=IPATH
       ELSE IF(IMOD(IPATH).EQ.6)THEN
C       model 6, PMR sideband transmission
        OUTPUT(IPATH,I)=(DPEXP(-TAUTMP(LAYINC(1,IPATH))*SCALE(1,IPATH))-
     1  DPEXP(-TAUTMP(LAYINC(2,IPATH))*SCALE(2,IPATH))) / 2.0
        LSTCEL=IPATH
       ELSE IF(IMOD(IPATH).EQ.7)THEN
C       model 7, last_cell*last_atmosphere
        OUTPUT(IPATH,I)=OUTPUT(LSTATM,I)*OUTPUT(LSTCEL,I)
       ELSE IF(IMOD(IPATH).EQ.8)THEN
C       model 8, product of two path outputs
        OUTPUT(IPATH,I)=OUTPUT(LAYINC(1,IPATH),I)*
     1  OUTPUT(LAYINC(2,IPATH),I)
       ELSE IF(IMOD(IPATH).EQ.9)THEN
C       model 9,product of one transmission and last output
        OUTPUT(IPATH,I)=OUTPUT(IPATH-1,I)*
     1  DPEXP(-TAUTMP(LAYINC(1,IPATH))*SCALE(1,IPATH))
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.10)THEN
        TROLD=1.
        OUTPUT(IPATH,I)=0.
        DO 383 LAYER=1,NLAYIN(IPATH)
        BB(LAYER,IPATH)=NLTE(LAYINC(LAYER,IPATH))*
     1   PLANCK(SNGL(V),EMTEMP(LAYER,IPATH))
        TAU=TAUTMP(LAYINC(LAYER,IPATH))*SCALE(LAYER,IPATH)
        TR=DPEXP(-TAU)
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+
     1  BB(LAYER,IPATH)*(TROLD-TR)
        TROLD=TR
383      CONTINUE
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.11)THEN
        TROLD=1.
        OUTPUT(IPATH,I)=0.
        DO 387 LAYER=1,NLAYIN(IPATH)
        TAU=TAUTMP(LAYINC(LAYER,IPATH))*SCALE(LAYER,IPATH)
        TR=DPEXP(-TAU)
        IOFF = BBOFF(NPATH,NLAYIN,NCALC,IPATH,LAYER,CURBIN,MAXBBBIN)
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+
     1  BBBIN(IOFF)*(TROLD-TR)
        TROLD=TR
387     CONTINUE
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.12)THEN
C       model 12, PMR wideband transmission
        OUTPUT(IPATH,I)=
     1  0.5*(DPEXP(-TAUTMP(LAYINC(1,IPATH))*SCALE(1,IPATH)) +
     1  DPEXP(-TAUTMP(LAYINC(2,IPATH))*SCALE(2,IPATH)))
        LSTCEL=IPATH
       ELSE IF(IMOD(IPATH).EQ.13)THEN
C       model 13, SCR sideband transmission (1-cell transmission)
        OUTPUT(IPATH,I)=
     1  1. - DPEXP(-TAUTMP(LAYINC(1,IPATH))*SCALE(1,IPATH))
        LSTCEL=IPATH
       ELSE IF(IMOD(IPATH).EQ.14)THEN
C       model 14 SCR wideband transmission
        OUTPUT(IPATH,I)=0.5*(1.0 + DPEXP(-TAUTMP(LAYINC(1,IPATH))*
     1   SCALE(1,IPATH)))
        LSTCEL=IPATH
       ELSE IF(IMOD(IPATH).EQ.15)THEN
C       model 15 Scattering.
        NLAYS=NLAYIN(IPATH)
        OUTPUT(IPATH,I)=0.
        DO 917 LAYER=1,NLAYIN(IPATH)
         DO 9181 IL=1,NCONT
          LFRAC1(IL,LAYER)=LFRAC(IL,LAYINC(LAYER,IPATH))
9181     CONTINUE
         TAUS(LAYER)=TAUTMP(LAYINC(LAYER,IPATH))*SCALE(LAYER,IPATH)
         BNU(LAYER)=NLTE(LAYINC(LAYER,IPATH))*
     1   PLANCK(SNGL(V),EMTEMP(LAYER,IPATH))
         IF(TAUSCAT(LAYINC(LAYER,IPATH)).GT.0.0)THEN
          EPS(LAYER)=1. - TAUSCAT(LAYINC(LAYER,IPATH))/
     1    TAUTMP(LAYINC(LAYER,IPATH))
         ELSE
          EPS(LAYER)=1.0
         END IF
917     CONTINUE

        DO 918 JL=1,NMU
         RADG(JL)=PLANCK(SNGL(V),EMTEMP(NLAYIN(IPATH),IPATH))
918     CONTINUE

        IF(ISOL.EQ.1)THEN
         SOLAR=PLANCK(SNGL(V),TSUN)*PI*(THETA0/DIST)**2
        ELSE
         SOLAR=0.
        END IF

        vc = sngl(v)
        nphi = 100

        call scloud8( rad, sol_ang, emiss_ang, aphi, radg, solar, 
     1   lowbc, galb, mu1, wt1, nmu, nf, igdist, vc, eps, bnu,
     2   taus, nlays, ncont, lfrac1, liscat, lcons,lncons, lnorm, nphi)


        OUTPUT(IPATH,I)=RAD
        LSTATM=IPATH

       ELSE IF(IMOD(IPATH).EQ.16)THEN
C       model 16. single scattering.
        WRITE(*,57)IMOD(IPATH),IPATH
57      FORMAT(' %no current code for model',I3,' Path: ',I3)
        STOP
       ELSE IF(IMOD(IPATH).EQ.17)THEN
C       model 17. Disc-averaged spectral irradiance. i.e. W cm-2 (cm-1)-1. 
C       To convert to radiance (i.e. W cm-2 sr-1 (cm-1)-1), the irradiance 
C       must be divided by pi.
        TAU=0.
        E3OLD=0.5
        OUTPUT(IPATH,I)=0.
        DO 93 LAYER=1,NLAYIN(IPATH)
        BB(LAYER,IPATH)=NLTE(LAYINC(LAYER,IPATH))*
     1   PLANCK(SNGL(V),EMTEMP(LAYER,IPATH))
        TAU=TAU+TAUTMP(LAYINC(LAYER,IPATH))*SCALE(LAYER,IPATH)
        E3NEW=E3(TAU)
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+
     1  2*PI*BB(LAYER,IPATH)*(E3OLD-E3NEW)
        E3OLD=E3NEW
93      CONTINUE
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.18)THEN
C       model 18. hemispherically integrated spectral irradiance. 
C            i.e. W cm-2 (cm-1)-1. 
C       To convert to radiance (i.e. W cm-2 sr-1 (cm-1)-1), the irradiance 
C       must be divided by pi.
C       As model 17 but uses Planck function at centre of bin.
        TAU=0.
        E3OLD=0.5
        OUTPUT(IPATH,I)=0.
        DO 97 LAYER=1,NLAYIN(IPATH)
        TAU=TAU+TAUTMP(LAYINC(LAYER,IPATH))*SCALE(LAYER,IPATH)
        E3NEW=E3(TAU)
        IOFF = BBOFF(NPATH,NLAYIN,NCALC,IPATH,LAYER,CURBIN,MAXBBBIN) 
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+
     1  2*PI*BBBIN(IOFF)*(E3OLD-E3NEW)
        E3OLD=E3NEW
97      CONTINUE
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.19)THEN
C       As model 17 but using Curtis-Godson
        E3OLD=0.5
        OUTPUT(IPATH,I)=0.
        DO 483 LAYER=1,NLAYIN(IPATH)
        BB(LAYER,IPATH)=NLTE(LAYINC(LAYER,IPATH))*
     1   PLANCK(SNGL(V),EMTEMP(LAYER,IPATH))
        TAU=TAUTMP(LAYINC(LAYER,IPATH))*SCALE(LAYER,IPATH)
        E3NEW=E3(TAU)
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+
     1  2*PI*BB(LAYER,IPATH)*(E3OLD-E3NEW)
        E3OLD=E3NEW
483     CONTINUE
        LSTATM=IPATH
       ELSE IF(IMOD(IPATH).EQ.20)THEN
C       As model 18 but using Curtis-Godson
        E3OLD=0.5
        OUTPUT(IPATH,I)=0.
        DO 487 LAYER=1,NLAYIN(IPATH)
        TAU=TAUTMP(LAYINC(LAYER,IPATH))*SCALE(LAYER,IPATH)
        E3NEW=E3(TAU)
        IOFF = BBOFF(NPATH,NLAYIN,NCALC,IPATH,LAYER,CURBIN,MAXBBBIN)
        OUTPUT(IPATH,I)=OUTPUT(IPATH,I)+
     1  2*PI*BBBIN(IOFF)*(E3OLD-E3NEW)
        E3OLD=E3NEW
487     CONTINUE
        LSTATM=IPATH
       ELSE
        WRITE(*,54)IMOD(IPATH),IPATH
54      FORMAT(' %no code for model encountered',I3,' Path: ',I3)
        STOP
       END IF

61    CONTINUE


C     here including filter profile if necessary
      IF(NOFILT)GOTO LABEL
      VTMP=SNGL(V-DBLE(VBIN(CURBIN)))
      FILCUR=FILCON(1,CURBIN)
      DO 82 ISUM=2,IORDP1
      FILCUR=FILCUR+FILCON(ISUM,CURBIN)*VTMP
      VTMP=VTMP*VTMP
82    CONTINUE
      DO 86 IPATH=1,NPATH
86    OUTPUT(IPATH,I)=OUTPUT(IPATH,I)*FILCUR

      GOTO LABEL
      END


CPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
