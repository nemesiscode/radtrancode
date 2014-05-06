C     $Id: laycom.f,v 1.5 2011-06-17 14:47:50 irwin Exp $
C***********************************************************************
C_TITL: LAYCOM
C
C_DESC:	Variables used by path-calculating routines but not by genlbl
C	routines.
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST: 
C***********************************************************************
C-----------------------------------------------------------------------
C
C	Overall control flags and variables. N.B. that this block assumes
C	that pathcom.f has already been included.
C
C-----------------------------------------------------------------------
      LOGICAL INTERV,LINED,COMB,CLRLAY
C INTERV: =.TRUE. if an interval has been defined.
C LINED: =.TRUE. if a linedata file has been defined.
C COMB: =.TRUE. if want to combine all cell and atmosphere paths.
C CLRLAY: =.TRUE. if want to remove unwanted layers

      REAL VMAX,ERRDEF
C VMAX: maximum wavenumber
C ERRDEF: default fractional error limit for integration over wavenumber

      COMMON /OVECOM/ INTERV,LINED,VMAX,ERRDEF,CLRLAY

C-----------------------------------------------------------------------
C
C	Variables to hold atmospheric model.
C
C-----------------------------------------------------------------------
      LOGICAL MODEL
C MODEL: =.TRUE. if an atmospheric model has been read in.

      INTEGER NPRO,NVMR,ID(MAXGAS),ISO(MAXGAS),ATMGAS(MAXGAS),AMFORM
C NPRO: number of points in each profile.
C NVMR: number of volume-mixing-ratios.
C ID: gas identifier for each gas.
C ISO: isotope identifier for each gas. 
C ATMGAS: index of each model gas in the main layer arrays.
C AMFORM: the file format for the model.

      INTEGER IPLANET,NDUST,JFP,JFC
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS)
C H: height [km] above some NOMINAL zero.
C P: pressure [atm (not bar!)].
C T: temperature [Kelvin].
C VMR: volume-mixing-ratio profiles for each of the gases.

      REAL MOLWT,DUST(MAXCON,MAXPRO),DUSTH(MAXPRO),FPH2(MAXPRO)
      REAL FPH2H(MAXPRO),FPH2I(MAXPRO)
      REAL FCLOUD(MAXPRO),FCLOUDH(MAXPRO),FCLOUDI(MAXPRO)
      INTEGER ICLOUD(MAXCON,MAXPRO),ICLOUDI(MAXCON,MAXPRO)
      REAL LATITUDE,RADIUS
C RADIUS: radius of the reference sphere (i.e. height=0), typically 
C the planetary radius).

      COMMON /MODCOM/MODEL,H,P,T,DUST,JFP,FPH2I,JFC,FCLOUDI,ICLOUDI,
     1 VMR,IPLANET,LATITUDE,MOLWT,NPRO,NVMR,ID,ISO,ATMGAS,NDUST,
     2 RADIUS,AMFORM

C-----------------------------------------------------------------------
C
C	Variables to describe dust/continuum model.
C
C-----------------------------------------------------------------------
      INTEGER DMOD
C DMOD: dust model: 0=none, 1=homogenous, 2= Conrath 74 martian dust.

      REAL DNORM,VINF
C DNORM: normalising factor so that optical depth per cm per cm2
C =DNORM*molecule mass density*model.
C VINF is the inflection for the Conrath Martian dust.

      COMMON /DUSCOM/ DMOD,DNORM,VINF

C-----------------------------------------------------------------------
C
C	Variables for layer calculations.
C
C-----------------------------------------------------------------------
      LOGICAL LAYERS,INFATM
C LAYERS: =true if any layers have been calculated for most recent model.
C INFATM: =true if top layer is to represent all molecules above base
C height.

      REAL LAYHT,LAYANG,LAYSF(MAXLAY)
C LAYHT: bottom height for calculating layers (eg tangent height for
C explicit limb paths).
C LAYANG: angle between the zentih and path for layer calculations.
C LAYSF: geometrical scale factors for each layer since the amounts are
C always scaled back to vertical paths.

      INTEGER LAYTYP,NLAY,FSTLAY,LSTLAY
C FSTLAY: first layer in the most recent atmospheric calc.
C LSTLAY: is the last in the most recent atmospheric calc.
C LAYTYP: type of layer splitting [0=dP, 1=dlnP, 2=dH, 3=dS].
C NLAY: number of layers to split the atmosphere into.

      COMMON /LAYCOM/ LAYERS,INFATM,NLAY,LAYHT,LAYANG,LAYSF,LAYTYP,
     1 FSTLAY,LSTLAY

C-----------------------------------------------------------------------
C
C	Variables to describe atmospheric calculations.
C
C-----------------------------------------------------------------------
      LOGICAL LIMB,NADLAY,WF,CG,THERM,BROAD,ABSORB,BINBB,SCATTER     
C LIMB: =.TRUE. if limb calculation.
C WF: =.TRUE. if weighting function calculation.
C CG: =.TRUE. if curtis-godson calculation.
C THERM: =.TRUE. if thermal emission calculation.
C BROAD: =.TRUE. if broad band calculation (i.e. emission outside lbl).
C ABSORB: =.TRUE. if calculate absorption instead of transmission.
C BINBB: =.TRUE. if calculate planck function at centre of bin in genlbl.
C SCATTER: = .TRUE. if full scattering calculation.

      LOGICAL SINGLE,HEMISP,NEARLIMB,NETFLUX
C SINGLE: =.TRUE. if single-scattering calculation.
C HEMISP: =.TRUE. if hemispherical calculation.
C NEARLIMB: =.TRUE. if calculating near limb in a spherical atmosphere
C NETFLUX: =.TRUE. if doing net flux calculation

      REAL HT,ANGLE
C HT: lowest height to use for the calculation, the tangent height
C for limb paths and normally the bottom of the model for nadir paths

      INTEGER BOTLAY,IPZEN
C BOTLAY is the bottom layer to use for atmosphere.
C IPZEN defines where the zenith angle is defined. 0 = at bottom of 
C       bottom layer; 1 = at the 0km altitude level; and 2 = at the very top
C       of the atmosphere.
      COMMON /ATMCOM/ LIMB,NADLAY,WF,CG,THERM,BROAD,ABSORB,BINBB,
     1 HT,ANGLE,BOTLAY,IPZEN

C-----------------------------------------------------------------------
C
C	Variables used to describe cell paths.
C
C-----------------------------------------------------------------------
