C     $Id: pathcom.f,v 1.8 2007-11-02 15:26:08 irwin Exp $
C***********************************************************************
C_TITL:	PATHCOM
C
C_DESC:	The variables defined here are those normally used when 
C	calculating atmospheric paths and will normally be explcitly
C	included in any code which computes these paths, performs final
C	calculations on the results or plots out the results. The
C	variables are divided into two groups:
C		1) variables which roughly correspond to the input
C		parameters of GENLBL which will be used by LBL/GENLBL or
C		other transmission calculation programs.
C
C		2) variables which can be used to pass additional
C		parameters to final calculation routines. These have
C		generic names so that they can be used for many different
C		purposes.
C
C	The variables are defined in advance so that a common driver file
C	format can be used for all transmission calculations, The type
C	number of a calculation can be used to differeniate between
C	different calculations.
C
C	Defined types are stored in [calcutt.radtran]types.txt and DO NOT
C	correspond to the types defined in the original version.
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST:	18nov92	SBC	Edded EMTEMP and BASET. removed WEIGHT
C	23nov02 PDP     Added to the header, minor formatting, etc.
C***********************************************************************
C-----------------------------------------------------------------------
C
C	Variables to be used by GENLBL.
C
C-----------------------------------------------------------------------
C MAXPAT: maximum number of paths.
C MAXINC: maximum number of the layers included in any one path i.e. each
C layer can be included more than once in each path or not at all.
C MAXLAY: maximum number of layers to be considered.
C MAXCON: maximum number of points in the continuum spectra.
C MAXFIL: maximum number of points in the filter profile.


      INTEGER NGAS,NFILT,NLAYER,NPATH,IDGAS(MAXGAS),IPROC(MAXGAS)
      INTEGER NPOINT,ICONV,IMOD(MAXPAT),NLAYIN(MAXPAT),FLAGH2P,FLAGC
      INTEGER LAYINC(MAXINC,MAXPAT),NCONT,ISOGAS(MAXGAS),INLTE
      REAL AMOUNT(MAXLAY,MAXGAS),TEMP(MAXLAY),PRESS(MAXLAY),VMIN,DELV
      REAL SCALE(MAXINC,MAXPAT),FWHM,WING,VREL,FILTER(MAXFIL)
      REAL VFILT(MAXFIL),CONT(MAXCON,MAXLAY),VCONT(MAXCON)
      REAL DOP(MAXLAY),PP(MAXLAY,MAXGAS),ERRLIM(MAXPAT),QH,QHe
      REAL BASEP(MAXLAY),BASEH(MAXLAY),DELH(MAXLAY),BASET(MAXLAY)
      REAL TOTAM(MAXLAY),EMTEMP(MAXINC,MAXPAT),HFP(MAXLAY)
      REAL HFC(MAXLAY)
      INTEGER IFC(MAXCON,MAXLAY)
      COMMON /PATHVA/ AMOUNT,TEMP,PRESS,VMIN,DELV,SCALE,FWHM,WING,VREL,
     1 FILTER,VFILT,CONT,VCONT,DOP,PP,ERRLIM,NGAS,NFILT,NLAYER,NPATH,
     2 IDGAS,IPROC,NPOINT,ICONV,IMOD,NLAYIN,LAYINC,NCONT,ISOGAS,
     3 BASEP,BASEH,DELH,BASET,TOTAM,EMTEMP,QH,QHe,INLTE,FLAGH2P,HFP,
     4 FLAGC,HFC,IFC

      CHARACTER*100 LINKEY
      COMMON /PATHCH/ LINKEY

C-----------------------------------------------------------------------
C
C	Calculation Record Arrays.
C
C-----------------------------------------------------------------------
      INTEGER MAXCAL,MAXCP,MAXIP,MAXRP
      PARAMETER (MAXCAL=300,MAXCP=2,MAXIP=10,MAXRP=50)
C MAXCAL: maximum number of calculation which can be performed.
C MAXCP: maximum number of character parameters allowed in each 
C calculation.
C MAXIP: maximum number of integer parameters allowed in each calculation.
C MAXRP: maximum number of real parameters allowed in each calculation.

      INTEGER NCALC,ITYPE(MAXCAL),ICALD(MAXIP,MAXCAL),NINTP(MAXCAL)
      INTEGER NREALP(MAXCAL),NCHP(MAXCAL)
      REAL RCALD(MAXRP,MAXCAL)
      COMMON /CALINT/ NCALC,ITYPE,ICALD,NINTP,NREALP,RCALD,NCHP

      CHARACTER*100 OPFILE,DRFILE,XFILE
      CHARACTER*30 CCALD(MAXCP,MAXCAL)
      COMMON /CALCH/ OPFILE,CCALD,DRFILE,XFILE


C-----------------------------------------------------------------------
C
C	Dust arrays.
C
C-----------------------------------------------------------------------
      INTEGER NSEC
      REAL XSEC(2,MAXCON,MAXSEC),VSEC(MAXSEC)
      COMMON /DUSTSPEC/ XSEC,VSEC,NSEC
