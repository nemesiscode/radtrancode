C       $Id:
C***********************************************************************
C_TITL: ARRDEF.F
C
C_DESC:
C
C_ARGS: See definitions below.
C
C_CALL: No calls.
C
C_HIST:	
C***********************************************************************

      INTEGER MAXLAY, MAXPAT, MAXGAS, MAXCON, MAXSEC, MAXBIN,
     1 MAXG, MAXK, MAXMU, MAXSCATPAR, MAXSCATLAY, MAXPHAS, MAXOUT, 
     2 MAXF, MAXV, MAXOUT3, MAXOUT4, MAXRANK, MAXINC,
     3 MAXFIL, MAXPRO, IDIM, IORDP1, IORDER, MPOINT
      PARAMETER (MAXPAT=200,MAXLAY=2*MAXPAT, MAXGAS=30, 
     1 MAXCON=10, MAXSCATPAR=10, MAXSEC=16000, MAXBIN=16000, MAXG=51,
     2 MAXK=30, MAXMU=21, MAXSCATLAY=40, MAXPHAS=50, MAXOUT=100000,
     3 MAXF=41, MAXV=401, MPOINT=16000, MAXOUT3=MAXOUT,
     4 MAXOUT4=MAXOUT*35, MAXRANK=10*MAXG*MAXG, MAXINC=2*MAXPAT, 
     5 MAXFIL=1000, MAXPRO=MAXPAT, IDIM=2048,IORDER=2,IORDP1=IORDER+1)
C     MAXLAY: maximum number of layers allowed.
C     MAXPAT: maximum number of paths allowed.
C     MAXGAS: maximum number of gases allowed.
C     MAXCON: maximum number of aerosol types
c       NB: MAXCON must equal MAXSCATPAR. These cannot be set to arbitrary values.
c       They are hard coded into nemesis/mod_scatter.f and scatter/phase1.f
c       Only 10 and 14 are currently defined. 10 is default/original. If you want 14
c       then need to also edit mod_scatter.f and phase1.f (comment/uncomment some lines)
C     MAXSEC: Maximum number of x-section wavelengths
C     MAXBIN: maximum number of bins allowed.
C     MAXG: maximum number of g-ordinates/weights allowed.
C     MAXRANK: Max size of arrays for overlapping gas k-distributions.
C     MAXK: maximum number of T/P-ordinates in k-tables
C     MAXMU: maximum number zenith angles for scattering calc.
C     MAXCON: MAXSCATPAR: Maximum number of parameters defining phase function
C     MAXSCATLAY: maximum number of layers for scattering calc.
C     MAXPHAS: Maximum number of phase angles in phase function
C     MAXOUT: maximum number of points output allowed (#1).
C     MAXF: maximum number of fourier azimuth components
C     MAXV: maximum number of variable parameters for gradiant calcs
C     MAXOUT3: maximum number of points output allowed (#3).
C     MAXOUT4: maximum number of points output allowed (#4).
C     MAXFIL: maxumum number of points in filter function.
C     MAXPRO: maximum number of levels in prf files.
C     IDIM: maximum array size in matrix inversion routines
C     IORDER: order of polynomial to fit to continuum
C     IORDP1: Number of elements of polynomial array
C     MPOINT: Maximum number of elements in k-table pre-spectra
