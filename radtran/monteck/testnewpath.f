      PROGRAM TESTNEWPATH
C     $Id:
C     ****************************************************************
C     New program, using monte-carlo multiple scattering routines, 
C     which calculates the thermal emission from a spherical atmosphere
C     at any zenith angle. Code may thus do nadir AND limb calculations.
C
C     Code uses correlated-k spectra model
C
C     Pat Irwin    13/5/99
C     Modified from nadirck	PGJI	27/1/00
C     Modified from mcck	PGJI	20/12/01
C
C     ****************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NPRO,NGAS,NCONT,NCONT1,I,MPHOT,J
      PARAMETER (MPHOT=100000)
      REAL P(MAXPRO),T(MAXPRO),H(MAXPRO)
      REAL VMR(MAXPRO,MAXGAS),VREF,ACC
      INTEGER ID(MAXGAS),ISO(MAXGAS),ISTEP,IDUMP
      REAL RADIUS,MOLWT,DUST(MAXPRO,MAXCON),XG,XCORR(MAXCON)
      REAL WCORR,GCORR
      REAL TOTAM,PRESS,TEMP,AMOUNT(MAXGAS),PP(MAXGAS)
      REAL CONT(MAXCON)
      REAL CTOTAM,CTAUG,CTAUA,CTAUS,DTR,SOLZEN,SOLPHI
      REAL AVEC(3),SVEC(3),PVEC(3),TAUA,GALB,TGROUND
      REAL THETA(MAXCON,100)
      REAL TAUS,DVEC1(3),RES(MPHOT,3),SOLVEC(3),SOLRAD,DEVSUN
      REAL OMEGA1,OMEGA2,XFAC,SOLAR
      REAL XOMEGA(MAXCON),XSEC(MAXCON),VV,PI,TAUG1,TAUA1,TAUS1
      REAL HTAN,TAUG,ZEN,ZANG,ALTITUDE,CALCALT
      REAL XSECNADIR(MAXCON),SDIST
      INTEGER NLAMBDA,ISPACE,IRAY
      PARAMETER (PI=3.1415927)
      REAL XLAMBDA(MAXSEC),PHASED(MAXCON,MAXSEC,5)
      REAL MEAN,SDEV,MSCAT
      CHARACTER*100 XHGFILE

      INTEGER NITER,IDUM,IGDIST,ITER
      REAL XHG(MAXCON,3),TAUREQ,TMEAN,FSCAT,TAUSCAT(MAXCON)

      INTEGER NG
      REAL DEL_G(MAXG),TABK(MAXG,MAXPRO)

      INTEGER K,NWAVE,MWAVE,IWAVE,NPHASE
      PARAMETER(MWAVE=500)
      REAL VMIN,VMAX,DELV,VWAVE(MWAVE)

      NPHASE=100
      DTR = PI/180.0

C     Read in T/P/vmr profile
      CALL READPROF(NPRO,NGAS,RADIUS,MOLWT,XG,P,T,H,VMR,ID,ISO)

C     Read in dust profile
      CALL READDUST(H,NPRO,NCONT,DUST)

C     Read in scattering properties of dust
      CALL PROMPT('Enter name of aerosl xsc file : ')
      READ(5,1)XHGFILE
1     FORMAT(A)

      CALL GETHG(XHGFILE,NCONT1,NLAMBDA,XLAMBDA,PHASED)

      IF(NCONT1.NE.NCONT)THEN
       PRINT*,'.pha file is incompatible with dust.prf file'
       PRINT*,'NCONT1,NCONT = ',NCONT1,NCONT
       STOP
      ENDIF

      CALL PROMPT('Enter starting altitude : ')
      READ*,HTAN
      PVEC(1)=0.0
      PVEC(2)=0.0
      PVEC(3)=RADIUS + HTAN

      print*,PVEC

      CALL PROMPT('Enter direction zenith and azimuth : ')
      READ*,SOLZEN,SOLPHI

      SOLVEC(1)=SIN(SOLZEN*DTR)*COS(SOLPHI*DTR)
      SOLVEC(2)=SIN(SOLZEN*DTR)*SIN(SOLPHI*DTR)
      SOLVEC(3)=COS(SOLZEN*DTR)

      print*,'SOLVEC : ',SOLVEC
       
      CALL PROMPT('Enter ISPACE, IRAY, VV : ')
      READ*,ISPACE,IRAY,VV
      NWAVE=1
      VWAVE(1)=VV
      IWAVE=1

C      Interpolate k-tables to VV
       CALL GENTABSCK(NPRO,NGAS,ID,ISO,P,T,VMR,NWAVE,VWAVE,
     1 VV,IWAVE,ISPACE,NG,DEL_G,TABK)
     

      CALL PROMPT('Enter required OD, IGDIST : ')
      READ*,TAUREQ,IGDIST

C     Interpolate scattering properties to VV
      CALL INTERPHG(VV,NCONT,NLAMBDA,XLAMBDA,PHASED,XSEC,XOMEGA,XHG)


      CALL NEWPATH(VV,IRAY,TAUREQ,PVEC,SOLVEC,NPRO,NGAS,NCONT,MOLWT,
     1 RADIUS,P,T,H,DUST,TABK,IGDIST,XSEC,XOMEGA,TMEAN,
     3 FSCAT,TAUSCAT)

      PRINT*,'TMEAN',TMEAN
      PRINT*,'PVEC',PVEC
      PRINT*,'FSCAT',FSCAT
      NCONT1=NCONT+IRAY
      PRINT*,'TAUSCAT : ',(TAUSCAT(I),I=1,NCONT1)

      PRINT*,'Altitude = ',CALCALT(PVEC,RADIUS)

      END
