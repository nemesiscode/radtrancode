************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C_TITLE:		SUBROUTINE CIRSRAD_WAVEMC
C
C_DESCR:
C
C_ARGS:	Input Variables
C
C       SDIST:REAL Distance from Sun (as prompted by CIRSDRV) in
C			units of AU.
C	NWAVE		The number of calculation wavenumbers.
C	VWAVE		Bin centres in wavenumber space.
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

      SUBROUTINE cirsrad_waveMC(opfile,IDUM,sol_ang,emiss_ang,
     1  aphi,sdist,INormal,Iray,ispace,nwave,vwave, nem, vem,
     2  emissivity,tsurf,flagh2p,hfp,output)

      IMPLICIT NONE

C		Internal dimensions

C     Defines the maximum values for a series of variables (layers,
C     bins, paths, etc.)
      INCLUDE '../includes/arrdef.f'

C		Passed variables

      INTEGER	Inormal,IRAY,flagh2p
      INTEGER	nwave
      REAL OUTPUT(NWAVE)

      INTEGER NPRO,NGAS,NCONT,NCONT1,I,MPHOT,J
      PARAMETER (MPHOT=100000)
      REAL P(MAXPAT),T(MAXPAT),H(MAXPAT),VMR(MAXPAT,MAXGAS),VREF,ACC
      INTEGER ID(MAXGAS),ISO(MAXGAS),ISTEP,IDUMP,NAB,NSOL,NGR
      REAL RADIUS,MOLWT,DUST(MAXPAT,maxcon),XG
      REAL TOTAM,PRESS,TEMP,AMOUNT(MAXGAS),PP(MAXGAS),CONT(maxcon)
      REAL DTR,SOLZEN,SOLPHI,GAMMA,SOLZEN1
      REAL AVEC(3),SVEC(3),PVEC(3),TAUA,GALB,TGROUND,THETA(maxcon,100)
      REAL TAUS,DVEC1(3),RES(MPHOT,3),SOLVEC(3),SOLRAD,DEVSUN
      REAL OMEGA1,OMEGA2,XFAC,SOLAR,GEMI
      REAL sol_ang,emiss_ang,aphi,PI
      INTEGER NLAMBDA,MLAMBDA,ISPACE,MCONT
      PARAMETER (PI=3.1415927,MLAMBDA=200,MCONT=10)
      REAL XOMEGA(MCONT),XSEC1(MCONT),VV
      REAL HTAN,TAUG,ZEN,ZANG,ALTITUDE,CALCALT,SDIST
      REAL XLAMBDA(MLAMBDA),PHASED(MCONT,MLAMBDA,5)
      REAL MEAN,SDEV,MSCAT
      CHARACTER*100 XHGFILE,IPFILE,OPFILE

      INTEGER NITER,IDUM,IGDIST,ITER
      REAL XHG(MCONT,3)

      INTEGER NG
      REAL DEL_G(MAXG),TABK(MAXG,MAXPAT)

      INTEGER K,IWAVE,NPHASE

C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL    vwave(nwave),hfp(maxpat), tsurf,dv,esurf
      REAL    vem(maxsec),emissivity(maxsec),interpem
      INTEGER nem

      REAL	XCOM,XNEXT,FPARA,XRAY
      double precision get_albedo,dpi,draddeg

C     Dust variables

      INTEGER	nsec
      REAL	vsec(maxsec), xsec(2,maxcon,maxsec)

C     Internal variables

      INTEGER	L, Ig

      common/dust/vsec,xsec,nsec,ncont

C-----------------------------------------------------------------------
C
C	Check input parameters for possible problems.
C
C-----------------------------------------------------------------------


	PRINT*,'CIRSRAD_WAVEMC calling input parameters'
C	PRINT*,'SDist = ',SDist
C	PRINT*,'nwave',nwave,(vwave(i),i=1,nwave)
        PRINT*,'sol_ang, emiss_ang',sol_ang,emiss_ang
        PRINT*,'NG = ',NG
      if (nwave.gt.maxbin) then
        write (*,*) ' CIRSRAD_WAVE: Too many bins'
        write (*,*) ' NWave = ',nwave,' Maxbin = ',maxbin
        stop
      endif

      if (nsec.gt.maxsec) then
        write (*,*) ' CIRSRAD_WAVE: Too many dust continua pts'
        write (*,*) ' Nsec = ',nsec,' Maxsec = ',maxsec
        stop
      endif


      dpi = 4.0d0 * datan(1.0d0)
      draddeg = 2.0d0 * dpi/360.0d0


      NPHASE=100
      DTR = PI/180.0

C      CALL PROMPT('Enter IDUMP : ')
C      READ*,IDUMP
      IDUMP=1

      IF(IDUMP.EQ.1)OPEN(34,FILE='monteck_dump.out',STATUS='UNKNOWN')

      CALL READPROF(OPFILE,NPRO,NGAS,RADIUS,MOLWT,XG,P,T,H,VMR,ID,ISO)

C     Read in dust profile
      IPFILE='aerosol.prf'
      CALL READDUST(IPFILE,H,NPRO,NCONT,DUST)

C      CALL PROMPT('Enter random -ve integer : ')
C      READ*,IDUM

      ZEN = ASIN(RADIUS/(RADIUS+H(NPRO)))
      PRINT*,'Zenith angle of limb = ',ZEN*180.0/PI
C      IF(IDUMP.EQ.1)WRITE(34,*)'Zenith angle of limb = ',ZEN*180.0/PI

C      PRINT*,'Enter required tangent height. ( -ve numbers indicate'
C      CALL PROMPT('Zenith angle)  : ')
C      READ*,HTAN 

C     ************* Need to change this to be more general. Precludes
C     doing limb scattering with sunlight on ***********************
 
      SOLZEN=180.0
      SOLPHI=0.0

      IF(EMISS_ANG.LT.0.0) THEN
         HTAN = SOL_ANG
         ZANG = (1.0/DTR)*ASIN((HTAN+RADIUS)/(H(NPRO)+RADIUS))
         print*,'HTAN,ZANG = ',HTAN,ZANG
C         PRINT*,'Enter zenith and aziumuth angle of Sun'
C         CALL PROMPT('at tangent point : ')
C         READ*,SOLZEN,SOLPHI        
C        N.B. SOLPHI = 0 indicates FORWARD scattering
         SOLVEC(1)=SIN(SOLZEN*DTR)*COS(SOLPHI*DTR)
         SOLVEC(2)=SIN(SOLZEN*DTR)*SIN(SOLPHI*DTR)
         SOLVEC(3)=COS(SOLZEN*DTR)
         GAMMA = 90.0-ZANG
         CALL YROTATE(SOLVEC,GAMMA)
         SOLZEN1 = ACOS(SOLVEC(3))/DTR
      ELSE
         ZANG = EMISS_ANG
         HTAN = (RADIUS+H(NPRO))*SIN(ZANG*DTR) - RADIUS
         print*,'A: HTAN,ZANG = ',HTAN,ZANG
C         PRINT*,'Enter zenith and aziumuth angle of Sun'
C         CALL PROMPT('at point where photons enter atmosphere : ')
C         READ*,SOLZEN,SOLPHI
C        N.B. SOLPHI = 0 indicates FORWARD scattering
         SOLVEC(1)=SIN(SOLZEN*DTR)*COS(SOLPHI*DTR)
         SOLVEC(2)=SIN(SOLZEN*DTR)*SIN(SOLPHI*DTR)
         SOLVEC(3)=COS(SOLZEN*DTR)
      ENDIF
        
C      IF(IDUMP.EQ.1)THEN
C        WRITE(34,*)'Tangent altitude, Zenith angle = ',HTAN,ZANG
C      ENDIF
 
C      CALL PROMPT('Enter distance from Sun (AU) : ')
C      READ*,SDIST

C     Calculate angular radius of Sun
      SOLRAD = 0.5*0.533128/SDIST

C      CALL PROMPT('Enter acceptable deviation from Sun (deg) : ')
C      READ*,DEVSUN
      DEVSUN=10.0

      OMEGA1 = 2*PI*(1.0-cos(SOLRAD*DTR))
      OMEGA2 = 2*PI*(1.0-cos(DEVSUN*DTR))

C     Calculate ratio of projected areas
      XFAC= OMEGA1/OMEGA2

      NITER=10000
C      CALL PROMPT('Enter max. number of photons : ')
C      READ*,NITER

      IF(NITER.GT.MPHOT)THEN
       Print*,'Monteck. NITER must be less than or equal to MPHOT'
       Print*,'MPHOT = ',MPHOT
       STOP
      ENDIF

C     Read in scattering properties of dust
C      CALL PROMPT('Enter name of aerosol H-G file : ')
C      READ(5,1)XHGFILE

      CALL GETHG(OPFILE,NCONT1,NLAMBDA,XLAMBDA,PHASED)

      IF(NCONT1.NE.NCONT)THEN
       PRINT*,'.pha file is incompatible with dust.prf file'
       PRINT*,'NCONT1,NCONT = ',NCONT1,NCONT
       STOP
      ENDIF


      TGROUND=TSURF

C      CALL PROMPT('Enter desired convergence accuracy (rad.units): ')
C      READ*,ACC

       ACC = 1e-32



C      IF(IDUMP.EQ.1)THEN
C       WRITE(34,*)NITER,'   ! Number of iterations'
C       WRITE(34,*)NWAVE,'    ! NWAVE'
C       WRITE(34,*)(VWAVE(I),I=1,NWAVE)
C      ENDIF

      DO 1001 IWAVE=1,NWAVE

       VV = VWAVE(IWAVE)

C      Get solar flux at this wavelength/wavenumber
       CALL GET_SOLAR_WAVE(VV,SDIST,SOLAR)
C      Output from get_solar_wave is W cm-2 um-1 or W cm-1 (cm-1)-1. Need
C      to convert this to surface radiance of sun.
       SOLAR=SOLAR/OMEGA1

C      Also need to correct for fact that we'll actually accept slightly
C      larger angles for calculation.
       SOLAR=SOLAR*XFAC
 
C      Interpolate k-tables and gas continua to VV
       CALL GENTABSCK1(OPFILE,NPRO,NGAS,ID,ISO,P,T,VMR,NWAVE,VWAVE,
     1 VV,ISPACE,NG,DEL_G,TABK)

       DO I=1,NPRO
C        print*,I,NPRO
C        WRITE(34,*)(TABK(J,I),J=1,NG)
C        print*,(TABK(J,I),J=1,NG)
       ENDDO

C      Interpolate scattering properties to VV
       CALL INTERPHG(VV,NCONT,NLAMBDA,XLAMBDA,PHASED,XSEC1,XOMEGA,XHG)

C      Interpolate emissivity and albedo to VV
       CALL VERINT(VEM,EMISSIVITY,NEM,GEMI,VV) 
       GALB = 1.0-GEMI
       print*,'GALB, GEMI  = ',GALB,GEMI

C      Regrid phase functions to equal probability steps
       NCONT1 = NCONT+1
       CALL PHASPROB(NCONT1,XHG,NPHASE,THETA)

       print*,'RADIUS,H(NPRO)',RADIUS,H(NPRO)
       print*,'ZANG',ZANG
       SVEC(1)=0
       SVEC(2)=0.0
       SVEC(3)=RADIUS+H(NPRO)

       DVEC1(1) = SIN(ZANG*PI/180.0)
       DVEC1(2) = 0.0
       DVEC1(3) = -COS(ZANG*PI/180.0)

       print*,'Initial position vector : ',SVEC
       print*,'Zenith angle : ',ZANG
       print*,'Initial direction vector : ',DVEC1     

       CALL MCPHOTONCK(NITER,IDUM,
     1    XSEC1,XOMEGA,NPHASE,THETA,
     2    SVEC,DVEC1,SOLVEC,DEVSUN,SOLAR,TABK,NG,DEL_G,
     3    NPRO,NGAS,NCONT,MOLWT,RADIUS,P,T,H,DUST,
     4    GALB,TGROUND,IRAY,RES,ACC,MEAN,SDEV,
     5    MSCAT,ITER,ISPACE,VV,NAB,NSOL,NGR)

       print*,'VV,NAB,NSOL,NGR',VV,NAB,NSOL,NGR
C       IF(IDUMP.EQ.1)THEN
C        WRITE(34,*)ITER
C        DO I=1,ITER
C         WRITE(34,*)(RES(I,J),J=1,3)
C        ENDDO
C       ENDIF

C       print*,'waveMC',VV,MEAN,SDEV,MSCAT,NAB,NSOL,NGR,ITER,NITER
       WRITE(34,*)VV,MEAN,SDEV,MSCAT,NAB,NSOL,NGR,ITER,NITER

       OUTPUT(IWAVE)=MEAN

1001  CONTINUE 

      IF(IDUMP.EQ.1)CLOSE(34)


      WRITE(*,*)'%cirsrad_waveMC.f :: calculation complete'

C-----------------------------------------------------------------------
C
C	Return and end.
C
C-----------------------------------------------------------------------


	RETURN

	END

************************************************************************
************************************************************************
