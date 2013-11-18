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
      INCLUDE '../includes/constdef.f'
C		Passed variables

      INTEGER	Inormal,IRAY,flagh2p
      INTEGER	nwave,IABSORB(5)
      REAL OUTPUT(NWAVE),DABSORB(7)

      INTEGER NPRO,NGAS,NCONT,NCONT1,I,MPHOT,J
      PARAMETER (MPHOT=100000)
      REAL P(MAXPAT),T(MAXPAT),H(MAXPAT),VMR(MAXPAT,MAXGAS),VREF,ACC
      INTEGER ID(MAXGAS),ISO(MAXGAS),ISTEP,IDUMP,NAB,NSOL,NGR
      REAL RADIUS,MOLWT,DUST(MAXPAT,maxcon),XG,TOTAM
      REAL AMOUNT(MAXGAS),PP(MAXGAS),AvgCONTMP
      REAL XLEN,PRESS,TEMP,CONT(maxcon)
      REAL DTR,SOLZEN,SOLPHI,GAMMA,SOLZEN1
      REAL AVEC(3),SVEC(3),PVEC(3),TAUA,GALB,TGROUND,THETA(maxcon,100)
      REAL TAUS,DVEC1(3),RES(MPHOT,3),SOLVEC(3),SOLRAD,DEVSUN
      REAL OMEGA1,OMEGA2,XFAC,SOLAR,GEMI
      REAL sol_ang,emiss_ang,aphi
      INTEGER NLAMBDA,ISPACE,id1,IGAS
      REAL XOMEGA(MAXCON),XSEC1(MAXCON),VV,X
      REAL HTAN,TAUG,ZEN,ZANG,ALTITUDE,CALCALT,SDIST
      REAL XLAMBDA(MAXSEC),PHASED(MAXCON,MAXSEC,5)
      REAL MEAN,SDEV,MSCAT
      CHARACTER*100 XHGFILE,IPFILE,OPFILE,ACCFILE
      LOGICAL ACCEXIST

      INTEGER NITER,IDUM,IGDIST,ITER
      REAL XHG(MAXCON,3)

      INTEGER NG
      REAL DEL_G(MAXG),TABK(MAXG,MAXPAT)
      REAL T1,TRANSREQ,TACTUAL

      INTEGER K,IWAVE,NPHASE

C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL    vwave(nwave),hfp(maxpat), tsurf,dv,esurf
      REAL    vem(maxsec),emissivity(maxsec),interpem
      INTEGER nem

      REAL	XCOM,XNEXT,FPARA,XRAY
      double precision dpi,draddeg

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


C	PRINT*,'CIRSRAD_WAVEMC calling input parameters'
C	PRINT*,'SDist = ',SDist
C	PRINT*,'nwave',nwave,(vwave(i),i=1,nwave)
C        PRINT*,'sol_ang, emiss_ang',sol_ang,emiss_ang
C        PRINT*,'NG = ',NG

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

      IF(IDUMP.EQ.1)THEN
        OPEN(34,FILE='monteck_dump.out',STATUS='UNKNOWN')
        WRITE(34,*)'Standard output on each line is : '
        WRITE(34,*)'X = Wavelength/wavenumber'
        WRITE(34,*)'MEAN = Mean radiance'
        WRITE(34,*)'SDEV = Standard deviation of radiance'
        WRITE(34,*)'MSCAT = Mean number of scattering events'    
        WRITE(34,*)'NAB = Number of photons absorbed in atm'
        WRITE(34,*)'NGR = Number of photons absorbed by ground'
        WRITE(34,*)'NSOL = Number of photons encountering Sun'
        WRITE(34,*)'X, MEAN, SDEV, MSCAT, NAB, NSOL, NGR, ITER, NITER,
     1   SDEV/MEAN'

        OPEN(35,FILE='cirsrad_waveMC.out',STATUS='UNKNOWN')
        WRITE(35,*)NWAVE
      ENDIF

      CALL READPROF(OPFILE,NPRO,NGAS,RADIUS,MOLWT,XG,P,T,H,VMR,ID,ISO)

C     Read in dust profile
      IPFILE='aerosol.prf'
      CALL READDUST(IPFILE,H,NPRO,NCONT,DUST)

C      CALL PROMPT('Enter random -ve integer : ')
C      READ*,IDUM

      IDUM=-7
      ZEN = ASIN(RADIUS/(RADIUS+H(NPRO)))
C      PRINT*,'Zenith angle of limb = ',ZEN*180.0/PI
C      IF(IDUMP.EQ.1)WRITE(34,*)'Zenith angle of limb = ',ZEN*180.0/PI

C      PRINT*,'Enter required tangent height. ( -ve numbers indicate'
C      CALL PROMPT('Zenith angle)  : ')
C      READ*,HTAN 

C     ************* Need to change this to be more general. Precludes
C     doing limb scattering with sunlight on ***********************
 

      IF(EMISS_ANG.LT.0.0) THEN
         HTAN = SOL_ANG
         ZANG = (1.0/DTR)*ASIN((HTAN+RADIUS)/(H(NPRO)+RADIUS))
C         print*,'HTAN,ZANG = ',HTAN,ZANG
C         PRINT*,'Enter zenith and aziumuth angle of Sun'
C         CALL PROMPT('at tangent point : ')
C         READ*,SOLZEN,SOLPHI        
C         N.B. SOLPHI = 0 indicates FORWARD scattering
          SOLZEN=180.0
         SOLPHI=0.0
         SOLVEC(1)=SIN(SOLZEN*DTR)*COS(SOLPHI*DTR)
         SOLVEC(2)=SIN(SOLZEN*DTR)*SIN(SOLPHI*DTR)
         SOLVEC(3)=COS(SOLZEN*DTR)
         GAMMA = 90.0-ZANG
         CALL YROTATE(SOLVEC,GAMMA)
         SOLZEN1 = ACOS(SOLVEC(3))/DTR
      ELSE
         ZANG = EMISS_ANG
         HTAN = (RADIUS+H(NPRO))*SIN(ZANG*DTR) - RADIUS
C         print*,'A: HTAN,ZANG = ',HTAN,ZANG
         SOLZEN=SOL_ANG
         SOLPHI=APHI
C         PRINT*,'Enter zenith and aziumuth angle of Sun'
C         CALL PROMPT('at point where photons enter atmosphere : ')
C         READ*,SOLZEN,SOLPHI
C        N.B. SOLPHI = 0 indicates FORWARD scattering
         SOLVEC(1)=SIN(SOLZEN*DTR)*COS(SOLPHI*DTR)
         SOLVEC(2)=SIN(SOLZEN*DTR)*SIN(SOLPHI*DTR)
         SOLVEC(3)=COS(SOLZEN*DTR)
      ENDIF
        
C      print*,'SOLZEN,SOLPHI = ',SOLZEN,SOLPHI
C      print*,'SOLVEC = ',SOLVEC

C      IF(IDUMP.EQ.1)THEN
C        WRITE(34,*)'Tangent altitude, Zenith angle = ',HTAN,ZANG
C      ENDIF
 
C      CALL PROMPT('Enter distance from Sun (AU) : ')
C      READ*,SDIST

C     Calculate angular radius of Sun in degrees
      SOLRAD = 0.5*0.533128/SDIST

C      CALL PROMPT('Enter acceptable deviation from Sun (deg) : ')
C      READ*,DEVSUN
      DEVSUN=10.0

C      OMEGA1 = 4*PI*(sin(SOLRAD*DTR*0.5)**2)

C     OMEGA2 is solid angle of Sun dor assumed inflated angular radius
      OMEGA2 = 4*PI*(sin(DEVSUN*DTR*0.5)**2)

C      NITER=20
      NITER=10000
C      NITER=1000
C      CALL PROMPT('Enter max. number of photons : ')
C      READ*,NITER

      IF(NITER.GT.MPHOT)THEN
       Print*,'cirsrad_waveMC. NITER must be less than or equal'
       print*,'to MPHOT. MPHOT = ',MPHOT
       STOP
      ENDIF


      CALL GETHG(OPFILE,NCONT1,NLAMBDA,XLAMBDA,PHASED)

      IF(NCONT1.NE.NCONT)THEN
       PRINT*,'.pha file is incompatible with dust.prf file'
       PRINT*,'NCONT1,NCONT = ',NCONT1,NCONT
       STOP
      ENDIF


      TGROUND=TSURF

C      CALL PROMPT('Enter desired convergence accuracy (rad.units): ')
C      READ*,ACC
      CALL FILE(OPFILE,ACCFILE,'acc')
      INQUIRE(FILE=ACCFILE,EXIST=ACCEXIST)

      IF(ACCEXIST)THEN
       OPEN(12,FILE=ACCFILE,STATUS='OLD')
        READ(12,*)ACC
       CLOSE(12)
      ELSE
       ACC = 1e-32
      ENDIF


C      IF(IDUMP.EQ.1)THEN
C       WRITE(34,*)NITER,'   ! Number of iterations'
C       WRITE(34,*)NWAVE,'    ! NWAVE'
C       WRITE(34,*)(VWAVE(I),I=1,NWAVE)
C      ENDIF
   

      DO 1001 IWAVE=1,NWAVE

       X = VWAVE(IWAVE)
C      Set vv to the current WAVENUMBER
       vv = x
       if(ispace.eq.1)then
          vv=1e4/x
       endif

C      Get solar flux at this wavelength/wavenumber     
       CALL GET_SOLAR_WAVE(X,SDIST,SOLAR)

C      Output from get_solar_wave is W cm-2 um-1 or W cm-1 (cm-1)-1. Need
C      to convert this to surface radiance of sun for the assumed solar size.
       SOLAR=SOLAR/OMEGA2

C      Need extra correction factor of 2.0. Don't know why!!
       SOLAR=SOLAR*0.5

C      Finally need to correct for the solar zenith angle
C      print*,'C',solar,solzen,dtr
C      Actually, I'm not sure we do!!
C      SOLAR=SOLAR*COS(SOLZEN*DTR)

C      Interpolate k-tables and gas continua to X
       CALL GENTABSCK1(OPFILE,NPRO,NGAS,ID,ISO,P,T,VMR,NWAVE,VWAVE,
     1 X,ISPACE,NG,DEL_G,TABK)

C       print*,'WT',(TABK(I,4),I=1,NG)

C      Interpolate scattering properties to X
       CALL INTERPHG(X,NCONT,NLAMBDA,XLAMBDA,PHASED,XSEC1,XOMEGA,XHG)

C      Regrid phase functions to equal probability steps
C      Also add in extra dust 'type' to deal with Rayleigh
C      scattering, by setting XHG(NCONT+1,1) to be negative
       NCONT1 = NCONT+1
       XHG(NCONT1,1)=-1
       CALL PHASPROB(NCONT1,XHG,NPHASE,THETA)

      
C       print*,'DD',(TABK(1,I),I=1,NPRO)
       DO 102 I=1,NPRO

C       Calculate CIA and continuum absorption in a 1km-long path


C       TOTAM is the number of molecules/cm2 along a path of 1km length.
        TOTAM=MODBOLTZ*P(I)/T(I)
        XLEN=1.0

C       NCIACON diagnostic
        id1=0

        DO IGAS=1,NGAS
         AMOUNT(IGAS)=TOTAM*VMR(I,IGAS)
         PP(IGAS)=P(I)*VMR(I,IGAS)
        ENDDO

C        print*,'FLAGH2P = ',FLAGH2P
        AvgCONTMP=0.
        IF(FLAGH2P.EQ.1)THEN
          FPARA = HFP(I)
          CALL NPARACON(VV,P(I),T(I),NGAS,id,iso,Amount,PP,FPARA,      
     1     XLen,AvgCONTMP,IABSORB,DABSORB,id1)
        ELSE
          CALL NCIACON(vv,P(I),T(I),INormal,NGas,id,iso,Amount,PP,
     1     XLen,AvgCONTMP,IABSORB,DABSORB,id1)
        ENDIF

C       AvgCONTMP contains opacity of a 1 km-long path. Need to convert to
C       same units as TABK, namely (molecule/cm2)-1 and add to each
C       element of the TABK array at this level.  

        DO IG=1,NG
C          IF(I.EQ.4)THEN
C           print*,'XX',X,VV,I,IG,TABK(IG,I),TOTAM,AvgCONTMP,
C     &		1e20*AvgCONTMP/TOTAM
C          ENDIF

          TABK(IG,I)=TABK(IG,I)+1e20*AvgCONTMP/TOTAM

        ENDDO


102    CONTINUE
C       print*,'DE',(TABK(1,I),I=1,NPRO)




C      Interpolate emissivity and albedo to X
C       print*,'GALB = ',GALB
       CALL VERINT(VEM,EMISSIVITY,NEM,GEMI,X) 
       GALB = 1.0-GEMI
C       print*,'GALB, GEMI  = ',GALB,GEMI




C      Use MonteCarlo optical depth code to compute transmission
C      to ground to check all is well.
       T1=0.
       TRANSREQ=0.0
       DO IG=1,NG

        SVEC(1)=0
        SVEC(2)=0.0
        SVEC(3)=RADIUS+H(NPRO)

        DVEC1(1) = SIN(ZANG*PI/180.0)
        DVEC1(2) = 0.0
        DVEC1(3) = -COS(ZANG*PI/180.0)

C        CALL TMPPATH(VV,IRAY,TRANSREQ,SVEC,DVEC1,NPRO,NGAS,NCONT,
C     1 MOLWT,RADIUS,P,T,H,DUST,TABK,IG,XSEC,XOMEGA,TACTUAL)
C        T1=T1+DEL_G(IG)*TACTUAL
C        print*,X,VV,IG,DEL_G(IG),TACTUAL
       ENDDO

C       print*,'Initial position vector : ',SVEC
C       print*,'Zenith angle : ',ZANG
C       print*,'Initial direction vector : ',DVEC1     
C       print*,'Initial solar vector : ',SOLVEC     
C       print*,'Niter, IDUM',NITER,IDUM
C       print*,'DEVSUN, SOLAR',DEVSUN,SOLAR
C       print*,'NPRO,NGAS,NCONT,MOLWT,RADIUS',NPRO,NGAS,
C     1   NCONT,MOLWT,RADIUS
C       print*,'GALB,TGROUND,IRAY',GALB,TGROUND,IRAY
C       print*,'XSEC,XOMEG',(XSEC1(J),J=1,NCONT),(XOMEGA(J),J=1,NCONT)

C       print*,'FF',(TABK(I,4),I=1,NG)

       SVEC(1)=0
       SVEC(2)=0.0
       SVEC(3)=RADIUS+H(NPRO)

       DVEC1(1) = SIN(ZANG*PI/180.0)
       DVEC1(2) = 0.0
       DVEC1(3) = -COS(ZANG*PI/180.0)


       CALL MCPHOTONCK(NITER,IDUM,
     1    XSEC1,XOMEGA,NPHASE,THETA,
     2    SVEC,DVEC1,SOLVEC,DEVSUN,SOLAR,TABK,NG,DEL_G,
     3    NPRO,NGAS,NCONT,MOLWT,RADIUS,P,T,H,DUST,
     4    GALB,TGROUND,IRAY,RES,ACC,MEAN,SDEV,
     5    MSCAT,ITER,ISPACE,X,NAB,NSOL,NGR)

C       print*,'X,NAB,NSOL,NGR',X,NAB,NSOL,NGR

       IF(IDUMP.EQ.1)THEN
C        WRITE(34,*)ITER
C        DO I=1,ITER
C         WRITE(34,*)(RES(I,J),J=1,3)
C        ENDDO

C        print*,'waveMC',X,MEAN,SDEV,MSCAT,NAB,NSOL,NGR,ITER,NITER
         IF(MEAN.NE.0.0)THEN
          WRITE(34,*)X,MEAN,SDEV,MSCAT,NAB,NSOL,NGR,ITER,NITER,
     1     SDEV/MEAN
         ELSE
          WRITE(34,*)X,MEAN,SDEV,MSCAT,NAB,NSOL,NGR,ITER,NITER
         ENDIF    

         WRITE(35,*)X,MEAN,SDEV,VV,T1

       ENDIF


       OUTPUT(IWAVE)=MEAN

1001  CONTINUE 

      CLOSE(12)
      IF(IDUMP.EQ.1)THEN
       CLOSE(34)
       CLOSE(35)
      ENDIF

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
