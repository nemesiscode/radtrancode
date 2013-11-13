      PROGRAM MONTECK
C     $Id:
C     ****************************************************************
C     Monte-carlo multiple scattering program which calculates the radiance 
C     observed at any angle in a spherical atmosphere, with and without
C     sunshine. Code may thus do nadir AND limb calculations.
C
C     Code uses correlated-k spectra model
C
C     Pat Irwin    			13/5/99
C     Modified from nadirck	PGJI	27/1/00
C     Modified from mcck	PGJI	20/12/01
C     Modified from Mcck	PGJI	1/8/05
C
C     ****************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NPRO,NGAS,NCONT,NCONT1,I,MPHOT,J
      PARAMETER (MPHOT=100000)
      REAL P(MAXPRO),T(MAXPRO),H(MAXPRO),VMR(MAXPRO,MAXGAS),VREF,ACC
      INTEGER ID(MAXGAS),ISO(MAXGAS),ISTEP,IDUMP,NAB,NSOL,NGR
      REAL RADIUS,MOLWT,DUST(MAXPRO,MAXCON),XG
      REAL TOTAM,PRESS,TEMP,AMOUNT(MAXGAS),PP(MAXGAS),CONT(MAXCON)
      REAL DTR,SOLZEN,SOLPHI,GAMMA,SOLZEN1
      REAL AVEC(3),SVEC(3),PVEC(3),TAUA,GALB,TGROUND,THETA(MAXCON,100)
      REAL TAUS,DVEC1(3),RES(MPHOT,3),SOLVEC(3),SOLRAD,DEVSUN
      REAL OMEGA1,OMEGA2,XFAC,SOLAR
      REAL XOMEGA(MAXCON),XSEC(MAXCON),VV,PI
      REAL HTAN,TAUG,ZEN,ZANG,ALTITUDE,CALCALT,SDIST
      INTEGER NLAMBDA,ISPACE,IRAY,ICHANNEL
      PARAMETER (PI=3.1415927)
      REAL XLAMBDA(MAXSEC),PHASED(MAXCON,MAXSEC,5)
      REAL MEAN,SDEV,MSCAT
      CHARACTER*100 XHGFILE,IPFILE,SOLNAME

      INTEGER NITER,IDUM,IGDIST,ITER
      REAL XHG(MAXCON,3)

      INTEGER NG
      REAL DEL_G(MAXG),TABK(MAXG,MAXPRO)

      INTEGER K,NWAVE,MWAVE,IWAVE,NPHASE
      PARAMETER(MWAVE=500)
      REAL VMIN,VMAX,DELV,VWAVE(MWAVE)

      REAL TOT_TIME
      DOUBLE PRECISION TIME,TIME1,TIME2

C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL GETTIME(TIME)
      TIME1= TIME

      NPHASE=100
      DTR = PI/180.0

      CALL PROMPT('Enter IDUMP : ')
      READ*,IDUMP

      IF(IDUMP.EQ.1)OPEN(34,FILE='monteck_dump.out',STATUS='UNKNOWN')

C     Read in T/P/vmr profile
      CALL PROMPT('Enter input profile name : ')
      READ(5,1)IPFILE
1     FORMAT(A)

      CALL READPROF(IPFILE,NPRO,NGAS,RADIUS,MOLWT,XG,P,T,H,VMR,ID,ISO)

C     Read in dust profile
      IPFILE='aerosol.prf'
      CALL READDUST(IPFILE,H,NPRO,NCONT,DUST)

      CALL PROMPT('Enter random -ve integer : ')
      READ*,IDUM

      ZEN = ASIN(RADIUS/(RADIUS+H(NPRO)))
      PRINT*,'Zenith angle of limb = ',ZEN*180.0/PI
      IF(IDUMP.EQ.1)WRITE(34,*)'Zenith angle of limb = ',ZEN*180.0/PI


      PRINT*,'Enter required tangent height. ( -ve numbers indicate'
      CALL PROMPT('Zenith angle)  : ')
      READ*,HTAN 
        
      IF(HTAN.LT.0.0)THEN
         ZANG = -HTAN
         HTAN = (RADIUS+H(NPRO))*SIN(ZANG*DTR) - RADIUS
         PRINT*,'Enter zenith and aziumuth angle of Sun'
         CALL PROMPT('at point where photons enter atmosphere : ')
         READ*,SOLZEN,SOLPHI
C        N.B. SOLPHI = 0 indicates FORWARD scattering
         SOLVEC(1)=SIN(SOLZEN*DTR)*COS(SOLPHI*DTR)
         SOLVEC(2)=SIN(SOLZEN*DTR)*SIN(SOLPHI*DTR)
         SOLVEC(3)=COS(SOLZEN*DTR)
      ELSE
         ZANG = (1.0/DTR)*ASIN((HTAN+RADIUS)/(H(NPRO)+RADIUS))
         print*,'ZANG = ',ZANG
         PRINT*,'Enter zenith and aziumuth angle of Sun'
         CALL PROMPT('at tangent point : ')
         READ*,SOLZEN,SOLPHI        
C        N.B. SOLPHI = 0 indicates FORWARD scattering
         SOLVEC(1)=SIN(SOLZEN*DTR)*COS(SOLPHI*DTR)
         SOLVEC(2)=SIN(SOLZEN*DTR)*SIN(SOLPHI*DTR)
         SOLVEC(3)=COS(SOLZEN*DTR)
         print*,'SOLVEC : ',SOLVEC
         GAMMA = 90.0-ZANG
         CALL YROTATE(SOLVEC,GAMMA)
         print*,'Rotated SOLVEC : ',SOLVEC
         SOLZEN1 = ACOS(SOLVEC(3))/DTR
         print*,'New zenith angle : ',SOLZEN1
      ENDIF
        
      IF(IDUMP.EQ.1)THEN
        WRITE(34,*)'Tangent altitude, Zenith angle = ',HTAN,ZANG
      ENDIF

C     Read in solar flux file
      CALL PROMPT('Enter solar file name : ')
      READ(5,1)IPFILE

      CALL PROMPT('Calculate wavenumber(0) or wavelength(1) spectra : ')
      READ*,ISPACE
 
      CALL PROMPT('Enter distance from Sun (AU) : ')
      READ*,SDIST

      CALL OPENSOL(IPFILE,ISPACE,SOLNAME)
      CALL INIT_SOLAR_WAVE(SOLNAME)


C     Calculate angular radius of Sun
      SOLRAD = 0.5*0.533128/SDIST

      CALL PROMPT('Enter acceptable deviation from Sun (deg) : ')
      READ*,DEVSUN

      OMEGA1 = 2*PI*(1.0-cos(SOLRAD*DTR))
      OMEGA2 = 2*PI*(1.0-cos(DEVSUN*DTR))
C     Calculate ratio of projected areas
      XFAC= OMEGA1/OMEGA2

      CALL PROMPT('Enter max. number of photons : ')
      READ*,NITER

      IF(NITER.GT.MPHOT)THEN
       Print*,'Monteck. NITER must be less than or equal to MPHOT'
       Print*,'MPHOT = ',MPHOT
       STOP
      ENDIF

C     Read in scattering properties of dust
      CALL PROMPT('Enter name of .xsc x-section file : ')
      READ(5,1)XHGFILE

      CALL GETHG(XHGFILE,NCONT1,NLAMBDA,XLAMBDA,PHASED)
      do i = 1,200
       print*,(phased(1,i,j),j=1,5)
      enddo

      IF(NCONT1.NE.NCONT)THEN
       PRINT*,'.pha file is incompatible with dust.prf file'
       PRINT*,'NCONT1,NCONT = ',NCONT1,NCONT
       STOP
      ENDIF

      CALL PROMPT('Enter ground albedo and temperature : ')
      READ*,GALB,TGROUND


      CALL PROMPT('Enter IRAY : ')
      READ*,IRAY


      PRINT*,'Channel integrated (0) or regularly spaced'
      CALL PROMPT('k-table(1) ?')
      READ*,ICHANNEL

      IF(ICHANNEL.EQ.1)THEN
       IF(ISPACE.EQ.1)THEN
        CALL PROMPT('Enter desired wavelength range : ')
        READ*,VMIN,VMAX
      
        PRINT*,'Enter reference wavelength and'
        CALL PROMPT('step in k-tables : ') 
        READ*,VREF,DELV
       ELSE
        CALL PROMPT('Enter desired wavenumber range : ')
        READ*,VMIN,VMAX

        PRINT*,'Enter reference wavenumber and'
        CALL PROMPT('step in k-tables : ')
        READ*,VREF,DELV
       ENDIF
     
       J = INT((VMIN-VREF)/DELV)
       VMIN = VREF+J*DELV

       J = 1 + INT((VMAX-VREF)/DELV)
       VMAX = VREF+J*DELV

       Print*,'Calculation range : ',VMIN,VMAX
       NWAVE = 1 + INT((VMAX-VMIN)/DELV)
       print*,'Delv, NWAVE : ',DELV,NWAVE     

       DO IWAVE=1,NWAVE
        VWAVE(IWAVE)=VMIN+(IWAVE-1)*DELV
       ENDDO

      ELSE

       CALL PROMPT('Enter number of required wavelengths : ')
       READ*,NWAVE
       CALL PROMPT('Enter calculation wavelengths/wavenumbers : ')
       READ*,(VWAVE(IWAVE),IWAVE=1,NWAVE)

      ENDIF


      CALL PROMPT('Enter desired convergence accuracy (rad.units): ')
      READ*,ACC

      OPEN(35,FILE='monteck_spec.dat',STATUS='UNKNOWN')

      WRITE(35,*)NWAVE

      IF(IDUMP.EQ.1)THEN
       WRITE(34,*)NITER,'   ! Number of iterations'
       WRITE(34,*)NWAVE,'    ! NWAVE'
       WRITE(34,*)(VWAVE(I),I=1,NWAVE)
      ENDIF


      DO 1001 IWAVE=1,NWAVE

       VV = VWAVE(IWAVE)
       PRINT*,'VV = ',VV
       IF(IDUMP.EQ.1)THEN
        IF(ISPACE.EQ.1)THEN
         WRITE(34,*)'Wavelength : ',VV
        ELSE
         WRITE(34,*)'Wavenumber : ',VV
        ENDIF
       ENDIF

C      Get solar flux at this wavelength/wavenumber
       CALL GET_SOLAR_WAVE(VV,IWAVE,SOLAR)
       print*,'SOLAR = ',SOLAR

       SOLAR=SOLAR*XFAC

C      Interpolate k-tables to VV
       CALL GENTABSCK(NPRO,NGAS,ID,ISO,P,T,VMR,NWAVE,VWAVE,
     1 VV,IWAVE,ISPACE,NG,DEL_G,TABK)

C      Interpolate scattering properties to VV
       CALL INTERPHG(VV,NCONT,NLAMBDA,XLAMBDA,PHASED,XSEC,XOMEGA,XHG)

       print*,'xsec,xomega =',xsec(1),xomega(1),ncont

C      Regrid phase functions to equal probability steps
       NCONT1 = NCONT+1
       CALL PHASPROB(NCONT1,XHG,NPHASE,THETA)

       SVEC(1)=0
       SVEC(2)=0.0
       SVEC(3)=RADIUS+H(NPRO)

       DVEC1(1) = SIN(ZANG*PI/180.0)
       DVEC1(2) = 0.0
       DVEC1(3) = -COS(ZANG*PI/180.0)

       print*,'Initial position vector : ',SVEC
       print*,'Zenith angle : ',ZANG
       print*,'Initial direction vector : ',DVEC1     
       print*,'NG = ',NG
       CALL MCPHOTONCK(NITER,IDUM,
     1    XSEC,XOMEGA,NPHASE,THETA,
     2    SVEC,DVEC1,SOLVEC,DEVSUN,SOLAR,TABK,NG,DEL_G,
     3    NPRO,NGAS,NCONT,MOLWT,RADIUS,P,T,H,DUST,
     4    GALB,TGROUND,IRAY,RES,ACC,MEAN,SDEV,
     5    MSCAT,ITER,ISPACE,VV,NAB,NSOL,NGR)

       IF(IDUMP.EQ.1)THEN
        WRITE(34,*)ITER
        DO I=1,ITER
         WRITE(34,*)(RES(I,J),J=1,3)
        ENDDO
       ENDIF

       WRITE(35,*)VV,MEAN,SDEV,MSCAT,NAB,NSOL,NGR,ITER

1001  CONTINUE 

      IF(IDUMP.EQ.1)CLOSE(34)
      CLOSE(35)


      WRITE(*,*)'%Monteck.f :: calculation complete'
      CALL GETTIME(TIME)
      TIME2= TIME
      TOT_TIME=SNGL(TIME2-TIME1)
      WRITE(*,200)TOT_TIME
200   FORMAT(' Elapsed time (sec)= ',F8.1)

      END
