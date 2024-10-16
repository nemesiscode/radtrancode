      SUBROUTINE MCPHOTONCK(NPHOT,IDUM,XSEC,XOMEGA,NPHASE,THETA,
     1 SVEC,DVEC1,SOLVEC,DEVSUN,SOLAR,TABK,NG,DEL_G,
     3 NPRO,NGAS,NCONT,MOLWT,RADIUS,P,T,H,DUST,GALB,TGROUND,IRAY,
     4 RES,ACC,MEAN,SDEVM,MSCAT,IPHOT,ISPACE,XX,NAB,NSOL,NGR,NGS)
C     $Id: 
C     ****************************************************************
C     Subroutine to perform monte-carlo multiple scattering 
C     for a spherically symmetric planet for solar light. Subroutine
C     allows individual photons to trickle through the atmosphere,
C     scattering as they go until eventually absorbed. The planck function
C     for the temperature of the absorbing region is then added to the
C     output sum.
C
C     Code used correlated-k model.
C      
C     Input variables
C	NPHOT	INTEGER		Number of photons to fire
C	IDUM	INTEGER		Random number seed
C	XSEC(NCONT) REAL	Particle cross-sections
C	XOMEGA(NCONT) REAL	Particle single scattering albedos
C	NPHASE	INTEGER		Number of angle points in phase function
C	THETA(MAXCON,100) REAL	Scattering probability angles (compiled
C				by PHASPROB).
C	SVEC(3)	REAL		Starting position vector of photons
C	DVEC1(3) REAL		Initial direction vector of photons
C	SOLVEC(3) REAL		Direction of centre of Sun
C	DEVSUN  REAL		Acceptable deviation from Sun in degrees
C	SOLAR	REAL		Radiance of solar surface
C	TABK(MAXG,MAXPRO) REAL	overlapped k-coefficents at each level in
C				the atmosphere. 
C	NG	INTEGER		Number of g-ordinates in k-distr.	
C	DEL_G(MAXG) INTEGER	k-distr. weights
C	NPRO	INTEGER		Number of levels in profile
C	NGAS	INTEGER		Number of gases 
C	NCONT	INTEGER		Number of aerosol types (terminates if
C				NCONT>1)
C	MOLWT	REAL		Molecular weight of atmosphere
C	RADIUS	REAL		Planetary radius (at 0km altitude reference level)
C	P(MAXPRO) REAL		Pressure profile
C	T(MAXPRO) REAL		Temperature profile
C	H(MAXPRO) REAL		height profile
C	DUST(MAXPRO,MAXCON) REAL	dust profiles
C	GALB	REAL		Ground albedo
C	TGROUND REAL		Ground temperature (K)
C	IRAY	INTEGER		0=Rayleigh scattering off, 1=on.
C	ACC	REAL		Required absolute accuracy
C	ISPACE 	INTEGER		0=wavelength, 1=wavenumber
C	XX	REAL		Required wavelength/wavenumber
C	
C		
C     Output variables
C	RES(MPHOT,3)	REAL	Outcome of each photon entering atmosphere
C				RES(*,1) holds final destination:
C				  1=ground
C				  2=space
C				  3=atmosphere
C				RES(*,2) holds radiance of destination
C				RES(*,3) records number of scattering events
C	MEAN	REAL		Mean scattered radiance
C	SDEVM	REAL		Standard deviation of radiance
C	MSCAT	REAL		Mean number of scattering events
C	IPHOT	INTEGER		Number of photons actually fired
C	NAB	INTEGER		Number of photons absorbed in atm
C	NSOL	INTEGER		Number of photons encountering Sun
C	NGR	INTEGER		Number of photons absorbed by ground
C	NGS	INTEGER		Number of photons reaching the ground
C
C     Pat Irwin    	Original	13/11/01
C			Revised		22/7/05
C
C     ****************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NPRO,NGAS,NCONT,I,J,MPHOT,ICONT
      PARAMETER (MPHOT = 100000)
      REAL P(MAXPRO),T(MAXPRO),H(MAXPRO),TRANTOT
      REAL RADIUS,MOLWT,DUST(MAXPRO,MAXCON)
      REAL SVEC(3),DVEC(3),PVEC(3)
      REAL RAN11,ALTITUDE,DVEC1(3),SOLVEC(3),DEVSUN,SOLAR
      REAL XOMEGA(MAXCON),XSEC(MAXCON),CALCALT,PI,TRANSREQ,FSCAT
      REAL ALT0,TABK(MAXG,MAXPRO),TAUSCAT(MAXCON),TAUTOT(MAXCON+1)
      REAL X,MSCAT,ACC,MEAN,MEAN2,ACC1,VV,PLANCK_WAVE
      REAL DEL_G(MAXG),G_ORD(MAXG),X1,GALB,TGROUND,TMEAN,SDEV
      REAL SUM,SDEVM,TNOW,XX,SCOS,ARCTAN,DTR,GALB1
      PARAMETER (PI=3.1415927,DTR=PI/180.)
      
      INTEGER IPHOT,NPHOT,NPHASE,IDUM,NSCAT,IHIT
      INTEGER IGDIST,NG,CIGDIST(MAXG),ISPACE,IRAY,NCONT1
      INTEGER NAB,NSOL,NGR,NGS
      REAL THETA(MAXCON,100),ALPHA,PHI,THETA1(100)
      REAL RES(MPHOT,3),THET
      CHARACTER*1 ANS

      MEAN=0.
      SDEV=0.
      MSCAT=0.
      IPHOT=0
      NAB=0
      NSOL=0
      NGR=0
      NGS=0
C      print*,NPHOT,IDUM,NPRO,NGAS,NCONT,MOLWT
C      print*,XSEC(1),XOMEGA(1),NPHASE,'NPHASE'
C      print*,(THETA(1,J),J=1,NPHASE)
      print*,'SVEC',SVEC
      print*,'DVEC = ',DVEC1
      print*,'SOLVEC = ',SOLVEC
      print*,'DEVSUN',DEVSUN
C      print*,SOLAR
C      print*,NG,'NG'
C      do i=1,npro
C       print*,(TABK(j,i),j=1,NG)
C      enddo
C      print*,(DEL_G(J),J=1,NG)
C      do i=1,npro
C       print*,H(i),P(i),T(i),(DUST(i,j),j=1,NCONT)
C      enddo
C      print*,RADIUS,TGROUND,IRAY,'IRAY'
C      print*,ACC,MEAN,SDEVM,MSCAT,IPHOT,ISPACE,XX,NAB,NSOL,NGR,NGS
C      print*,GALB

      IF(NPHOT.GT.MPHOT)THEN
       PRINT*,'Error in MCPHOTONCK NPHOT > MPHOT'
       PRINT*,NPHOT,MPHOT
       STOP
      ENDIF

      SCOS = COS(DEVSUN*PI/180.0)
C      print*,'SCOS, DEVSUN',SCOS,DEVSUN
      IF(IRAY.GT.0)THEN
       NCONT1 = NCONT+1
      ELSE
       NCONT1 = NCONT
      ENDIF

      G_ORD(1)=0.0
      DO I=1,NG  
       G_ORD(I+1)=G_ORD(I)+DEL_G(I)
       CIGDIST(I)=0
      ENDDO
      ALT0 = H(NPRO)

      MEAN = 0.0
      MSCAT = 0.0

      open(37,file='test.dat',status='unknown')

      GALB1=1.0
      DO 1000 IPHOT=1,NPHOT 
C       PRINT*,'IPHOT,NPHOT',IPHOT,NPHOT
C      Select IGDIST by size of g_interval (amounts to implicit
C               integration in g-space)
       X1 = RAN11(IDUM)
C       print*,'X1,IDUM = ',X1,IDUM
       DO I=1,NG
        IF(X1.GE.G_ORD(I).AND.X1.LT.G_ORD(I+1))IGDIST=I
       ENDDO
       IF(X1.EQ.1.0)IGDIST=NG

C       print*,'IGDIST = ',IGDIST
C      Add one to count of IGDIST distribution (to make sure g-space is uniformly
C      sampled.
       CIGDIST(IGDIST)=CIGDIST(IGDIST)+1

       NSCAT = 0
       DO I=1,3
         PVEC(I) = SVEC(I)
         DVEC(I) = DVEC1(I)
         RES(IPHOT,I)=0.0
       ENDDO
       ALTITUDE = CALCALT(PVEC,RADIUS)
C       PRINT*,'Initial altitude, solar = ',ALTITUDE,SOLAR
C       PRINT*,'Initial theta,phi = ',ACOS(DVEC(3))/DTR,
C     1  ARCTAN(DVEC(2),DVEC(1))/DTR

C      Calculate random transmission to pass
101    TRANSREQ = RAN11(IDUM)

C      Set VV to current WAVENUMBER and pass to intpath
       IF(ISPACE.EQ.1)THEN
        VV = 1E4/XX
       ELSE
        VV=XX
       ENDIF

C       print*,'PVEC',(PVEC(I),I=1,3)
C       print*,'DVEC',(DVEC(I),I=1,3)

       CALL INTPATH(VV,IRAY,TRANSREQ,PVEC,DVEC,NPRO,NGAS,NCONT,MOLWT,
     1 RADIUS,P,T,H,DUST,TABK,IGDIST,XSEC,XOMEGA,TRANTOT,TMEAN,FSCAT,
     2 TAUSCAT)
  
       ALTITUDE = CALCALT(PVEC,RADIUS)

C       print*,'ALTITUDE, H(1), TRANTOT',ALTITUDE,H(1),TRANTOT
    
C       PRINT*,'ALTITUDE,TRANSREQ, FSCAT = ',ALTITUDE,TRANSREQ,FSCAT
C       PRINT*,'Current theta,phi = ',ACOS(DVEC(3))/DTR,
C     1  ARCTAN(DVEC(2),DVEC(1))/DTR


       IF(ALTITUDE.LE.H(1))THEN
C          PRINT*,'Photon hits surface of reflectance',GALB
          NGS=NGS+1

          IF(RAN11(IDUM).LE.GALB)THEN
C           PRINT*,'Photon is reflected (LAMBERT) from surface',GALB
           SUM=0.0
           DO I=1,3
            SUM=SUM+PVEC(I)**2
           ENDDO
           SUM=SQRT(SUM)


C          Calculate the angle that the photon will go off into if
C          reflected from a Lambertian surface
C          First, PHI can be anywhere between 0 and 360 degrees
           PHI = 2*PI*RAN11(IDUM)
           
C          Second, the zenith angle. Flux should be equal in terms
C          of solid angle, which depends as 
C          dOmeg=sin(thet).dthet.dphi=-d(cos(thet)).dphi
           THET = ACOS(RAN11(IDUM))

C          Less photons are arriving per unit area at higher reflection 
C          angles, so incorporate this into effective reflectivity correction
           GALB1=COS(THET)

C          Also Correct for surface albedo radiance/irradiance
           GALB1=PI*GALB1

C          Reset length of PVEC (if HEIGHT < H(1))
           DO I=1,3
            PVEC(I)=PVEC(I)*(RADIUS+H(1))/SUM
           ENDDO
 

C          Reset the photon direction vector to be straight up
C          from surface 
           DO I=1,3
            DVEC(I)=PVEC(I)/(RADIUS+H(1))
           ENDDO

C          Modify the direction vector with the new THET and PHI
C           print*,'THET, PHI = ',THET,PHI
           CALL MODVEC(DVEC,THET,PHI)

           GOTO 101

          ENDIF 
C          print*,'Photon is absorbed'
          RES(IPHOT,1)=1
          RES(IPHOT,2)=PLANCK_WAVE(ISPACE,XX,TGROUND)
          NGR=NGR+1
          GOTO 999

       ENDIF
  
       IF(ALTITUDE.GT.ALT0)THEN
C          PRINT*,'Photon leaves atmosphere',DVEC

          RES(IPHOT,1)=2
          CALL HITSUN(SOLVEC,DVEC,SCOS,IHIT)
C          print*,'----------------'
C          print*,(SOLVEC(I),I=1,3)
C          print*,(DVEC(I),I=1,3)
C          print*,SCOS,IHIT
C          PRINT*,SOLVEC,DVEC,IHIT
           write(37,*)(SOLVEC(I),I=1,3),(DVEC(J),J=1,3),IHIT
          IF(IHIT.EQ.0)THEN
           RES(IPHOT,2)=0.0
          ELSE
C           PRINT*,'photon comes near the Sun'

C

           RES(IPHOT,2)=SOLAR*GALB1/PI
           NSOL=NSOL+1
          ENDIF
          GOTO 999

       ENDIF

C      COMPUTE WHETHER BEAM WILL BE ABSORBED       
       IF(RAN11(IDUM).GT.FSCAT)THEN
C         print*,'photon absorbed in atmosphere'

         RES(IPHOT,1)=3
         RES(IPHOT,2)=PLANCK_WAVE(ISPACE,XX,TMEAN)
         NAB=NAB+1
         GOTO 999

       ELSE

C        print*,'Photon scattered in atmosphere'
        NSCAT = NSCAT+1

        TAUTOT(1)=TAUSCAT(1)
        DO I=2,NCONT1
         TAUTOT(I)=TAUTOT(I-1)+TAUSCAT(I)
        ENDDO
        X1=RAN11(IDUM)
        ICONT=1
        J=1
        DO I=1,NCONT1
         TAUTOT(I)=TAUTOT(I)/TAUTOT(NCONT1)
         IF(X1.LT.TAUTOT(I).AND.J.GT.0)THEN
           ICONT=I
           J=-1
         ENDIF
        ENDDO

C        print*,'Particle type : ',ICONT,' of ',NCONT1
        DO I=1,NPHASE
         THETA1(I)=THETA(ICONT,I)
        ENDDO

        CALL DEFLECT(IDUM,THETA1,NPHASE,ALPHA,PHI)
      
        CALL MODVEC(DVEC,ALPHA,PHI)
  
       ENDIF

       GOTO 101

999    CONTINUE

       RES(IPHOT,3)=NSCAT
       X = RES(IPHOT,2)

C       print*,'RES : ',(RES(IPHOT,J),J=1,3)  
       MEAN = ((IPHOT-1)*MEAN + X)/FLOAT(IPHOT)
       MEAN2 = ((IPHOT-1)*MEAN2 + X*X)/FLOAT(IPHOT)
       MSCAT = ((IPHOT-1)*MSCAT + RES(IPHOT,3))/FLOAT(IPHOT)

       SDEV = SQRT(ABS(MEAN2 - MEAN**2))       

       SDEVM = SDEV/SQRT(FLOAT(IPHOT))

C      Abort if already sufficiently converged after 200 photons
c       at least 50 of which have hit the Sun.
C       IF(IPHOT.GT.200.AND.NSOL.GT.50.AND.SDEVM.LE.ACC) GOTO 1001      

1000  CONTINUE

1001  CONTINUE

C     IPHOT is now actually a loop, rather than a counter, but reset 
C     to NPHOT if necessary
      IF(IPHOT.GT.NPHOT)IPHOT=NPHOT

      close(37)

      RETURN

      END
