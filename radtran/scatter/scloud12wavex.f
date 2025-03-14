      subroutine SCLOUD12WAVEX( rad, sol_ang, emiss_ang, aphi, radg, 
     1 solar, lowbc, galb, iray, mu1, wt1, nmu, nf, igdist, vwave, vv,
     2 bnu, tauc, taus, icloud, fcover, taug, tauray, nlay, ncont)
C     $Id:
C     ***********************************************************************
C     Compute emergent intensity at top of multilayer cloud using the
C     matrix operator algorithm.  Diffuse incident radiation is allowed
C     at the bottom and single-beam incident radiation (sunlight) at
C     the top. 
C
C     The layer numbers here are 1 at the top increasing to NLAY at the 
C     bottom. If a Lambertian reflector is added at the bottom then the 
C     number of layers is increased by 1.
C
C     NOTE:  the angle arrays passed from the calling pgm. are assumed to be 
C     in order of increasing MU.  However, the CLOUD subroutines use the 
C     opposite convention, so we must reverse the order.  (The order has no 
C     effect within this routine, other than in using the supplied boundary 
C     conditions, but it does affect pre-/post-processing in the calling s/w.)
C
C     Optimised for maximum speed by cutting out any code from scloud6
C     which is not actually used for NIMS retrieval runs.
C
C     Modified from scloud10 to now include rayleigh air scattering
C
C     Modified from scloud11wave to include broken clouds.
C
C	Input variables
C       radg(maxmu)     real  Incident intensity at the bottom of the atm
C	sol_ang		real  solar zenith angle (degrees)
C	emiss_ang	real  emission zenith angle (degrees)
C       solar           real  Incident solar flux at the top of the atm
C       aphi            real  Azimuth angle between Sun and observer
C       lowbc           integer Lower boundary condition:
C                               0 = thermal, 1 = Lambert reflection.
C       galb            double  Ground albedo at the bottom.
C       mu1(maxmu)      double  zenith angle point quadrature
C       wt1(maxmu)      double  zenith angle point quadrature
C       nmu             integer Number of points in quadrature
C       nf              integer Number of terms in azimuth Fourier expansion   
C	igdist          integer Switch to recalculate phase functions. If > 1
C                               then phase functions left alone
C       bnu(MAXSCATLAY)     real  Mean Planck function in each layer
C       tauc(MAXCON,MAXSCATLAY) real  Total extinction optical thickness of
C				each layer of each cloud type
C       taus(MAXCON,MAXSCATLAY) real  Total scattering optical thickness of
C				each layer of each cloud type
C	icloud(MAXCON,MAXSCATLAY) integer Array with 0 if particle is part
C				of haze and 1 if particle part of
C				broken cloud
C	taug(MAXSCATLAY)	real	Total optical thickness of each layer of 
C				gas continuum, line absorption
C	tauray(MAXSCATLAY)  real  Total optical thickness of each layer due
C                               to Rayleigh scattering
C       fcover(MAXSCATLAY)  real  Fractional cloud cover of each layer
C       nlay            integer Number of layers
C       ncont           integer Number of aerosol types/modes
C       liscat(MAXCON)   integer Phase function type identifier
C       lncons(MAXCON)   integer Number of phase function parameters
C       lnorm(MAXCON)    integer Normalisation identifier
C       nphi            integer Number of azimuth integration ordinates
C
C	Output variables
C	rad		real	Emergent radiance 
C
C	New History
C
C	3/11/94	PGJI	Adapted for RADTRAN
C	7/1/97  PGJI   	Revised again to pass out all fourier coefficients
C	7/5/97	PGJI	Optimised from Scloud6 for maximum speed
C
C     ***********************************************************************

      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'

      INTEGER LTOT,NLAY,LT1,LOWBC,NMU,I,J,K
      INTEGER IC,I1,J1,NCONT,NCONS,NORM,ISCAT,KL,NPHI,IP,JP
      INTEGER IPOW0,L,IL,IMU0,IDUMP,NCONT1
      double precision mu1(maxmu), wt1(maxmu), PI
      real rad, radf(0:MAXF-1), radg(maxmu), bnu(MAXSCATLAY), 
     1 tauc(MAXCON,MAXLAY), taus(MAXCON,MAXLAY), fcover(MAXLAY),
     2 vwave, solar, vv, xf, lfrac(MAXCON,MAXLAY),aphi,drad,solar1,
     3 taug(MAXLAY),conv1,tauray(MAXLAY)
      integer icloud(MAXCON,MAXLAY)
      real defconv,conv
      integer nf,imu,iray
      integer igdist,imie
      logical lrep,checkconv
      common /imiescat/imie


      DOUBLE PRECISION MU(MAXMU), WTMU(MAXMU), MM(MAXMU,MAXMU),
     1 MMINV(MAXMU,MAXMU), RL(MAXMU,MAXMU,MAXSCATLAY),EPS(MAXSCATLAY),
     2 TL(MAXMU,MAXMU,MAXSCATLAY), JL(MAXMU,1,MAXSCATLAY),
     3 RBASE(MAXMU,MAXMU,MAXSCATLAY),F,G1,G2,
     4 TBASE(MAXMU,MAXMU,MAXSCATLAY), JBASE(MAXMU,1,MAXSCATLAY),
     5 UMI(MAXMU,1,MAXSCATLAY), 
     6 ACOM(MAXMU,1), BCOM(MAXMU,1),
     7 PPLPL(MAXMU,MAXMU), PPLMI(MAXMU,MAXMU), CONS8(MAXSCATPAR),
     8 PPLPLS(MAXMU,MAXMU),PPLMIS(MAXMU,MAXMU),
     9 PPLPLR(MAXMU,MAXMU),PPLMIR(MAXMU,MAXMU)

      DOUBLE PRECISION RLF(MAXMU,MAXMU,MAXSCATLAY),
     1 TLF(MAXMU,MAXMU,MAXSCATLAY), JLF(MAXMU,1,MAXSCATLAY)

      DOUBLE PRECISION PPLN(MAXCON,0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION PMIN(MAXCON,0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION PPLR(0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION PMIR(0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION E(MAXMU,MAXMU),CCINV(MAXMU,MAXMU),CC(MAXMU,MAXMU)
      DOUBLE PRECISION U0PL(MAXMU,1),UTMI(MAXMU,1),TAUT,BC,OMEGA
      DOUBLE PRECISION FRAC,SFRAC
      DOUBLE PRECISION STEST,GALB,ZMU0,ZMU,TEX, TAUSCAT, TAUR, TAU1 
      DOUBLE PRECISION STEST1
      DOUBLE PRECISION RTOP(MAXMU,MAXMU,MAXSCATLAY),
     1 UPL(MAXMU,1,MAXSCATLAY),JTOP(MAXMU,1,MAXSCATLAY),
     2 UTPL(MAXMU,1,MAXSCATLAY),TTOP(MAXMU,MAXMU,MAXSCATLAY),
     3 U0MI(MAXMU,1,MAXSCATLAY)
      REAL UMIF(MAXMU,MAXSCATLAY),UPLF(MAXMU,MAXSCATLAY)
      REAL FPL,FMI,IRAD(MAXSCATLAY),F1,F2

      DOUBLE PRECISION OMEGA1,TAUTG,TAUTF,OMEGAF,TAUSCATF,OMEGAT

      REAL SOL_ANG,EMISS_ANG
      REAL YX(4),T,U,FEMM,FSOL,XFAC
      INTEGER ICO,ISOL,IEMM,ISCL(MAXSCATLAY)
      LOGICAL RAMAN
      REAL DENS(MAXSCATLAY),FPRAMAN(MAXSCATLAY)
      REAL JRAMAN(1000,MAXSCATLAY),VRAM0,VRAMST
      REAL LDENS1,LDENS2,SUM1,SUM2,LD,FF,XD,FMEAN(MAXSCATLAY)
      REAL H2ABUND(MAXSCATLAY),KWT
      INTEGER IRAMAN,ILAY,NRAMAN

C     Common blocks
      COMMON/UNIT/ E
      COMMON/AREA1/CCINV,MM
      COMMON/AREA2/CC,MMINV,MU,PI
      COMMON/INPUT/U0PL,UTMI
      COMMON/HOMOG/ TAUT, BC, OMEGA, IPOW0
      COMMON/PHMAT/ PPLPL, PPLMI
      COMMON/NEWPH/PPLN,PMIN,PPLR,PMIR
      COMMON /SCATDUMP/ IDUMP
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet
 
      COMMON/RAMAN/RAMAN,DENS,JRAMAN,VRAM0,VRAMST,NRAMAN,FPRAMAN,
     1 H2ABUND,IRAMAN,KWT

C--------------------------------------------------------------------------

C     Define additional debugging 'reporting variable', lrep
C      lrep=.true.
      lrep=.false.
      if(idump.gt.0.or.(lrep.and.igdist.eq.1))then
      print*,'SCLOUD12WAVEX : sol_ang',sol_ang, 'emiss_ang',emiss_ang,
     1 'aphi', aphi, 'solar', solar, 'lowbc',lowbc, 'galb',galb,
     2 'nmu',nmu, 'nf',nf, 'nlay',nlay, 'ncont',ncont,
     3 'igdist',igdist
      print*,'vwave : ',vwave,'vv : ',vv
      print*,(mu1(i),i=1,nmu)
      print*,(wt1(i),i=1,nmu)
      print*, (radg(i),i=1,nmu) 
      print*,'i, tauc(1:ncont),taus(1:ncont),taug,bnu'
      do i=1,nlay
       print*,i,(tauc(j,i),j=1,ncont),(taus(j,i),j=1,ncont),taug(i),
     &  bnu(i)
      end do
      print*,'i,tauray,fcover,icloud(1:ncont)'
      do i=1,nlay
       print*,i,tauray(i),fcover(i),(icloud(j,i),j=1,ncont)
      end do
      endif

      IF(RAMAN)THEN
C      If Raman scattering turned on, find location in Raman scattered
C      source function array, JRAMAN, of any scattered radiation to add to
C      the scattered radiation at this wavelength.
C      Raman-scattered light is assumed to be isotropic and is added to
C      the calculation via double1.f as an additional emission term.
       IRAMAN=1+INT((1E4/VV-VRAM0)/VRAMST)
      ENDIF

C     Find correction for any quadrature errors
      xfac=0.
      do i=1,nmu
       xfac=xfac+sngl(mu1(i)*wt1(i))
      enddo
      xfac=0.5/xfac

      LTOT = NLAY	! Set internal number of layers
      LT1 = LTOT

C     In case of Lambertian reflection, add extra dummy layer at bottom,
C     whose transmission matrix = (1-A)*Unit-Matrix. This layer must be
C     omitted from computation by doubling

      IF (LOWBC.EQ.1) LTOT = LTOT+1
      IF (LTOT.GT.MAXSCATLAY) THEN
       print*, 'ltot ', ltot, ' MAXSCATLAY ', MAXSCATLAY
       CALL ABEND(' SCLOUD12WAVEX: TOO MANY LAYERS')
      ENDIF
      IF (NMU.GT.MAXMU) THEN
        CALL ABEND(' SCLOUD12WAVEX: TOO MANY ANGLE POINTS')
      ENDIF


C     Reset the order of angles:
      DO I=1,NMU
	MU(I) = MU1(NMU+1-I)
	WTMU(I) = WT1(NMU+1-I)
      ENDDO

      PI = 4.0D0*DATAN(1.0D0)

C ********************************
C     SET UP CONSTANT MATRICES
C ********************************

      DO J = 1,NMU
 	DO I = 1,NMU
	  E(I,J) = 0.0D0
	  IF(I.EQ.J) E(I,J) = 1.0D0
	  MM(I,J) = 0.0D0
 	  IF(I.EQ.J) MM(I,J) = MU(I)
	  MMINV(I,J) = 0.0D0
  	  IF(I.EQ.J) MMINV(I,J) = 1.0D0/MU(I)
	  CC(I,J) = 0.0D0
 	  IF(I.EQ.J) CC(I,J) = WTMU(I)
 	  CCINV(I,J) = 0.0D0
 	  IF(I.EQ.J) CCINV(I,J) = 1.0D0/WTMU(I)
	ENDDO
      ENDDO

C     PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C     Precalculate phase function arrays if new wavelength
      IF(IGDIST.EQ.1) THEN
       DO J1=1,NCONT
        CALL READ_HG(VWAVE,J1,NCONT,F,G1,G2)

C       If imie=1 then read in phase function from PHASEN.DAT files
C       otherwise read in from hgphaseN.dat files

        if(imie.eq.1)then
          ISCAT=4
          NCONS=0
        else
          CALL READ_HG(VWAVE,J1,NCONT,F,G1,G2)
          if(idump.gt.0.or.lrep)print*,'f,g1,g2 = ',
     1      f,g1,g2

          ISCAT =2
          NCONS = 3
          CONS8(1)=F
          CONS8(2)=G1
          CONS8(3)=G2
        endif

        NORM = 1
        NPHI = 100

        DO 900 IC=0,NF
         CALL CALC_PMAT6(NF, IC, PPLPL, PPLMI, MU, WTMU,
     1    NMU, ISCAT, CONS8, NCONS, NORM, J1, NCONT, VWAVE, NPHI)

C        Transfer matrices to those for each scattering particle
         DO JP = 1,NMU
          DO IP = 1,NMU
           PPLN(J1,IC,IP,JP)=PPLPL(IP,JP)
           PMIN(J1,IC,IP,JP)=PPLMI(IP,JP)
          END DO
         END DO
900     CONTINUE
       ENDDO



       IF(IRAY.GT.0)THEN
C       Rayleigh scattering included so we need to calculate
C       the Rayleigh scattering phase function

        ISCAT = 0
        NCONS = 0
        NPHI = 100
        NCONT1=NCONT+1
        J1=NCONT+1
        NORM=1
        DO 800 IC=0,NF
         CALL CALC_PMAT6(NF, IC, PPLPL, PPLMI, MU, WTMU,
     1    NMU, ISCAT, CONS8, NCONS, NORM, J1, NCONT1, VWAVE, NPHI)

C        Transfer matrices to storage
         DO JP = 1,NMU
          DO IP = 1,NMU
           PPLR(IC,IP,JP)=PPLPL(IP,JP)
           PMIR(IC,IP,JP)=PPLMI(IP,JP)
          END DO
         END DO
C         print*,IC,PPLR(IC,1,1),PMIR(IC,1,1)
800     CONTINUE
       ENDIF
C     PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
      ENDIF



      RAD=0.0
      CONV = 100.0
      CONV1 = 100.0

C     ******* Define convergence (max % change in last two iterations
C     Set very small now to make sure we never apply this.
      DEFCONV=1e-20

C     Set convergence to be a little less draconian: 0.0001%
C     Should be fine for back-scattering conditions.
C      DEFCONV=0.0001

C     Flag for printing out convergence information
      checkconv=.false.
C      checkconv=.true.

      DO 1000 IC=0,NF


C **********************************************************************
C     CALCULATE REFLECTION,TRANSMISSION AND SOURCE MATRICES FOR EACH    
C	  HOMOGENEOUS LAYER L; INSERT INTO ARRAY X(I,J,L).	  	      
C **********************************************************************    
C
      IPOW0 = 16
      DO 2000 L = 1,LT1
        if(idump.ne.0)print*,'L,IGDIST,PPLN(1,0,1,1) =',L,IGDIST,
     &		PPLN(1,0,1,1)
        TAUT=TAUG(L)
        TAUSCAT=0.0
        TAUTF=TAUG(L)
        TAUSCATF=0.0
C       For particles which are part of extended haze icloud = 0.
C       For particles which are part of broken cloud, icloud = 1.

        I1=0
        DO 2005 J = 1, NCONT
         IF(ICLOUD(J,L).EQ.0)THEN
          TAUT = TAUT + 1.0D0*TAUC(J,L)
          TAUSCAT = TAUSCAT + 1.0D0*TAUS(J,L)
          TAUTF = TAUTF + 1.0D0*TAUC(J,L)
          TAUSCATF = TAUSCATF + 1.0D0*TAUS(J,L)
         ELSE
          TAUTF = TAUTF + 1.0D0*TAUC(J,L)
          TAUSCATF = TAUSCATF + 1.0D0*TAUS(J,L)
          I1 = I1+1
         ENDIF
2005    CONTINUE
	BC = 1.0D0*BNU(L)
        XF=FCOVER(L)
        if(idump.ne.0)print*,'L, FCOVER',L,XF

        IF(IRAY.GT.0)THEN
		TAUR = TAURAY(L)
		TAUSCAT=TAUSCAT+TAUR
		TAUT=TAUT+TAUR
		TAUSCATF=TAUSCATF+TAUR
		TAUTF=TAUTF+TAUR
        ENDIF
        OMEGA=TAUSCAT/TAUT
        OMEGAF=TAUSCATF/TAUTF



        if(idump.ne.0)then
         print*,'OMEGA, OMEGAF = ',OMEGA,OMEGAF
         print*,TAUT,TAUSCAT
         print*,TAUTF,TAUSCATF
        endif
        IF(OMEGA.GT.1.0)THEN
         if(idiag.gt.0)print*,'Omega too big! Reducing to 1.0',OMEGA
         OMEGA = 1.D0
        ELSE IF(OMEGA.LT.0.0)THEN
         if(idiag.gt.0)print*,'Omega too small! Setting to 0.0',OMEGA
         OMEGA = 0.D0
        ENDIF
        IF(OMEGAF.GT.1.0)THEN
         if(idiag.gt.0)print*,'OmegaF too big! Reducing to 1.0',OMEGAF
         OMEGAF = 1.D0
        ELSE IF(OMEGAF.LT.0.0)THEN
         if(idiag.gt.0)print*,'OmegaF too small! Setting to 0.0',OMEGAF
         OMEGAF = 0.D0
        ENDIF


C       P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P
C
C       Assign phase function coefficients for each atmospheric layer

        DO KL=1,NMU
          DO J1=1,NMU
           PPLPLS(J1,KL)=0.
           PPLMIS(J1,KL)=0.
           PPLPLR(J1,KL)=0.
           PPLMIR(J1,KL)=0.
          END DO
        END DO
        SFRAC=0.D0

        DO J1 = 1,NCONT

C         Transfer phase-function array for current particle type
          DO JP = 1,NMU
             DO IP = 1,NMU
              PPLPL(IP,JP) = PPLN(J1,IC,IP,JP)
              PPLMI(IP,JP) = PMIN(J1,IC,IP,JP)
             END DO
          END DO

          IF(ICLOUD(J1,L).EQ.0)THEN
C          Horizontally homogeneous haze. Add proportion of
C          phase function to both 'clear' and 'cloudy' areas

           if(idump.ne.0)print*,'homogeneous',J1
   
           FRAC=TAUS(J1,L)/TAUSCAT

C           IF(L.EQ.15)THEN
C            print*,J1,PPLPL(1,1),PPLMI(1,1)
C           ENDIF

           IF(FRAC.GT.0.0)THEN
            if(idump.gt.0)print*,'Layer ',L,' Adding fraction ',
     &		FRAC,' of mode ',J1
            CALL MADD(FRAC,PPLPLR,PPLPL,PPLPLR,NMU,NMU,
     1       MAXMU,MAXMU)
            CALL MADD(FRAC,PPLMIR,PPLMI,PPLMIR,NMU,NMU,
     1       MAXMU,MAXMU)
           ENDIF

C           IF(L.EQ.15)THEN
C            print*,J1,FRAC,PPLPLR(1,1),PPLMIR(1,1)
C           ENDIF

           FRAC=TAUS(J1,L)/TAUSCATF

           IF(FRAC.GT.0.0)THEN
            if(idump.gt.0)print*,'Layer ',L,' Adding fraction ',
     &		FRAC,' of mode ',J1
            CALL MADD(FRAC,PPLPLS,PPLPL,PPLPLS,NMU,NMU,
     1       MAXMU,MAXMU)
            CALL MADD(FRAC,PPLMIS,PPLMI,PPLMIS,NMU,NMU,
     1       MAXMU,MAXMU)

           ENDIF
C           IF(L.EQ.15)THEN
C            print*,J1,FRAC,PPLPLS(1,1),PPLMIS(1,1)
C           ENDIF

          ELSE
C          Broken cloud. Add proportion of scattering to
C          'cloudy' areas only
           if(idump.ne.0)print*,'Broken ',J1
           FRAC=TAUS(J1,L)/TAUSCATF
           IF(FRAC.GT.0.0)THEN
            if(idump.gt.0)print*,'Layer ',L,' Adding fraction ',
     &		FRAC,' of mode ',J1
            CALL MADD(FRAC,PPLPLS,PPLPL,PPLPLS,NMU,NMU,
     1       MAXMU,MAXMU)
            CALL MADD(FRAC,PPLMIS,PPLMI,PPLMIS,NMU,NMU,
     1       MAXMU,MAXMU)
           ENDIF
          ENDIF
        ENDDO


        IF(IRAY.GT.0)THEN
C         Add Rayleigh scattering from air itself to both
C         'cloudy' and 'clear' areas
           DO JP = 1,NMU
             DO IP = 1,NMU
              PPLPL(IP,JP) = PPLR(IC,IP,JP)
              PPLMI(IP,JP) = PMIR(IC,IP,JP)
             END DO
           END DO

           FRAC = TAUR/TAUSCAT

C           IF(L.EQ.15)THEN
C            print*,'RAY',TAUR,TAUSCAT,FRAC,PPLPL(1,1),
C     &       PPLMI(1,1)
C           ENDIF

           if(idump.ne.0)print*,'Rayleigh fraction ',FRAC
           IF(FRAC.GT.0.0)THEN
            CALL MADD(FRAC,PPLPLR,PPLPL,PPLPLR,NMU,NMU,
     1       MAXMU,MAXMU)
            CALL MADD(FRAC,PPLMIR,PPLMI,PPLMIR,NMU,NMU,
     1       MAXMU,MAXMU)
           ENDIF
C           IF(L.EQ.15)THEN
C            print*,J1,FRAC,PPLPLR(1,1),PPLMIR(1,1)
C           ENDIF


           FRAC = TAUR/TAUSCATF
           if(idump.ne.0)print*,'Rayleigh fraction ',FRAC
           IF(FRAC.GT.0.0)THEN
            CALL MADD(FRAC,PPLPLS,PPLPL,PPLPLS,NMU,NMU,
     1       MAXMU,MAXMU)
            CALL MADD(FRAC,PPLMIS,PPLMI,PPLMIS,NMU,NMU,
     1       MAXMU,MAXMU)
           ENDIF
C           IF(L.EQ.15)THEN
C            print*,J1,FRAC,PPLPLS(1,1),PPLMIS(1,1)
C           ENDIF

        ENDIF 



C       P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P

C       Need to calculate 2 sets of R,T and J matrices, one with
C       just the extended hazes and Rayleigh scattering and one
C       with the broken clouds as well.

C       Set layer scattering flag to zero (default)
        ISCL(L)=0

C       First do the 'clear areas'
        IF(TAUT.EQ.0)THEN
         do i1=1,nmu
          do j1=1,nmu
           RL(I1,J1,L)=0.0
           TL(I1,J1,L)=0.0
          end do
          JL(i1,1,L) =0.0
          TL(i1,i1,L)=1.0
         end do
         if(idump.gt.0)print*,'zero layer ',l,iscl(l)
        ELSE IF(OMEGA.EQ.0)THEN
         do i1=1,nmu
          do j1=1,nmu
           RL(I1,J1,L)=0.0
           TL(I1,J1,L)=0.0
          end do
          TEX = -MMINV(i1,i1)*TAUT
          if(idump.gt.0)print*,i1,MMINV(i1,i1),TAUT
          if(idump.gt.0)print*,TEX
          if(TEX.GT.-200.0D0)THEN
            TL(i1,i1,L)=DEXP(TEX)
          ELSE
            TL(i1,i1,L)=0.D0
          ENDIF
          if(idump.gt.0)print*,TL(i1,i1,L)
          JL(i1,1,L)= BC*(1.0 - TL(i1,i1,L))
         end do
         if(idump.gt.0)print*,'gas layer ',l,iscl(l)
        ELSE
         ISCL(L) = 1
         if(idump.gt.0)print*,'cloud layer ',l,iscl(l)

         DO KL=1,NMU
           DO J1=1,NMU
            PPLPL(J1,KL)=PPLPLR(J1,KL)
            PPLMI(J1,KL)=PPLMIR(J1,KL)
           END DO
         END DO

 	 CALL DOUBLE1(IC,L,RL(1,1,L),TL(1,1,L),JL(1,1,L),NMU,MAXMU)

C         IF(L.EQ.15)THEN
C          print*,L,TAUT,OMEGA,RL(1,1,L),TL(1,1,L),JL(1,1,L)
C         ENDIF

        ENDIF

C       Now do the 'cloudy' areas
        OMEGAT = OMEGAF
        IF(TAUTF.EQ.0)THEN
         do i1=1,nmu
          do j1=1,nmu
           RLF(I1,J1,L)=0.0
           TLF(I1,J1,L)=0.0
          end do
          JLF(i1,1,L) =0.0
          TLF(i1,i1,L)=1.0
         end do
         if(idump.gt.0)print*,'zero layer ',l,iscl(l)
        ELSE IF(OMEGAT.EQ.0)THEN
         do i1=1,nmu
          do j1=1,nmu
           RLF(I1,J1,L)=0.0
           TLF(I1,J1,L)=0.0
          end do
          TEX = -MMINV(i1,i1)*(TAUTF+TAUT)
          if(idump.gt.0)print*,i1,MMINV(i1,i1),TAUTF+TAUT
          if(idump.gt.0)print*,TEX
          if(TEX.GT.-200.0D0)THEN
            TLF(i1,i1,L)=DEXP(TEX)
          ELSE
            TLF(i1,i1,L)=0.D0
          ENDIF
          if(idump.gt.0)print*,TLF(i1,i1,L)
          JLF(i1,1,L)= BC*(1.0 - TLF(i1,i1,L))
         end do
         if(idump.gt.0)print*,'gas layer ',l,iscl(l)
        ELSE
         ISCL(L) = 1
         if(idump.gt.0)print*,'cloud layer ',l,iscl(l)
         DO KL=1,NMU
           DO J1=1,NMU
            PPLPL(J1,KL)=PPLPLS(J1,KL)
            PPLMI(J1,KL)=PPLMIS(J1,KL)
           END DO
         END DO

         TAU1 = TAUT
         OMEGA1 = OMEGA
       
         TAUT = TAUTF
         OMEGA = OMEGAF
 	 CALL DOUBLE1(IC,L,RLF(1,1,L),TLF(1,1,L),JLF(1,1,L),
     1      NMU,MAXMU)

C         IF(L.EQ.15)THEN
C          print*,L,TAUT,OMEGA,RL(1,1,L),TL(1,1,L),JL(1,1,L)
C          STOP
C         ENDIF

         TAUT = TAU1
         OMEGA = OMEGA1

        ENDIF

C       Now add R,L and J matrices together with XF as the weight
        DO KL=1,NMU
           DO J1=1,NMU
            TL(KL,J1,L)=XF*TLF(KL,J1,L)+(1.0-XF)*TL(KL,J1,L)
            RL(KL,J1,L)=XF*RLF(KL,J1,L)+(1.0-XF)*RL(KL,J1,L)
           ENDDO
           JL(KL,1,L)=XF*JLF(KL,1,L)+(1.0-XF)*JL(KL,1,L)
        ENDDO

2000  CONTINUE

C     SPECIAL MATRICES FOR LAMBERTIAN REFLECTION AT BOTTOM:
      IF (LOWBC.EQ.1) THEN
	DO J=1,NMU
	  JL(J,1,LTOT) = 0.D0
          IF(IC.EQ.0)THEN
 	   DO I=1,NMU
	    TL(I,J,LTOT) = 0.D0
	    RL(I,J,LTOT) = 2.0D0*GALB*MU(J)*WTMU(J)  !Sum of MU*WTMU = 0.5
C           Make any necessary quadrature correction.
            RL(I,J,LTOT) = RL(I,J,LTOT)*XFAC
	   ENDDO
C	   TL(J,J,LTOT) = 1.D0-GALB
	  ELSE
 	   DO I=1,NMU
	    TL(I,J,LTOT) = 0.0D0
	    RL(I,J,LTOT) = 0.0D0
	   ENDDO
          ENDIF
	ENDDO
      ENDIF


C ***********************************************************************
C  CALCULATE UPWARD MATRICES FOR COMPOSITE OF L LAYERS FROM BASE OF CLOUD.
C    XBASE(I,J,L) IS THE X MATRIX FOR THE BOTTOM L LAYERS OF THE CLOUD.
C  AS FOR "TOP", R01 = R10 & T01 = T10 IS VALID FOR LAYER BEING ADDED ONLY.
C ***********************************************************************
      DO J = 1,NMU
	JBASE(J,1,1) = JL(J,1,LTOT)
	DO I = 1,NMU
	  RBASE(I,J,1) = RL(I,J,LTOT)
	  TBASE(I,J,1) = TL(I,J,LTOT)
	ENDDO
      ENDDO
      DO L = 1,LTOT-1
        if(idump.gt.0)print*,'Multiply ',L,' of ',LTOT-1
	K = LTOT-L
	CALL ADDP( RL(1,1,K), TL(1,1,K), JL(1,1,K), ISCL(K),
     1  RBASE(1,1,L), TBASE(1,1,L), JBASE(1,1,L), RBASE(1,1,L+1), 
     2  TBASE(1,1,L+1), JBASE(1,1,L+1), NMU, MAXMU)
      ENDDO


      IF(IC.NE.0)THEN
       DO IL=1,NMU
        DO J1=1,LTOT
         JBASE(IL,1,J1)=0.
        END DO
       END DO
      END IF
     

      IF(SOL_ANG.GT.90.0)THEN
       ZMU0 = DCOS((180 - SOL_ANG)*PI/180.0)
       SOLAR1 = 0.0
      ELSE
       ZMU0 = DCOS(SOL_ANG*PI/180.0)
       SOLAR1 = SOLAR
      ENDIF
      ZMU = DCOS(EMISS_ANG*PI/180.0)

      ISOL=1
      DO J=1,NMU-1
       IF(ZMU0.LE.MU(J).AND.ZMU0.GT.MU(J+1))ISOL = J
      END DO
      IF(ZMU0.LE.MU(NMU))ISOL=NMU-1 
       
      
      IEMM=1
      DO J=1,NMU-1
       IF(ZMU.LE.MU(J).AND.ZMU.GT.MU(J+1))IEMM = J
      END DO
      IF(ZMU.LE.MU(NMU))IEMM=NMU-1

      FSOL = SNGL((MU(ISOL)-ZMU0)/(MU(ISOL)-MU(ISOL+1)))
      FEMM = SNGL((MU(IEMM)-ZMU)/(MU(IEMM)-MU(IEMM+1)))

      if(idump.gt.0)then
       print*,'isol,fsol',isol,fsol
       print*,'iemm,femm',iemm,femm
      endif

      DO J=1,NMU
       U0PL(J,1) = 0.D0
       IF(IC.EQ.0)THEN
        UTMI(J,1) = RADG(NMU+1-J)
       ELSE
        UTMI(J,1) = 0.
       END IF
      END DO  

      ICO=0
      DO IMU0=ISOL,ISOL+1
       U0PL(IMU0,1) = SOLAR1/(2.0D0*PI*WTMU(IMU0))


      CALL MMUL(1.0D0,RBASE(1,1,LTOT),U0PL,ACOM,NMU,NMU,1,MAXMU,MAXMU,1)
      CALL MMUL(1.0D0,TBASE(1,1,LTOT),UTMI,BCOM,NMU,NMU,1,MAXMU,MAXMU,1)
      CALL MADD(1.0D0,ACOM,BCOM,ACOM,NMU,1,MAXMU,1)
      CALL MADD(1.0D0,JBASE(1,1,LTOT),ACOM,UMI,NMU,1,MAXMU,1)
       

C      IF(IGDIST.EQ.1.AND.VWAVE.EQ.200.0)THEN
C       DO IMU=1,NMU
C        print*,IMU,MU(IMU),SNGL(UMI(IMU,1,1))
C       ENDDO
C      ENDIF

       DO IMU = IEMM,IEMM+1
        ICO=ICO+1
        YX(ICO)=SNGL(UMI(IMU,1,1))
       END DO

       U0PL(IMU0,1)=0.D0

      END DO

      T = FEMM
      U = FSOL

      RADF(IC)=(1-T)*(1-U)*YX(1) + T*(1-U)*YX(2) + T*U*YX(4) + 
     &(1-T)*U*YX(3)
   
      if(idump.gt.0)then
        print*,'radf,y',radf(ic),(yx(imu),imu=1,4)       
      end if

      CONV1 = CONV

      DRAD = RADF(IC)*COS(IC*APHI*sngl(PI)/180.0)
      IF(IC.GT.0)DRAD=DRAD*2

      RAD=RAD+DRAD
      CONV = ABS(100*DRAD/RAD)

      IF(CONV.LT.DEFCONV.AND.CONV1.LT.DEFCONV)THEN
       GOTO 2001
      ENDIF


C     For azimuthal case and RAMAN scattering compute internal radiances
      IF(IC.EQ.0.AND.RAMAN)THEN

C       print*,'Calculating Raman source'


C **********************************************************************
C  CALCULATE DOWNWARD MATRICES FOR COMPOSITE OF L LAYERS FROM TOP OF CLOUD.
C        XTOP(I,J,L) IS THE X MATRIX FOR THE TOP L LAYERS OF CLOUD.
C  NOTE THAT R21 = R12 & T21 = T12 VALID FOR THE HOMOGENEOUS LAYER BEING ADD$
C  BUT NOT FOR THE INHOMOGENEOUS RESULTING "TOP" LAYER.
C
C  i.e. XTOP(I,J,L) is the effective reflectivity, transmission and emission
C  of the top L layers of the atmosphere (i.e. layers 1-L)
C  Specifically
C    RTOP(I,J,L) is RL0
C    TTOP(I,J,L) is T0L
C    JTOP(J,1,L) is JP0L
C **********************************************************************
      DO J = 1,NMU
        JTOP(J,1,1) = JL(J,1,1)
        DO I = 1,NMU
          RTOP(I,J,1) = RL(I,J,1)
          TTOP(I,J,1) = TL(I,J,1)
        ENDDO
      ENDDO
      DO L = 1,LTOT-1
        if(idump.gt.0)print*,'Multiply ',L,' of ',LTOT-1
        CALL ADDP( RL(1,1,L+1), TL(1,1,L+1), JL(1,1,L+1), ISCL(L+1),
     1   RTOP(1,1,L),TTOP(1,1,L), JTOP(1,1,L), RTOP(1,1,L+1),
     2   TTOP(1,1,L+1),JTOP(1,1,L+1), NMU, MAXMU)
      ENDDO

C ***********************************************************************
C  CALCULATE UPWARD MATRICES FOR COMPOSITE OF L LAYERS FROM BASE OF CLOUD.
C    XBASE(I,J,L) IS THE X MATRIX FOR THE BOTTOM L LAYERS OF THE CLOUD.
C  AS FOR "TOP", R01 = R10 & T01 = T10 IS VALID FOR LAYER BEING ADDED ONLY.
C
C  i.e. XBASE(I,J,L) is the effective reflectivity, transmission and emission
C  of the bottom L layers of the atmosphere (i.e. layers LTOT-L+1 to LTOT)
C  Specifically (M = LTOT)
C    RBASE(I,J,L) is RLM
C    TBASE(I,J,L) is TML
C    JBASE(J,1,L) is JMLM
C ***********************************************************************
      DO J = 1,NMU
        JBASE(J,1,1) = JL(J,1,LTOT)
        DO I = 1,NMU
          RBASE(I,J,1) = RL(I,J,LTOT)
          TBASE(I,J,1) = TL(I,J,LTOT)
        ENDDO
      ENDDO


      DO L = 1,LTOT-1
        if(idump.gt.0)print*,'Multiply ',L,' of ',LTOT-1
        K = LTOT-L
        CALL ADDP( RL(1,1,K), TL(1,1,K), JL(1,1,K), ISCL(K),
     1  RBASE(1,1,L), TBASE(1,1,L), JBASE(1,1,L), RBASE(1,1,L+1),
     2  TBASE(1,1,L+1), JBASE(1,1,L+1), NMU, MAXMU)
      ENDDO


C ****************************************************
C      CALCULATING INTERIOR INTENSITIES
C         UPL(J,1,L) GOES DOWN OUT OF LAYER L.
C          UMI(J,1,L) GOES UP OUT OF LAYER L.
C ****************************************************

C      Assume solar zenith angle is aligned with one of the quadrature zenith a>
       DO I=1,NMU
        U0PL(I,1)=0.0
       ENDDO
       U0PL(ISOL,1) = SOLAR1/(2.0D0*PI*WTMU(ISOL))

       DO I=1,NMU
        UMI(I,1,1)=JBASE(I,1,1)
       ENDDO
       DO L = 1,LTOT-1
        K = LTOT-L
C       Calculate I(L)-
        CALL IUP( RTOP(1,1,L), TTOP(1,1,L), JTOP(1,1,L),
     1    RBASE(1,1,K), TBASE(1,1,K), JBASE(1,1,K), UMI(1,1,L),
     2    NMU, MAXMU)

C       Calculate I(L)+
        CALL IDOWN( RTOP(1,1,L), TTOP(1,1,L), JTOP(1,1,L),
     1    RBASE(1,1,K), TBASE(1,1,K), JBASE(1,1,K), UPL(1,1,L),
     2    NMU, MAXMU)

       ENDDO
C       print*,ISOL,SOLAR1,U0PL(ISOL,1)*WTMU(ISOL)
       do L=1,LTOT-1
        SUM1=0.0
        SUM2=0.0
        DO IMU=1,NMU
         SUM1=SUM1+SNGL(UPL(IMU,1,L)*WTMU(IMU))
         SUM2=SUM2+SNGL(UMI(IMU,1,L)*WTMU(IMU))
        ENDDO
C        print*,'aa',L,SUM1,SUM2,SUM1+SUM2
       enddo

C *******************************************************
C     CALCULATING EXTERIOR INTENSITIES UTPL AND U0MI
C *******************************************************

C      Basically upward radiation at top of atmosphere:
C      U0MI = U0PL*RBASE(LTOT) + TBASE*UTMI(LTOT) + JBASE(LTOT)

       CALL MMUL(1.0D0,RBASE(1,1,LTOT),U0PL,ACOM,NMU,NMU,1,
     1   MAXMU,MAXMU,1)
       CALL MMUL(1.0D0,TBASE(1,1,LTOT),UTMI,BCOM,NMU,NMU,1,
     1   MAXMU,MAXMU,1)
       CALL MADD(1.0D0,ACOM,BCOM,ACOM,NMU,1,MAXMU,1)
       CALL MADD(1.0D0,JBASE(1,1,LTOT),ACOM,U0MI,NMU,1,MAXMU,1)

C      Basically downward radiation at bottom of atmosphere:
C      UTPL = U0PL*TTOP(NLAY) + RTOP(NLAY)*UTMI + JTOP(NLAY)

       CALL MMUL(1.0D0,TTOP(1,1,LT1),U0PL,ACOM,NMU,NMU,1,
     1  MAXMU,MAXMU,1)
       CALL MMUL(1.0D0,RTOP(1,1,LT1),UTMI,BCOM,NMU,NMU,1,
     1  MAXMU,MAXMU,1)
       CALL MADD(1.0D0,ACOM,BCOM,ACOM,NMU,1,MAXMU,1)
       CALL MADD(1.0D0,JTOP(1,1,LTOT),ACOM,UTPL(1,1,LT1),
     1  NMU,1,MAXMU,1)

       DO 302 IMU = 1, NMU
C       Keep order of the zenith angles in the arrays
C       such that they go from large to small mu, opposite to way listed in the
C       <runname>.sca files
        DO 301 L=1,LT1
          IF(L.NE.LT1)THEN
            UPLF(IMU,L)=SNGL(UPL(IMU,1,L))
          ELSE
            UPLF(IMU,L)=SNGL(UTPL(IMU,1,L))
          ENDIF
          IF(L.NE.1)THEN
            UMIF(IMU,L)=SNGL(UMI(IMU,1,L))
          ELSE
            UMIF(IMU,L)=SNGL(U0MI(IMU,1,L))
          ENDIF
301     CONTINUE

302    CONTINUE


C      Now compute total irradiance at each level in atmosphere
C      by integrating over solid angle at each point in atmosphere.
C      Remember that UPLF is radiance DOWN out of layer L and UMIF is
C      radiance UP out of layer L. 2*PI term of azimuth integration is
C      already included fluxes computed and so we don't need to add anything
C      here. 

       DO 303 I=1,LT1
        IF(I.EQ.1)THEN
         FPL = SNGL(U0PL(ISOL,1)*WTMU(ISOL))
        ELSE
         FPL=0.0
         DO IMU=1,NMU
          FPL=FPL+UPLF(IMU,I-1)*SNGL(WTMU(IMU))
         ENDDO
        ENDIF
        FMI = 0.0
        DO IMU=1,NMU
         FMI=FMI+UMIF(IMU,I)*SNGL(WTMU(IMU))
        ENDDO
        IRAD(I)=(FMI+FPL)
303    CONTINUE


C      Now find mass-weighted values of irradiance
       DO 304 ILAY=1,LT1
        IF(ILAY.EQ.1)then
         LDENS1=ALOG(0.001*DENS(1))
         F1=SNGL(U0PL(ISOL,1)*2*PI*WTMU(ISOL))
        ELSE
         LDENS1=ALOG(DENS(ILAY-1))
         F1=IRAD(ILAY-1)
        ENDIF
        LDENS2=ALOG(DENS(ILAY))
        F2=IRAD(ILAY)
        SUM1=0.0
        SUM2=0.0
        DO J=1,50
         LD = LDENS1 + (J-1)*(LDENS2-LDENS1)/49.0
         FF = F1 + (J-1)*(F2-F1)/49.0
         XD = EXP(LD)
         SUM1=SUM1+XD
         SUM2=SUM2+XD*FF
        ENDDO
        FMEAN(ILAY)=SUM2/SUM1
304    CONTINUE

C      Calculate RAMAN x-section and add to source array, JRAMAN, which
C      is held in a COMMON block.
       CALL RAMANJSOURCE(VV,LT1,FMEAN)

      ENDIF

1000  CONTINUE

2001  CONTINUE
      if(idump.gt.0)print*,'rad = ',rad

      RETURN
      END

