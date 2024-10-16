      subroutine scloud11wave( rad, sol_ang, emiss_ang, aphi, radg, 
     1 solar, lowbc, galb, iray, mu1, wt1, nmu, nf, igdist, vwave, vv,
     2 eps, omegas, bnu, tau, tauray, nlay, ncont, lfrac)
C     $Id: scloud11wave.f,v 1.11 2011-06-17 15:57:54 irwin Exp $
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
C       eps(MAXSCATLAY)     double  Fraction of thermal absorption in each layer
C	omegas(MAXSCATLAY)  double  Single scattering albedo of each layer
C       bnu(MAXSCATLAY)     real  Mean Planck function in each layer
C       tau(MAXSCATLAY)     real  Total optical thickness of each layer
C       tauray(MAXSCATLAY)  real  Rayleigh optical thickness of each layer
C       nlay            integer Number of layers
C       ncont           integer Number of aerosol types/modes
C       lfrac(MAXCON,MAXSCATLAY) real Fraction of scattering contributed by each
C                                  type in each layer
C       liscat(MAXCON)   integer Phase function type identifier
C       lncons(MAXCON)   integer Number of phase function parameters
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
      INTEGER IC,I1,J1,NCONT,NCONT1,NCONS,NORM,ISCAT,KL,NPHI,IP,JP,MF
      INTEGER IPOW0,L,IL,IMU0,IDUMP
      double precision mu1(maxmu), wt1(maxmu), PI
      real rad, radf(0:MAXF-1), radg(maxmu), bnu(MAXSCATLAY), 
     1 tau(MAXLAY),vwave, solar, vv, tauray(MAXLAY),
     2 lfrac(MAXCON,MAXLAY),aphi,drad,solar1,sum
      integer nf,imu,iray
      integer igdist,imie
      logical lrep,checkconv
      common /imiescat/imie

      DOUBLE PRECISION MU(MAXMU), WTMU(MAXMU), MM(MAXMU,MAXMU),
     1 MMINV(MAXMU,MAXMU), RL(MAXMU,MAXMU,MAXSCATLAY),
     2 TL(MAXMU,MAXMU,MAXSCATLAY), JL(MAXMU,1,MAXSCATLAY),
     3 RBASE(MAXMU,MAXMU,MAXSCATLAY),F,G1,G2,EPS(MAXSCATLAY),
     4 TBASE(MAXMU,MAXMU,MAXSCATLAY), JBASE(MAXMU,1,MAXSCATLAY),
     5 UMI(MAXMU,1,MAXSCATLAY),OMEGAS(MAXSCATLAY),
     6 ACOM(MAXMU,1), BCOM(MAXMU,1),
     7 PPLPL(MAXMU,MAXMU), PPLMI(MAXMU,MAXMU), CONS8(MAXSCATPAR),
     8 PPLPLS(MAXMU,MAXMU),PPLMIS(MAXMU,MAXMU)
      DOUBLE PRECISION PPLN(MAXCON,0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION PMIN(MAXCON,0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION PPLR(0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION PMIR(0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION E(MAXMU,MAXMU),CCINV(MAXMU,MAXMU),
     1  CC(MAXMU,MAXMU),U0PL(MAXMU,1),UTMI(MAXMU,1),TAUT,
     2  BC,OMEGA,FRAC,SFRAC
      DOUBLE PRECISION STEST,GALB,ZMU0,ZMU,TEX, TAUSCAT, TAUR, TAU1, 
     1 STEST1
      DOUBLE PRECISION RTOP(MAXMU,MAXMU,MAXSCATLAY),
     1 UPL(MAXMU,1,MAXSCATLAY),JTOP(MAXMU,1,MAXSCATLAY),
     2 UTPL(MAXMU,1,MAXSCATLAY),TTOP(MAXMU,MAXMU,MAXSCATLAY),
     3 U0MI(MAXMU,1,MAXSCATLAY)
      REAL UMIF(MAXMU,MAXSCATLAY),UPLF(MAXMU,MAXSCATLAY)
      REAL FPL,FMI,IRAD(MAXSCATLAY),F1,F2

      REAL SOL_ANG,EMISS_ANG,CONV,DEFCONV
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

1     FORMAT(A)

C     Define additional debugging 'reporting variable', lrep
C      lrep=.true.
      lrep=.false.
      if(idump.gt.0.or.(lrep.and.igdist.eq.1))then
      print*,'scloud11wave : sol_ang',sol_ang, 'emiss_ang',emiss_ang,
     1 'aphi', aphi, 'solar', solar, 'lowbc',lowbc, 'galb',galb,
     2 'nmu',nmu, 'nf',nf, 'nlay',nlay, 'ncont',ncont,
     3 'igdist',igdist
      print*,'vwave : ',vwave,'vv : ',vv
      print*,(mu1(i),i=1,nmu)
      print*,(wt1(i),i=1,nmu)
      print*, (radg(i),i=1,nmu) 
      do i=1,nlay
       print*,i,tau(i),(1.0 - eps(i)),omegas(i),bnu(i),tauray(i)
      end do
      do i=1,nlay
       print*,i,(1.0 - eps(i)),omegas(i),(lfrac(j,i),j=1,ncont)
       sum=0.0
       do j=1,ncont
        sum=sum+lfrac(j,i)
       enddo
       print*,'sum = ',sum
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

      LTOT = NLAY     ! Set internal number of layers
      LT1 = LTOT

C     In case of Lambertian reflection, add extra dummy layer at bottom,
C     whose transmission matrix = (1-A)*Unit-Matrix. This layer must be
C     omitted from computation by doubling

      IF (LOWBC.EQ.1) LTOT = LTOT+1
      IF (LTOT.GT.MAXSCATLAY) THEN
        if(idiag.gt.0)print*, 'ltot ', ltot,' MAXSCATLAY ', MAXSCATLAY
        CALL ABEND(' SCLOUD11WAVE: TOO MANY LAYERS')
      ENDIF
      IF (NMU.GT.MAXMU) THEN
        CALL ABEND(' SCLOUD11WAVE: TOO MANY ANGLE POINTS')
      ENDIF
      if(idump.gt.0.or.lrep)print*,'scloud11wave: NLAY,LTOT,IGDIST',
     1  NLAY,LTOT,IGDIST

C     Reset the order of angles:
      DO I=1,NMU
        MU(I) = MU1(NMU+1-I)
        WTMU(I) = WT1(NMU+1-I)
        IF(IDUMP.GT.0)PRINT*,I,MU(I),WTMU(I)
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
         IF(IDUMP.GT.0)THEN
          print*,I,J,CC(I,J),CCINV(I,J)
         ENDIF
        ENDDO
      ENDDO

C      lrep=.true.
      lrep=.false.
C     PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C     Precalculate phase function arrays if new wavelength
      IF(IGDIST.EQ.1) THEN
       DO J1=1,NCONT
        if(idump.gt.0.or.lrep)print*,'j1,vwave',j1,
     1   vwave


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
         if(idump.gt.0.or.lrep)print*,
     1    'Calculating matrix from scratch IC = ',IC
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
        IF(IDUMP.GT.0)PRINT*,'Calculating Rayleigh scattering phase'
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

800     CONTINUE
       ENDIF
C     PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
      ENDIF




      RAD=0.0
      CONV = 100.0
      
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
        LREP=.FALSE.
C        if(L.EQ.33.and.IGDIST.EQ.1.AND.IC.EQ.0)LREP=.TRUE.
        if(idump.ne.0)print*,'L,IGDIST =',L,IGDIST
        TAUT = 1.0D0*TAU(L)
        BC = 1.0D0*BNU(L)
C        OMEGA = 1.0D0*(1. - EPS(L))
        OMEGA = OMEGAS(L)

        TAUSCAT = TAUT*OMEGA
        TAUR = TAURAY(L)

        if(lrep)then
         print*,'l',L
         print*,'tau',TAU(L),TAUT
         print*,'bnu',BNU(L),BC
         print*,'eps',EPS(L),OMEGA
         print*,'tauscat,taur,tauray',TAUSCAT,TAUR,TAURAY(L)
        endif

        if(idump.ne.0)then
         print*,L
         print*,TAU(L),TAUT
         print*,BNU(L),BC
         print*,EPS(L),OMEGA
         print*,TAUSCAT,TAUR,TAURAY(L)
        endif

C       Calling codes now already include Rayleigh optical depth in 
C       tauscat if IRAY>0, so we need to subtract it first here
        TAUSCAT = TAUSCAT-TAUR

C       Add in an error trap to counter single-double subtraction overflows 
        IF(TAUSCAT.LT.0.)TAUSCAT=0.

        if(idump.gt.0)then
         print*,'L,TAUT,BC,OMEGA = ',L,TAUT,BC,OMEGA
         print*,'TAUSCAT,TAUR',TAUSCAT,TAUR
        endif
        if(lrep)then
         print*,'L,TAUT,BC,OMEGA = ',L,TAUT,BC,OMEGA
         print*,'TAUSCAT,TAUR',TAUSCAT,TAUR
        endif

        IF(TAUT.LT.0.0)THEN
         IF(IDIAG.GT.0)THEN
          PRINT*,'Error in scloud11wave TAUT < 0. Setting to zero'
          PRINT*,'L,TAUT,BC,OMEGA'
          PRINT*,L,TAUT,BC,OMEGA
         ENDIF
         TAUT = 0.0
        END IF

        IF(OMEGA.GT.1.0)THEN
         OMEGA = 1.D0
        ELSE IF(OMEGA.LT.0.0)THEN
         OMEGA = 0.D0
        ENDIF

C       P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P
C
C       Assign phase function coefficients for each atmospheric layer

        DO KL=1,NMU
          DO J1=1,NMU
           PPLPLS(J1,KL)=0.
           PPLMIS(J1,KL)=0.
          END DO
        END DO
        SFRAC=0.D0

        DO J1 = 1,NCONT

          FRAC = LFRAC(J1,L)
          TAU1 = FRAC*TAUSCAT

          IF((TAUSCAT+TAUR).GT.0.0)THEN
            FRAC = TAU1/(TAUSCAT+TAUR)
          ELSE
            FRAC = 0.0
          ENDIF


          SFRAC = SFRAC+FRAC
          if(idump.gt.0)print*,J1,FRAC,SFRAC
          if(lrep)print*,'j1,frac,sfrac',J1,FRAC,SFRAC

          IF(FRAC.GT.0.0)THEN
           DO JP = 1,NMU
            DO IP = 1,NMU
             PPLPL(IP,JP) = PPLN(J1,IC,IP,JP)
             PPLMI(IP,JP) = PMIN(J1,IC,IP,JP)
            END DO
           END DO
           if(idump.gt.0)print*,'Layer ',L,' Adding fraction ',
     &	    FRAC,' of mode ',J1
           if(lrep)print*,'Layer ',L,' Adding fraction ',
     &	    FRAC,' of mode ',J1
           CALL MADD(FRAC,PPLPLS,PPLPL,PPLPLS,NMU,NMU,
     1      MAXMU,MAXMU)
           CALL MADD(FRAC,PPLMIS,PPLMI,PPLMIS,NMU,NMU,
     1      MAXMU,MAXMU)
           if(idump.gt.0)then
            print*,'PPLPL'
            DO JP=1,NMU
             PRINT*,(PPLPLS(JP,IP),IP=1,NMU)
            ENDDO
            PRINT*,'PPLMI'
            DO JP=1,NMU
             PRINT*,(PPLMIS(JP,IP),IP=1,NMU)
            ENDDO
           endif
          ENDIF

         END DO
      
         IF((TAUSCAT+TAUR).GT.0.0)THEN
            FRAC = TAUR/(TAUSCAT+TAUR)
         ELSE
            FRAC = 0.0
         ENDIF

         SFRAC = SFRAC+FRAC

         if(lrep)print*,'Ray ',FRAC,SFRAC      
         if(idump.gt.0)print*,'Ray ',FRAC,SFRAC
         IF(FRAC.GT.0)THEN
          DO JP = 1,NMU
            DO IP = 1,NMU
             PPLPL(IP,JP) = PPLR(IC,IP,JP)
             PPLMI(IP,JP) = PMIR(IC,IP,JP)
            END DO
          END DO
          if(idump.gt.0)print*,'Layer ',L,' Adding fraction ',
     &		FRAC,' of Rayleigh scattering'
          if(lrep)print*,'Layer ',L,' Adding fraction ',
     &		FRAC,' of Rayleigh scattering'
          CALL MADD(FRAC,PPLPLS,PPLPL,PPLPLS,NMU,NMU,
     1      MAXMU,MAXMU)
          CALL MADD(FRAC,PPLMIS,PPLMI,PPLMIS,NMU,NMU,
     1      MAXMU,MAXMU)

           if(idump.gt.0)then
            print*,'PPLPL'
            DO JP=1,NMU
             PRINT*,(PPLPLS(JP,IP),IP=1,NMU)
            ENDDO
            PRINT*,'PPLMI'
            DO JP=1,NMU
             PRINT*,(PPLMIS(JP,IP),IP=1,NMU)
            ENDDO
           endif


         ENDIF

C         TAUT = TAUT + TAUR
         OMEGA = (TAUSCAT+TAUR)/TAUT
         
         if(idump.gt.0)then
          print*,'Rayleigh: TAUT,TAUR,TAUSCAT,TAUT-TAUR,OMEGA',
     &       TAUT,TAUR,TAUSCAT,TAUT-TAUR,OMEGA
         endif
         if(lrep)then
          print*,'Rayleigh: TAUR,TAUSCAT,TAUT-TAUR,OMEGA',TAUR,
     &          TAUSCAT,TAUT-TAUR,OMEGA
         endif

         STEST=ABS(SFRAC-1.D0)
	 STEST1 = ABS(SFRAC-0.D0)
         IF(STEST.GT.0.02D0.AND.STEST1.GT.0.02D0.AND.IDIAG.GT.0)THEN
          PRINT*,'Scloud11wave. WARNING.  SUM(FRAC) must = 0 OR 1.'
          print*,sfrac
         END IF

C         IF(IC.EQ.0.AND.SFRAC.GT.0.0)THEN
C           CALL HANSEN( IC, PPLPLS, PPLMIS, MAXMU, WTMU, NMU)
C         ENDIF

         DO KL=1,NMU
          DO J1=1,NMU
           PPLPL(J1,KL)=PPLPLS(J1,KL)
           PPLMI(J1,KL)=PPLMIS(J1,KL)
          END DO
         END DO


C       P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P


        IF(TAUT.EQ.0)THEN
         do i1=1,nmu
          do j1=1,nmu
           RL(I1,J1,L)=0.0
           TL(I1,J1,L)=0.0
          end do
          JL(i1,1,L) =0.0
          TL(i1,i1,L)=1.0
         end do
         iscl(l) = 0
         if(idump.gt.0.or.lrep)print*,'zero layer ',l,iscl(l)
        ELSE IF(OMEGA.EQ.0)THEN
         do i1=1,nmu
          do j1=1,nmu
           RL(I1,J1,L)=0.0
           TL(I1,J1,L)=0.0
          end do
          TEX = -MMINV(i1,i1)*TAUT
          if(idump.gt.0.or.lrep)print*,i1,MMINV(i1,i1),TAUT
          if(idump.gt.0.or.lrep)print*,TEX
          if(TEX.GT.-200.0D0)THEN
            TL(i1,i1,L)=DEXP(TEX)
          ELSE
            TL(i1,i1,L)=0.D0
          ENDIF
          if(idump.gt.0.or.lrep)print*,TL(i1,i1,L)
          JL(i1,1,L)= BC*(1.0 - TL(i1,i1,L))
         end do
         iscl(l) = 0
         if(idump.gt.0.or.lrep)print*,'gas layer ',l,iscl(l)
        ELSE
         iscl(l) = 1
         if(idump.gt.0.or.lrep)print*,'cloud layer ',l,iscl(l)
C         OMEGA=0.3
 	 CALL DOUBLE1(IC,L,RL(1,1,L),TL(1,1,L),JL(1,1,L),NMU,MAXMU)
        END IF

2000  CONTINUE

C  SPECIAL MATRICES FOR LAMBERTIAN REFLECTION AT BOTTOM:
      IF (LOWBC.EQ.1) THEN
	DO J=1,NMU
	  JL(J,1,LTOT) = (1.0-GALB)*RADG(NMU+1-J)
          IF(IC.EQ.0)THEN
 	   DO I=1,NMU
C           Transmission of surface layer set to zero following
C           Plass et al.,1973.
	    TL(I,J,LTOT) = 0.D0
	    RL(I,J,LTOT) = 2.0D0*GALB*MU(J)*WTMU(J)  !Sum of MU*WTMU = 0.5
C           Make any necessary quadrature correction.
            RL(I,J,LTOT) = RL(I,J,LTOT)*XFAC
	   ENDDO
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
C
C  i.e. XBASE(I,J,L) is the effective reflectivity, transmission and emission
C  of the bottom L layers of the atmosphere (i.e. layers LTOT-L+1 to LTOT)
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

      DRAD = RADF(IC)*COS(IC*APHI*SNGL(PI)/180.0)
      IF(IC.GT.0)DRAD=DRAD*2

      RAD=RAD+DRAD
      CONV = ABS(100*DRAD/RAD)
      if(igdist.eq.1.and.checkconv.and.idiag.gt.0)then
       print*,'Conv: IMU0,IMU,NF,IC,RAD,DRAD,100*DRAD/RAD',
     1  INT((1-FSOL)*ISOL+FSOL*(ISOL+1)+0.5),
     2  INT((1-FEMM)*IEMM+FEMM*(IEMM+1)+0.5),NF,IC,RAD,DRAD,
     3  100*DRAD/RAD
      endif

      IF(CONV.LT.DEFCONV)THEN
       if(idiag.gt.0)then
        PRINT*,'scloud11wave: Converged after ',IC,
     1  ' fourier components'
       endif
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

C      Assume solar zenith angle is aligned with one of the quadrature zenith angles
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

