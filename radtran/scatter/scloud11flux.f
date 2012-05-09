      subroutine scloud11flux(radg, solar, sol_ang, lowbc, galb,
     1 iray,mu1, wt1, nmu, nf, igdist, vwave, vv,eps, bnu, tau,
     2 tauray, nlay, ncont, lfrac, umif, uplf)
C     $Id: scloud11flux.f,v 1.10 2011-06-17 15:57:54 irwin Exp $
C     ***********************************************************************
C     Compute and return interbal radiation fields in a scattering atmosphere
C     Code uses matrix operator algorithm.  Diffuse incident radiation 
C     is allowed at the bottom and single-beam incident radiation 
C     (sunlight) at the top. 
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
C     Modified from scloud11wave.
C
C	Input variables
C       radg(maxmu)     real  Incident intensity at the bottom of the atm
C       solar           real  Incident solar flux at the top of the atm
C	sol_ang		real  solar zenith angle (degrees)
C       lowbc           integer Lower boundary condition:
C                               0 = thermal, 1 = Lambert reflection.
C       galb            double precision  Ground albedo at the bottom.
C	iray		integer	Rayleigh scattering on/off
C       mu1(maxmu)      double precision  zenith angle point quadrature
C       wt1(maxmu)      double precision  zenith angle point quadrature
C       nmu             integer Number of points in quadrature
C       nf              integer Number of terms in azimuth Fourier expansion   
C	igdist          integer Switch to recalculate phase functions. If > 1
C                               then phase functions left alone
C       eps(MAXSCATLAY)     real*8  Fraction of thermal absorption in each layer
C       bnu(MAXSCATLAY)     real  Mean Planck function in each layer
C       tau(MAXSCATLAY)     real  Total optical thickness of each layer
C       tauray(MAXSCATLAY)  real*4  Rayleigh optical thickness of each layer
C       nlay            integer Number of layers
C       ncont           integer Number of aerosol types/modes
C       lfrac(MAXCON,MAXSCATLAY) real*4 Fraction of scattering contributed by each
C                                  type in each layer
C
C	Output variables
C       uplf(maxmu,MAXSCATLAY,MAXF) real      Internal radiances (downwards)
C       umif(maxmu,MAXSCATLAY,MAXF) real      Internal radiances (upwards)
C
C	New History
C
C	3/11/94	PGJI	Adapted for RADTRAN
C	7/1/97  PGJI   	Revised again to pass out all fourier coefficients
C	7/5/97	PGJI	Optimised from Scloud6 for maximum speed
C	3/8/05  PGJI	Adapted to output internal radiances
C
C     ***********************************************************************

      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER LTOT,NLAY,LT1,LOWBC,NMU,I,J,K
      INTEGER IC,I1,J1,NCONT,NCONS,NORM,ISCAT,KL,NPHI,IP,JP
      INTEGER IPOW0,L,IL,IMU0,IDUMP,NCONT1
      double precision mu1(maxmu), wt1(maxmu), PI
      real radg(maxmu), bnu(MAXSCATLAY), 
     1 tau(MAXLAY),vwave, solar, vv, tauray(MAXLAY),
     2 lfrac(MAXCON,MAXLAY),drad,solar1
      integer nf,imu,iray,jmu
      integer igdist

      DOUBLE PRECISION MU(MAXMU), WTMU(MAXMU), MM(MAXMU,MAXMU),
     1 MMINV(MAXMU,MAXMU), RL(MAXMU,MAXMU,MAXSCATLAY),EPS(MAXSCATLAY),
     2 TL(MAXMU,MAXMU,MAXSCATLAY), JL(MAXMU,1,MAXSCATLAY),
     3 RBASE(MAXMU,MAXMU,MAXSCATLAY),F,G1,G2,
     9 RTOP(MAXMU,MAXMU,MAXSCATLAY),
     4 TBASE(MAXMU,MAXMU,MAXSCATLAY), JBASE(MAXMU,1,MAXSCATLAY),
     5 TTOP(MAXMU,MAXMU,MAXSCATLAY),
     6 ACOM(MAXMU,1), BCOM(MAXMU,1),JTOP(MAXMU,1,MAXSCATLAY),
     7 PPLPL(MAXMU,MAXMU), PPLMI(MAXMU,MAXMU), CONS8(MAXSCATPAR),
     8 PPLPLS(MAXMU,MAXMU),PPLMIS(MAXMU,MAXMU),
     9 UMI(MAXMU,1,MAXSCATLAY), UPL(MAXMU,1,MAXSCATLAY),
     & U0MI(MAXMU,1,MAXSCATLAY), UTPL(MAXMU,1,MAXSCATLAY)    
      REAL UMIF(MAXMU,MAXSCATLAY,MAXF),UPLF(MAXMU,MAXSCATLAY,MAXF)

      DOUBLE PRECISION PPLN(MAXCON,0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION PMIN(MAXCON,0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION PPLR(0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION PMIR(0:MAXF-1,MAXMU,MAXMU)
      DOUBLE PRECISION E(MAXMU,MAXMU),CCINV(MAXMU,MAXMU),CC(MAXMU,MAXMU)
      DOUBLE PRECISION U0PL(MAXMU,1),UTMI(MAXMU,1),TAUT,BC,
     1 OMEGA,FRAC,SFRAC,STEST,GALB,ZMU0,ZMU,TEX, TAUSCAT, 
     2 TAUR, TAU1, STEST1

      REAL SOL_ANG,FSOL
      INTEGER ISOL,ISCL(MAXSCATLAY)

C     Common blocks
      COMMON/UNIT/ E
      COMMON/AREA1/CCINV,MM
      COMMON/AREA2/CC,MMINV,MU,PI
      COMMON/INPUT/U0PL,UTMI
      COMMON/HOMOG/ TAUT, BC, OMEGA, IPOW0
      COMMON/PHMAT/ PPLPL, PPLMI
      COMMON/NEWPH/PPLN,PMIN,PPLR,PMIR
      COMMON /SCATDUMP/ IDUMP


C--------------------------------------------------------------------------


      if(idump.gt.0)then
      print*,'scloud11flux : sol_ang',sol_ang,'solar', solar,
     1  'lowbc',lowbc,'galb',galb,'nmu',nmu, 'nf',nf, 'nlay',nlay,
     2  'ncont',ncont,'igdist',igdist
      print*,'vwave : ',vwave,'vv : ',vv
      print*,(mu1(i),i=1,nmu)
      print*,(wt1(i),i=1,nmu)
      print*, (radg(i),i=1,nmu) 
      do i=1,nlay
       print*,tau(i),(1.0 - eps(i)),bnu(i),tauray(i)
      end do
      do i=1,nlay
       print*,(1.0 - eps(i)),(lfrac(j,i),j=1,ncont)
      end do
      endif

      LTOT = NLAY      ! Set internal number of layers
      LT1 = LTOT

C     In case of Lambertian reflection, add extra dummy layer at bottom,
C     whose transmission matrix = (1-A)*Unit-Matrix. This layer must be
C     omitted from computation by doubling

      IF (LOWBC.EQ.1) LTOT = LTOT+1
      IF (LTOT.GT.MAXSCATLAY) THEN
       print*, 'ltot ', ltot, ' MAXSCATLAY ', MAXSCATLAY
       CALL ABEND(' SCLOUD11FLUX: TOO MANY LAYERS')
      ENDIF
      IF (NMU.GT.MAXMU) THEN
       CALL ABEND(' SCLOUD11FLUX: TOO MANY ANGLE POINTS')
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

        ISCAT =2
        NCONS = 3
        CONS8(1)=F
        CONS8(2)=G1
        CONS8(3)=G2
        NPHI = 100
        NORM = 1

        DO 900 IC=0,NF
C         print*,'Calculating matrix from scratch IC = ',IC
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
        NCONT1=NCONT+11
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

800     CONTINUE
       ENDIF
C     PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
      ENDIF


      DO 1000 IC=0,NF


C **********************************************************************
C     CALCULATE REFLECTION,TRANSMISSION AND SOURCE MATRICES FOR EACH    
C	  HOMOGENEOUS LAYER L; INSERT INTO ARRAY X(I,J,L).	  	      
C **********************************************************************    
C
      IPOW0 = 16
      DO 2000 L = 1,LT1
        if(idump.ne.0)print*,'L,IGDIST =',L,IGDIST
        TAUT = 1.0D0*TAU(L)
        BC = 1.0D0*BNU(L)
        OMEGA = 1.0D0*(1. - EPS(L))

        TAUSCAT = TAUT*OMEGA
        TAUR = TAURAY(L)

C       Calling codes now already include Rayleigh optical depth in
C       tauscat if IRAY>1, so we need to subtract it first here
        TAUSCAT = TAUSCAT-TAUR

        if(idump.gt.0)then
         print*,'L,TAUT,BC,OMEGA = ',L,TAUT,BC,OMEGA
         print*,'TAUSCAT,TAUR',TAUSCAT,TAUR
        endif

        IF(TAUT.LT.0.0)THEN
         PRINT*,'Error in scloud11flux TAUT < 0. Setting to zero'
         PRINT*,'L,TAUT,BC,OMEGA'
         PRINT*,L,TAUT,BC,OMEGA
         TAUT = 0.0
        END IF

        IF(OMEGA.GT.1.0)THEN
C         print*,'Omega too big! Reducing to 1.0',OMEGA
         OMEGA = 1.D0
        ELSE IF(OMEGA.LT.0.0)THEN
C         print*,'Omega too small! Setting to 0.0',OMEGA
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

          IF(FRAC.GT.0.0)THEN
           DO JP = 1,NMU
            DO IP = 1,NMU
             PPLPL(IP,JP) = PPLN(J1,IC,IP,JP)
             PPLMI(IP,JP) = PMIN(J1,IC,IP,JP)
            END DO
           END DO
           if(idump.gt.0)print*,'Layer ',L,' Adding fraction ',
     &       FRAC,' of mode ',J1
           CALL MADD(FRAC,PPLPLS,PPLPL,PPLPLS,NMU,NMU,
     1      MAXMU,MAXMU)
           CALL MADD(FRAC,PPLMIS,PPLMI,PPLMIS,NMU,NMU,
     1      MAXMU,MAXMU)

          ENDIF

         END DO
      
         IF((TAUSCAT+TAUR).GT.0.0)THEN
            FRAC = TAUR/(TAUSCAT+TAUR)
         ELSE
            FRAC = 0.0
         ENDIF

         SFRAC = SFRAC+FRAC
      
         IF(FRAC.GT.0)THEN
          DO JP = 1,NMU
            DO IP = 1,NMU
             PPLPL(IP,JP) = PPLR(IC,IP,JP)
             PPLMI(IP,JP) = PMIR(IC,IP,JP)
            END DO
          END DO

          if(idump.gt.0)print*,'Layer ',L,' Adding fraction ',
     &     FRAC,' of Rayleigh scattering '
          CALL MADD(FRAC,PPLPLS,PPLPL,PPLPLS,NMU,NMU,
     1      MAXMU,MAXMU)
          CALL MADD(FRAC,PPLMIS,PPLMI,PPLMIS,NMU,NMU,
     1      MAXMU,MAXMU)

         ENDIF

         TAUT = TAUT + TAUR
         OMEGA = (TAUSCAT+TAUR)/TAUT
         
         if(idump.gt.0)then
          print*,'Rayleigh: TAUR,TAUSCAT,TAUT-TAUR,OMEGA',TAUR,
     &          TAUSCAT,TAUT-TAUR,OMEGA
         endif

         STEST=ABS(SFRAC-1.D0)
         STEST1 = ABS(SFRAC-0.D0)
         IF(STEST.GT.0.02D0.AND.STEST1.GT.0.02D0)THEN
          PRINT*,'Scloud11flux. ERROR.  SUM(FRAC) must = 0 OR 1.'
          print*,sfrac
          STOP
         END IF

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
         iscl(l) = 0
         if(idump.gt.0)print*,'gas layer ',l,iscl(l)
        ELSE
         iscl(l) = 1
         if(idump.gt.0)print*,'cloud layer ',l,iscl(l)
         CALL DOUBLE1(IC,L,RL(1,1,L),TL(1,1,L),JL(1,1,L),NMU,MAXMU)

        END IF


2000  CONTINUE

C  SPECIAL MATRICES FOR LAMBERTIAN REFLECTION AT BOTTOM:
      IF (LOWBC.EQ.1) THEN
        DO J=1,NMU
          JL(J,1,LTOT) = 0.D0
          IF(IC.EQ.0)THEN
            DO I=1,NMU
             TL(I,J,LTOT) = 0.D0
             RL(I,J,LTOT) = 2.0D0*GALB*MU(J)*WTMU(J)  !Sum of MU*WTMU
C                                                           = 0.5
            ENDDO
            TL(J,J,LTOT) = 1.D0-GALB
          ELSE
            DO I=1,NMU
             TL(I,J,LTOT) = 0.0D0
             RL(I,J,LTOT) = 0.0D0
            ENDDO
          ENDIF
        ENDDO
      ENDIF


C **********************************************************************
C  CALCULATE DOWNWARD MATRICES FOR COMPOSITE OF L LAYERS FROM TOP OF CLOUD.
C        XTOP(I,J,L) IS THE X MATRIX FOR THE TOP L LAYERS OF CLOUD.
C  NOTE THAT R21 = R12 & T21 = T12 VALID FOR THE HOMOGENEOUS LAYER BEING ADD$
C  BUT NOT FOR THE INHOMOGENEOUS RESULTING "TOP" LAYER.
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
         JTOP(IL,1,J1)=0.0
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

      ISOL=1
      DO J=1,NMU-1
       IF(ZMU0.LE.MU(J).AND.ZMU0.GT.MU(J+1))ISOL = J
      END DO
      IF(ZMU0.LE.MU(NMU))ISOL=NMU-1 
       
      
      FSOL = SNGL((MU(ISOL)-ZMU0)/(MU(ISOL)-MU(ISOL+1)))

      if(idump.gt.0)then
       print*,'isol,fsol',isol,fsol
      endif

      DO J=1,NMU
       U0PL(J,1) = 0.D0
       IF(IC.EQ.0)THEN
        UTMI(J,1) = RADG(NMU+1-J)
       ELSE
        UTMI(J,1) = 0.
       END IF
      END DO  

      DO IMU0=ISOL,ISOL+1
       U0PL(IMU0,1) = SOLAR1/(2.0D0*PI*WTMU(IMU0))


C ****************************************************
C      CALCULATING INTERIOR INTENSITIES FOR CLOUD.
C         UPL(J,1,L) GOES DOWN OUT OF LAYER L.
C          UMI(J,1,L) GOES UP OUT OF LAYER L.
C ****************************************************
       DO L = 1,LTOT-1
        K = LTOT-L
        CALL IUP( RTOP(1,1,L), TTOP(1,1,L), JTOP(1,1,L),
     1    RBASE(1,1,K), TBASE(1,1,K), JBASE(1,1,K), UMI(1,1,L+1),
     2    NMU, MAXMU)
        CALL IDOWN( RTOP(1,1,L), TTOP(1,1,L), JTOP(1,1,L),
     1    RBASE(1,1,K), TBASE(1,1,K), JBASE(1,1,K), UPL(1,1,L),
     2    NMU, MAXMU)
       ENDDO


C *******************************************************
C     CALCULATING EXTERIOR INTENSITIES UTPL AND U0MI
C *******************************************************
       CALL MMUL(1.0D0,RBASE(1,1,LTOT),U0PL,ACOM,NMU,NMU,1,
     1   MAXMU,MAXMU,1)
       CALL MMUL(1.0D0,TBASE(1,1,LTOT),UTMI,BCOM,NMU,NMU,1,
     1   MAXMU,MAXMU,1)
       CALL MADD(1.0D0,ACOM,BCOM,ACOM,NMU,1,MAXMU,1)
       CALL MADD(1.0D0,JBASE(1,1,LTOT),ACOM,U0MI,NMU,1,MAXMU,1)

       CALL MMUL(1.0D0,TTOP(1,1,LTOT),U0PL,ACOM,NMU,NMU,1,
     1  MAXMU,MAXMU,1)
       CALL MMUL(1.0D0,RTOP(1,1,LTOT),UTMI,BCOM,NMU,NMU,1,
     1  MAXMU,MAXMU,1)
       CALL MADD(1.0D0,ACOM,BCOM,ACOM,NMU,1,MAXMU,1)
       CALL MADD(1.0D0,JTOP(1,1,LTOT),ACOM,UTPL(1,1,LT1),
     1  NMU,1,MAXMU,1)


       DO 302 IMU = 1, NMU

        JMU  = NMU+1-IMU

        DO 301 L=1,LT1
         IF(L.EQ.1)THEN 
          UMIF(JMU,L,IC+1)=U0MI(IMU,1,1)
         ENDIF
         IF(L.EQ.LT1)THEN 
          UPLF(JMU,L,IC+1)=UTPL(IMU,1,LT1)
         ENDIF
         IF(IMU0.EQ.ISOL)THEN
          UPLF(JMU,L,IC+1)=(1-FSOL)*SNGL(UPL(IMU,1,L))
          UMIF(JMU,L,IC+1)=(1-FSOL)*SNGL(UMI(IMU,1,L))
         ELSE
          UPLF(JMU,L,IC+1)=UPLF(JMU,L,IC+1)+FSOL*SNGL(UPL(IMU,1,L))
          UMIF(JMU,L,IC+1)=UMIF(JMU,L,IC+1)+FSOL*SNGL(UMI(IMU,1,L))
         ENDIF
301     CONTINUE

302    CONTINUE

       U0PL(IMU0,1)=0.D0

      END DO

1000  CONTINUE

      RETURN
      END

