      subroutine scloud8( rad, sol_ang, emiss_ang, aphi, radg, solar, 
     1lowbc, galb, mu1, wt1, nmu, nf, igdist, vwave, eps, bnu, tau,
     2nlay, ncont, lfrac, liscat, lcons, lncons, lnorm, nphi)
C     $Id: scloud8.f,v 1.5 2011-06-17 15:57:55 irwin Exp $
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
C     Optimised for maximum speed by cutting out any code from scloud6 which
C     is not actually used for NIMS retrieval runs.
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
C       eps(MAXSCATLAY)     real  Fraction of thermal absorption in each layer
C       bnu(MAXSCATLAY)     real  Mean Planck function in each layer
C       tau(MAXSCATLAY)     real  Total optical thickness of each layer
C       nlay            integer Number of layers
C       ncont           integer Number of aerosol types/modes
C       lfrac(MAXCON,MAXSCATLAY) real Fraction of scattering contributed by each
C                                  type in each layer
C       liscat(MAXCON)   integer Phase function type identifier
C       lncons(MAXCON)   integer Number of phase function parameters
C       lcons(MAXCON,10) real  Phase function parameters
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
      INTEGER IC,I1,J1,NCONT,NCONS,NORM,ISCAT,KL,NPHI,IP,JP,MF
      INTEGER IPOW0,L,IL,JN,KN,IMU0,IMU1,IDUMP
      DOUBLE PRECISION mu1(maxmu), wt1(maxmu), PI
      real rad, radf(0:MAXF-1), radg(maxmu), eps(MAXSCATLAY), 
     1 bnu(MAXSCATLAY), 
     2 tau(MAXSCATLAY),vwave,flux, dff, solar, lcons(MAXCON,MAXSCATPAR),
     3 lfrac(MAXCON,MAXSCATLAY),aphi,fin1,fout1,drad
      integer liscat(MAXCON),lnorm(MAXCON),lncons(MAXCON),nf,imu
      integer igdist

      DOUBLE PRECISION MU(MAXMU), WTMU(MAXMU), MM(MAXMU,MAXMU),
     1 MMINV(MAXMU,MAXMU), RL(MAXMU,MAXMU,MAXSCATLAY),
     2 TL(MAXMU,MAXMU,MAXSCATLAY), JL(MAXMU,1,MAXSCATLAY),
     3 RBASE(MAXMU,MAXMU,MAXSCATLAY),
     4 TBASE(MAXMU,MAXMU,MAXSCATLAY), JBASE(MAXMU,1,MAXSCATLAY),
     5 UMI(MAXMU,1,MAXSCATLAY), 
     6 ACOM(MAXMU,1), BCOM(MAXMU,1),
     7 PPLPL(MAXMU,MAXMU), PPLMI(MAXMU,MAXMU), CONS8(MAXSCATPAR),
     8 PPLPLS(MAXMU,MAXMU),PPLMIS(MAXMU,MAXMU)
      DOUBLE PRECISION PPLSTO(0:MAXSCATLAY,MAXMU,MAXMU,MAXSCATLAY)
      DOUBLE PRECISION PMISTO(0:MAXSCATLAY,MAXMU,MAXMU,MAXSCATLAY)
      DOUBLE PRECISION PPLN(MAXCON,MAXMU,MAXMU)
      DOUBLE PRECISION PMIN(MAXCON,MAXMU,MAXMU)
      DOUBLE PRECISION E(MAXMU,MAXMU),CCINV(MAXMU,MAXMU)
      DOUBLE PRECISION CC(MAXMU,MAXMU)
      DOUBLE PRECISION U0PL(MAXMU,1),UTMI(MAXMU,1),TAUT,BC,OMEGA
      DOUBLE PRECISION FRAC,SFRAC
      DOUBLE PRECISION STEST,GALB,ZMU0,ZMU

      real SOL_ANG,EMISS_ANG
      real YX(4),T,U,FEMM,FSOL
      INTEGER ICO,ISOL,IEMM

C     Common blocks
      COMMON/UNIT/ E
      COMMON/AREA1/CCINV,MM
      COMMON/AREA2/CC,MMINV,MU,PI
      COMMON/INPUT/U0PL,UTMI
      COMMON/HOMOG/ TAUT, BC, OMEGA, IPOW0
      COMMON/PHMAT/ PPLPL, PPLMI
      COMMON/PHASESTO/PPLSTO,PMISTO
      COMMON /SCATDUMP/ IDUMP
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet


C--------------------------------------------------------------------------


C	print*,'SCLOUD8 called OK'
C                        print*,sol_ang,emiss_ang,aphi
C                        print*,'radg',(radg(j),j=1,nmu)
C                        print*,solar,lowbc,galb
C                        print*,(mu1(j),j=1,nmu)
C                        print*,(wt1(j),j=1,nmu)
C                        print*,nf,Igdist,vwave
C                        print*,(eps(j),j=1,nlay)
C                        print*,(bnu(j),j=1,nlay)
C                        print*,(tau(j),j=1,nlay)
C                        print*,ncont
C                        print*,(liscat(j),j=1,ncont)
C                        print*,(lncons(j),j=1,ncont)
C                        print*,(lnorm(j),j=1,ncont)
C			do j=1,ncont
C                         print*,(lcons(j,i),i=1,lncons(j))
C                        enddo
C                        print*,nphi


      LTOT = NLAY	! Set internal number of layers
      LT1 = LTOT

C     In case of Lambertian reflection, add extra dummy layer at bottom,
C     whose transmission matrix = (1-A)*Unit-Matrix. This layer must be
C     omitted from computation by doubling
      if(igdist.eq.5.and.idump.gt.0)then
       print*,'tau(i),eps(i),bnu(i)'
       do i=1,nlay
        print*,i,tau(i),eps(i),bnu(i)
       enddo
      endif
      IF (LOWBC.EQ.1) LTOT = LTOT+1
      IF (LTOT.GT.MAXSCATLAY) CALL ABEND(' SCLOUD8: TOO MANY LAYERS')
      IF (NMU.GT.MAXMU) CALL ABEND(' SCLOUD8: TOO MANY ANGLE POINTS')

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


      DO 1000 IC=0,NF

C     PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
C     Calculate the Phase matrices for the aerosol types defined

      IF(IGDIST.EQ.1)THEN		! IGDIST=1 for the first k-coeff
C					  at a given wavelength

      DO J1 = 1,NCONT

        NCONS=LNCONS(J1)
        NORM=LNORM(J1)
        ISCAT=LISCAT(J1)


        DO KL=1,NCONS
          CONS8(KL)=LCONS(J1,KL)
        END DO
C       Calculate the phase function fourier coefficients for each particle
C       type and each illumination/viewing zenith angles and for either
C	light going through layer (PLPL) or being reflected from it (PLMI)
        CALL CALC_PMAT6(NF, IC, PPLPL, PPLMI, MU, WTMU, 
     1    NMU, ISCAT, CONS8, NCONS, NORM, J1, NCONT, VWAVE, NPHI)
C       Transfer matrices to those for wach scattering particle 
        DO IP = 1,NMU
         DO JP = 1,NMU
          PPLN(J1,IP,JP)=PPLPL(IP,JP)
          PMIN(J1,IP,JP)=PPLMI(IP,JP)
         END DO
        END DO

      END DO

      END IF

C     PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP


C **********************************************************************
C     CALCULATE REFLECTION,TRANSMISSION AND SOURCE MATRICES FOR EACH    
C	  HOMOGENEOUS LAYER L; INSERT INTO ARRAY X(I,J,L).	  	      
C **********************************************************************    
C
      IPOW0 = 16
      DO 2000 L = 1,LT1
	TAUT = 1.0D0*TAU(L)
	BC = 1.0D0*BNU(L)
	OMEGA = 1.0D0*(1. - EPS(L))

        IF(TAUT.LT.0.0)THEN
         IF(IDIAG.GT.0)THEN
          PRINT*,'Error in scloud8. TAUT < 0. Setting to zero'
          PRINT*,'L,TAUT,BC,OMEGA'
          PRINT*,L,TAUT,BC,OMEGA
         ENDIF
         TAUT = 0.0
C         STOP
        END IF


C       P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P*P
C
C       Calculate phase function coefficients for each atmospheric layer
C       averaged over particle types and relative concentrations

        IF(IGDIST.EQ.1) THEN
C       For first kcoef, calculate the averaged phase function fourier
C	coefficients.
         DO J1=1,NMU
          DO KL=1,NMU
           PPLPLS(J1,KL)=0.
           PPLMIS(J1,KL)=0.
          END DO
         END DO

         SFRAC=0.D0
         DO J1 = 1,NCONT

         FRAC = LFRAC(J1,L)
         SFRAC = SFRAC+FRAC

         DO IP = 1,NMU
          DO JP = 1,NMU
           PPLPL(IP,JP) = PPLN(J1,IP,JP)
           PPLMI(IP,JP) = PMIN(J1,IP,JP)
          END DO
         END DO


         CALL MADD(FRAC,PPLPLS,PPLPL,PPLPLS,NMU,NMU,
     1    MAXMU,MAXMU)
         CALL MADD(FRAC,PPLMIS,PPLMI,PPLMIS,NMU,NMU,
     1    MAXMU,MAXMU)

         END DO

         STEST=ABS(SFRAC-1.D0)
         IF(STEST.GT.0.02D0)THEN
          PRINT*,'Scloud8. ERROR.  SUM(FRAC) must = 1.'
          print*,sfrac
          STOP
         END IF

C        Transfer averaged phase fourier coefficient matrices for each layer 
C	 to 'holding' matrices to prevent repeated calculation for IGDIST <> 1
         DO I1=1,NMU
          DO J1=1,NMU
           PPLSTO(IC,I1,J1,L)=PPLPLS(I1,J1)
           PMISTO(IC,I1,J1,L)=PPLMIS(I1,J1)
          END DO
         END DO
         CALL MEQU(PPLPL,NMU,MAXMU,PPLPLS)
         CALL MEQU(PPLMI,NMU,MAXMU,PPLMIS)

        ELSE

C        If IGDIST <> 1, then simply read phase fourier coefficient matrices
C	 from 'holding matrices'
         DO I1=1,NMU
          DO J1=1,NMU
           PPLPL(I1,J1)=PPLSTO(IC,I1,J1,L)
           PPLMI(I1,J1)=PMISTO(IC,I1,J1,L)
          END DO
         END DO

        END IF

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

        ELSE IF(OMEGA.EQ.0)THEN
         do i1=1,nmu
          do j1=1,nmu
           RL(I1,J1,L)=0.0
           TL(I1,J1,L)=0.0
          end do
          TL(i1,i1,L)=EXP(-TAUT*MMINV(i1,i1))
          JL(i1,1,L)= BC*(1.0 - TL(i1,i1,L))
         end do

        ELSE
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
	    RL(I,J,LTOT) = 2.0D0*GALB*MU(J)*WTMU(J)  !Sum of MU*WTMU = 0.5
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
	K = LTOT-L
	CALL ADD( RL(1,1,K), TL(1,1,K), JL(1,1,K), RBASE(1,1,L), 
     1 TBASE(1,1,L), JBASE(1,1,L), RBASE(1,1,L+1), 
     2 TBASE(1,1,L+1), JBASE(1,1,L+1), NMU, MAXMU)
      ENDDO


      IF(IC.NE.0)THEN
       DO IL=1,NMU
        DO J1=1,LTOT
         JBASE(IL,1,J1)=0.
        END DO
       END DO
      END IF
     

      ZMU0 = DCOS(SOL_ANG*PI/180.0)
      ZMU = DCOS(EMISS_ANG*PI/180.0)

      ISOL=-1
      IEMM=-1
      DO J=1,NMU-1
       IF(ZMU0.LE.MU(J).AND.ZMU0.GT.MU(J+1))ISOL = J
       IF(ZMU.LE.MU(J).AND.ZMU.GT.MU(J+1))IEMM = J
      END DO

      IF(ISOL.LE.0)THEN
       IF(ZMU0.GT.MU(1))THEN
        ISOL=1
        FSOL=0.0
       ELSEIF(ZMU0.LT.MU(NMU))THEN
        ISOL=NMU-1
        FSOL=1.0
       ENDIF
      ELSE
       FSOL = SNGL((MU(ISOL)-ZMU0)/(MU(ISOL)-MU(ISOL+1)))
      ENDIF


      IF(IEMM.LE.0)THEN
       IF(ZMU.GT.MU(1))THEN
        IEMM=1
        FEMM=0.0
       ELSEIF(ZMU.LT.MU(NMU))THEN
        IEMM=NMU-1
        FEMM=1.0
       ENDIF
      ELSE
       FEMM = SNGL((MU(IEMM)-ZMU)/(MU(IEMM)-MU(IEMM+1)))
      ENDIF


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
       U0PL(IMU0,1) = SOLAR/(2.0D0*PI*WTMU(IMU0))


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
   
       
1000  CONTINUE

      IF(IGDIST.EQ.5.AND.IDUMP.GT.0)THEN
        PRINT*,'IGDIST=5: Fcoeff, rad, drad, %change'
      ENDIF
      RAD=0.0
      DO IC=0,NF
       DRAD = RADF(IC)*COS(IC*APHI*SNGL(PI)/180.0)
       RAD=RAD+DRAD
       IF(IGDIST.EQ.5.AND.IDUMP.GT.0)PRINT*,IC,RAD,DRAD,100*DRAD/RAD
      END DO

C
      RETURN
      END

