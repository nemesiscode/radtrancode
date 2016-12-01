      program test_pmat6

C ***********************************************************************
C     CALCULATE P++ AND P+- (Plass et al 1973)
C ***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '../includes/arrdef.f'
      INTEGER ICONT,NCONT,NPHI
      REAL*8 MU(MAXMU), WTMU(MAXMU),RL(MAXMU,MAXMU,MAXSCATLAY),
     1 TL(MAXMU,MAXMU,MAXSCATLAY),RTOP(MAXMU,MAXMU,MAXSCATLAY),
     2 PPLPL(MAXMU,MAXMU), PPLMI(MAXMU,MAXMU), CONS8(MAXSCATPAR)
      REAL*8 PTPL(MAXCON,MAXF,MAXMU,MAXMU),PTMI(MAXCON,MAXF,MAXMU,MAXMU)
      REAL*4 VWAVE

      NMU=5
      MU(1)=0.165278957666387
      MU(2)=0.477924949810444
      MU(3)=0.738773865105505
      MU(4)=0.919533908166459
      MU(5)=1.00000000000000 
      WTMU(1)=0.327539761183898
      WTMU(2)=0.292042683679684
      WTMU(3)=0.224889342063117
      WTMU(4)=0.133305990851069
      WTMU(5)=2.222222222222220E-002

      print*,'Enter IC, NF'
      READ*,IC,NF

      ISCAT=2
      NCONS=3
      CONS8(1)=0.5
      CONS8(2)=0.3
      CONS8(3)=-0.3
      NORM=1
      NPHI=100
      NCONT=5
      ICONT=1

      CALL CALC_PMAT6(NF, IC, PPLPL, PPLMI, MU, WTMU,
     1 NMU, ISCAT, CONS8, NCONS, NORM, ICONT, NCONT, VWAVE, NPHI)

      PRINT*,'PPLPL: HG RAY'
      DO I=1,NMU
       PRINT*,(PPLPL(I,J),J=1,NMU)
      ENDDO

      PRINT*,'PPLMI: HG RAY'
      DO I=1,NMU
       PRINT*,(PPLMI(I,J),J=1,NMU)
      ENDDO


      CONS8(1)=1.0
      CONS8(2)=0.0
      CONS8(3)=0.0
      ICONT=2

      CALL CALC_PMAT6(NF, IC, PPLPL, PPLMI, MU, WTMU,
     1 NMU, ISCAT, CONS8, NCONS, NORM, ICONT, NCONT, VWAVE, NPHI)

      PRINT*,'PPLPL: HG ISO 1'
      DO I=1,NMU
       PRINT*,(PPLPL(I,J),J=1,NMU)
      ENDDO

      PRINT*,'PPLMI: HG ISO 1'
      DO I=1,NMU
       PRINT*,(PPLMI(I,J),J=1,NMU)
      ENDDO


      CONS8(1)=0.0
      CONS8(2)=0.0
      CONS8(3)=0.0
      ICONT=3
      CALL CALC_PMAT6(NF, IC, PPLPL, PPLMI, MU, WTMU,
     1 NMU, ISCAT, CONS8, NCONS, NORM, ICONT, NCONT, VWAVE, NPHI)

      PRINT*,'PPLPL: HG ISO 2'
      DO I=1,NMU
       PRINT*,(PPLPL(I,J),J=1,NMU)
      ENDDO

      PRINT*,'PPLMI: HG ISO 2'
      DO I=1,NMU
       PRINT*,(PPLMI(I,J),J=1,NMU)
      ENDDO


      ISCAT = 0
      NCONS = 0
      NPHI = 100
      ICONT=4

      CALL CALC_PMAT6(NF, IC, PPLPL, PPLMI, MU, WTMU,
     1 NMU, ISCAT, CONS8, NCONS, NORM, ICONT, NCONT, VWAVE, NPHI)

      PRINT*,'PPLPL: Rayleigh'
      DO I=1,NMU
       PRINT*,(PPLPL(I,J),J=1,NMU)
      ENDDO

      PRINT*,'PPLMI: Rayleigh'
      DO I=1,NMU
       PRINT*,(PPLMI(I,J),J=1,NMU)
      ENDDO

      ISCAT = 1
      NCONS = 0
      NPHI = 100
      ICONT=5

      CALL CALC_PMAT6(NF, IC, PPLPL, PPLMI, MU, WTMU,
     1 NMU, ISCAT, CONS8, NCONS, NORM, ICONT, NCONT, VWAVE, NPHI)

      PRINT*,'PPLPL: ISO'
      DO I=1,NMU
       PRINT*,(PPLPL(I,J),J=1,NMU)
      ENDDO

      PRINT*,'PPLMI: ISO'
      DO I=1,NMU
       PRINT*,(PPLMI(I,J),J=1,NMU)
      ENDDO

      END
