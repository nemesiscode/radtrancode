      SUBROUTINE CALC_PMAT6(NF, IC, PPLPL, PPLMI, MU, WTMU,
     1 NMU, ISCAT, CONS8, NCONS, NORM, ICONT, NCONT, VWAVE, NPHI)
C     $Id: calc_pmat6.f,v 1.4 2011-06-17 15:57:52 irwin Exp $
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
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C     Common block PSHARE to store all the fourier components
      COMMON/PSHARE/ PTPL,PTMI
      COMMON/UNIT/ E(MAXMU,MAXMU)

      PI = 4.0D0*DATAN(1.0D0)

      IF(IC.EQ.0)THEN

       CALL PHASINT2( NF, MU, NMU, NPHI, ISCAT,
     1 CONS8, NCONS, ICONT, NCONT, VWAVE)


C       do i=1,nmu
C        tran = 0.0
C        refl = 0.0
C        do j=1,nmu
C         tran = tran + wtmu(j)*ptpl(icont,1,i,j)*2.0*pi
C         refl = refl + wtmu(j)*ptmi(icont,1,i,j)*2.0*pi
C        end do
C        err = 100*abs(1.0 - (tran+refl))
C        if(err.gt.10.0)then    
C         print*,'Calc_pmat6: prenormalisation zenith quadrature
C     & > 10% off:'
C         write(*,1000)icont,i,tran,refl,tran+refl
C1000  format('      Scat : ',i2,' imu : ',i2,' T = ',f8.3,' R = ',
C     &f8.3,' T+R = ',f8.3)
C        endif
C       end do

      ENDIF
 
      
      DO I=1,NMU
       DO J=1,NMU
        PPLPL(I,J)=PTPL(ICONT,IC+1,I,J)
        PPLMI(I,J)=PTMI(ICONT,IC+1,I,J)
       END DO
      END DO


      IF (NORM.EQ.1) THEN
	CALL HANSEN( IC, PPLPL, PPLMI, MAXMU, WTMU, NMU)
      ELSEIF (NORM.EQ.2) THEN
	if(idiag.gt.0)PRINT*,'CALC_PMAT6. NORM=2 Option disabled'
      ENDIF


      RETURN

      END

