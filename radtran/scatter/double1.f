      SUBROUTINE DOUBLE1(IC,L,RNEXT,TNEXT,JNEXT,NMU,JDIM)
C     $Id: double1.f,v 1.2 2008-04-09 11:35:12 irwin Exp $
C **********************************************************************
C
C     Subroutine to calculate the reflection, transmission and source 
C     matrices for a single homogeneous layer using the doubling method.
C
C     The formalism used here is that of Matrix Operators as defined by the 
C     reference: Plass et al.,  Appl. Optics 12, pp 314 (1973) 
C				(NIMS/PAPER/SCAT/1973/2)
C     Input Variables:
C	IC	INTEGER		Coefficient of azimuth expansion
C	L	INTEGER		Layer ID( not used)
C	NMU	INTEGER		Number of cosine angles
C	JDIM    INTEGER		Principal size of output arrays
C       Other inputs come through COMMON Blocks
C
C     Output variables
C	RNEXT(JDIM,JDIM)	DOUBLE	Reflection matrix of layer
C	TNEXT(JDIM,JDIM)	DOUBLE	Transmission matrix of layer
C	JNEXT(JDIM,1)		DOUBLE	Source matrix of layer
C
C      Pat Irwin	17/9/96
C
C **********************************************************************    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../includes/arrdef.f'
      DOUBLE PRECISION OMEGA, PPLPL(MAXMU,MAXMU), PPLMI(MAXMU,MAXMU),
     1 ACOM(MAXMU,MAXMU), BCOM(MAXMU,MAXMU), GPLPL(MAXMU,MAXMU),
     2 GPLMI(MAXMU,MAXMU), TINIT(MAXMU,MAXMU), RINIT(MAXMU,MAXMU),
     3 T1(MAXMU,MAXMU), R1(MAXMU,MAXMU), RNEXT(JDIM,JDIM),
     4 TNEXT(JDIM,JDIM), JNEXT(JDIM,1), MM(MAXMU,MAXMU),
     5 MMINV(MAXMU,MAXMU), MU(MAXMU), J1(MAXMU,1)
      INTEGER IC
      COMMON/UNIT/ E(MAXMU,MAXMU)
      COMMON/AREA1/ CCINV(MAXMU,MAXMU), MM
      COMMON/AREA2/  CC(MAXMU,MAXMU), MMINV, MU, PI
      COMMON/HOMOG/ TAUT, BC, OMEGA, IPOW0
      COMMON/PHMAT/ PPLPL, PPLMI
      integer idiag,iquiet,IRAMAN,NRAMAN
      common/diagnostic/idiag,iquiet
C
      LOGICAL RAMAN
      REAL DENS(MAXSCATLAY),FPRAMAN(MAXSCATLAY)
      REAL JRAMAN(1000,MAXSCATLAY),VRAM0,VRAMST,VTMP
      REAL H2ABUND(MAXSCATLAY),KWT
      COMMON/RAMAN/RAMAN,DENS,JRAMAN,VRAM0,VRAMST,NRAMAN,FPRAMAN,
     1 H2ABUND,IRAMAN,KWT

      IF (JDIM.NE.MAXMU) CALL ABEND(' DOUBLE: DIMENSION ERROR')



      sum=0
      do I=1,nmu
       sum=sum+cc(i,i)
      end do
      dsum = abs(sum - 1.0)
      if(dsum.gt.0.02.and.idiag.gt.0)then
       print*,'Double1: Warning - Sum of weights <> 1.0'
       print*,sum
      endif
      
      IF(IC.EQ.0)THEN

        DO J=1,NMU
         SUM= 0
         DO I=1,NMU
           SUM=SUM+(PPLPL(I,J)+PPLMI(I,J))*CC(I,I)
         END DO
         DSUM = ABS(SUM*2.0*PI - 1.0)
         IF(DSUM.GT.0.02)THEN
          if(idiag.gt.0)then
           PRINT*,'Double1: Warning - Sum of phase function <> 1'
 	  PRINT*,'IC,L  = ',IC,L
           PRINT*,'J,SUM = ',J,SUM*2.0*PI
           PRINT*,'PPLPL'
           DO I=1,NMU
            PRINT*,(PPLPL(I,K),K=1,NMU)
           ENDDO
           PRINT*,'PPLMI'
           DO I=1,NMU
            PRINT*,(PPLMI(I,K),K=1,NMU)
           ENDDO

           PRINT*,'Doing a brutal renormalisation'
          endif

          XFAC=1./(2.0*PI*SUM)
          DO I=1,NMU
           PPLPL(I,J)=PPLPL(I,J)*XFAC
           PPLMI(I,J)=PPLMI(I,J)*XFAC
          END DO

         ENDIF
        END DO
      ENDIF 



C     ****************************************************************
C     First calculate properties of initial very thin layer.
C
C     ******************************
C     COMPUTATION OF GAMMA++	(Plass et al. 1973)
C     GPLPL = MMINV*(E - CON*PPLPL*CC)
C     ******************************
      CON=OMEGA*PI

      DEL01 = 0.0D0
      IF(IC.EQ.0)DEL01 = 1.0D0

      CON=CON*(1.0D0 + DEL01)

      CALL MMUL(CON,PPLPL,CC,ACOM,NMU,NMU,NMU,MAXMU,MAXMU,MAXMU)
      CALL MADD(-1.0D0,E,ACOM,BCOM,NMU,NMU,MAXMU,MAXMU)
      CALL MMUL(1.0D0,MMINV,BCOM,GPLPL,NMU,NMU,NMU,MAXMU,MAXMU,MAXMU)

C      print*,'GPLPL'
C      DO I=1,NMU
C       print*,(GPLPL(I,J),J=1,NMU)
C      ENDDO   
C101   format(' ',5(e12.3))

C
C     ******************************
C     COMPUTATION OF GAMMA+-   
C     GPLMI = MMINV*CON*PPLMI*CC
C     ******************************
      CALL MMUL(CON,PPLMI,CC,ACOM,NMU,NMU,NMU,MAXMU,MAXMU,MAXMU)
      CALL MMUL(1.0D0,MMINV,ACOM,GPLMI,NMU,NMU,NMU,MAXMU,MAXMU,MAXMU)


C      print*,'GPLMI'
C      DO I=1,NMU
C       print*,(GPLMI(I,J),J=1,NMU)
C      ENDDO

C     N.b: GAMMA++ = GAMMA--
C	   GAMMA+- = GAMMA-+

C      print*,'TAUT,IPOW0 = ',TAUT,IPOW0
      NN=INT(DLOG(TAUT)/DLOG(2.0D0))+IPOW0

      IF(NN.GE.1)THEN
       XFAC=1.0/(2.0D0**NN)
      ELSE
       XFAC=1.0
      ENDIF
      TAU0 = TAUT*XFAC

C      PRINT*,'TAU0, NN',TAU0,NN

C **********************************************************************
C     COMPUTATION OF R, T AND J FOR INITIAL LAYER (SINGLE SCATTERING)
C **********************************************************************
      CALL MADD(-TAU0,E,GPLPL,TINIT,NMU,NMU,MAXMU,MAXMU)
      CALL MMUL(TAU0,E,GPLMI,RINIT,NMU,NMU,NMU,MAXMU,MAXMU,MAXMU)
      DO J=1,NMU
       J1(J,1)=0.0
      ENDDO
      IF(IC.EQ.0)THEN
       DO 80 J=1,NMU
C       If Raman scattering is on, add radiation scattered from shorter 
C       wavelengths here. Raman-scattered radiation is assumed to be isotropic 
C       and so is added in the same way as thermal emission.
        IF(RAMAN.AND.(IRAMAN.GE.1).AND.(IRAMAN.LE.NRAMAN))THEN
         J1(J,1)=(1.0D0-OMEGA)*BC*TAU0*MMINV(J,J)+
     1     XFAC*JRAMAN(IRAMAN,L)*MMINV(J,J)
        ELSE
 	 J1(J,1)=(1.0D0-OMEGA)*BC*TAU0*MMINV(J,J)
        ENDIF
   80  ENDDO
      ENDIF
C      print*,'Tinit' 
C      DO I=1,NMU
C       print*,(TINIT(I,J),J=1,NMU)
C      ENDDO

C      print*,'Rinit'
C      DO I=1,NMU
C       print*,(RINIT(I,J),J=1,NMU)
C      ENDDO

      IF(NN.LT.1)THEN

C       WRITE(*,*)'DOUBLE1 : TAUT,NN,TAU0',TAUT,NN,TAU0
C       WRITE(*,*)'Optical thickness < 2**(-(IPOW0-1)) => single scat'

       CALL MEQU(TNEXT,NMU,MAXMU,TINIT)
       CALL MEQU(RNEXT,NMU,MAXMU,RINIT)
       DO 81 J=1,NMU
	JNEXT(J,1)=J1(J,1)
81     ENDDO
       RETURN 
      END IF

      CALL MEQU(T1,NMU,MAXMU,TINIT)
      CALL MEQU(R1,NMU,MAXMU,RINIT)


C
C **************************************************************
C     COMPUTATION OF R AND T FOR SUBSEQUENT LAYERS: DOUBLING    
C **************************************************************
      DO 100 N=1,NN
 	CALL ADD(R1,T1,J1,R1,T1,J1,RNEXT,TNEXT,JNEXT,NMU,JDIM)
	TAUL=TAU0*(2.0D0**N)

C        PRINT*,'TAUL',TAUL
C        print*,'T1'
C        DO I=1,NMU
C         print*,(T1(I,J),J=1,NMU)
C        ENDDO
C        print*,'R1'
C        DO I=1,NMU  
C         print*,(T1(I,J),J=1,NMU)
C        ENDDO

	IF(TAUL.EQ.TAUT) GOTO 100
	CALL MEQU(R1,NMU,MAXMU,RNEXT)
	CALL MEQU(T1,NMU,MAXMU,TNEXT)
	DO 85 J=1,NMU
	  J1(J,1)=JNEXT(J,1)
   85 	ENDDO
  100 ENDDO
C
      RETURN
      END
