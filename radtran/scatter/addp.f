      SUBROUTINE ADDP(R1,T1,J1,ISCAT1,RSUB,TSUB,JSUB,RANS,TANS,JANS,
     &NMU,JDIM)
C     $Id: addp.f,v 1.2 2011-06-17 15:57:52 irwin Exp $
C     *****************************************************************
C
C     Subroutine to add the diffuse reflection, transmission and reflection
C     matrices for two adjacent atmospheric layers
C
C     Input variables:
C	R1(JDIM,JDIM)	DOUBLE	Diffuse reflection operator for 1st layer
C	T1(JDIM,JDIM)	DOUBLE	Diffuse transmission operator for 1st layer
C	J1(JDIM,1)	DOUBLE	Diffuse source function for 1st layer
C	ISCAT1		INTEGER Flag to indicate if 2nd layer is scattering
C	RSUB(JDIM,JDIM)	DOUBLE	Diffuse reflection operator for 2nd layer
C	TSUB(JDIM,JDIM)	DOUBLE	Diffuse transmission operator for 2nd layer
C	JSUB(JDIM,1)	DOUBLE	Diffuse source function for 2nd layer
C	JDIM		INTEGER	Array size
C	NMU		INTEGER	Number of elements used
C
C     Output variables
C	RANS(JDIM,JDIM)	DOUBLE	Combined diffuse reflection operator
C	TANS(JDIM,JDIM)	DOUBLE	Combined diffuse transmission operator
C	JANS(JDIM,1)	DOUBLE	Combined diffuse source function
C
C     Pat Irwin		17/9/96
C
C**********************************************************************


C---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../includes/arrdef.f'
      COMMON/UNIT/ E(MAXMU,MAXMU)
      DOUBLE PRECISION CCOM(MAXMU,MAXMU), JCOM(MAXMU,1), WKSPCE(MAXMU),
     1 AA(MAXMU,MAXMU), BB(MAXMU,MAXMU), ACOM(MAXMU,MAXMU),
     2 BCOM(MAXMU,MAXMU)
      INTEGER L(MAXMU), M(MAXMU)
      DOUBLE PRECISION R1(JDIM,JDIM), T1(JDIM,JDIM), J1(JDIM,1), 
     1 RSUB(JDIM,JDIM), TSUB(JDIM,JDIM), JSUB(JDIM,1), 
     2 RANS(JDIM,JDIM), TANS(JDIM,JDIM), JANS(JDIM,1)
      INTEGER ISCAT1
C
      IF (JDIM.NE.MAXMU) CALL ABEND(' ADD: DIMENSION ERROR')
C
      IF(ISCAT1.EQ.1)THEN
C      2nd layer is scattering.
       CALL MMUL(-1.0D0,RSUB,R1,BCOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)
       CALL MADD(1.0D0,E,BCOM,BCOM,NMU,NMU,JDIM,JDIM)
       CALL MATINV8(BCOM,NMU,JDIM,ACOM)
       CALL MEQU(BCOM,NMU,JDIM,ACOM)
       CALL MMUL(1.0D0,T1,BCOM,CCOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)
       CALL MMUL(1.0D0,CCOM,RSUB,RANS,NMU,NMU,NMU,JDIM,JDIM,JDIM)
       CALL MMUL(1.0D0,RANS,T1,ACOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)
       CALL MADD(1.0D0,R1,ACOM,RANS,NMU,NMU,JDIM,JDIM)
       CALL MMUL(1.0D0,CCOM,TSUB,TANS,NMU,NMU,NMU,JDIM,JDIM,JDIM)
       CALL MMUL(1.0D0,RSUB,J1,JCOM,NMU,NMU,1,JDIM,JDIM,1)
       CALL MADD(1.0D0,JSUB,JCOM,JCOM,NMU,1,JDIM,1)
       CALL MMUL(1.0D0,CCOM,JCOM,JANS,NMU,NMU,1,JDIM,JDIM,1)
       CALL MADD(1.0D0,J1,JANS,JANS,NMU,1,JDIM,1)
  
      ELSE
C       2nd layer is non-scattering
       CALL MMUL(1.0D0,RSUB,J1,JCOM,NMU,NMU,1,JDIM,JDIM,1)
       CALL MADD(1.0D0,JSUB,JCOM,JCOM,NMU,1,JDIM,1)

       DO I=1,NMU
        TA = T1(I,I)
        DO J=1,NMU
         TB = T1(J,J)
         TANS(I,J) = TSUB(I,J)*TA
         RANS(I,J) = RSUB(I,J)*TA*TB
        END DO
        JANS(I,1) = J1(I,1) + TA*JCOM(I,1)
       END DO

      ENDIF


      RETURN
      END
