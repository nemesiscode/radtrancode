      SUBROUTINE ADDP(R1,T1,J1,ISCAT1,RSUB,TSUB,JSUB,RANS,TANS,JANS,
     &NMU,JDIM)
C     $Id: addp.f,v 1.2 2011-06-17 15:57:52 irwin Exp $
C     *****************************************************************
C
C     Subroutine to add the diffuse reflection, transmission and reflection
C     matrices for two adjacent atmospheric layers
C
C     Input variables:
C	R1(JDIM,JDIM)	DOUBLE	Diffuse reflection operator for 2nd layer
C	T1(JDIM,JDIM)	DOUBLE	Diffuse transmission operator for 2nd layer
C	J1(JDIM,1)	DOUBLE	Diffuse source function for 2nd layer
C	ISCAT1		INTEGER Flag to indicate if 2nd layer is scattering
C	RSUB(JDIM,JDIM)	DOUBLE	Diffuse reflection operator for existing layer
C	TSUB(JDIM,JDIM)	DOUBLE	Diffuse transmission operator for existing layer
C	JSUB(JDIM,1)	DOUBLE	Diffuse source function for existing layer
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
C     Subroutine solves Eq. 7b,8b,9b of Plass et al. (1973) 
C     Here :
C        R1 is R10 for the homogenous layer 1 being added 
C                    (and thus equal to R01)
C        T1 is T10 for the homogenous layer 1 being added
C		     (and thus equal to T10)
C        J1 is the source function for the homegenous layer 1
C                    (and thus same in +ve and -ve directions)
C                    +ve is going in direction from layer 1 to 2nd layer
C 		     -ve is going in direction from 2nd layer to layer 1
C	 RSUB is R12 for the 2nd layer (homegenous or composite)
C	 TSUB is T21 for the 2nd layer (homegenous or composite)
C        JSUB is source function for 2nd layer (homegenous or composite)
C		     going in -ve direction (i.e. JM21)
C        RANS is R02 for combined layers
C	 TANS is T20 for combineds layers
C        JANS is JM20 for combined layers (i.e. in -ve direction)

      IF(ISCAT1.EQ.1)THEN
C      2nd layer is scattering. Solve Eq. 7b,8b,9b of Plass et al. (1973)
       CALL MMUL(-1.0D0,RSUB,R1,BCOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)
C       BCOM=-RSUB*R1
       CALL MADD(1.0D0,E,BCOM,BCOM,NMU,NMU,JDIM,JDIM)
C       BCOM=E-RSUB*R1
       CALL MATINV8(BCOM,NMU,JDIM,ACOM)
C       ACOM = INV(E-RSUB*R1)
       CALL MEQU(BCOM,NMU,JDIM,ACOM)
C       BCOM=INV(E-RSUB*R1)
       CALL MMUL(1.0D0,T1,BCOM,CCOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)
C       CCOM=T1*INV(E-RSUB*R1)
       CALL MMUL(1.0D0,CCOM,RSUB,RANS,NMU,NMU,NMU,JDIM,JDIM,JDIM)
C       RANS=T1*INV(E-RSUB*R1)*RSUB
       CALL MMUL(1.0D0,RANS,T1,ACOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)
C       ACOM=T1*INV(E-RSUB*R1)*RSUB*T1
       CALL MADD(1.0D0,R1,ACOM,RANS,NMU,NMU,JDIM,JDIM)
C       RANS = R1+T1*INV(E-RSUB*R1)*RSUB*T1
       CALL MMUL(1.0D0,CCOM,TSUB,TANS,NMU,NMU,NMU,JDIM,JDIM,JDIM)
C       TANS=T1*INV(E-RSUB*R1)*TSUB
       CALL MMUL(1.0D0,RSUB,J1,JCOM,NMU,NMU,1,JDIM,JDIM,1)
C       JCOM=RSUB*J1
       CALL MADD(1.0D0,JSUB,JCOM,JCOM,NMU,1,JDIM,1)
C       JCOM=JSUB+RSUB*J1
       CALL MMUL(1.0D0,CCOM,JCOM,JANS,NMU,NMU,1,JDIM,JDIM,1)
C       JANS=T1*INV(E-RSUB*R1)*(JSUB+RSUB*J1)
       CALL MADD(1.0D0,J1,JANS,JANS,NMU,1,JDIM,1)
C       JANS = J1+T1*INV(E-RSUB*R1)*(JSUB+RSUB*J1)  
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
        JANS(I,1) = J1(I,1) + TA*JANS(I,1)
       END DO

      ENDIF


      RETURN
      END
