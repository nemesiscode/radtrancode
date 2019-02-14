C**********************************************************************
C
      SUBROUTINE IUP(RA,TA,JA,RB,TB,JB,UMI,NMU,JDIM)
C----------------------------------------------------
C     Original		Lucas Kamp	1990s
C     Documented	Pat Irwin	2/7/07
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../includes/arrdef.f'
      DOUBLE PRECISION RA(JDIM,JDIM), TA(JDIM,JDIM), JA(JDIM,1),
     1 RB(JDIM,JDIM), TB(JDIM,JDIM), JB(JDIM,1), UMI(JDIM,1), 
     2 ACOM(MAXMU,MAXMU), BCOM(MAXMU,MAXMU), XCOM(MAXMU,1),
     3 YCOM(MAXMU,1), WKSPCE(MAXMU), AA(MAXMU,MAXMU), BB(MAXMU,MAXMU)
      DIMENSION L(MAXMU), M(MAXMU)
      COMMON/UNIT/ E(MAXMU,MAXMU)
      COMMON/INPUT/U0PL(MAXMU,1), UTMI(MAXMU,1)
C
      IF (JDIM.NE.MAXMU) CALL ABEND(' IUP: DIMENSION ERROR')
C     See Plass et al.(1993), Apl. Opt. 12, pp 314-329.
C
C     This is equation 5
C	RA = R10
C	RB = R12
C	TA = T01
C       TB = T21
C	JA = JP01
C       JB = JM21
C	UMI = I2-
C       U0PL(common) is I0+
C	UTMI(common) is I2-
C
C       Output UMI is I1-
C
C     calculate r12*r10 -> ACOM
      CALL MMUL(1.0D0,RB,RA,ACOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)

C     Calculate (E - r12*r10) -> BCOM
      CALL MADD(-1.0D0,E,ACOM,BCOM,NMU,NMU,JDIM,JDIM)

C     Calculate (E - r12*r10)^-1 -> ACOM
      CALL MATINV8(BCOM,NMU,JDIM,ACOM)

C     Transfer result to BCOM
      CALL MEQU(BCOM,NMU,JDIM,ACOM)

C     Calculate t21*I2- -> XCOM
      CALL MMUL(1.0D0,TB,UTMI,XCOM,NMU,NMU,1,JDIM,JDIM,1)

C     Calculate r12*t01 -> ACOM
      CALL MMUL(1.0D0,RB,TA,ACOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)

C     Calculate r12*t01*I0+ -> YCOM
      CALL MMUL(1.0D0,ACOM,U0PL,YCOM,NMU,NMU,1,JDIM,JDIM,1)

C     Add: t21*I2- + r12*t01*I0+ -> XCOM
      CALL MADD(1.0D0,XCOM,YCOM,XCOM,NMU,1,JDIM,1)

C     Calculate r12*J01+ -> YCOM
      CALL MMUL(1.0D0,RB,JA,YCOM,NMU,NMU,1,JDIM,JDIM,1)

C     Add total and put in UMI
      CALL MADD(1.0D0,XCOM,YCOM,UMI,NMU,1,JDIM,1)

C     Add J21- to UMI
      CALL MADD(1.0D0,UMI,JB,XCOM,NMU,1,JDIM,1)

C     Multiply and put result in UMI
      CALL MMUL(1.0D0,BCOM,XCOM,UMI,NMU,NMU,1,JDIM,JDIM,1)
      RETURN
      END
