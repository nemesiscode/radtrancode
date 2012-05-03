C**********************************************************************
C
      SUBROUTINE IDOWN(RA,TA,JA,RB,TB,JB,UPL,NMU,JDIM)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../includes/arrdef.f'
      DOUBLE PRECISION RA(JDIM,JDIM), RB(JDIM,JDIM), TA(JDIM,JDIM),
     1 TB(JDIM,JDIM),JA(JDIM,1), JB(JDIM,1), ACOM(MAXMU,MAXMU),
     2  BCOM(MAXMU,MAXMU), UPL(JDIM,1), XCOM(MAXMU,1), YCOM(MAXMU,1),
     3  WKSPCE(MAXMU), AA(MAXMU,MAXMU), BB(MAXMU,MAXMU)
      DIMENSION L(MAXMU), M(MAXMU)
      COMMON/UNIT/ E(MAXMU,MAXMU)
      COMMON/INPUT/U0PL(MAXMU,1), UTMI(MAXMU,1)
C
      IF (JDIM.NE.MAXMU) CALL ABEND(' IDOWN: DIMENSION ERROR')
C
C     See Plass et al.(1993), Apl. Opt. 12, pp 314-329.
C
C     This is equation 6
C
C
C     Calculate r10*r12
      CALL MMUL(1.0D0,RA,RB,ACOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)

C     Calculate E-r10*r12
      CALL MADD(-1.0D0,E,ACOM,BCOM,NMU,NMU,JDIM,JDIM)

C     Calculate (E-r10*r12)^-1 -> ACOM
      CALL MATINV8(BCOM,NMU,JDIM,ACOM)

C     Transfer to BCOM
      CALL MEQU(BCOM,NMU,JDIM,ACOM)

C     Calculate t01*I0+
      CALL MMUL(1.0D0,TA,U0PL,XCOM,NMU,NMU,1,JDIM,JDIM,1)

C     Calculate r10*t21
      CALL MMUL(1.0D0,RA,TB,ACOM,NMU,NMU,NMU,JDIM,JDIM,JDIM)

C     Calculate r10*t21*I2-
      CALL MMUL(1.0D0,ACOM,UTMI,YCOM,NMU,NMU,1,JDIM,JDIM,1)

C     Add previous two results
      CALL MADD(1.0D0,XCOM,YCOM,XCOM,NMU,1,JDIM,1)

C     calculate r10*J21-
      CALL MMUL(1.0D0,RA,JB,YCOM,NMU,NMU,1,JDIM,JDIM,1)

C     Add to total
      CALL MADD(1.0D0,XCOM,YCOM,UPL,NMU,1,JDIM,1)

C     Add J01+ to total and put in UPL
      CALL MADD(1.0D0,UPL,JA,XCOM,NMU,1,JDIM,1)

C     Multiply by (E-r10*r12)^-1 for result in UPL
      CALL MMUL(1.0D0,BCOM,XCOM,UPL,NMU,NMU,1,JDIM,JDIM,1)
      RETURN
      END

