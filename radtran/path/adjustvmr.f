      SUBROUTINE ADJUSTVMR(NGAS,VMR,ISCALE)
C     *****************************************************************
C     Subroutine to adjust the vmrs at a particular level to add up to 1.0.
C
C     Input variables
C	NGAS		INTEGER	Number of gases
C	VMR(NGAS)	REAL	Gas vmrs
C	ISCALE(NGAS)	INTEGER	Flags to indicate if gas vmr can be scaled(1)
C				or not (0).
C
C     Output variables
C	VMR(NGAS)	REAL	Scaled vmrs.
C
C     Pat Irwin		Original	9/5/12
C
C     *****************************************************************
      IMPLICIT NONE
      INTEGER NGAS,ISCALE(NGAS),MAXGAS,IGAS,J
      REAL VMR(NGAS),SUM,SUM1
      PARAMETER(MAXGAS=20)
      REAL RATIO(MAXGAS),XFAC

      SUM=0.0
      SUM1=0.0
      DO 10 IGAS=1,NGAS

       RATIO(IGAS)=VMR(IGAS)/VMR(1)

       IF((IGAS.GT.1).AND.(VMR(IGAS).GT.VMR(1)))THEN
        PRINT*,'Warning from ADJUSTVMR.F - first gas does not have'
        PRINT*,'highest abundance'
        PRINT*,(VMR(J),J=1,NGAS)
       ENDIF

       SUM=SUM+VMR(IGAS)
       IF(ISCALE(IGAS).EQ.0)THEN
        SUM1=SUM1+VMR(IGAS)
       ENDIF

10    CONTINUE
      
      IF(SUM.NE.1.0)THEN
C      Need to adjust the VMRs of those gases that can be scaled to 
C      bring the total sum to 1.0. We assume here that the scaleable gases have
C      fixed ratios with respect to the first gas, which is assumed to have the
C      highest abundance.

       XFAC = (1.0-SUM1)/(SUM-SUM1)

       DO 20 IGAS=1,NGAS
C       Apply scaling factor to gases that can be scaled
        IF(ISCALE(IGAS).EQ.1)THEN
         VMR(IGAS)=VMR(IGAS)*XFAC
        ENDIF
20     CONTINUE

      ENDIF

      RETURN

      END
