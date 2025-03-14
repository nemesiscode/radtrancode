      SUBROUTINE ADJUSTVMR(NPRO,NGAS,VMR,ISCALE,IERR)
C     *****************************************************************
C     Subroutine to adjust the vmrs at a particular level to add up to 1.0.
C
C     Input variables
C	NPRO			INTEGER Number of vertical layers
C	NGAS			INTEGER	Number of gases
C	VMR(MAXPRO,MAXGAS)	REAL	Gas vmrs
C	ISCALE(MAXGAS)		INTEGER	Flags to indicate if gas vmr can be 
C				scaled(1) or not (0).
C
C     Output variables
C	VMR(MAXPRO,MAXGAS)	REAL	Scaled vmrs.
C	IERR			INTEGER Error flag. Set to 0 if all OK,
C					 but set to 1 if any vmr has gone
C					 negative
C
C     Pat Irwin		Original	9/5/12
C     Pat Irwin		Added IERR	1/9/14
C
C     *****************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NGAS,NPRO,ISCALE(MAXGAS),IGAS,J,IPRO,IERR
      REAL VMR(MAXPRO,MAXGAS),SUM,SUM1,XFAC

      IERR=0
      DO 100 IPRO=1,NPRO

       SUM=0.0
       SUM1=0.0
       DO 10 IGAS=1,NGAS

        SUM=SUM+VMR(IPRO,IGAS)
        IF(ISCALE(IGAS).EQ.0)THEN
         SUM1=SUM1+VMR(IPRO,IGAS)
        ENDIF

10     CONTINUE

      
       IF(SUM.NE.1.0)THEN
C       Need to adjust the VMRs of those gases that can be scaled to 
C       bring the total sum to 1.0.

        XFAC = (1.0-SUM1)/(SUM-SUM1)

        DO 20 IGAS=1,NGAS
C        Apply scaling factor to gases that can be scaled
         IF(ISCALE(IGAS).EQ.1)THEN
          VMR(IPRO,IGAS)=VMR(IPRO,IGAS)*XFAC
         ENDIF
20      CONTINUE

       ENDIF

       DO 30 IGAS=1,NGAS
        IF(VMR(IPRO,IGAS).LT.0.0)THEN
          IERR=1
        ENDIF
30     CONTINUE

100   CONTINUE

      RETURN

      END
