      SUBROUTINE WRITE_K(LUN1,NG,K_G,IREC)
C     ****************************************************************
C     Subroutine to write a single k-distribution to an output .kta 
C     lookup file
C
C     Input variables
C	LUN1		INTEGER	Logical unit number
C	NG		INTEGER	Number of k-coefficients in distribution
C	K_G(MAXG)		REAL	k-distribution
C	IREC		INTEGER index of first output record
C
C     Output variables
C    	IREC		INTEGER Index of next record
C
C     Pat Irwin		23/2/96
C	
C     ****************************************************************
      INCLUDE '../includes/arrdef.f'
      REAL K_G(MAXG)
      INTEGER LUN1,NG,IREC

      DO 40 LOOP=1,NG
          WRITE(LUN1,REC=IREC)K_G(LOOP)
          IREC=IREC+1
40    CONTINUE

      RETURN

      END

