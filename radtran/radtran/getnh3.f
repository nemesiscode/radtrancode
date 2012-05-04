      REAL FUNCTION GETNH3(P,T,FNH3,FH2,LLQ,TR)
C     *******************************************************************
C     Function to compute the Rosenkrantz-Ben-Reuven lineshape correction
C     parameter for the NH3 inversion  doublets.
C
C     Input parameters:
C	P	REAL	Pressure (atm)
C	T	REAL	Temperature (K)
C	FNH3	REAL	NH3 vmr
C	FH2	REAL	H2 vmr
C	LLQ	CHARACTER*15	Lower state local quanta
C	TR	REAL	Distance to other doublet member (cm-1).
c	output
c	GETNH3	real	y_k in the birnbaum et al 2000
C 
C     Pat Irwin		22/9/00
c
c	n teanby	11/03/04	re-organised so clearer
C
C     *******************************************************************
      IMPLICIT NONE

      REAL P,T,FNH3,FH2,BRNH3,BRH2,TR,deltav_jk
      CHARACTER*15 LLQ
      INTEGER J,K
      PARAMETER (BRNH3 = -2.5E-2, BRH2 = -4.5E-3)
 
      READ(LLQ(1:2),*)J
      READ(LLQ(3:4),*)K

      if (tr.eq.0.) then
         GETNH3 = 0.0
      else
c	** eqn 13 from birnbaum et al 2000 **
         deltav_jk = (Brnh3*Fnh3 + Brh2*Fh2)*(K**2)*(2*J+1)/((J+1)**2)
c	** apply temperature dependence **
c	** NB because P is in atmospheres P=P/P0 **
	   deltav_jk = deltav_jk * P * (296.0/T)**0.5
c	** divide by Tr to get y_k from birnbaum et al 2000 **
	   GETNH3 = deltav_jk/Tr
      endif

      RETURN

      END
