      REAL FUNCTION CALC_NLTE(P,INLTE)
C     $Id: calc_nlte.f,v 1.2 2011-06-17 15:40:25 irwin Exp $
C     *****************************************************************
C     Function to calculate the Non Local Thermodynamic Equilibrium
C     correction to the Planck function.
C
C     Hack at will.
C
C     Pat Irwin		16/8/94
C     *****************************************************************
      REAL P
      INTEGER INLTE

      IF(INLTE.EQ.1)THEN
C     Assume NLTE is as that estimated by Glenn Orton for Saturn 
C     and Titan NLTE_Sat
      
       CALC_NLTE = 1.E6*P/(1.+1.E6*P) 

      ELSE
C     Assume no correction

       CALC_NLTE = 1.
 
      END IF

      RETURN
      END
