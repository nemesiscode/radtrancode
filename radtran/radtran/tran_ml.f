      real function tran_ml(SS,BB,U)
C     $Id: tran_ml.f,v 1.2 2011-06-17 15:40:27 irwin Exp $
C     **********************************************************************
C
C     Calculates the transmission of a Malkmus-Lorentz band using Lacis and
C     Oinas conventions
C
C     Pat Irwin		26/10/94
C
C     **********************************************************************
      implicit none
      double precision SS,BB
      real U,pi,arg
      parameter (pi=3.1415927)

      arg = sngl(sqrt(1 + 4*SS*U/(pi*BB)))
      arg = sngl(pi*BB*0.5*(arg - 1.))

      tran_ml = exp(-arg)

      return
 
      end
