      real function tau_mg_lor(mg,knu0,yv,EL,SFB,qrot,P,T,U,q)
C     $Id:
C     *************************************************************
C     Calculates the optical depth of a path using malkmus-lorentz or 
C     goody-lorentz band approximations.
C
C     Pat Irwin	?/?/??	Original version
C     Pat Irwin 26/4/12	Commented
C
C     *************************************************************

      implicit none
      real knu,U,T,EL,knu0,C1,T0,yv,SFB,P,q,qrot,pi
      double precision A,B,denom_goodylor,denom_malklor
      parameter (T0=296.0, pi=3.1415927, C1=1.439)
      logical mg


      if(knu0.eq.0.or.U.eq.0.)then
       tau_mg_lor=0.
      else


      knu=(knu0*(T0/T)**qrot)*exp(C1*EL*(1/T0 - 1/T))


      A=dble(knu*U)
      B=dble(pi*A*yv*P*(q+(1-q)/SFB)*sqrt(T0/T))
     
      if(mg)then
       denom_malklor = 1/(2*A) + sqrt( 1/((2*A)**2) +1/(4*B) )
       tau_mg_lor = sngl(1./denom_malklor)
      else
       denom_goodylor = sqrt( 1/(A**2) + 1/B ) 
       tau_mg_lor = sngl(1./denom_goodylor)
      end if


      end if

      return
      end
