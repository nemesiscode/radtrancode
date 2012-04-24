      real function tau_mg_lor(mg,knu0,yv,EL,SFB,qrot,P,T,U,q)
C     $Id: tau_mg_lor.f,v 1.2 2011-06-17 15:40:27 irwin Exp $
      implicit none
      real knu,U,T,EL,knu0,C1,T0,yv,SFB,P,q,qrot,pi
      double precision A,B,denom_goodylor,denom_malklor
      parameter (T0=296.0, pi=3.1415927, C1=1.439)
      logical mg


      if(knu0.eq.0.or.U.eq.0.)then
       tau_mg_lor=0.
      else


      knu=(knu0*(T0/T)**qrot)*exp(C1*EL*(1/T0 - 1/T))

C      print*,'knu0,yv,EL,SFB,qrot'
C      print*,knu0,yv,EL,SFB,qrot
C      print*,'P,T,U,q'
C      print*,P,T,U,q

C      print*,'mg = ',mg

      A=dble(knu*U)
      B=dble(pi*A*yv*P*(q+(1-q)/SFB)*sqrt(T0/T))
     
      if(mg)then
       denom_malklor = 1/(2*A) + sqrt( 1/((2*A)**2) +1/(4*B) )
       tau_mg_lor = sngl(1./denom_malklor)
      else
       denom_goodylor = sqrt( 1/(A**2) + 1/B ) 
       tau_mg_lor = sngl(1./denom_goodylor)
      end if

C      print*,'knu,A,B,tau = ',knu,A,B,tau_mg_lor


      end if

      return
      end
