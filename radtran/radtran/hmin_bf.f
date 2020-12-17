C***********************************************************************
C Program for H- opacity
C Purpose: For temperatures > 2000 K , H- begins to become a source of
C opacity in the near-IR (mainly WFC3 observing band)
C Version: 1.0 J Taylor 27/04/18
C Descriptions: Information was taken from Bell & Berrington (1987)
C and John (1988), further work was conducted by Lenzuni (1991) which
C provided analytical formulae for H- opacity (but in terms for density).
C (Thanks to Vivien Parmentier for giving my guidance on this)
C The most information was taken from Bell 1988, it contains formulae
C and the coefficients needed to calculate the absorption.
C Eq 6 gives formula needed for free free, Eq 3/4/5 gives for bound free
C overall H- opacity is ktot = kff + kbf
C***********************************************************************
       real function HMIN_BF(V0)
       implicit none
       integer i
       real wv, V0, lambda0, coeff(6) , sig, summ
C Now begin calculation for bound-free
       DATA coeff/152.519,49.534,-118.858,92.536,-34.194,4.982/
       DATA lambda0/1.6419/
   
       wv = 1E4/V0
       if (wv.lt.lambda0) then
C	 WRITE(*,*)'hmin wv',wv
         do i=1,6
 	        summ= summ+coeff(i)*((1.0/wv)-(1.0/lambda0)) 
     &		**((real(i)-1.0)/2.0)
         enddo
    	 WRITE(*,*), 'hmin_bf',summ
       else
      		 summ = 0.0
       endif
         sig = 1E-18*(wv**3.0)*((1.0/wv)-(1.0/lambda0))**(3.0/2.0)
C       WRITE(*,*) 'This is hmin_bf',sig
       
       if (wv.lt.lambda0) then
	HMIN_BF = sig*summ
       else
	HMIN_BF =0.0
       endif
       end

