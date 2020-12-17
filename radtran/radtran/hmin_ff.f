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
       real function HMIN_FF(V0,TEMP)
       implicit none
       real  V0, TEMP

       integer  K, i

       real  hn(6),An1(6),Bn1(6),Cn1(6),Dn1(6),En1(6),Fn1(6)
       real An2(4),Bn2(4),Cn2(4),Dn2(4),En2(4),Fn2(4)

       real  wv, t_theta, summ 
C table 3a of Bell
       DATA An1 / 0.0, 2483.346, -3449.889, 2200.040, -696.271, 88.283 /
       DATA Bn1 / 0.0, 285.827, -1158.382, 2427.719, -1841.400, 
     & 444.517 /  
       DATA Cn1 / 0.0, -2054.291, 8746.523, -13651.105, 8624.970,
     &  -1863.864 /
       DATA Dn1 / 0.0, 2827.776, -11485.632, 16755.524, -10051.530,
     &  2095.288 /
       DATA En1 / 0.0, -1341.537, 5303.609, -7510.494, 4400.067,
     &  -901.788 /
       DATA Fn1 / 0.0, 208.952, -812.939, 1132.738, -655.020, 132.985 /
C table 3b of Bell     
       DATA An2 / 518.1021, 473.2636, -482.2089, 115.5291 /
       DATA Bn2 / -734.8666, 1443.4137, -737.1616, 169.6374 /
       DATA Cn2 / 1021.1775, -1977.3395, 1096.8827, -245.649 / 
       DATA Dn2 / -479.0721, 922.3575, -521.1341, 114.243 /
       DATA En2 / 93.1373, -178.9275, 101.7963, -21.9972 /
       DATA Fn2 / -6.4285, 12.3600, -7.0571, 1.5097 /

       wv = 1E4/V0 !convert to wavelength
       t_theta = 5040.0/TEMP 
C  convertion needed for calculation, see papers
C  now need to do the two wavelength conditions for the 6 different
C  coefficients
      if (wv.gt.0.3645) then
       do K=1,6
        summ = summ + t_theta**((real(K)+1.0)/2.0) * (wv**2*An1(K) +
     &  Bn1(K) + Cn1(K)/wv + Dn1(K)/wv**2 + En1(K)/wv**3 + Fn1(K)/wv**4)
       enddo
C       WRITE(*,*)'hmin_ff',summ
      elseif ((wv.lt.0.3645) .and. (wv.gt.0.1823)) then
       do K = 1,4
        summ = summ + t_theta**((real(K)+1.0)/2.0) * (wv**2*An2(K) + 
     &  Bn2(K) + Cn2(K)/wv + Dn2(K)/wv**2 + En2(K)/wv**3 + Fn2(K)/wv**4)       
       enddo
C       WRITE(*,*) 'hmin_ff', summ
      else
        summ = 0.0
      endif
C note: this k_ff has units cm^4dyne^-1, need to convert..
       HMIN_FF = 1E-29 * summ
      
       if(summ.lt.0.0)	HMIN_FF=0.0
      
       WRITE(*,*) 'hmin_ff output', summ, HMIN_FF

       end
