      real function vvweiss(gamma,nu,nu0)
C     ****************************************************************
C     Function to calculate the VanVleck-Weisskopf lineshape which is
C     more accurate as nu --> 0.0. 
C
C     Doppler broadening is ignored since the places where VVW is most
C     important are toward the far, far ir and submillimeter through
C     microwave, where (for the outer planets, at least) one
C     tends to be working with opacities that are significant only toward
C     the troposphere, with minimal concern about Doppler effects.
C
C     Developed for CIRS data analysis
C
C     Current implementation taken from Microwave Remote Sensing, by
C     Ulaby, Moore and Fung. BR13D
C
C     Pat Irwin  	18/9/00 
C     ****************************************************************
      implicit none
      real pi,gamma,nu,nu0,dnu1,dnu2,sum
      parameter (pi=3.1415927)

      dnu1 = (nu-nu0)**2
      dnu2 = (nu+nu0)**2

      sum = gamma/(dnu1 + gamma**2) + gamma/(dnu2 + gamma**2)

      vvweiss = sum*(nu/nu0)/pi

      return

      end
