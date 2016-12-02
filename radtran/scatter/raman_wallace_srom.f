      subroutine raman_wallace_srom(v,T,inormal,fH2, xsec, omega)
C     ************************************************************
C     Subroutine to estimate Raman-scattering x-section and single-
C     scattering albedo using modified version of Wallace 
C     approximation recommended by Sromovsky.
C
C     Taken from Sromovsky (2005) Icarus 173, 254-283, Eq. 3
C
C     Input variables
C       v       real    Wavenumber (cm-1)
C       T       real    Temperature (K)
C       inormal integer 0 = equilibrium ortho/para, 1= 3:1 ortho/para
C	fH2	real	H2 vmr
C
C     Output variables
C	xsec	real	Extinction x-section (cm2)
C	omega	real	Single scattering albedo
C
C     Pat Irwin	2/12/16	Original version.
C
C     ************************************************************
      real v,T,fH2,xsec,xsca,xrot,xvib,beta
      parameter (beta = 0.433)

      call calc_raman(v,T,inormal,xrot,xvib)


      xsec = xrot+xvib
      xsca = xrot+beta*xvib

      omega=xsca/xsec
      xsec=xsec*fH2

      return

      end
