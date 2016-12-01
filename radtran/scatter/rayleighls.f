      real function rayleighls(v,heoverh2, ch4overh2, NH3MIX)
C     **************************************************************
C     Function to evaluate the Rayleigh scattering cross section for
C     Jovian air adapted from IDL code developed by Larry Sromovsky
C     Computes Rayleigh scattering cross section per molecule
C     considering only H2, He, CH4, and NH3 with only NH3 expressed 
C     as a volume mixing ratio
C
C     Pat Irwin   30/11/16
C
C     Input variables
C       v       	real    wavenumber (cm-1)
C       heoverh2       	real    He/H2
C       ch4overh2       real    CH4/H2
C       NH3MIX		real	NH3 vmr
C
C     Output variable
C       rayleighj       real    scattering cross section (cm2)
C
C     **************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/constdef.f'
      REAL v,heoverh2, ch4overh2, NH3MIX
      REAL h2overtot
      INTEGER ngas,J
      PARAMETER(ngas=4)
      REAL COMP(ngas),A(ngas),B(ngas),D(ngas),loschpm3,wl
      DATA A/13.58e-5, 3.48e-5, 37.0e-5, 37.0e-5/
      DATA B/7.52e-3,  2.3e-3, 12.0e-3, 12.0e-3/     
      DATA D/ 0.0221,   0.025,    .0922, .0922/
C     refr index equation coefficients from Allen, Astrophys. Quant., p 87 (1964)
C     where n-1=A(1+B/wl^2), where wl is wavelength
C     and n is the refractive index at STP (0C, 1 Atm=1.01325bar)
C     used NH3 value as a guess for CH4 which is not listed
C     depol. factors from Penndorf, J. Opt. Soc. of Amer., 47, 176-182 (1957)
C     used Parthasarathy (1951) values from Table II.
C     used CO2 value as a guess for CH4 which is not listed
      REAL sumwt,xc1,nr,fact


      h2overtot=(1.0-NH3MIX)/(1+heoverh2+ch4overh2)
      comp(1)=h2overtot
      comp(2)=heoverh2*h2overtot
      comp(3)=ch4overh2*h2overtot
      comp(4)=NH3MIX

C     loschpm3 is molecules per cubic micron at STP
      loschpm3=2.687E19*1E-12
C     wl is wavelength in microns
      wl=1e4/v


C     Compute summation over molecule-dependent scattering properties
C     Cross section formula also given in van de Hulst (1957)
C     Here 
      xc1=0.
      sumwt=0.
      DO J=1,NGAS
       nr=1.0+A(J)*(1.0+B(J)/wl**2)
       xc1=xc1+(nr**2-1.0)**2*comp(J)*(6+3*D(J))/(6-7*D(J))
       sumwt=sumwt+comp(j)
      ENDDO

      fact=8.0*(PI**3)/(3.0*(wl**4)*(loschpm3**2))

C     average cross section in cm^2 per molecule 
      rayleighls=fact*1e-8*xc1/sumwt 

      return

      end
