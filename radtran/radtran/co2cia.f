      real function co2cia(v0)
C     *************************************************************
C     Subroutine to return CIA absorption coefficients for CO2
C
C     Pat Irwin		21/10/05
C     *************************************************************
      real v0,xl,aco2

C     Compute wavelength for CO2-CO2 CIA approximation
      xl = 1e4/v0
      aco2 =0.0

C     2.3 micron window. Assume de Bergh 1995 a = 4e-8 cm-1/amagat^2
      if(xl.ge.2.15.and.xl.le.2.55)aco2 = 4e-8

C     1.73 micron window. Assume mean a = 6e-9 cm-1/amagat^2
      if(xl.ge.1.7.and.xl.le.1.76)aco2 = 6e-9

C     1.18 micron window. Assume a mean a = 1.5e-9 cm-1/amagat^2
      if(xl.ge.1.05.and.xl.le.1.35)aco2 = 1.5e-9

      co2cia = aco2

      return

      end
