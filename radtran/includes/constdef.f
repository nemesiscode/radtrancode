C     ***************************************************************
C     Useful constants for planetary RT calculations
C
C     Pat Irwin	29/2/12	Original
C
C     ***************************************************************
      REAL AVOGAD,MODBOLTZ,RGAS,AMAGAT,KBOLTZMANN,AU,PI
      PARAMETER (KBOLTZMANN=1.37947E-23)
      PARAMETER (MODBOLTZ=1.013E4/KBOLTZMANN)
      PARAMETER (AVOGAD=6.022045E23,RGAS=8310.0,AMAGAT=2.68675E19)
      PARAMETER (AU=1.49598E13,PI=3.1415927)
C     MODBOLTZ=(1/KBOLTZ)*1.013E5*1E-6*1E5. When MODBOLTZ is multiplied by the
C              pressure of a path in atm and divided by the temperature in K
C              the result is the number of molecules/cm2 per km.
C     AVOGAD is Avogadro's number
C     RGAS is the Molar gas constant in units of mJ mol-1 K-1 
C     AMAGAT is the number of molecules per unit volume at STP (1 atm. 273.15K),
C               here in units of mol cm-3 
C     AU  is an Astronomical Unit in units of cm. The units of cm are used 
C         since the solar flux has to be calculated in units of W cm-2 um-1
C         or W cm-2 (cm-1)-1
C     PI is pi!
