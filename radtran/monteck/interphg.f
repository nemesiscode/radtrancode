      SUBROUTINE INTERPHG(VV,NCONT,NLAMBDA,XLAMBDA,PHASED,XSEC,
     1 XOMEGA,XHG)
C     $Id: interphg.f,v 1.1 2005/08/01 11:11:57 irwin Exp $
C     *******************************************************************
C     Subroutine to interpolate the H-G phase functions to the required
C     wavelength/wavenumber
C
C     Input variables
C	VV	REAL	Required wavenumber/wavelength
C	NCONT	INTEGER	Number of particle types
C	NLAMBDA INTEGER Number of wavenumbers/wavelengths in scattering
C			spectra
C	XLAMBDA(NLAMBDA) REAL Wavenumbers/wavelengths of scat. spectra
C	PHASED(MAXCON,MAXSEC,5) REAL Tabulated array of scatterimg 
C                                    properties: xsec,omega,f,g1,g2
C
C     Output variables
C	XSEC(NCONT)	REAL	Particle cross-sections
C	XOMEGA(NCONT)	REAL	Particle S.S. albedos
C	XHG(NCONT,3)	REAL	Particle f,g1,g2 values
C
C     Pat Irwin		Original	27/1/00
C			Revised		1/8/05
C
C     *******************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      INTEGER NLAMBDA,I,J,K,NCONT
      REAL XLAMBDA(MAXSEC),PHASED(MAXCON,MAXSEC,5)
      REAL XSEC(MAXCON),XOMEGA(MAXCON),XHG(MAXCON,3),F,VV

      J=-1
      DO 10 I=1,NLAMBDA-1
      
       IF(VV.GE.XLAMBDA(I).AND.VV.LT.XLAMBDA(I+1))THEN
        J=I
        F=(VV-XLAMBDA(I))/(XLAMBDA(I+1)-XLAMBDA(I))
       ENDIF

10    CONTINUE

      IF(J.LT.0)THEN
       IF(VV.GE.XLAMBDA(NLAMBDA))THEN
        J=NLAMBDA-1
        F=1.0
        PRINT*,'InterpHG. VV > VMAX'
        PRINT*,VV,(XLAMBDA(I),I=1,NLAMBDA)
       ENDIF
       IF(VV.LT.XLAMBDA(1))THEN
        J=1
        F=0.0
        PRINT*,'InterpHG. VV < VMIN'
        PRINT*,VV,(XLAMBDA(I),I=1,NLAMBDA)
       ENDIF
       
      ENDIF

      DO 15 I=1,NCONT

       XSEC(I) = (1.0-F)*PHASED(I,J,1) + F*PHASED(I,J+1,1)
       XOMEGA(I) = (1.0-F)*PHASED(I,J,2) + F*PHASED(I,J+1,2)
     
       DO 20 K=1,3
        XHG(I,K)  = (1.0-F)*PHASED(I,J,K+2) + F*PHASED(I,J+1,K+2)
20     CONTINUE

15    CONTINUE

      RETURN

      END
