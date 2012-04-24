	SUBROUTINE GASCON(V0,DV,ID,ISO,IORDP1,AMOUNT,PPRESS,
     &  PRESS,TEMP,POLY)
C     $Id: gascon.f,v 1.2 2010-10-29 12:38:43 irwin Exp $
C----------------------------------------------------------------------------
C_TITLE:  GASCON: to compute gaseous continuum spectra
C      
C     Input variables                             
C         V0:REAL         lowest wavenumber
C         DV:REAL         wavenumber range
C         ID:INTEGER      local gas identifier
C         ISO:INTEGER     local isotope identifier
C         ORDER1:INTEGER  order of polynomial required +1
C         AMOUNT:REAL     absorber amount
C         PPRESS:REAL     partial pressure of gas in atm
C         PRESS:REAL      total pressure in atm
C         TEMP:REAL       temp in K
C
C     Output variable
C         POLY(IORDP1):REAL  on exit holds polynomial fit to optical depth
C
C_DESCR:  computes a polynomial approximation to any known continuum spectra
C         for a particular gas over a defined wavenumber region.
C
C_HIST:   26nov86 SBC ORIGINAL VERSION
C          3feb88 SBC modified to return only one polynomial rather than
C                 for all layers in GENLBL
C	  27mar96 ILH modified to return a H20 continuum (based on GENLN2 
C		  sub-routine CONTUM.F)
C          3feb12 Stripped down to just call NCIACON and fit polynomial 
C----------------------------------------------------------------------------
C
	INTEGER ID,ISO,IORDP1,I
	REAL V0,DV,AMOUNT,PPRESS,PRESS,TEMP,POLY(IORDP1)
        REAL XX(3),YY(3),ABSORB



	DO 10 IP=1,IORDP1

	 FF = V0 + FLOAT(IP-1)*0.5*DV

         CALL NGASCON(FF,ID,ISO,AMOUNT,PPRESS,
     &  PRESS,TEMP,ABSORB)

         POLY(IP)=ABSORB

   10	CONTINUE

C       Convert polynomial calculations to a polynomial fit.

        DO I=1,IORDP1
         XX(I)=0.5*DV*(I-1)
         YY(I)=POLY(I)
        END DO
        XMIN=XX(1)
        XMAX=XX(IORDP1)
        NX=IORDP1

        CALL CALC_PCOEFF(NX,YY,XX,XMIN,XMAX,IORDP1,POLY)

C        print*,'Fitted quadratic : ',(POLY(I),I=1,IORDP1)

C        DO I=1,IORDP1
C         XX(I)=0.5*DV*(I-1)
C         YY(I)=POLY(1) + POLY(2)*XX(I) + POLY(3)*XX(I)*XX(I)
C         print*,xx(i),yy(i)
C        END DO       
C        print*,'Reconstruct : ',(YY(I),I=1,IORDP1)

	RETURN
	END

