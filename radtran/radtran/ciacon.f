      SUBROUTINE CIACON(V0,DV,P,T,INORMAL,NGAS,IDGAS,
     &     ISOGAS, AMOUNT, PP, XLEN, POLY, IDUMP)
C----------------------------------------------------------------------------
C_TITLE:  CIACON: to compute gaseous continuum spectra from
C         a variety of gas pairs.
C
C_ARGS:   V0:REAL         lowest wavenumber
C         DV:REAL         wavenumber range
C         P:REAL          total pressure in atm
C         T:REAL          temp in K
C         INORMAL:INT     flag for ortho:para ratio (0=equilib. =1:1)
C                                                   (1=normal   =3:1)
C         NGAS:INT         no. of gases
C         IDGAS:INT(NGAS)  array of gas identifiers
C         ISOGAS:INT(NGAS) array of gas isotope identifiers
C         AMOUNT:REAL(NGAS)  amount of each gas (no. molecules/cm2)
C         PP:REAL(NGAS)    partial pressure of gas (atm)
C	  XLEN:REAL 	   path length in km
C         POLY(IORDP1):REAL   on exit holds poynomial fit to optical depth
C	  IDUMP: INT	  Flag for additional print statements.
C
C_KEYS:   SUBR,ATMO,SPEC,VMS
C
C_DESCR:  computes a polynomial approximation to any known continuum spectra
C         for a particular gas over a defined wavenumber region.
C
C_FILES:  none
C
C_CALLS   none
C
C_BUGS:
C
C_HIST:   26nov86 SBC ORIGINAL VERSION (HYDCON.F)
C          3feb88 SBC modified to return only one polynomial rather than
C                 for all layers in GENLBL
C         10feb97 This program adapted from HYDCON.F (H2-H2 and H2-He 
C                 continuum code) and its subroutine OD_HYDTAB.F
C                 to create CIACON.F, which calculates general gas pairs
C                 (mainly for Titan work) (C.Nixon)
C  	  26apr12 Updated for consistency between Nemesis and Radtrans
C  		  method of calculating CIA absorption
C----------------------------------------------------------------------------
        IMPLICIT NONE
        INCLUDE '../includes/arrdef.f'
	INTEGER I,J,NX, INORMAL, NGAS
        INTEGER IDGAS(MAXGAS), ISOGAS(MAXGAS),IDUMP
	REAL V0,DV,AMOUNT(MAXGAS),P,T,POLY(IORDP1)
        real XLEN,XW
        REAL X(IORDP1),Y(IORDP1),XMIN,XMAX,PP(MAXGAS)
        REAL ABSORB,DABSORB(7)
        INTEGER IABSORB(5)
   

C the number of points across the bin

        NX=IORDP1

C ******************** loop for each point across the bin **************
        DO 20 I=1,NX
         X(I)=V0 + (I-1)*DV/FLOAT(IORDER)
         if(IDUMP.EQ.1)print*,'I,X(I) : ',I,X(I)
         XW=X(I)

         CALL NCIACON(XW,P,T,INORMAL,NGAS,IDGAS,ISOGAS,AMOUNT,PP,
     1    XLEN,ABSORB,IABSORB,DABSORB,IDUMP)

         Y(I)= ABSORB

         if(IDUMP.eq.1)then
          print*, X(i), 'cm-1: total (alpha*tau) is ',y(i)
         endif

20      CONTINUE
C ************************ END WAVENUMBER LOOP *********************
  
        XMIN=V0
        XMAX=V0+DV

c fit the points with a polynomial, return co-efficients
        CALL CALC_PCOEFF(NX,Y,X,XMIN,XMAX,POLY)
        if(IDUMP.EQ.1)THEN
         print*,XMIN,XMAX
         do i=1,IORDP1
           print*,X(I),Y(I),POLY(I)
         enddo
        endif

	RETURN
	END

