      SUBROUTINE GENTABSCK1(KEYFIL,NPRO,NGAS,ID,ISO,P,T,VMR,NWAVE,VWAVE,
     1 VV,ISPACE,NG,DEL_G,TABK)
C     $Id: gentabsck.f,v 1.1 2005/08/01 11:11:57 irwin Exp $
C     ****************************************************************
C     Subroutine to open and read in k-distribution tables for each gas in
C     the atmosphere and interpolate their values at each point in the
C     atmospheric profile.
C
C     k-distributions are then combined assuming uncorrelated lines and
C     using the vmr's as the weight. 
C
C     Finally, gas contunuum absorption and CIA is added on.
C
C     Input variables
C	KEYFIL		character*100  	Run name
C	NPRO		integer		Number of levels in atm. profile
C	NGAS		integer		Number of gases in profile.
C	ID(MGAS)	integer		Gas Ids
C	ISO(MGAS)	integer		Gas isotope Ids
C	P(NPRO)		real		Profile pressures
C	T(NPRO)		real		Profile temperatures
C	VMR(MPRO,MGAS)	real		Gas v.m.r. profiles
C	VV		real		Calc wavenumber or wavelength
C	ISPACE		integer		1=wavelength,0=wavenumber
C
C     Output variables
C	NG		integer		Number of points in k-quadrature
C	DEL_G(MAXG)	real		k-space weights
C	TABK(MAXG,MPRO)	real		interpolated k-distributions 
C					including continuum absorptions. 
C					Units of absorption are (molecule/cm2)-1
C	
C     Pat Irwin    Original	13/5/99
C    		   Revised	1/1/05
C
C     ****************************************************************
      IMPLICIT NONE
      include '../includes/arrdef.f' 
      integer nwave

      INTEGER NPRO,NGAS,I,J,K
      REAL P(MAXPAT),T(MAXPAT),VMR(MAXPAT,MAXGAS)
      INTEGER ID(MAXGAS),ISO(MAXGAS),ISPACE
      REAL VV,X,QH,QHe,TAU1,TAU2,QH_He,PP,TOTAMH,TAUTMP      
      
      INTEGER IWAVE

      CHARACTER*100 KEYFIL,KLIST
      INTEGER NG
 
      REAL K_G1(MAXG),K_G2(MAXG),K_G(MAXG),Q1,Q2
      REAL TABK(MAXG,MAXPAT),DV
      REAL KOUT(MAXLAY,MAXGAS,MAXG)
      
      integer lun(maxbin,maxgas),ireck(maxbin,maxgas),npk,ntk,
     1 ngk
      real  xmink(maxbin,maxgas),delk(maxbin,maxgas),pk(maxk), 
     1 tk(maxk),g_ord(maxg),delg(maxg),vwave(nwave),del_g(maxg)
      real delvk, fwhmk

      common/interpk/lun, ireck, xmink, delk, pk, npk, tk, ntk,
     1   	ngk,delvk, fwhmk, g_ord, delg, kout

      NG=NGK


      DV=VWAVE(2)-VWAVE(1)
      IWAVE = 1+ INT((VV-VWAVE(1))/DV)

      CALL GET_K(NPRO,P,T,NGAS,IWAVE,VV)

      DO 1000 I=1,NPRO

       QH=0.0
       QHe=0.0
       TAU1=0.0
       TAU2=0.0
       DO 243 J=1,NGAS

          IF(ID(J).EQ.39) QH=VMR(I,J)
          IF(ID(J).EQ.40) QHe=VMR(I,J)
          IF(J.EQ.1)THEN
           DO K=1,NG
            K_G(K) = KOUT(I,J,K)
           END DO
           Q1=VMR(I,J)
          ELSE
           DO K=1,NG
            K_G1(K) = K_G(K)
            K_G2(K) = KOUT(I,J,K)
           END DO
           Q2=VMR(I,J)
           CALL OVERLAP(DELG,NG,K_G1,Q1,K_G2,Q2,K_G)
           Q1 = Q1 + Q2
          END IF
243    CONTINUE

C      Make sure we pass wavenumbers to NGASCON
       X = VV
       IF(ISPACE.EQ.1)THEN
        X=1E4/VV
       ENDIF

       TAU1=0.0
       DO K=1,NGAS
        PP = VMR(I,K)*P(I)
        CALL NGASCON(X,ID(K),ISO(K),1.0,PP,P(I),T(I),TAUTMP)
        TAU1=TAU1+TAUTMP
       ENDDO

       TAU2=0.0

       TAU1=TAU1*1E20
       TAU2=TAU2*1E20

       DO K=1,NG
        TABK(K,I)= K_G(K) + TAU1 + TAU2
       ENDDO

1000  CONTINUE

      DO K=1,NG
       DEL_G(K)=DELG(K)
      ENDDO

      RETURN

      END
