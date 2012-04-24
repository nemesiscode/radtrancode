      SUBROUTINE CALCTEMPPARAM(NLAYER,NGAS,PRESS,TEMP,AAMOUNT,
     1 IDGAS,ISOGAS)
C     $Id:
C***********************************************************************
C_TITL:	CALCTEMPPARAM.f
C
C_DESC:	Calculate the line contributions to the continuum bins
C
C_ARGS:	Input variables:
C	WING		REAL	Bin width.
C	MAXDV		REAL	Line wing cut-off
C	NLAYER		INTEGER	Number of Layers
C	NGAS		INTEGER	Number of gases
C	PRESS(NLAYER)	REAL	Total pressure [atm].
C	TEMP(NLAYER)	REAL	Temperature [Kelvin].
C	IDGAS(NGAS)	INTEGER	Gas ID
C       ISOGAS(NGAS)	INTEGER Isotope ID
C	IPROC(NGAS)	INTEGER	Line wing processing parameter.
C	IB		INTEGER Buffer 1 or 2
C
C	../includes/*.f variables:
C	VLIN(2,MAXLIN)	REAL	Line position [cm^-1].
C	SLIN(2,MAXLIN)	REAL	Line strength [cm^-1 molecule^-1 cm^-2] at
C				STP.
C	ALIN(2,MAXLIN)	REAL	Air-broadened halfwidth [cm^-1/atm] @ STP.
C	ELIN(2,MAXLIN)	REAL	Lower state energy line position [cm^-1].
C	IDLIN(2,MAXLIN)	REAL	Air Force Geospace Lab. identifier.
C	SBLIN(2,MAXLIN)	REAL	Self broadening coefficient. NOTE: the
C				self-broadening coefficient used in this
C				program is the 'air'-broadened halfwidth
C				minus the self-broadened halfwidth.
C	TDW(2,MAXLIN)	REAL	Temperature coefficient of air-broadened
C				halfwidth.
C	TDWS(2,MAXLIN)	REAL	Temperature coefficient of self-broademed
C				halfwidth.
C	DOUBV(2,MAXLIN)	REAL	The inversion doublet separation (only
C				applies to longwave NH3 lines at present.   
C	LLQ(2,MAXLIN)	CHARA*9	The lower state local quanta index.
C	NLINE(2)	INTEGER	Number of lines stored in each buffer
C
C_HIST:	15apr11	PGJI	Modified from lbl_kcont.f
C
C***************************** VARIABLES *******************************

      IMPLICIT NONE


      INTEGER NGAS,IDGAS(NGAS),ISOGAS(NGAS),IPROC(NGAS),NLAYER
      REAL PRESS(NLAYER),TEMP(NLAYER)

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/dbcom.f'
C ../includes/dbcom.f stores the line database variables (e.g. MASSNO and
C RELABU).

C     General variables ...
      INTEGER I,J,K,LAYER,IPTF
      REAL XMASS,PARTF,GETMASS
       

C Continuum variables ...
      INCLUDE '../includes/contdef.f'

      REAL AAMOUNT(MAXLAY,MAXGAS)

C******************************** CODE *********************************


      DO 16 I=1,NGAS
C     Check isotopes included in model for each gas and set mass for doppler
C     width calculation

       XMASS = GETMASS(IDGAS(I),ISOGAS(I))

       IPTF=0

       DO 16 J=1,NLAYER
C      Calculate temperature tempendance parameters for the gas lines
C      (Note: TCORS1 includes factor of 1.E-47 for scaling of stored line
C       strengths. Scaling is applied in two stages to avoid numerical
C       overflow)
        IF(ISOGAS(I).EQ.0)THEN
         K=1
        ELSE
         K=ISOGAS(I)
        END IF
        TCORS1(J,I)=PARTF(IDGAS(I),K,TEMP(J),IPTF)*AAMOUNT(J,I)*1e-20
        TCORS1(J,I)=1.E-27*TCORS1(J,I)
        TCORDW(J,I)=4.301E-7*SQRT(TEMP(J)/XMASS)
16    CONTINUE

      DO 17 J=1,NLAYER
        TCORS2(J)=1.439*(TEMP(J)-296.)/(296.*TEMP(J))
17    CONTINUE

C      print*,'CALCTEMPPARAM'
C      DO 18 J=1,NLAYER
C       PRINT*,'LAYER = ',J
C       PRINT*,(TCORS1(J,I),I=1,NGAS)
C       PRINT*,(TCORDW(J,I),I=1,NGAS)
C       PRINT*,TCORS2(J)
C18    CONTINUE

      RETURN

      END
