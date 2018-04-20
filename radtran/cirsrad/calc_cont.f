      SUBROUTINE CALC_CONT(WING,MAXDV,NLAYER,NGAS,PRESS,TEMP,
     1 FRAC,IDGAS,ISOGAS,IPROC,IB,IX)
C     $Id:
C***********************************************************************
C_TITL:	CALC_CONT.f
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
C	FRAC(NLAYER,NGAS) REAL Fractional abundance of each gas
C	IDGAS(NGAS)	INTEGER	Gas ID
C       ISOGAS(NGAS)	INTEGER Isotope ID
C	IPROC(NGAS)	INTEGER	Line wing processing parameter.
C	IB		INTEGER Buffer 1 or 2
C	IX		INTEGER Call number to routine
C
C	../includes/*.f variables:
C	VLIN(2,MAXLIN)	REAL*8	Line position [cm^-1].
C	SLIN(2,MAXLIN)	REAL*8	Line strength [cm^-1 molecule^-1 cm^-2] at
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
C	29feb12	PGJI	Updated for Radtrans2.0
C
C***************************** VARIABLES *******************************

      IMPLICIT NONE


      INTEGER NGAS,IDGAS(NGAS),ISOGAS(NGAS),IPROC(NGAS),NLAYER
      REAL VMIN,VMAX,WING,VREL,MAXDV

      INCLUDE '../includes/arrdef.f'
C     Continuum and temperature parameter variables ...
      INCLUDE '../includes/contdef.f'
      REAL TRATIO,TSTIM,VOUT(MAXBIN),YOUT(MAXBIN)
      REAL FRAC(MAXLAY,MAXGAS),PRESS(NLAYER),TEMP(NLAYER)

      INCLUDE '../includes/lincomc.f'
C ../includes/lincom.f stores the linedata variables (including MAXLIN,
C VLIN, SLIN, ALIN, ELIN, SBLIN, TDW, TDWS and that lot).


C     General variables ...
      REAL DV,LINECONTRIB,X,FNH3,FH2
      INTEGER I,J,K,L,LINE,IBIN,CURBIN,IGAS,IB
      INTEGER LAYER,IX
      DOUBLE PRECISION VV,VJ,STR
      REAL LSE,XWIDA,XWIDS,XTDEP

      REAL DPEXP,ABSCO,Y,CONVAL
       
C     GASCON switches
      INCLUDE '../includes/gascom.f'
      INCLUDE '../includes/lcocom.f'

C******************************** CODE *********************************


      print*,'CALC_CONT. IX,IB,NLAYER = ',IX,IB, NLAYER      
      DO 13 LINE=1,NLINE(IB)
       CURBIN = 1 + INT((VLIN(IB,LINE)-VBIN(1))/WING)
       IGAS=IDLIN(IB,LINE)
       DO 15 J=1,NBIN
	
C       Computing continuum for all except adjacent bins
        IF(ABS(CURBIN-J).LE.1)GOTO 15

        DO 21 K=1,IORDP1
          VV = DBLE(CONWAV(K) + VBIN(J)) 
          DV = SNGL(VV - VLIN(IB,LINE))
C         Don't calculate at wavenumbers more than MAXDV away
          IF(ABS(DV).LE.MAXDV)THEN

           DO 101 LAYER=1,NLAYER


             IF(IDGAS(IGAS).EQ.1.AND.IH2O.GT.0.AND.
     &         ABS(DV).GT.25.0)THEN
C              Don't calc continuum more than 25cm-1 from H2O lines
C              if H2O continuum is turned on
               CONVAL = 0.0

             ELSE
              FH2=-1.
              FNH3=-1.
              CONVAL = LINECONTRIB(IPROC(IGAS),IDGAS(IGAS),VV,
     1       TCORDW(LAYER,IGAS),TCORS1(LAYER,IGAS),TCORS2(LAYER),
     2       PRESS(LAYER),TEMP(LAYER),FRAC(LAYER,IGAS),VLIN(IB,LINE),
     3       SLIN(IB,LINE),ELIN(IB,LINE),ALIN(IB,LINE),SBLIN(IB,LINE),
     4       TDW(IB,LINE),TDWS(IB,LINE),LLQ(IB,LINE),DOUBV(IB,LINE),
     5       FNH3,FH2)
C              IF(LINE.EQ.1)THEN
C               PRINT*,IPROC(IGAS),IDGAS(IGAS),VV
C               PRINT*,TCORDW(LAYER,IGAS),TCORS1(LAYER,IGAS),
C     1  TCORS2(LAYER)
C               PRINT*,PRESS(LAYER),TEMP(LAYER),FRAC(LAYER,IGAS),
C     1  VLIN(IB,LINE)
C               PRINT*,SLIN(IB,LINE),ELIN(IB,LINE),ALIN(IB,LINE),
C     1  SBLIN(IB,LINE)
C               PRINT*,TDW(IB,LINE),TDWS(IB,LINE),LLQ(IB,LINE),
C     1  DOUBV(IB,LINE)
C               PRINT*,FNH3,FH2
C              ENDIF
             ENDIF

            CONVALS(K,LAYER,J)=CONVALS(K,LAYER,J)+CONVAL

101        CONTINUE

          ENDIF

21       CONTINUE

15     CONTINUE

13    CONTINUE


C     Extra continuum for linedata files that are so large that they
C     have been stripped of weaker lines.
      IF(IJLCO.GT.0.AND.IX.EQ.1)THEN
       DO 299 IGAS=1,NGAS
        IF(IDGAS(IGAS).EQ.IDLCO.AND.ISOGAS(IGAS).EQ.ISOLCO)THEN
         print*,'Adding LCO for gas : ',IGAS,IDLCO,ISOLCO
         FNH3=-1.0
         FH2=-1.0
         DO 300 LAYER=1,NLAYER
          CALL CALC_LCO(VBIN,NBIN,WING,MAXDV,IPROC(IGAS),
     1     IDGAS(IGAS),TCORDW(LAYER,IGAS),TCORS1(LAYER,IGAS),
     2     TCORS2(LAYER),PRESS(LAYER),TEMP(LAYER),FRAC(LAYER,IGAS),
     3     FNH3,FH2,VOUT,YOUT)
          DO 301 I=1,NBIN
           DO 302 K=1,IORDP1
            VV=VBIN(I)+CONWAV(K)
            CALL VERINT(VOUT,YOUT,NBIN+1,CONVAL,SNGL(VV))
            print*,K,LAYER,I,CONVALS(K,LAYER,I),CONVAL
            CONVALS(K,LAYER,I)=CONVALS(K,LAYER,I)+CONVAL
302        CONTINUE
301       CONTINUE
300      CONTINUE
        ENDIF
299    CONTINUE
      ENDIF

C     Convert continuum values to polynomial coefficients
      DO 200 LAYER=1,NLAYER
       DO 202 J=1,NBIN
        DO 205 K=1,IORDP1
         CONTINK(K,LAYER,J)=0.0
         DO 210 L=1,IORDP1
          CONTINK(K,LAYER,J) = CONTINK(K,LAYER,J) + 
     &         UNIT(L,K)*CONVALS(L,LAYER,J)
210      CONTINUE
205     CONTINUE
202    CONTINUE
200   CONTINUE



      RETURN

      END
************************************************************************
************************************************************************
