      SUBROUTINE CALC_FINE(VV,WING,MAXDV,LAYER,NLAYER,NGAS,PRESS,TEMP,
     1 FRAC,IDGAS,ISOGAS,IPROC,IBS,IFCONT,XK)
C     $Id:
C***********************************************************************
C_TITL:	CALC_FINE.f
C
C_DESC:	Calculate the line contributions to the continuum bins
C
C_ARGS:	Input variables:
C	VV		REAL 	Calculation wavenumber
C	WING		REAL	Bin width.
C	MAXDV		REAL	Line wing cut-off
C	LAYER		INTEGER	Layer number
C	NLAYER		INTEGER Number of layers
C	NGAS		INTEGER	Number of gases
C	PRESS(NLAYER)	REAL	Total pressure [atm].
C	TEMP(NLAYER)	REAL	Temperature [Kelvin].
C	FRAC(NLAYER,NGAS) REAL Fractional abundance of each gas
C	IDGAS(NGAS)	INTEGER	Gas ID
C       ISOGAS(NGAS)	INTEGER Isotope ID
C	IPROC(NGAS)	INTEGER	Line wing processing parameter.
C	IBS(2)		INTEGER Buffers 1 and 2
C
C	../includes/*.f variables:
C	VLIN(2,MAXLIN)	REAL	Line position [cm^-1].
C	SLIN(2,MAXLIN)	REAL	Line strength [cm^-1 molecule^-1 cm^-2] at
C				STP.
C	ALIN(2,MAXLIN)	REAL	Air-broadened halfwidth [cm^-1/atm] @ STP.
C	ELIN(2,MAXLIN)	REAL	Lower state energy line position [cm^-1].
C	IDLIN(2,MAXLIN)	INTEGER	Air Force Geospace Lab. identifier.
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
C	LLQ(2,MAXLIN)	CHARA*15The lower state local quanta index.
C	NLINE(2)	INTEGER	Number of lines stored in each buffer
C
C_HIST:	15apr11	PGJI	Modified from lbl_kcont.f
C_HIST:	29feb12	PGJI	Updated for Radtrans2.0
C
C***************************** VARIABLES *******************************

      IMPLICIT NONE


C     Continuum and temperature parameter variables ...
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/contdef.f'
      REAL TRATIO,TSTIM
      LOGICAL NANTEST,ISNAN

      INTEGER NGAS,IDGAS(NGAS),ISOGAS(NGAS),IPROC(NGAS)
      INTEGER LAYER,NLAYER
      REAL VMIN,VMAX,WING,VREL,MAXDV,XK
      REAL FRAC(MAXLAY,MAXGAS),PRESS(NLAYER),TEMP(NLAYER)


      INCLUDE '../includes/lincomc.f'
C ../includes/lincom.f stores the linedata variables (including MAXLIN,
C VLIN, SLIN, ALIN, ELIN, SBLIN, TDW, TDWS and that lot).


C     General variables ...
      REAL DV,LINECONTRIB,VV,X,FNH3,FH2
      INTEGER I,J,K,L,LINE,IBIN,CURBIN,IGAS,IB,IBS(2),IBX
      INTEGER IFCONT,F1,L1,B1
      REAL CONVAL,VTMP
       



C******************************** CODE *********************************

      NANTEST=.FALSE.

      CURBIN = 1+INT((VV-VBIN(1))/WING)
      IFCONT = 0
      XK=0.0


      DO 12 IBIN=CURBIN-1,CURBIN+1

       DO 200 IBX=1,2

        IB=IBS(IBX)

        F1=FSTLIN(IB,IBIN)
        L1=LSTLIN(IB,IBIN)
        B1=LASTBIN(IB)
        IF(IBIN.EQ.CURBIN-1.AND.IBX.EQ.1.AND.B1.LT.IBIN) THEN
C        We have now run past the end of buffer 1 and need to read in a 
C        new load of lines
C        Pass back flag  
         print*,IBIN,B1,F1,L1
         IFCONT = 1
        ENDIF

        IF(F1.GT.0) THEN
         DO 13 LINE=F1,L1

          IGAS=IDLIN(IB,LINE)

          DV = (VV - VLIN(IB,LINE))

C         Ignore lines more than MAXDV widths away
          IF(ABS(DV).LE.MAXDV)THEN

            FNH3=-1.0
            FH2=-1.0 
            CONVAL = LINECONTRIB(IPROC(IGAS),IDGAS(IGAS),VV,
     1       TCORDW(LAYER,IGAS),TCORS1(LAYER,IGAS),TCORS2(LAYER),
     2       PRESS(LAYER),TEMP(LAYER),FRAC(LAYER,IGAS),VLIN(IB,LINE),
     3       SLIN(IB,LINE),ELIN(IB,LINE),ALIN(IB,LINE),SBLIN(IB,LINE),
     4       TDW(IB,LINE),TDWS(IB,LINE),LLQ(IB,LINE),DOUBV(IB,LINE),
     5       FNH3,FH2)

             NANTEST=ISNAN(CONVAL)

             IF(NANTEST)THEN
              print*,'Calc_fine.f: Code reports that calculated'
              print*,'absorption is Not A Number'
              print*,'CONVAL = ',CONVAL
              print*,'LAYER,IGAS = ',LAYER,IGAS
              print*,'IPROC,IDGAS,VV = ',IPROC(IGAS),IDGAS(IGAS),VV
              print*,'TCORDW = ',TCORDW(LAYER,IGAS)
              print*,'TCORS1 = ',TCORS1(LAYER,IGAS)
              print*,'TCORS2 = ',TCORS2(LAYER)
              print*,'PRESS, TEMP, FRAC = ',PRESS(LAYER),TEMP(LAYER),
     1         FRAC(LAYER,IGAS)
              print*,'VLIN, SLIN, ELIN, ALIN, SBLIN = ',VLIN(IB,LINE),
     3         SLIN(IB,LINE),ELIN(IB,LINE),ALIN(IB,LINE),SBLIN(IB,LINE)
              print*,'TDW, TDWS, LLQ, DOUBV = ',
     4       TDW(IB,LINE),TDWS(IB,LINE),LLQ(IB,LINE),DOUBV(IB,LINE)
              print*,'FNH3, FH2 = ',FNH3,FH2
               
              STOP
             ENDIF

            XK=XK+CONVAL

          ENDIF

13       CONTINUE

        ENDIF

200    CONTINUE

12    CONTINUE

C     Now add on the line wing continuum contribution
      DV = VV-VBIN(CURBIN)
      VTMP=1.0

      DO I=1,IORDP1
       XK=XK+CONTINK(I,LAYER,CURBIN)*VTMP
       VTMP=VTMP*DV
      ENDDO

      NANTEST=ISNAN(XK)
      IF(NANTEST)THEN
       PRINT*,'CALC_FINE: XK is NAN'       
      ENDIF
      RETURN

      END
************************************************************************
************************************************************************
