      SUBROUTINE LBL_KCONT(VMIN,VMAX,WING,VREL,PRESS,TEMP,IDGAS,
     1 ISOGAS,FRAC,IPROC,IP,IT,MAXDV,IPTF,IEXO)
C     $Id:
C***********************************************************************
C_TITL:	LBL_KCONT.f
C
C_DESC:	General purpose absorption coefficient and cumulative 
C	K-distribution calculation routine. Subroutine to calculate the
C	line continuum in the range VMIN-VMAX in bins of width WING.
C
C_ARGS:	Input variables:
C	VMIN		REAL	Wavenumber of start of complete range.
C	VMAX		REAL	Wavenumber at end of range.
C	WING		REAL	Bin width.
C	VREL		REAL	Extra wavenumber range to consider
C	PRESS		REAL	Total pressure [atm].
C	TEMP		REAL	Temperature [Kelvin].
C	IDGAS		INTEGER	The LOCAL gas identifier for each gas.
C				The local gas identifier agrees with
C				HITRAN id's as far as possible (Jan
C				1988) with extra id's for gases not in the
C				HITRAN compilation. eg. those in the GEISA
C				compilation.
C        ISOGAS		INTEGER The LOCAL isotopic identifier for each
C				gas (if zero, all the isotopes of the gas
C				are included). Isotope ID's also agree as
C				far as possible with HITRAN ID's. Similar
C				(1,2,3...n) ID's have been defined for
C				additional gases. If zero, then line
C				strengths are used as tabulated
C				(i.e. corrected for normal terrestrial
C				distribution). If a specific isotope is
C				selected then its strength is scaled to
C				the pure isotope using the HITRAN
C				linedatabase values.
C	 FRAC		REAL	Broadening fraction (0 = air-broadened,
C				1 = self-broadened).
C	 IPROC		INTEGER	Line wing processing parameter.
C	 IP		INTEGER	Pressure ordinate in calculated CONTINK
C				continuum array.
C	 IT		INTEGER	Temperature ordinate in calculated
C				CONTINK continuum array.
C	 MAXDV		REAL	Line wing cut-off parameter: The maximum  
C				wavenumber-distance away within which to
C				consider the contribution of the line 
C				wings.
C	 IPTF		INTEGER	Partition function flag for CH4.
C	 IEXO		INTEGER IEXO=0, same line data for all temperatures
C				IEXO=1, different line data for different
C				  temperatures.
C	 IEXTRA		INTEGER	IEXTRA=1 
C
C	../includes/*.f variables:
C	VLIN(MAXLIN)	REAL*8	Line position [cm^-1].
C	SLIN(MAXLIN)	REAL*8	Line strength [cm^-1 molecule^-1 cm^-2] at
C				STP.
C	ALIN(MAXLIN)	REAL	Air-broadened halfwidth [cm^-1/atm] @ STP.
C	ELIN(MAXLIN)	REAL	Lower state energy line position [cm^-1].
C	IDLIN(MAXLIN)	REAL	Air Force Geospace Lab. identifier.
C	SBLIN(MAXLIN)	REAL	Self broadening coefficient. NOTE: the
C				self-broadening coefficient used in this
C				program is the 'air'-broadened halfwidth
C				minus the self-broadened halfwidth.
C	TDW(MAXLIN)	REAL	Temperature coefficient of air-broadened
C				halfwidth.
C	TDWS(MAXLIN)	REAL	Temperature coefficient of self-broademed
C				halfwidth.
C	DOUBV(MAXLIN)	REAL	The inversion doublet separation (only
C				applies to longwave NH3 lines at present.   
C	LLQ(MAXLIN)	CHARA*9	The lower state local quanta index.
C
C_FILE:	unit=72		contink.out
C
C_CALL:	loadbins	Loads linedata "bins" for GENLBL.
C	matinv		Inverts an input matrix.
C	dpexp		Performs exponentiation in double precision.
C	partf		Uses four-term polynomial fits to compute ratio
C			of total partition functions for use by line by
C			line programs.
C
C_HIST:	17feb93	PGJI	Modified from genlbl.f
C	02mar93	PGJI	Reduced and simplified to be of more general use
C	10may94	PGJI	Revised and further debugged
C	1Sept95	ALW	Revised and further debugged
C	3Nov95	ALW	Introduced fracb variable
C	30Mar00	PGJI	Modified from lbl_kdists to compute the continuum
C			contribution.
C	8/3/04	NT added iproc=6 option, Rosenkrantz-Ben-Reuven lineshape
C			with voigt lineshape as default (NH3 only)
C***************************** VARIABLES *******************************

      IMPLICIT NONE

C Definition of input parameters ...
      INTEGER IP,IT,IDGAS,ISOGAS,IPROC,IREAD,IEXO
      REAL VMIN,VMAX,WING,VREL,MAXDV
      REAL FRAC,PRESS,TEMP
      DOUBLE PRECISION VV,VJ,STR
      REAL LSE,XWIDA,XWIDS,XTDEP

C The include files ...
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
C ../includes/bincom.f stores the line bin variables (e.g. NBIN, FSTBIN,
C LSTBIN, NLINES, FSTLIN, LSTLIN) and band parameters.
      INCLUDE '../includes/dbcom.f'
C ../includes/dbcom.f stores the line database variables (e.g. MASSNO and
C RELABU).
      INCLUDE '../includes/lincom.f'
C ../includes/lincom.f stores the linedata variables (including MAXLIN,
C VLIN, SLIN, ALIN, ELIN, SBLIN, TDW, TDWS and that lot).
      INCLUDE '../includes/parcom.f'
      INCLUDE '../includes/lcocom.f'
C ../includes/parcom.f stores the parameter values such as MAXLAY,


C General variables ...
      INTEGER NGAS
      REAL DV,VBOT,VTOP,RANGE
C NGAS: Number of gases.
C DV: Distance from central wavenumber [cm^-1].
C VBOT: =VMIN - VREL, bottom of calculation range [cm^-1].
C VTOP: =VMAX + VREL, top of calculation range [cm^-1].
C RANGE: =VTOP-VBOT, the calculation range [cm^-1].


C Continuum variables ...
      REAL CONTINK(IORDP1,MAXK,MAXK,MAXBIN)
      REAL CONVAL(NWAV),CONWAV(NWAV)
      REAL MATRIX(IORDP1,IORDP1),UNIT(IORDP1,IORDP1)
      DOUBLE PRECISION DMATRIX(IORDP1,IORDP1),DUNIT(IORDP1,IORDP1)
C IORDER is the order of the continuum polynomial.
C IORDP1 is the number of parameters in polynomial fit.
C CONTINK holds continuum polynomial for each bin.
C CONVAL holds the absorption coefficients at wavenumbers CONWAV prior to
C fitting continuum.
C MATRIX and UNIT are both used for matrix inversion when computing
C polynomials where insufficient points for least squares.

      COMMON /CONCOMK/CONTINK,CONVAL,CONWAV,MATRIX,UNIT


C Misc (and UNDOCUMENTED!!) variables or variables defined in the code ...
      INTEGER I,J,K,L,LINE,IPTF
      REAL XMASS,GETMASS,PARTF
      REAL DPEXP,LINECONTRIB
C XMASS: Molar mass [g mole^-1].
C DPEXP: Returns a REAL value after performing exponentiation in double
C precision.


C Partition function variables ...
      REAL TCORS1,TCORS2,TCORDW,VOUT(MAXBIN),YOUT(MAXBIN)
      REAL TCORS1LCO,TCORS2LCO
C TCORS1: the partition function temperature dependence times absorber
C fraction.
C TCORS2: the temperature dependence for the Boltzman distribution.
C TCORDW: the temperature dependence of the doppler width.


C Variables needed for the new lineshapes ...
      REAL PI
      PARAMETER (PI=3.1415927)

      REAL FNH3,FH2
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C******************************** CODE *********************************

cc      WRITE(*,*)'LBL_KCONT :: VMIN, VMAX = ',VMIN,VMAX
cc      WRITE(*,*)'LBL_KCONT :: WING, VREL = ',WING,VREL
cc      WRITE(*,*)'LBL_KCONT :: PRESS, TEMP = ',PRESS,TEMP
cc      WRITE(*,*)'LBL_KCONT :: IDGAS, ISOGAS = ',IDGAS,ISOGAS
cc      WRITE(*,*)'LBL_KCONT :: FRAC, IPROC = ',FRAC,IPROC
cc      WRITE(*,*)'LBL_KCONT :: IP, IT = ',IP,IT
cc      WRITE(*,*)'LBL_KCONT :: MAXDV = ',MAXDV

C Checking isotope included in model and setting mass for doppler width
C calculation
      XMASS=GETMASS(IDGAS,ISOGAS)

cc      WRITE(*,420)GASNAM(IDGAS),IDGAS,ISOGAS,XMASS
cc420   FORMAT('LBL_KCONT.f :: ',1A6,'(',I2,',',I2,') has mass = ',F7.2)

C NOTE: TCORS1 includes factor of 1.e-27 for scaling of stored lines


      TCORS1 = PARTF(IDGAS,ISOGAS,TEMP,IPTF)*1.E-27
      TCORDW = 4.301E-7*SQRT(TEMP/XMASS)
      TCORS2 = 1.439*(TEMP - 296.)/(296.*TEMP)

      IF(IJLCO.GT.0)THEN
       IF(ABS(TCALCLCO-296.0).GT.0.5)THEN
        TCORS1LCO = 1E-27*PARTF(IDLCO,ISOLCO,TEMP,IPTFLCO)/
     1       PARTF(IDLCO,ISOLCO,TCALCLCO,IPTFLCO)
        TCORS2LCO = 1.439*(TEMP - TCALCLCO)/(TCALCLCO*TEMP)
C        print*,TEMP,TCALCLCO,TCORS1LCO
       ELSE
        TCORS1LCO = TCORS1
        TCORS2LCO = TCORS2
       ENDIF
      ENDIF

      VBOT = VMIN - VREL
      VTOP = VMAX + VREL
      IF(VBOT.LT.0.0)VBOT = 0.0
      RANGE = VTOP - VBOT
      NBIN = INT(RANGE/WING) + 1
      if(idiag.gt.0)then
       print*,'LBL_KCONT: NBIN, RANGE, WING = ',NBIN, RANGE, WING
       print*,'LBL_KCONT: VMIN, VMAX, VBOT, VTOP = ',
     1 VMIN,VMAX,VBOT,VTOP
      endif
      IF(NBIN.GT.MAXBIN)THEN
        WRITE(*,*)'LBL_KCONT.f :: *ERROR* NBIN > MAXBIN'
        WRITE(*,*)'Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'NBIN, MAXBIN = ',NBIN,MAXBIN
        WRITE(*,*)'VBOT,VTOP,RANGE,WING = ',VBOT,VTOP,RANGE,WING
	STOP
      ENDIF
	
      IF(IP.GT.MAXK)THEN
        WRITE(*,*)'LBL_KCONT.f :: *ERROR* IP > MAXK'
        WRITE(*,*)'Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'IP, MAXK = ',IP,MAXK
	STOP
      ENDIF

      IF(IT.GT.MAXK)THEN
        WRITE(*,*)'LBL_KCONT.f :: *ERROR* IT > MAXK'
        WRITE(*,*)'Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'IT, MAXK = ',IT,MAXK
	STOP
      ENDIF

      DO 127 I=1,NBIN
        VBIN(I) = VBOT + (I - 1)*WING
127   CONTINUE

C-----------------------------------------------------------------------
C Initialising continuum polynomial
      DO 16 I=1,IORDP1
        DO 14 J=1,NBIN
          CONTINK(I,IP,IT,J) = 0.0
14      CONTINUE
16    CONTINUE

C-----------------------------------------------------------------------
C Read in the lines if necessary
      IREAD=0
      IF(IEXO.EQ.0)THEN
C      Same line database for all temperatures. Only need to read in lines
C      for IP=1 and IT=1
       IF(IP.EQ.1.AND.IT.EQ.1)THEN
        IREAD=1
       ENDIF
      ELSE
C      Different line database for different temperatures. 
C      Assume loop over pressure takes place within loop over temperature in
C      calling program)
       IF(IP.EQ.1)THEN
        IREAD=1
       ENDIF
      ENDIF

      IF(IREAD.EQ.1)THEN
	NGAS = 1
        CALL LOADBINS(WING,NGAS,IDGAS,ISOGAS)

        WRITE(*,*)'LBL_KCONT :: computing the continuum polynomial from'
        WRITE(*,*)'NBIN bins where NBIN = ',NBIN
        WRITE(*,*)'Temperature = ',TEMP
        WRITE(*,*)'Pressure = ',PRESS
      ENDIF

C=======================================================================
C
C	Now step through each bin computing continuum polynomial of line
C	wings from other bins if LBL calculation. Only need to compute for
C	central bins if monochromatic calc.
C      
C	For other models which use direct linedata, read in the lines
C	for all the bins.
C
C=======================================================================

      FSTBIN = INT((VMIN - VBOT)/WING)
      LSTBIN = INT((VMAX - VBOT)/WING) + 2
      IF(FSTBIN.LT.1)THEN
       FSTBIN=1
      ENDIF

      DO 444 K=1,IORDP1
C Wavenumbers to compute continuum at
        CONWAV(K) = FLOAT(K - 1)*WING/FLOAT(IORDER)
444   CONTINUE
C Setting up matrix of wavenumbers for use in polynomial calculation
      DO 18 K=1,IORDP1
        MATRIX(1,K) = 1.0
        DMATRIX(1,K) = 1.0
        DO 19 J=2,IORDP1
          MATRIX(J,K) = MATRIX(J-1,K)*CONWAV(K)
          DMATRIX(J,K) = DMATRIX(J-1,K)*CONWAV(K)
19      CONTINUE
18    CONTINUE


      L = IORDP1
C Find the inverse of the matrix
      CALL DMATINV(DMATRIX,L,IORDP1,DUNIT)

      DO K=1,IORDP1
       DO J=1,IORDP1
        UNIT(J,K)=SNGL(DUNIT(J,K))
       ENDDO
      ENDDO

      DO 13 I=1,NBIN
C Computing continuum for all except adjacent bins
        DO 15 J=FSTBIN,LSTBIN
          IF(ABS(I-J).LE.1)GOTO 15
C For each line
            DO 22 K=1,IORDP1
              CONVAL(K) = 0.0
22          CONTINUE
C Compute absorption coefficient for normal incidence
            DO 20 LINE=FSTLIN(I),LSTLIN(I)

              DO 21 K=1,IORDP1
                VV=DBLE(VBIN(J)+CONWAV(K))
                DV=SNGL(VV-VLIN(LINE))
                IF(ABS(DV).LE.MAXDV)THEN
                 FNH3=-1.0
                 FH2=-1.0
                 CONVAL(K)=CONVAL(K)+LINECONTRIB(IPROC,IDGAS,VV,TCORDW,
     1 TCORS1,TCORS2,PRESS,TEMP,FRAC,VLIN(LINE),SLIN(LINE),ELIN(LINE),
     2 ALIN(LINE),SBLIN(LINE),TDW(LINE),TDWS(LINE),LLQ(LINE),
     3 DOUBV(LINE),FNH3,FH2)
C                IF(K.EQ.1)print*,'AA',VV,K,IP,IT,J,CONVAL(K)

                ENDIF
21           CONTINUE
20         CONTINUE


           DO 23 K=1,IORDP1
             DO 314 L=1,IORDP1
               CONTINK(K,IP,IT,J) = CONTINK(K,IP,IT,J) + 
     1         UNIT(L,K)*CONVAL(L)
314          CONTINUE
23         CONTINUE
15      CONTINUE
13    CONTINUE  

C     Extra continuum for linedata files that are so large that they
C     have been stripped of weaker lines.
      IF(IJLCO.GT.0)THEN
       IF(IDGAS.EQ.IDLCO.AND.ISOGAS.EQ.ISOLCO)THEN
        if(idiag.gt.0)print*,'Adding LCO for gas : ',IDLCO,ISOLCO
        FNH3=-1.0
        FH2=-1.0

        CALL CALC_LCO(VBIN,NBIN,WING,MAXDV,IPROC,IDGAS,TCORDW,
     1   TCORS1LCO,TCORS2LCO,PRESS,TEMP,FRAC,FNH3,FH2,VOUT,YOUT)
        
        DO 301 I=1,NBIN
C         print*,I,VBIN(I),PRESS,TEMP,VOUT(I),YOUT(I)
          IF (ISNAN(YOUT(I)))THEN
           print*,'NAN gas come out of CALC_LCO'
           print*,I,VBIN(I)
           stop
          ENDIF
         DO 302 K=1,IORDP1
          VV=DBLE(VBIN(I)+CONWAV(K))
          CALL VERINT(VOUT,YOUT,NBIN+1,CONVAL(K),SNGL(VV))
C          print*,'A',VV,CONVAL(K),CONTINK(K,IP,IT,I)
302      CONTINUE
         DO 203 K=1,IORDP1

C            IF(K.EQ.1)print*,'AA',VV,IP,IT,I,CONTINK(K,IP,IT,I),
C     1        CONVAL(K)

             DO 304 L=1,IORDP1
               CONTINK(K,IP,IT,I) = CONTINK(K,IP,IT,I) + 
     1         UNIT(L,K)*CONVAL(L)
304          CONTINUE
C          print*,'B',VV,CONVAL(K),CONTINK(K,IP,IT,I)          
203      CONTINUE

301     CONTINUE
       ENDIF
      ENDIF

      RETURN

      END
************************************************************************
************************************************************************
