      SUBROUTINE LBL_KCONT_EXO3(INORM,FIRST,NLIN,VBOT,VMIN,
     1 VMAX,WING,IP,PRESS,TEMP,IDGAS,ISOGAS,FRAC,IPROC,MAXDV)
C     $Id:
C***********************************************************************
C_TITL:	LBL_KCONT_EXO3.f
C
C_DESC:	General purpose absorption coefficient and cumulative 
C	K-distribution calculation routine. Subroutine to calculate the
C	line continuum in the range VMIN-VMAX in bins of width WING.
C
C_ARGS:	Input variables:
C	VMIN		REAL	Wavenumber of start of complete range.
C	VMAX		REAL	Wavenumber at end of range.
C	WING		REAL	Bin width.
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
C
C	../includes/*.f variables:
C	VLIN(MAXLIN)	REAL	Line position [cm^-1].
C	SLIN(MAXLIN)	REAL	Line strength [cm^-1 molecule^-1 cm^-2] at
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
      INTEGER IP,IT,IDGAS,ISOGAS,IPROC,INORM,IBIN,NLIN,LINE1
      REAL VMIN,VMAX,WING,MAXDV,VBOT
      REAL FRAC,PRESS,TEMP


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
C ../includes/parcom.f stores the parameter values such as MAXLAY,


C General variables ...
      INTEGER NGAS,FIRST,ISUM
      REAL DV,VTOP,RANGE,VTMP1,TAUTMP1,VTMP2,TAUTMP2,VTMP3,TAUTMP3


C Continuum variables ...
      INTEGER IREAD
      INTEGER MP,MT
      PARAMETER (MP=20,MT=20)
      REAL CONTINK1(IORDP1,MP,MAXBIN),TAUSTORE(IORDP1,MAXBIN)
      REAL CONVAL(NWAV),CONWAV(NWAV)
      REAL MATRIX(IORDP1,IORDP1),UNIT(IORDP1,IORDP1)
      DOUBLE PRECISION DMATRIX(IORDP1,IORDP1),DUNIT(IORDP1,IORDP1)
C IORDER is the order of the continuum polynomial.
C IORDP1 is the number of parameters in polynomial fit.
C CONTINK1 holds continuum polynomial for each bin.
C CONVAL holds the absorption coefficients at wavenumbers CONWAV prior to
C fitting continuum.
C MATRIX and UNIT are both used for matrix inversion when computing
C polynomials where insufficient points for least squares.

      COMMON /CONCOMK1/CONTINK1,CONVAL,CONWAV,MATRIX,UNIT,IREAD


C Misc (and UNDOCUMENTED!!) variables or variables defined in the code ...
      INTEGER I,J,K,L,LINE
      REAL ABSCO,AD,X,Y
      REAL DPEXP,SUBLINE
C ABSCO is, I think, the absorption coefficient strength
C AD is, I think, the Doppler-broadened line width
C X is the wavenumber-distance away from the unperturbed line frequency
C multiplied by the Doppler-broadened line width (for some reason)
C Y is, I think, the Lorentz/collision-broadened line width (divided by
C the Doppler-broadened line width for some reason).
C DPEXP: Returns a REAL value after performing exponentiation in double
C precision.


C Partition function variables ...
      REAL TCORS1,TCORS2,TRATIO,TCORDW,TSTIM
C TCORS1: the partition function temperature dependence times absorber
C fraction.
C TCORS2: the temperature dependence for the Boltzman distribution.
C TCORDW: the temperature dependence of the doppler width.
C TRATIO: =296/T for use in calculating T dependence of line width.
C TSTIM: the correction for stimulated emission.


C Variables needed for the new lineshapes ...
      REAL PI
      PARAMETER (PI=3.1415927)

      REAL VV,FNH3,FH2

C******************************** CODE *********************************

c      WRITE(*,*)'LBL_KCONT_EXO1 :: VMIN, VMAX = ',VMIN,VMAX
c      WRITE(*,*)'LBL_KCONT_EXO1 :: WING = ',WING
c      WRITE(*,*)'LBL_KCONT_EXO1 :: PRESS, TEMP = ',PRESS,TEMP
c      WRITE(*,*)'LBL_KCONT_EXO1 :: IDGAS, ISOGAS = ',IDGAS,ISOGAS
c      WRITE(*,*)'LBL_KCONT_EXO1 :: FRAC, IPROC = ',FRAC,IPROC
c      WRITE(*,*)'LBL_KCONT_EXO1 :: MAXDV = ',MAXDV

C-----------------------------------------------------------------------
C Initialising continuum polynomial
      IF(INORM.EQ.1)THEN 
       Print*,'Initialising continuum arrays'
       DO 16 I=1,IORDP1
        DO 17 K=1,MP
         DO 14 J=1,NBIN        
          CONTINK1(I,K,J) = 0.0
14       CONTINUE
17      CONTINUE
16     CONTINUE


       DO 444 K=1,IORDP1
C      Wavenumbers to compute continuum at
        CONWAV(K) = FLOAT(K - 1)*WING/FLOAT(IORDER)
444    CONTINUE
C      Setting up matrix of wavenumbers for use in polynomial calculation
       DO 18 K=1,IORDP1
        MATRIX(1,K) = 1.0
        DMATRIX(1,K) = 1.0
        DO 19 J=2,IORDP1
          MATRIX(J,K) = MATRIX(J-1,K)*CONWAV(K)
          DMATRIX(J,K) = DMATRIX(J-1,K)*CONWAV(K)
19      CONTINUE
18     CONTINUE


       L = IORDP1
C      Find the inverse of the matrix
       CALL DMATINV(DMATRIX,L,IORDP1,DUNIT)

       DO K=1,IORDP1
        DO J=1,IORDP1
         UNIT(J,K)=SNGL(DUNIT(J,K))
        ENDDO
       ENDDO
      ENDIF

C-----------------------------------------------------------------------


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

C      PRINT*,'VMIN,VBOT,WING',VMIN,VBOT,WING
      FSTBIN = (VMIN - VBOT)/WING
      IF(FSTBIN.LT.1)FSTBIN=1
      LSTBIN = (VMAX - VBOT)/WING + 2

      print*,'LBL_KCONT_EXO3: FSTBIN,LSTBIN = ',FSTBIN,LSTBIN

      TRATIO = 296./TEMP

      DO K=1,IORDP1
       DO J=FSTBIN,LSTBIN
        TAUSTORE(K,J)=0.0
       ENDDO
      ENDDO

      DO 20 LINE1=1,NLIN
       LINE=FIRST+LINE1-1
       IF(LINE.GT.MAXLIN)LINE=LINE-MAXLIN
       IBIN = 1+ (VLIN(LINE)-VBOT)/WING

       ABSCO = ABSCO_ARR(LINE)
       AD = AD_ARR(LINE)
       Y = Y_ARR(LINE)

C Computing continuum for all except adjacent bins
       DO 15 J=FSTBIN,LSTBIN
          IF(ABS(IBIN-J).LE.1)GOTO 15
C For each line

C Compute absorption coefficient for normal incidence

            DO 21 K=1,IORDP1
                VV = VBIN(J)+CONWAV(K)
                X = (VV - VLIN(LINE))/AD
                DV = X*AD
                FNH3=-1.0
                FH2=-1.0
                CONVAL(K)=SUBLINE(IDGAS,PRESS,TEMP,IPROC,VV,VLIN(LINE),
     1           ABSCO,X,Y,AD,FNH3,FH2,LLQ(LINE),DOUBV(LINE))

               TAUSTORE(K,J)=TAUSTORE(K,J)+CONVAL(K)

21          CONTINUE

15     CONTINUE

20    CONTINUE

      DO 606 J=FSTBIN,LSTBIN
             DO 23 K=1,IORDP1
              DO 314 L=1,IORDP1
               CONTINK1(K,IP,J) = CONTINK1(K,IP,J) + 
     1         UNIT(L,K)*TAUSTORE(L,J)
314           CONTINUE
23           CONTINUE
606   CONTINUE

      RETURN

      END
************************************************************************
************************************************************************
