      SUBROUTINE LBL_FKNEW(IWAVE,VSTART,VEND,PRESS,TEMP,IDGAS,ISOGAS,
     1 IPROC,IP,IT,FRAC,MAXDV,IPTF,NPOINT)
C     $Id: lbl_fknew.f,v 1.7 2011-09-06 15:32:20 irwin Exp $
C***********************************************************************
C_TITL:	LBL_FKNEW.f
C
C_DESC: General purpose absorption coefficient and cumulative
C	k-distribution calculation routine. Calculates the cumulative
C	k-distribution for a spectral interval VSTART to VEND for a single
C	gas. This is done by first generating the lbl absorption
C	coefficient spectra and then analysing these according to the
C	equation given by Lacis and Oinas (1991):
C	         f(k) = (1/(V2-V1))*SUM(ABS(dV/DK))
C	 f(k) is then summed to give the cumulative k distribution.
C
C_ARGS:	Input variables:
C	IWAVE		INTEGER	Wavespace of final k-table
C	VSTART		REAL	Wavenumber of beginning of range.
C	VEND		REAL	Wavenumber of end of range.
C	PRESS		REAL	Total pressure [atm].
C	TEMP		REAL	Temperature [Kelvin].
C	IDGAS		INTEGER	The LOCAL gas identifier for each gas.
C				The local gas identifier agrees with
C				HITRAN id's as far as possible (Jan 
C				1988) with extra ID's for gases not in the
C				HITRAN compilation. eg. those in the GEISA
C				compilation.
C	ISOGAS		INTEGER	The LOCAL isotopic identifier for each 
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
C	FRAC		REAL	Broadening fraction (0 = air-broadened,
C				1 = self-broadened).
C	IPROC		INTEGER	Line wing identifier.
C	IP		INTEGER	ID of corresponding pressure in
C				accompanying CONTINK continuum array.
C	IT		INTEGER	ID of corresponding temperature in
C				accompanying CONTINK continuum array.
C	MAXDV		REAL	Line wing cut-off parameter: The maximum
C				wavenumber-distance away within which to
C				consider the contribution of the line
C				wings.
C	IPTF		INTEGER Partition function flag for CH4
C
C	../includes/*.f variables:
C	VLIN(MAXLIN)	REAL	Line position [cm^-1].
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
C	Output variable:
C	K_G(MAXG)		REAL	Calculated k-distribution.
C
C_FILE: No files openned.
C
C_CALL:	calc_kdist	Subroutine to calculate the cumulative 
C			k-distribution in a band by binning a
C			pre-calculated lbl absorption coefficient
C			spectrum.
C	dpexp		Performs exponentiation in double precision.
C	partf		Uses four-term polynomial fits to compute ratio
C			of total partition functions for use by line by
C			line programs.
C
C_HIST:	17feb93	PGJI	Modified from genlbl.f
C	02mar93	PGJI	Reduced and simplified to be of more general use.
C	10may94	PGJI	Revised and further debugged.
C	1Sept95	ALW	Revised and further debugged.
C	3Nov95	ALW	Introduced fracb variable.
C	30Mar00	PGJI	Modified to read in continuum arrays.
C	8/3/04	NT added iproc=6 option, Rosenkrantz-Ben-Reuven lineshape
C			with voigt lineshape as default (NH3 only)
C	4/4/06	NT added in some pre-tabluation of variables to speed
C			up the code.
C***************************** VARIABLES *******************************

      IMPLICIT NONE

C Definition of input parameters ...
      INTEGER IDGAS,ISOGAS,IPROC,IP,IT,IWAVE
      REAL PRESS,TEMP,VSTART,VEND,MAXDV,FRAC,FPOINT


C The include files ...
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
C ../includes/bincom.f stores the line bin variables (including NLINES, 
C FSTLIN, LSTLIN) and band parameters.
      INCLUDE '../includes/dbcom.f' 
C ../includes/dbcom.f stores the linedata base variables (e.g. MASSNO and
C RELABU).
      INCLUDE '../includes/lincom.f'
C ../includes/lincom.f stores the linedata variables including MAXLIN,
C VLIN, SLIN, ALIN, ELIN, SBLIN, TDW, TDWS and that lot).
      INCLUDE '../includes/parcom.f'
C ../includes/parcom.f stores the parameter values such as MAXLAY,


      REAL OUTPUT(MPOINT),TX(MPOINT)
      CHARACTER*1 ANS

      COMMON /SPECTRUM/ OUTPUT


C Continuum variables ...
      INTEGER ISUM,IREAD,IPTF
      REAL TAUTMP,PARTF
      REAL CONTINK(IORDP1,MAXK,MAXK,MAXBIN)
      REAL CONVAL(NWAV),CONWAV(NWAV)
      REAL MATRIX(IORDP1,IORDP1),UNIT(IORDP1,IORDP1)
C IORDER: order of the continuum polynomial.
C IORDP1: number of parameters in polynomial fit
C TAUTMP: holds the normal incidence case optical depth
C CONTINK: holds continuum polynomial for each bin.
C CONVAL holds the absorption coefficients at wavenumbers CONWAV prior to
C fitting continuum
C MATRIX and UNIT are both used for matrix inversion when computing
C polynomials where insufficient points for least squares.

      COMMON /CONCOMK/ CONTINK,CONVAL,CONWAV,MATRIX,UNIT,IREAD


C General variables ...
      REAL K_BOLTZ,C_LIGHT,R_GAS
      PARAMETER (K_BOLTZ=1.38E-23,C_LIGHT=2.998E8,R_GAS=8.3143)
C K_BOLTZ: Boltzman constant [J K^-1 molecule^-1].
C C_LIGHT: Speed of light [m sec^-1].
C R_GAS: Universal gas constant [J K^-1 mol^-1].

      INTEGER I,J,K,JBIN,IBIN1,IBIN2,NPOINT
      REAL VMIN,DELV,WING,DV
      REAL XPW,V0,XPD,MXWID
C I,J,K,JBIN: Incremators.
C IBIN1: =CURBIN - 1 
C IBIN2: =CURBIN + 1
C NPOINT: =3*(VEND-VSTART)/MXWID, the number of points within the
C interval. Must be a minimum of three.
C VMIN: Wavenumber [cm^-1] minimum.
C DELV: Spectral point-spacing [cm^-1].
C WING: Bin size [cm^-1] used store linedata.
C V0: Most probable speed assuming a Maxwellian distribution [m sec^-1].
C XPW: Average Lorentz half width [cm^-1].
C XPD: Average Doppler half width [cm^-1].
C MXWID: Root mean square of XPW, XPD.


C Misc (and UNDOCUMENTED!!) variables defined in the code ...
      INTEGER LINE,LABEL
      REAL XMASS,GETMASS,VTMP
      REAL V
c--------------------------------
cc NT*** these  are now arrays **
ccc      REAL ABSCO,AD,X,Y,DV
c--------------------------------
      REAL DPEXP,SUBLINE
C XMASS: Molar mass [g mole^-1].
C ABSCO is, I think, the absorption coefficient strength.
C AD is, I think, the Doppler-broadened line width.
C X is the wavenumber-distance away from the unperturbed line frequency
C multiplied by the Doppler-broadened line width (for some reason).
C Y is, I think, the Lorentz/collision-broadened line width (divided by
C the Doppler-broadened line width for some reason).
C DV is the wavenumber-distance away from the unperturbed line frequency.
C DPEXP: Returns a REAL value after performing exponentiation in double  
C precision.
C HUMLIC, GVOICO2 and PARTF are all programs that return a REAL value for
C use in calculating the lineshapes.


C Partition function variables ...
      REAL TCORS1,TCORS2,TRATIO,TCORDW,TSTIM
C TCORS1: Partition-function temperature-dependence times absorber 
C fraction.
C TCORS2: Temperature dependence for the Boltzman distribution.
C TRATIO: =296/TEMP, used in calculating T-dependence of line width.
C TCORDW: Temperature dependence of the doppler width.
C TSTIM: Correction for stimulated emission.


C Variables needed for the new lineshapes ...
      REAL PI
      PARAMETER (PI=3.1415927)

      REAL FNH3,FH2,X

c arrays to pre-tabulate some stuff so speed it up
      integer itest
c maxlin is read in from lincom.f
      integer line_done(maxlin)

C******************************** CODE *********************************

C Checking isotope included in model and setting mass for doppler width
C calculation
      XMASS=GETMASS(IDGAS,ISOGAS)

1     FORMAT(A)

cc      WRITE(*,420)GASNAM(IDGAS),IDGAS,ISOGAS,XMASS
cc420   FORMAT('LBL_FKNEW.f :: ',1A6,'(',I2,',',I2,') has mass = ',F7.2)

C NOTE: TCORS1 includes a factor of 1.e-27 for scaling of stored lines

      TCORS1 = PARTF(IDGAS,ISOGAS,TEMP,IPTF)*1.E-27
      TCORS2 = 1.439*(TEMP - 296.)/(296.*TEMP)
      TCORDW = 4.301E-7*SQRT(TEMP/XMASS)
      TRATIO = 296./TEMP

C=======================================================================
C
C	Determining NPOINT (the number of points within the interval
C	VEND - VSTART).
C
C	An adjustment of the fine-scale structure within the calculation
C	bin has been made from a fixed value to the root-mean-square of
C	the average Lorentz- (XPW) and Doppler- (XDW) half widths
C	(whichever is larger at that calculation pressure and 
C	temperature).
C
C	Previously, the fine-scale contribution of the lines
C	over the range of each FWHM-wide bin was calculated at a
C	user-inputted value (usually 0.001 cm^-1) which limited the number
C	of lines included in the calculation at low pressures, but was
C	excessive over-calculation at high pressures when the lines are
C	all very fat.
C
C	The Lorentz half width is given by ...
C	            P        Ts               where a_Ls, Ps and Ts are 
C	a_L = a_Ls --- SQRT(---)              the half-width, pressure and
C	            Ps       T                temperature at STP.
C
C	The Doppler half width is given by ...
C	       v0        2*Pi*R*T             where R = Uni. gas constant,
C	a_D = --- SQRT( ----------)                 M = Molecular mass.
C	       c            M
C
C
C	As I was unable to successfully code a method so as to include the
C	the actual Lorentz half-width at STP, I have substituted 0.1 in
C	its place -- after so much effort to get the calculatin correct
C	in the first place, this option seems a poor finish. However, that
C	was what I was able to do with the ability and time I have; better
C	luck to you. -Paul
C
C	For reference, the previous values were:
C		XPW = 0.1*PRESS
C		XPD = 4.301E-7*SQRT(TEMP/XMASS)
C
C=======================================================================
C
C      XPW = 0.1*PRESS*SQRT(TRATIO)
C      V0 = SQRT(2*K_BOLTZ*TEMP/XMASS)
C      XPD = (V0/C_LIGHT)*SQRT(2*PI*R_GAS*TEMP/XMASS)
C      MXWID = SQRT(XPD**2 + XPW**2)
      XPW = 0.1*PRESS*SQRT(TRATIO)
      V0=0.5*(VSTART+VEND)
      XPD = 4.301E-7*V0*SQRT(TEMP/XMASS)
      MXWID = SQRT(XPD**2 + XPW**2)
      print*,'MXWID=',MXWID

C Multiply by three so as to have the center point, and one at both VSTART
C (= VMIN - 0.5*FWHM) and VEND (= VSTART + FWHM).
C      NPOINT = 3*INT((VEND - VSTART)/MXWID)
C Actually a factor of 3 still gives quite a jumpy spectrum. Increase
C to 6.
      FPOINT = 6.*(VEND - VSTART)/MXWID
cc      WRITE(*,*)'LBL_FKNEW.f :: NPOINT = ',NPOINT

      IF(FPOINT.GE.FLOAT(MPOINT))THEN
cc        WRITE(*,*)'LBL_FKNEW.f :: *WARNING* NPOINT > MPOINT'
cc        WRITE(*,*)'NPOINT, MPOINT = ',NPOINT,MPOINT
cc        WRITE(*,*)'Setting NPOINT equal to MPOINT.'
        NPOINT = MPOINT - 1
      ELSE
        NPOINT=INT(FPOINT)
      ENDIF

      IF(FPOINT.LE.4.)THEN
cc        WRITE(*,*)'LBL_FKNEW.f :: *WARNING* NPOINT < 4 (= ',NPOINT,').'
cc        WRITE(*,*)'Setting NPOINT equal to 3.'
        NPOINT = 3
      ENDIF

      DELV = (VEND - VSTART)/NPOINT
      NPOINT = NPOINT + 1

C Initialising pointers and bins
      WING = VBIN(2) - VBIN(1)

c###################################################
c     pretabulate some line calculation parameters
      do itest=1,maxlin
         line_done(itest) = 0
      enddo

      do i=1,npoint
        V = VSTART + FLOAT(I - 1)*DELV
        CURBIN = INT((V - VBIN(1))/WING) + 1
        IBIN1 = CURBIN - 1
        IBIN2 = CURBIN + 1
        IF(IBIN1.LT.1)IBIN1 = 1
        IF(IBIN2.GT.NBIN)IBIN2 = NBIN

        DO JBIN=IBIN1,IBIN2
          DO LINE=FSTLIN(JBIN),LSTLIN(JBIN)
            if (line_done(line).eq.0) then
             line_done(line)=1
             tstim_arr(line) = (1.0 - dpexp(-1.439*vlin(line)/temp))/
     >          (1.0 - dpexp(-1.439*vlin(line)/296.0))
             absco_arr(line) = 
     >        sngl(slin(line)*tcors1*dpexp(tcors2*elin(line))*
     >        tstim_arr(line))
             ad_arr(line) = tcordw*vlin(line)
             y_arr(line) = (alin(line)*(1 - frac)*tratio**tdw(line) +
     >           (alin(line) - sblin(line))*frac*tratio**tdws(line))*
     >           press/ad_arr(line)
            endif
          enddo
        enddo
      enddo
c###################################################


C Infinite resolution
      DO 103 I=1,NPOINT
        V = VSTART + FLOAT(I - 1)*DELV
        ASSIGN 2001 TO LABEL
        GOTO 2000

2001    CONTINUE

103   CONTINUE

      GOTO 3000

C=======================================================================
C
C	Now compute the absorption coefficient.
C
C	This section calculates the transmission (etc.) at wavenumber
C	V. The output values are stored in the array OUTPUT(I). This would
C	be better as a subroutine, but there are so many parameters to
C	pass that it would take forever so both entry and exit is via a
C	GOTO statement.
C
C	First computing normal path optical depth (TAUTMP).
C
C=======================================================================

2000  CONTINUE

C First calculate continuum contribution. CURBIN is the current BIN for
C which calculations are being made.
      CURBIN = INT((V - VBIN(1))/WING) + 1
C Calculate the continuum absorption via the IORDP1 polynomial
C coefficients held in CONTINK.
      TAUTMP = CONTINK(1,IP,IT,CURBIN)
      VTMP = V - VBIN(CURBIN)
      DO 51 ISUM=2,IORDP1
        TAUTMP = TAUTMP + CONTINK(ISUM,IP,IT,CURBIN)*VTMP
        VTMP = VTMP*VTMP
51    CONTINUE

C The text within the IF statement below should probably be kept
C commented-out unless you are specifically debugging for errors. I make
C this suggestion because, I've found, it is not uncommon for the
C continuum TAU to be less than zero due to numerical inaccuracy (at least
C that is what I think the reason for it is). Otherwise, these diagnostic
C statements only serve to blow-up the size of your .prc file.
      IF(TAUTMP.LT.0.0)THEN
cc        WRITE(*,*)'LBL_FKNEW.f :: *ERROR* Continuum TAUTMP < 0.0!'
cc        WRITE(*,*)' '
cc        WRITE(*,*)'V, VBIN, VTMP = ',V,VBIN(CURBIN),VTMP
cc        WRITE(*,*)'I, CURBIN, TAUTMP = ',I,CURBIN,TAUTMP
cc        WRITE(*,*)'CONTINK = ',(CONTINK(K,IP,IT,CURBIN),K=1,3)
cc        WRITE(*,*)' '
cc        WRITE(*,*)'Resetting TAUTMP to zero.'
        TAUTMP = 0.0
      ENDIF


C=======================================================================
C
C	Step through and select the spectral model required
C
C=======================================================================
      IBIN1 = CURBIN - 1
      IBIN2 = CURBIN + 1

      IF(IBIN1.LT.1.OR.IBIN2.GT.NBIN)THEN
        WRITE(*,*)'LBL_FKNEW.f :: *ERROR* Either IBIN1 or IBIN2 is out'
        WRITE(*,*)'of range.'
        WRITE(*,*)' '
        WRITE(*,*)'IBIN1, IBIN2, NBIN = ',IBIN1,IBIN2,NBIN
        WRITE(*,*)' '
        WRITE(*,*)'Resetting according to ...'
        WRITE(*,*)'  IF (IBIN1.LT.1) IBIN1 = 1'
        WRITE(*,*)'  IF (IBIN2.GT.NBIN) IBIN2 = NBIN',NBIN
      ENDIF
      IF(IBIN1.LT.1)IBIN1 = 1
      IF(IBIN2.GT.NBIN)IBIN2 = NBIN

      DO 507 JBIN=IBIN1,IBIN2

C Compute absorption coefficient for normal incidence
        DO 52 LINE=FSTLIN(JBIN),LSTLIN(JBIN)
          DV=vlin(line)-v
          IF(ABS(DV).LE.MAXDV)THEN
           X  = abs(vlin(line)-v)/ad_arr(line)
           FNH3=-1.0
           FH2=-1.0
           TAUTMP=TAUTMP+SUBLINE(IDGAS,PRESS,TEMP,IPROC,V,
     1  VLIN(LINE),ABSCO_arr(line),X,Y_arr(line),ad_arr(line),
     1  FNH3,FH2,LLQ(line),DOUBV(line))
          ENDIF

52      CONTINUE
        IF(TAUTMP.LT.0.0)THEN
C          WRITE(*,*)' LBL_FKNEW :: *ERROR* TAUTMP < 0.0!'
C          WRITE(*,*)' V, VBIN, VTMP = ',V,VBIN(CURBIN),VTMP
C          WRITE(*,*)' I, CURBIN, TAUTMP = ',I,CURBIN,TAUTMP
C          WRITE(*,*)' '
C          WRITE(*,*)' Stopping program.'
C          STOP
           TAUTMP=0.0
        ENDIF
507   CONTINUE

C Now scaling for each path
      OUTPUT(I) = TAUTMP
      TX(I) = V
      GOTO LABEL

3000  CONTINUE

C Calculate cumulative k-distribution from LBL spectrum. NOTE: OUTPUT is
C passed via the /SPECTRUM/ common block for speed.

C      open(12,file='spectrum.out',status='unknown')
C      write(12,*)npoint
C      do i=1,npoint
C       write(12,*)tx(i),output(i)
C      enddo
C      close(12)

C      print*,'press a key to continue'
C      read(5,1)ans


      RETURN

      END
