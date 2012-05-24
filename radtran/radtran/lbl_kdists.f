      SUBROUTINE LBL_KDISTS(PP,PRESS,TEMP,NGAS,IDGAS,ISOGAS,IPROC,VMIN,
     1DELV,NPOINT,G_ORD,DELG,K_G,NG,LAYER1,IPTF)
C     $Id: lbl_kdists.f,v 1.5 2011-06-17 15:40:26 irwin Exp $
C-------------------------------------------------------------------------
C_TITLE: GENLBL_KDIST: calculates the cumulative K-Distribution for a spectral
C  	 interval for a mixture of gases at a number of pressures and 
C        temperatures. This is done by first generating the lbl absorption
C        coefficient spectra and then analysing these according to the 
C  	 equation given by Lacis and Oinas (1991):
C	         f(k) = (1/(V2-V1))*SUM(ABS(dV/DK))
C	 f(k) is then summed to give the cumulative k distribution.
C
C-------------------------------------------------------------------------
C	Input Variables
C
C        PP:REAL	        partial pressure of each gas (atm)
C        PRESS:REAL             total pressure for each layer (atm)
C        TEMP:REAL              temperature (Kelvin) of each layer
C        NGAS:INTEGER     	number of gases to consider
C        IDGAS(NGAS):INTEGER   	the LOCAL gas identifier for each gas
C                         	The local gas identifier agrees with HITRAN id's
C                         	as far as possible (Jan 1988) with extra id's
C                         	for gases not in the HITRAN compilation. eg.
C                         	those in the GEISA compilation
C        ISOGAS(NGAS):INTEGER   the local isotopic identifier, if zero all
C                         	isotopes of the gas are included.
C                         	Isotope id's also agree as far as possible with
C                         	HITRAN id's. Similar (1,2,3...n) id's have been
C                         	defined for additional gases.               
C                         	If zero then line strengths are used as 
C				tabulated (i.e. corrected for normal 
C				terrestrial distribution)
C                         	If a specific isotope is selected then its 
C				strength is scaled to the pure isotope using 
C				the HITRAN data base values
C        VMIN:REAL        	minimum wavenumber for lbl spectra
C        DELV:REAL        	wavenumber spacing for lbl spectra
C        NPOINT:INTEGER   	the number of points in lbl spectra
C   	                      	i.e. each wavenumber, v=vmin+(i-1)*delv for
C                         	i=1 to npoint
C	 IPTF:INTEGER		Partition function flag
C-------------------------------------------------------------------------
C	Output Variable
C
C
C-------------------------------------------------------------------------
C
C_DESCR: General purpose absorption coefficient and cumulative K-Distribution
C	 calculation routine.
C        The spectral region is divided into bins WING wavenumbers wide.
C        Line data is stored in blocks corresponding to the WING sized
C        bins. Each calculation first checks that the appropriate
C        bins are in memory then swaps bins in and out if necessary. Lines 
C        from the bin at that wavenumber and the two adjacent bins are treated 
C        explicitly.
C        Absorption by lines from other bins is calculated from
C        precomputed continuum polynomials based upon fits to the line
C        wings.
C
C	 Resulting internal lbl spectra are analysed to output the cumulative
C	 K-distribution for each layer by calling the subroutine KDIST_LBL
C
C-------------------------------------------------------------------------
C
C_CALLS: 
C
C_BUGS:
C
C_HIST:  17feb93 PGJI modified from genlbl.f
C	 02mar93 PGJI reduced and simplified to be of more general use
C        10may94 PGJI revised and further debugged
C        1Sept95 ALW revised and further debugged
C         3Nov95 ALW introduced fracb variable
C	18/3/04	NT added iproc=6 option, Rosenkrantz-Ben-Reuven lineshape
C			with voigt lineshape as default (NH3 only)
C------------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f' 

      INTEGER NGAS,IDGAS(NGAS),LAYER1,IWAVE
      INTEGER ISOGAS(NGAS),IPROC(NGAS),NPOINT,NG,JBIN
      REAL FNH3,FH2
      INTEGER IGAS,JH2,JNH3,J1
      REAL PI,LINECONTRIB
      PARAMETER (PI=3.1415927)
      REAL PP(NGAS),PRESS,DPEXP,NU,DV1,WY
      REAL G_ORD(MAXG),K_G(MAXG),DELG(MAXG)
      REAL TEMP,VMIN,DELV,WING
C------------------------------------------------------------------------------
C     general parameters
      integer mpoint,MFIL,NFIL
      PARAMETER (MPOINT=50000,MFIL=1000)
      REAL OUTPUT(MPOINT),VFIL(MFIL),TFIL(MFIL)

C------------------------------------------------------------------------------
      COMMON /SPECTRUM/ OUTPUT
C------------------------------------------------------------------------------


C------------------------------------------------------------------------------
      INCLUDE '../includes/lincom.f'
      INCLUDE '../includes/parcom.f'
      INCLUDE '../includes/bincom.f'
C------------------------------------------------------------------------------
C     continuum variables
      INTEGER ISUM
C     IORDER is the order of the continuum polynomial
      REAL CONTINK(IORDP1,MAXLAY,MAXBIN)
      REAL CONTIN(IORDP1,MAXLAY,MAXBIN),CONTMP(IORDP1)
      REAL CONVAL(NWAV),CONWAV(NWAV),FILCON(IORDP1,MAXBIN)
      REAL MATRIX(IORDP1,IORDP1),UNIT(IORDP1,IORDP1)
C     CONTIN holds continuum polynomial for each bin, CONTMP is a
C     temporary store for these, CONVAL holds the absorption coefficients
C     at wavenumbers CONWAV prior to fitting continuum
C     MATRIX and UNIT are both used for matrix inversion when computing
C     polynomials where insufficient points for least squares
      COMMON /CONCOM/CONTIN,CONTINK,CONTMP,CONVAL,CONWAV,FILCON,
     &MATRIX,UNIT
C------------------------------------------------------------------------------
C     general variables
      REAL TAUTMP,FRAC(MAXGAS),FRACB(MAXGAS),QTOT
C     TAUTMP holding the normal incidence case optical depth
C     FRAC holds volume fraction of each gas for use in combining gases
C     FRACB holds volume fraction of each gas for use in broadening calculations
C     (FRAC and FRACB can differ in calls from CALC_TABLE, in that only 1 gas
C     may be present, but it can be completely foreign broadened).
C
C     misc variables or variables defined in the code
      INTEGER I,J,K,LINE,IPTF
      INTEGER IORD1,ISOMAX
      REAL DOPMAX,XMASS
      REAL VTMP, V
      REAL ABSCO,AD,X,Y
      REAL PARTF
      INTEGER LABEL
C------------------------------------------------------------------------------
C     note dbcom defines the linedata base variables. it is not normally stored
C     in the same directory as the rest of the code
      INCLUDE '../includes/dbcom.f' 
C------------------------------------------------------------------------------
C     partition function variables
      REAL TCORS1(MAXGAS),TCORS2,TRATIO
      REAL TCORDW(MAXGAS)
      REAL TSTIM
C     TCORS1 is the partition function temperature dependence times absorber
C     fraction
C     TCORS2 is the temperature dependence for the Boltzman distribution
C     TCORDW is the temperature dependence of the doppler width
C     TRATIO is just 296/T for use in calculating T dependence of line width
C     TSTIM is the correction for stimulated emission
C------------------------------------------------------------------------------

C     code section
C----------


C      print*,PP,PRESS,TEMP,NGAS,IDGAS,ISOGAS,IPROC,VMIN,
C     1DELV,NPOINT,G_ORD,K_G,NG,VBOT,LAYER1


C     Identify which gas is NH3 and H2 (if applicable) so that the
C     weird new NH3 lineshape can be used.
      JH2 = -1
      JNH3 = -1
      DO IGAS=1,NGAS
       IF(IDGAS(IGAS).EQ.39.AND.
     1	(ISOGAS(IGAS).EQ.0.OR.ISOGAS(IGAS).EQ.1))THEN
        JH2 = IGAS
       ENDIF
       IF(IDGAS(IGAS).EQ.11.AND.
     1	(ISOGAS(IGAS).EQ.0.OR.ISOGAS(IGAS).EQ.1))THEN
        JNH3 = IGAS
       ENDIF
      ENDDO

      IF(NPOINT.GT.MPOINT)THEN
 	PRINT*,'Error in lbl_kdists: NPOINT > MPOINT'
        PRINT*,'NPOINT,MPOINT = ',NPOINT,MPOINT
        STOP
      ENDIF


C     defining some useful variables
      IORD1=IORDP1
C     IORD1 is the number of parameters in polynomial fit
      ISOMAX=MAXISO
C     ISOMAX used in subroutine calls - maximum number of isotopes that
C     can be used (per gas)
      LINMAX=MAXLIN
C----------
C     computing volume fractions


C     now precomputing temperature dependences etc.
      QTOT=0.
      DO 16 I=1,NGAS
   	FRACB(I)=PP(I)/PRESS
   	FRAC(I)=FRACB(I)
        QTOT=QTOT+FRAC(I)
C     checking isotope included in model and setting mass for doppler
C     width calculation
      IF(ISOGAS(I).EQ.0)THEN
C       have to compute mass for normal combination if using all isotopes
        XMASS=0.
        DO 255 J=1,DBNISO(IDGAS(I))
        XMASS=XMASS+MASSNO(J,IDGAS(I))*RELABU(J,IDGAS(I))
255     CONTINUE
C       it is possible that the relative abundance for one of the isotopes is
C       wrong (eg set to 1 because unknown) so checking that the final
C       mass is within 20% of the mass of the main isotope
        IF(MASSNO(1,IDGAS(I)).GT.1.E-32.AND.
     1  ABS((XMASS-MASSNO(1,IDGAS(I)))/MASSNO(1,IDGAS(I)))
     1  .GT.0.2)THEN
         XMASS=MASSNO(1,IDGAS(I))
         WRITE(*,287)
287      FORMAT(' %warning - using main isotope mass number')
         END IF
      ELSE
        IF(ISOGAS(I).LE.DBNISO(IDGAS(I)))GOTO 418
C       isotopes are included in arrays in same order as new HITRAN 86
C       definitions 
        WRITE(*,417)ISOGAS(I),IDGAS(I)
417     FORMAT(' %model doesn"t include isotope ',I5,' for gas ',I2)
        STOP
418     CONTINUE
        XMASS=MASSNO(ISOGAS(I),IDGAS(I))
C
      END IF
C     WRITE (*,420) IDGAS(I), ISOGAS(I), GASNAM(IDGAS(I)), XMASS
420   FORMAT(' %gas:',I2,'/',I2,2X,1A6,' mass=',F7.2)

C     note TCORS1 includes factor of 1.e-27 for scaling of stored lines
      IF(ISOGAS(I).EQ.0)THEN
        K=1
       ELSE
        K=ISOGAS(I)
      END IF
      TCORS1(I)=PARTF(IDGAS(I),K,TEMP,IPTF)*1.E-27
      TCORDW(I)=4.301E-7*SQRT(TEMP/XMASS)
16    CONTINUE

c      IF (IFLAG.EQ.2) THEN	! Required for Calc_table
c	QTOT=1
c     ENDIF

      TRATIO=296./TEMP
      TCORS2=1.439*(TEMP-296.)/(296.*TEMP)


C----------
C     initialising pointers and bins
C     increasing WING by maximum doppler shift so that even for shifted
C     layers enough lines are read in for accurate calculation
      DOPMAX=0.
      WING=VBIN(2)-VBIN(1)


C----------
C     infinite resolution
      DO 103 I=1,NPOINT
      V=VMIN+FLOAT(I-1)*DELV
      ASSIGN 2001 TO LABEL
      GOTO 2000
2001  CONTINUE
103   CONTINUE
C------------------
      GOTO 3000

C=============================================================================
C     this section calculates the transmission etc at wavenumber V
C     The output values are stored in the array OUTPUT(,I)
C     This would be better as a subroutine, but there are so many parameters
C     to pass that it would take forever so entry is via a GOTO and exit
C     via assigned GOTO
2000  CONTINUE
C     first computing normal path optical depth through each layer (TAUTMP)

C     now compute absorption coefficient
C     first continuum contribution
      CURBIN=INT((V-VBIN(1))/WING)+1
      TAUTMP=CONTINK(1,LAYER1,CURBIN)
      VTMP=V-VBIN(CURBIN)
      DO 51 ISUM=2,IORDP1
      TAUTMP=TAUTMP+CONTINK(ISUM,LAYER1,CURBIN)*VTMP
      VTMP=VTMP*VTMP
51    CONTINUE


      DO 507 JBIN=CURBIN-1,CURBIN+1 

      DO 52 LINE=FSTLIN(JBIN),LSTLIN(JBIN)
C     compute absorption coefficient for normal incidence

        J1=IDLIN(LINE)
        FNH3=-1.0
        FH2=-1.0
        TAUTMP=TAUTMP+LINECONTRIB(IPROC(J1),
     1       IDGAS(J1),V,TCORDW(J1),TCORS1(J1),
     2       TCORS2,PRESS,TEMP,FRAC(J1),
     3       VLIN(LINE),SLIN(LINE),ELIN(LINE),ALIN(LINE),SBLIN(LINE),
     4       TDW(LINE),TDWS(LINE),LLQ(LINE),DOUBV(LINE),FNH3,FH2)
     
52    CONTINUE
507   CONTINUE
59    CONTINUE
C     now scaling for each path

      OUTPUT(I)=TAUTMP

      GOTO LABEL


3000  CONTINUE


      OPEN(14,FILE='OUTPUT.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
      WRITE(14)OUTPUT
      CLOSE(14)


C      print*,'OK to calc_kdist'

C     CALCULATE CUMULATIVE K-DISTRIBUTION FROM LBL SPECTRUM
C     Note SPECTRUM is passed in a common block

C     Set NFIL to -1 to stop CALC_FKDIST_WAVE trying to average over
C     a filter function
      NFIL=-1

C     Radtrans calculations are in wavenumber space
      IWAVE=0

      CALL CALC_FKDIST_WAVEC(IWAVE,VMIN,DELV,NPOINT,NFIL,VFIL,
     1 TFIL,G_ORD,DELG,K_G,MAXG,NG)

      RETURN
      END


  
