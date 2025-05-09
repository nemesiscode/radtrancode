C***********************************************************************
C***********************************************************************
C_TITL:	SUBPATHG.f
C
C_DESC:	Subroutine version of PATH: produces path driver files
C	for use by Radtran (or its associated derivatives).
C
C	Computes atmospheric absorber paths and outputs in a form
C	suitable for reading by RADTRAN software.
C	Reads commands from an input path file with extension .pat and
C	produces an output driver file of the same name but extension .drv
C	The .pat files are block structured, each block starting with a
C	keyword and ending with a blank line. Available keywords are
C	described in the code.
C
C	The program calculates gaseous LAYERS based upon an input
C	atmospheric model or gas cell details. It then calculates PATHS
C	through these layers, thus allowing multiple use of any individual
C	optical depth calculation for a layer. It also records CALCULATION
C	types which may correspond to a single path or to complex
C	combinations of paths. For example a single cell transmission or a
C	single atmospheric limb emission require only one path but a
C	weighting function consists of multiple atmospheric paths. Both
C	paths and calculations can be defined as combinations of previous
C	ones.
C
C	Within the code calculation types are defined by logical variables
C	but these correspond to a single integer code ITYP which is passed
C	via the driver file to the genlbl software. In general one bit in
C	ITYP corresponds to one of the logical variables where possible
C	to simplify evaluation of ITYP in the genlbl software. Calculation
C	codes defined are:
C	* Atmospheric Codes have ...
C			bit 7 (128) zero
C			bit 6  (64) 1 if limb path
C			bit 5  (32) not used
C			bit 4  (16) not used
C			bit 3   (8) 1 if curtis godson
C			bit 2   (4) 1 if weighting function
C			bit 1   (2) 1 if emission calculation
C			bit 0   (1) 1 if calculate 1-transmission
C
C	Non atmosphere codes have ...
C			bit 7 (128) 1 128-159 are reserved for cell types,
C				defined in cell.f. 160 is a combined cell
C				and atmosphere
C
C	A path type IMOD is also defined for each path. These are defined
C	in genlbl.f
C
C_ARGS: Input variables:
C	runname		CHARACTER*100	Operational filename.
C
C_FILE:	unit=2		The path (.pat) file.
C
C_CALL: file		Forces file extension.
C	remsp		Removes leading spaces from a character string.
C	upcase		Capitolises text string.
C	addgas		Adds a new gas to the gas arrays.
C	rdmod		Reads in a model atmosphere.
C	rddmod		Reads in a dust/cloud model atmosphere.
C	rfpmod		Reads in a H2 ortho/para fraction profile.
C	locase		Forces the text string to lower case.
C	layerg		Calculates atmospheric layers.
C	newlayer	Calculates atmospheric layers.
C	muselayer	Calculates atmospheric layers.
C	mewlayer	Calculates atmospheric layers.
C	cloudlayer	Calculates atmospheric layers.
C	complayer	Calculates atmospheric layers.
C	atmg		Reads in atmospheric calculation details,
C			calculates layers and paths.
C	reflatm		Reads atmospheric calculation and calculates
C			layers.
C	sngatm		Reads in a single transmission path and treats it
C			as an atmospheric path so that it can be combined
C			with cells.
C	cell		Reads in cell details, calculates layers and
C			paths and combines cell calculations with any
C			atmospheric paths defined so far.
C	wrlbld		Writes out path and calculation details to driver
C			file.
C
C_HIST:	20feb93	SBC	Original version based upon limb.f
C	16aug94	PGJI	Updated for new dust.
C	29jul96	PGJI	Modified to run as subroutine
C	04sep01 PGJI	Upgraded to implicit gradient calculation
C       29feb12 PGJI	Updated for Radtrans2.0
C
C***********************************************************************

      SUBROUTINE SUBPATHG(runname)

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C ../includes/pathcom.f holds the variables used by the software when
C calculating atmospheric paths (e.g. NPATH, IMOD and LINKEY).
      INCLUDE '../includes/laycom.f'
C ../includes/laycom.f holds variables used only by the path software
C parameters are passed between routines mostly using common blocks
C because of the extensive use of large arrays. NOTE: laycom uses
C parameters defined in pathcom.
      INCLUDE '../includes/laygrad.f'
C ../includes/laygrad.f holds the variables for use in gradient
C calculation.

C Miscellaneous variables used in main code.
      INTEGER i,j,k,l,n,npath1
C NPATH1: = NPATH + 1 for combined cell and atmosphere calculations.
      INTEGER kr,kr1,ke,ke1

      CHARACTER*100 runname,text,buffer
      CHARACTER*100 proffile,aerofile,fparfile,fcfile
C PROFFILE: Atmospheric profile filename.
C AEROFILE: Dust/Aerosol profile filename.

      LOGICAL OK,DEF,BIT,REFL
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C********************************* CODE ********************************

C Initialise flags.
      INTERV = .FALSE.
      LINED = .FALSE.
      MODEL = .FALSE.
      LAYERS = .FALSE.
      CLRLAY = .TRUE.
      COMB = .TRUE.
      JFP = 0

      DO I=1,MAXPRO
       DO J=1,MAXCON
        ICLOUD(J,I)=0
       ENDDO
      ENDDO


      CALL FILE(runname,runname,'pat')
      OPEN(UNIT=2,FILE=runname,STATUS='OLD')

      CALL FILE(runname,drfile,'drv')

C Initialise layer counters and defaults.
      NLAYER = 0
      NCALC = 0
      NPATH = 0
      NFILT = 1
      NCONT = 0
      ERRDEF = 0.01
      NGAS = 0

C Initialising layer variables not set for all layer types.
      DO 40 I=1,maxlay
        BASEP(I) = 0.0
        BASEH(I) = 0.0
        DELH(I) = 0.0
        BASET(I) = 0.0
        DOP(I) = 0.0
        DO 45 J=1,MAXGAS
          AMOUNT(I,J) = 0.0
          PP(I,J) = 0.0
45      CONTINUE
40    CONTINUE


C Main keyword loop.
2     READ(2,1,END=9)text
1     FORMAT(A)
      CALL remsp(text)
      CALL upcase(text)

C Checking for negated keywords. NOTE: the code doesn't bother to check
C for meaningless negations. eg noatm has the same effect as atm.
      IF (text(1:2).EQ.'NO') THEN
        def = .FALSE.
        text(1:2) = '  '
        CALL remsp(text)
      ELSE
        def = .TRUE.
      ENDIF

      IF (text(1:1).EQ.' ') THEN
C Skipping blank lines.
        GOTO 2
      ELSE IF (text(1:8).EQ.'INTERVAL') THEN
        READ(2,*)vmin,vmax,delv,fwhm
        if(delv.gt.0.0)then
         NPOINT = NINT((vmax - vmin)/delv) + 1
        else
         NPOINT = 1
        endif
        READ(2,*)iconv,wing,vrel
        interv = .TRUE.
        if(idiag.gt.0)then
         WRITE(*,*)' '
         WRITE(*,*)' SUBPATHG.f :: iconv = ', iconv
         WRITE(*,*)' SUBPATHG.f :: vmin = ',vmin,', vmax = ',vmax
         WRITE(*,*)' SUBPATHG.f :: delv = ',delv,', fwhm = ',fwhm
         WRITE(*,*)' SUBPATHG.f :: wing = ',wing,', vrel = ',vrel
         IF (iconv.EQ.0) WRITE(*,*)' iconv = 0 ==> infinite resolution.'
        endif
      ELSE IF (text(1:9).EQ.'SPEC DATA') THEN
        linkey(1:) = text(10:)
        if(idiag.gt.0)then
         WRITE(*,*)' SUBPATHG.f :: reading spectral data from: ',linkey
        endif
        lined = .TRUE.
      ELSE IF (text(1:7).EQ.'PROCESS') THEN
        READ(text(8:),*)i,j,k
        CALL ADDGAS(i,j,l)
        iproc(l) = k
      ELSE IF (text(1:5).EQ.'MODEL') THEN
        proffile(1:94) = text(6:100)
        proffile(95:100) = '     '            ! must be this way for IFC
        CALL RDMOD(proffile)
      ELSE IF (TEXT(1:10).EQ.'DUST MODEL') THEN
        aerofile(1:89) = text(11:100)
        aerofile(90:100) = '          '       ! must be this way for IFC
        CALL RDDMOD(aerofile)
      ELSE IF (TEXT(1:13).EQ.'FPARAH2 MODEL')THEN
        fparfile(1:87) = text(14:100)
        fparfile(88:100) = '             '    ! must be this way for IFC
        CALL RFPMOD(fparfile)   
      ELSE IF (TEXT(1:12).EQ.'FCLOUD MODEL')THEN
        fcfile(1:88) = text(13:100)
        fcfile(89:100) = '            '    ! must be this way for IFC
        CALL RFCMOD(fcfile)   
      ELSE IF (TEXT(1:12).EQ.'DUST SPECTRA') THEN
        XFILE(1:) = TEXT(14:)
        CALL LOCASE(XFILE)
      ELSE IF (TEXT(1:5).EQ.'LAYER') THEN
        CALL LAYERG(TEXT(6:))
      ELSE IF (TEXT(1:6).EQ.'NLAYER') THEN
        CALL NEWLAYER(TEXT(7:))
      ELSE IF (TEXT(1:6).EQ.'ULAYER') THEN
        CALL URANLAYER(TEXT(7:))
      ELSE IF (TEXT(1:9).EQ.'MUSELAYER') THEN
        CALL MUSELAYER(TEXT(10:))
      ELSE IF (TEXT(1:6).EQ.'MLAYER') THEN
        CALL MEWLAYER(TEXT(7:))
      ELSE IF (TEXT(1:6).EQ.'CLAYER') THEN
        CALL CLOUDLAYER(TEXT(7:))
      ELSE IF(TEXT(1:9).EQ.'COMPLAYER')THEN
        CALL COMPLAYER(TEXT(10:))
      ELSE IF (TEXT(1:3).EQ.'ATM') THEN
        CALL ATMG(TEXT(4:))
      ELSE IF (TEXT(1:7).EQ.'REFLATM') THEN
        CALL REFLATM(TEXT(8:))
      ELSE IF (TEXT(1:6).EQ.'SNGATM') THEN
        CALL SNGATM(TEXT(4:))
      ELSE IF (TEXT(1:4).EQ.'CELL') THEN
        CALL CELL(TEXT(4:))
      ELSE IF (TEXT(1:5).EQ.'ERROR') THEN
        READ (TEXT(6:),*) ERRDEF
        ERRDEF = ERRDEF/100.
      ELSE IF (TEXT(1:6).EQ.'CLRLAY') THEN
        CLRLAY = DEF
      ELSE IF (TEXT(1:7).EQ.'COMBINE') THEN
        COMB = DEF
      ELSE
        WRITE(*,*)' SUBPATHG.f :: unrecognised main keyword.'
        WRITE(*,*)' Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)' Unknown keyword: ',text(1:40)
        STOP
      ENDIF
      GOTO 2

9     CONTINUE

      IF(COMB)THEN

C       combining cell and atmosphere paths
C       checking each calculation to see if it's a cell calculation
        N=NCALC
        DO 200 I=1,N
        IF(BIT(7,ITYPE(I)))THEN
          DO 202 J=1,N
          IF(.NOT.BIT(7,ITYPE(J)))THEN
C           combining cell calc I with atmosphere calc J
C           for each path in the cell calculation, combining with all
C           atmosphere paths
            NPATH1=NPATH+1
            DO 203 K=ICALD(1,I),ICALD(2,I)
            DO 201 L=ICALD(1,J),ICALD(2,J)
            NPATH=NPATH+1
            ERRLIM(NPATH)=ERRDEF
            IMOD(NPATH)=8
            LAYINC(1,NPATH)=K
            LAYINC(2,NPATH)=L
            SCALE(1,NPATH)=1.
            SCALE(2,NPATH)=1.
            NLAYIN(NPATH)=2
201         CONTINUE
203         CONTINUE
C           now adding the calculation type
            NCALC=NCALC+1
            ITYPE(NCALC)=160
            IF(BIT(5,ITYPE(J)))ITYPE(NCALC)=161
C           note that these calculation types also have bit 7 set
            NINTP(NCALC)=4
            NREALP(NCALC)=0
            ICALD(1,NCALC)=NPATH1
            ICALD(2,NCALC)=NPATH
            ICALD(3,NCALC)=I
            ICALD(4,NCALC)=J
            END IF
202       CONTINUE
          END IF
200     CONTINUE


      ENDIF


C Defining Reflecting Layer Calculations: first do combined atmosphere
C cell paths.
      REFL = .FALSE.
      DO 433 I=1,NCALC
        IF (ITYPE(I).EQ.161) REFL = .TRUE.
433   CONTINUE
      IF(REFL)THEN
        NCALC = NCALC + 1
        ITYPE(NCALC) = 201
        KR = 0
        KE = 0
        DO 434 I=1,NCALC-1
          IF (ITYPE(I).EQ.160) KE = KE + 1
          IF (ITYPE(I).EQ.161) KR = KR + 1
434     CONTINUE
        NINTP(NCALC) = 2 + KE + KR
        ICALD(1,NCALC) = KE
        ICALD(2,NCALC) = KR
        KE1 = 3
        KR1 = KE1 + KE
        DO 435 I=1,NCALC-1
          IF (ITYPE(I).EQ.160) THEN
            ICALD(KE1,NCALC) = I
            KE1 = KE1 + 1
          ENDIF
          IF (ITYPE(I).EQ.161) THEN
            ICALD(KR1,NCALC) = I
            KR1 = KR1 + 1
          ENDIF
435     CONTINUE
      ENDIF


      REFL = .FALSE.
      DO 453 I=1,NCALC
        IF (BIT(5,ITYPE(I))) REFL=.TRUE.
453   CONTINUE
      IF (REFL) THEN
        NCALC = NCALC + 1
        ITYPE(NCALC) = 200
        KR = 0
        KE = 0
        DO 454 I=1,NCALC-1
          IF (ITYPE(I).LT.128) THEN
            IF (BIT(5,ITYPE(I))) THEN
              KR = KR + 1
            ELSE
              KE = KE + 1
            ENDIF
          ENDIF
454     CONTINUE
        NINTP(NCALC) = 2 + KE + KR
        ICALD(1,NCALC) = KE
        ICALD(2,NCALC) = KR
        KE1 = 3
        KR1 = KE1+KE
        DO 455 I=1,NCALC-1
          IF (ITYPE(I).LT.128) THEN
            IF (BIT(5,ITYPE(I))) THEN
              ICALD(KR1,NCALC) = I
              KR1 = KR1 + 1
            ELSE
              ICALD(KE1,NCALC) = I
              KE1 = KE1 + 1
            ENDIF
          ENDIF
455     CONTINUE
      ENDIF

      OK = .TRUE.
      IF (.NOT.LINED) THEN
        WRITE(*,*)' SUBPATHG.f :: Error: no line data key defined.'
        OK = .FALSE.
      ENDIF
      IF (.NOT.INTERV) THEN
        WRITE(*,*)' SUBPATHG.f :: Error: no wavenum-interval defined.'
        OK = .FALSE.
      ENDIF
      IF (OK) THEN
C       Write driver file (.drv)
	if(iquiet.eq.0)CALL WRLBLD
      ELSE
	WRITE(*,*)' SUBPATHG.f :: no driver file (.drv) written.'
        STOP
      ENDIF

      CLOSE(UNIT=2)

      RETURN

      END
