      PROGRAM PATH
C     $Id: path.f,v 1.10 2011-06-23 09:14:32 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  BANDPATH: produces path driver files for use by genband, gencorrk
C
C_KEYS:   RADTRAN,PROG,VMS
C
C_DESCR:  Computes atmospheric absorber paths and outputs in a form
C         suitable for reading by RADTRAN software.
C         Reads commands from an input path file with extension .pat and
C         produces an output driver file of the same name but extension .drv
C         The .pat files are block structured, each block starting with a
C         keyword and ending with a blank line. Available keywords are
C         described in the code.
C
C         The program calculates gaseous LAYERS based upon an input
C         atmospheric model or gas cell details. It then calculates PATHS
C         through these layers, thus allowing multiple use of any individual
C         optical depth calculation for a layer. It also records CALCULATION
C         types which may correspond to a single path or to complex
C         combinations of paths. For example a single cell transmission or a
C         single atmospheric limb emission require only one path but a
C         weighting function consists of multiple atmospheric paths. Both
C         paths and calculations can be defined as combinations of previous
C         ones.
C
C         Within the code calculation types are defined by logical variables
C         but these correspond to a single integer code ITYP which is passed
C         via the driver file to the genlbl software. In general one bit in
C         ITYP corresponds to one of the logical variables where possible
C         to simplify evaluation of ITYP in the genlbl software. Calculation
C         codes defined are:
C           Atmospheric Codes have bit 7 (128) zero
C                                  bit 6  (64) 1 if limb path
C                                  bit 5  (32) not used
C                                  bit 4  (16) not used
C                                  bit 3   (8) 1 if curtis godson
C                                  bit 2   (4) 1 if weighting function
C                                  bit 1   (2) 1 if emission calculation
C                                  bit 0   (1) 1 if calculate 1-transmission
C
C           non atmosphere codes have bit 7 (128) 1
C           128-159 are reserved for cell types, defined in cell.f
C           160 is a combined cell and atmosphere
C
C
C         A path type IMOD is also defined for each path. These are defined
C         in genlbl.f
C
C_ARGS:   
C
C_FILES : unit 1 - atmospheric model, dust and cell input files
C         unit 2 - the path file [.pat]
C         unit 4 - the driver file [.drv] opened by wrlbld.f
C
C_CALLS:  PROMPT     prompts for user input
C         FILE       forces file extension
C         ASKYN      prompts user for logical Yes/No input
C         WTEXT      writes text to screen
C         REMSP      removes leading spaces from a character string
C         ADDGAS     adds a new gas to the gas arrays
C         RDMOD      reads in a model atmosphere
C         ATM        reads in atmospheric calculation details, calculates
C                    layers and paths.
C         SNGATM     reads in a single transmission path and treats it
C                    as an atmospheric path so that it can be combined with
C                    cells
C         CELL       reads in cell details, calculates layers and paths and
C                    combines cell calculations with any atmospheric paths
C                    defined so far
C         WRLBLD     writes out path and calculation details to driver file
C
C_BUGS:
C
C_HIST:   20feb93 SBC  Original version based upon limb.f
C	  16aug94 PGJI Updated for new dust.
C
C_END:
C--------------------------------------------------------------
C     Variables to hold calculated layers and the details of each paths and
C     calculation requested.
C     pathcom holds the variables used bu the genlbl software too
C     laycom holds variables used only by the path software
C     parameters are passed between routines mostly using common blocks
C     because of the extensive use of large arrays
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/laycom.f'
C     note that laycom uses parameters defined in pathcom
C--------------------------------------------------------------
C     miscellaneous variables used in main code
      REAL MAXDV
      INTEGER I,J,K,L,N,ICHECK
      CHARACTER*100 IPFILE,NMFILE
      CHARACTER*80 TEXT,TEXT1
      LOGICAL ASKYN,OK,DEF,BIT,REFL
C--------------------------------------------------------------
C
C     initialise flags
      INTERV=.FALSE.
      LINED=.FALSE.
      MODEL=.FALSE.
      LAYERS=.FALSE.
      CLRLAY=.TRUE.
      COMB=.TRUE.
      JFP = 0

      DO I=1,MAXPRO
       DO J=1,MAXDUST
        ICLOUD(J,I)=0
       ENDDO
      ENDDO

C
      CALL PROMPT('path file')
      READ(*,1)IPFILE
1     FORMAT(A)
      CALL FILE(IPFILE,IPFILE,'pat')
      OPEN(UNIT=2,FILE=IPFILE,STATUS='OLD')
      CALL FILE(IPFILE,DRFILE,'drv')
C
C
C     initialise layer counters and defaults
      NCALC=0
      NPATH=0
      NFILT=1
      NCONT=0
      ERRDEF=0.01
      NGAS=0
C
C     initialising layer variables not set for all layer types
      DO 40 I=1,MAXLAY
      BASEP(I)=0.
      BASEH(I)=0.
      DELH(I)=0.
      BASET(I)=0.
      DOP(I)=0.
      DO 45 J=1,LIMGAS
      AMOUNT(I,J)=0.
      PP(I,J)=0.
45    CONTINUE
40    CONTINUE
C
C     main keyword loop
2     READ(2,1,END=9)TEXT
      CALL REMSP(TEXT)
      TEXT1=TEXT
      CALL UPCASE(TEXT)
C     checking for negated keywords - note that I don't bother to check for
C     meaningless negations. eg noatm has the same effect as atm
      IF(TEXT(1:2).EQ.'NO')THEN
        DEF=.FALSE.
        TEXT(1:2)='  '
        CALL REMSP(TEXT)
       ELSE
        DEF=.TRUE.
        END IF 
C      
      IF(TEXT(1:1).EQ.' ')THEN
C       skipping blank lines
	GOTO 2
       ELSE IF(TEXT(1:8).EQ.'INTERVAL')THEN
	READ(2,*)VMIN,VMAX,DELV,FWHM
	NPOINT=NINT((VMAX-VMIN)/DELV)+1
	READ(2,*)ICONV,WING,VREL
	INTERV=.TRUE.
        WRITE(*,*)' ICONV = ',ICONV
        WRITE(*,42)VMIN,VMAX,DELV
        WRITE(*,41)FWHM,WING,VREL
42      FORMAT(' VMIN = ',F11.3,' VMAX = ',F11.3,' DELV = ',F11.3)
41      FORMAT(' FWHM = ',F11.3,' WING = ',F11.3,' VREL = ',F11.3)
        IF (ICONV.EQ.1) THEN
           WRITE(*,*)'*** Warning (path): Radtrans will ignore FWHM.'
           WRITE(*,*)'      Radtrans will assume FWHM = DELV to '
           WRITE(*,*)'      avoid calculating bits of the'
           WRITE(*,*)'      spectrum more than once'
           END IF
        IF(ICONV.EQ.0)CALL WTEXT('"infinite resolution"')
       ELSE IF(TEXT(1:9).EQ.'SPEC DATA')THEN
	 READ(2,1)LINKEY
         WRITE(*,43)LINKEY
43       FORMAT(' reading spectral data from:',A)
	 LINED=.TRUE.
       ELSE IF(TEXT(1:7).EQ.'PROCESS')THEN
        READ(TEXT(8:),*)I,J,K
        CALL ADDGAS(I,J,L)
        IPROC(L)=K
       ELSE IF(TEXT(1:5).EQ.'MODEL')THEN
        CALL RDMOD(TEXT1(6:))
       ELSE IF(TEXT(1:10).EQ.'DUST MODEL')THEN
        CALL RDDMOD(TEXT1(11:))
       ELSE IF(TEXT(1:13).EQ.'FPARAH2 MODEL')THEN
        CALL RFPMOD(TEXT1(14:))
       ELSE IF(TEXT(1:13).EQ.'FCLOUD MODEL')THEN
        CALL RFCMOD(TEXT1(14:))
       ELSE IF(TEXT(1:12).EQ.'DUST SPECTRA')THEN
        XFILE(1:)=TEXT1(14:)
       ELSE IF(TEXT(1:5).EQ.'LAYER')THEN
        CALL LAYER(TEXT(6:))
       ELSE IF(TEXT(1:6).EQ.'NLAYER')THEN
        CALL NEWLAYER(TEXT(7:))
       ELSE IF(TEXT(1:6).EQ.'ULAYER')THEN
        CALL URANLAYER(TEXT(7:))
       ELSE IF(TEXT(1:9).EQ.'MUSELAYER')THEN
        CALL MUSELAYER(TEXT(10:))
       ELSE IF(TEXT(1:6).EQ.'CLAYER')THEN
        CALL CLOUDLAYER(TEXT(7:))
       ELSE IF(TEXT(1:9).EQ.'COMPLAYER')THEN
        CALL COMPLAYER(TEXT(10:))
       ELSE IF(TEXT(1:6).EQ.'MLAYER')THEN
        CALL MEWLAYER(TEXT(7:))
       ELSE IF(TEXT(1:3).EQ.'ATM')THEN
        CALL ATM(TEXT(4:))
       ELSE IF(TEXT(1:7).EQ.'REFLATM')THEN
        CALL REFLATM(TEXT(8:))
       ELSE IF(TEXT(1:6).EQ.'SNGATM')THEN
        CALL SNGATM(TEXT(4:))
       ELSE IF(TEXT(1:4).EQ.'CELL')THEN
        CALL CELL(TEXT(4:))
       ELSE IF(TEXT(1:5).EQ.'ERROR')THEN
        READ(TEXT(6:),*)ERRDEF
        ERRDEF=ERRDEF/100.
       ELSE IF(TEXT(1:6).EQ.'CLRLAY')THEN
        CLRLAY=DEF
       ELSE IF(TEXT(1:7).EQ.'COMBINE')THEN
        COMB=DEF
       ELSE
        WRITE(*,44)TEXT(1:40)
44	FORMAT(' unrecognised main keyword:',A)
	STOP
	END IF
      GOTO 2
C
9     CONTINUE
C
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
      END IF


C     Defining Reflecting Layer Calculations
C     First do combined atmosphere cell paths
      REFL=.FALSE.
      DO 433 I=1,NCALC
       IF(ITYPE(I).EQ.161)REFL=.TRUE.
433   CONTINUE
      IF(REFL)THEN
       NCALC=NCALC+1
       ITYPE(NCALC)=201
       KR=0
       KE=0
       DO 434 I=1,NCALC-1
        IF(ITYPE(I).EQ.160)KE=KE+1
        IF(ITYPE(I).EQ.161)KR=KR+1
434    CONTINUE
       NINTP(NCALC)=2+KE+KR
       ICALD(1,NCALC)=KE
       ICALD(2,NCALC)=KR
       KE1=3
       KR1=KE1+KE
       DO 435 I=1,NCALC-1
        IF(ITYPE(I).EQ.160)THEN
         ICALD(KE1,NCALC)=I
         KE1=KE1+1
        END IF
        IF(ITYPE(I).EQ.161)THEN
         ICALD(KR1,NCALC)=I
         KR1=KR1+1
        END IF
435    CONTINUE
      END IF



      REFL=.FALSE.
      DO 453 I=1,NCALC
       IF((ITYPE(I).LT.160).AND.(BIT(5,ITYPE(I))))REFL=.TRUE.
453   CONTINUE
      IF(REFL)THEN
       NCALC=NCALC+1
       ITYPE(NCALC)=200
       KR=0
       KE=0
       DO 454 I=1,NCALC-1
        IF(ITYPE(I).LT.128)THEN
         IF(BIT(5,ITYPE(I)))THEN
          KR=KR+1
         ELSE
          KE=KE+1
         END IF
        END IF
454    CONTINUE
       NINTP(NCALC)=2+KE+KR
       ICALD(1,NCALC)=KE
       ICALD(2,NCALC)=KR
       KE1=3
       KR1=KE1+KE
       DO 455 I=1,NCALC-1
        IF(ITYPE(I).LT.128)THEN
        IF(BIT(5,ITYPE(I)))THEN
         ICALD(KR1,NCALC)=I
         KR1=KR1+1
        ELSE
         ICALD(KE1,NCALC)=I
         KE1=KE1+1
        END IF
        END IF
455    CONTINUE
      END IF

      DO 224 IPATH=1,NPATH
	ERRLIM(IPATH)=ERRDEF
224   CONTINUE
C
      OK=.TRUE.
      IF(.NOT.LINED)THEN
	CALL WTEXT('no line data key defined')
	OK=.FALSE.
      END IF
      IF(.NOT.INTERV)THEN
	CALL WTEXT('no wavenumber interval defined')
	OK=.FALSE.
      END IF
      IF(OK)THEN
	CALL WRLBLD
C       creating a new file containing just the name of the driver file
C       this is useful for running lbl in background under Unix
        CALL FILE(IPFILE,NMFILE,'nam')
        OPEN(UNIT=1,FILE=NMFILE,STATUS='UNKNOWN')
         WRITE(1,3)DRFILE
3       FORMAT(1X,A)
        CLOSE(UNIT=1)
      ELSE
	CALL WTEXT('no driver file written')
      END IF

      STOP

      END
