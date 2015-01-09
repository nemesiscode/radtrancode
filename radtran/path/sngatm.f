      SUBROUTINE SNGATM(TEXT)
C     $Id: sngatm.f,v 1.2 2002-07-12 09:05:08 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  SNGATM: reads in single atmospheric transmission paths
C
C_KEYS:   RADTRAN,SUBR
C
C_DESCR:  
C
C_ARGS:   
C
C_FILES : unit 2 - the path file [.pat]
C         unit 3 - the (optional) file containing the paths
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:   16jul93 SBC Original version
C
C_END:
C--------------------------------------------------------------
      CHARACTER TEXT*(*)
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
      INCLUDE '../includes/constdef.f'
C     note that laycom uses parameters defined in pathcom
C--------------------------------------------------------------
C     miscellaneous variables used in code
      LOGICAL DEF
      INTEGER LUNIT,IDIP,ISOIP,LOCID,INCGAS
      REAL AMIP,VMRIP,Q1,MODBOLTZA,AMIP1
      PARAMETER (MODBOLTZA=KBOLTZMANN/1.013)
C     MODBOLTZA = KBOLTZ/1.013 (where KBOLTZ = 1.381E-23) and multiplied
C     by 10. Then AMOUNT*MODBOLTZA*T(K)/P(ATM) gives path length in cm.

      CHARACTER*100 IPFILE
C--------------------------------------------------------------
C     note that TEXT(1:6) must contain "SNGATM"
      TEXT(1:6)='      '
      CALL REMSP(TEXT)
      IF(TEXT(1:1).EQ.' ')THEN
C       reading data directly from path file
        LUNIT=2
       ELSE
C       reading data from separate file
        OPEN(UNIT=3,FILE=TEXT,STATUS='OLD')
        LUNIT=3
        END IF
10    CONTINUE
      READ(UNIT=LUNIT,FMT=12,END=11,ERR=11)TEXT
12    FORMAT(A)
      CALL REMSP(TEXT)
      CALL UPCASE(TEXT)
      ABSORB = .TRUE.
      IF(TEXT(1:1).EQ.' ')GOTO 11
      IF(TEXT(1:2).EQ.'NO')THEN
        TEXT(1:2)='  '
        CALL REMSP(TEXT)
        DEF=.FALSE.
       ELSE
        DEF=.TRUE.
        END IF
      IF(TEXT(1:6).EQ.'ABSORB')ABSORB=DEF
      IF(TEXT(1:6).EQ.'THERM')THERM=DEF
      NPATH=NPATH+1
      ERRLIM(NPATH)=ERRDEF
      IMOD(NPATH)=0
      IF(ABSORB)IMOD(NPATH)=1
      IF(THERM)IMOD(NPATH)=11
      NLAYER=NLAYER+1
      READ(LUNIT,*)INCGAS,PRESS(NLAYER),TEMP(NLAYER)
      WRITE(*,101)PRESS(NLAYER),TEMP(NLAYER),INCGAS
101   FORMAT(' reading single atmospheric path. P=',E12.5,', T=',F6.1,I3,
     1'gases')
      AMTOT=0.
      DO 100 I=1,INCGAS
      READ(LUNIT,*)IDIP,ISOIP,VMRIP,AMIP
      CALL ADDGAS(IDIP,ISOIP,LOCID)
      PP(NLAYER,LOCID)=VMRIP*PRESS(NLAYER)
      IF(AMIP.LT.1.E10)THEN
C       note assuming that if amount <1.e10 then value is actual path length
C       in Km
        AMOUNT(NLAYER,LOCID)=1.E5*7.3398E21*VMRIP*
     1  PRESS(NLAYER)*AMIP/TEMP(NLAYER)
	DELH(NLAYER)=AMIP
       ELSE
        AMOUNT(NLAYER,LOCID)=AMIP
        AMIP1 = AMIP*MODBOLTZA*TEMP(NLAYER)/PRESS(NLAYER)
	DELH(NLAYER) = AMIP1/(VMRIP*1e5)		! Path length in km
        IF(DELH(NLAYER).GT.99999.9) THEN
         print*,'Layer Height too large in sngatm : ',DELH(NLAYER)
         print*,'Capping to 99999.9'
         DELH(NLAYER)=99999.9
        ENDIF
      END IF
100   CONTINUE
      READ(LUNIT,*)NCONT
      IF(NCONT.GT.0)THEN
       READ(LUNIT,12)XFILE
       DO 230 I=1,NCONT
        READ(LUNIT,*)CONT(I,NLAYER)
230    CONTINUE
      END IF
      NLAYIN(NPATH)=1
      LAYINC(1,NPATH)=NLAYER
      EMTEMP(1,NPATH)=TEMP(NLAYER)
      DOP(NLAYER)=0.
      SCALE(1,NPATH)=1.
52    CONTINUE
C     now defining the calculation
      NCALC=NCALC+1
      ITYPE(NCALC)=16
      IF(ABSORB)ITYPE(NCALC)=17
      IF(THERM)ITYPE(NCALC)=10
      NINTP(NCALC)=2
      NREALP(NCALC)=0
      ICALD(1,NCALC)=NPATH
      ICALD(2,NCALC)=NPATH
      GOTO 10
C
11    IF(LUNIT.EQ.3)CLOSE(3)
      RETURN
      END
