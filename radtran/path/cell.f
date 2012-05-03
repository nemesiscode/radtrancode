      SUBROUTINE CELL(TEXT)
C     $Id: cell.f,v 1.1.1.1 2000-08-17 09:26:56 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  CELL: reads cell calculation for path.f and calculates layers
C
C_KEYS:   RADTRAN,SUBR
C
C_DESCR:  
C
C_ARGS:   
C
C_FILES : unit 2 - the path file [.pat]
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:   26feb93 SBC Original version
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
C     note that laycom uses parameters defined in pathcom
C--------------------------------------------------------------
C     miscellaneous variables used in code
      INTEGER INCGAS,LOCID(MAXGAS),ISNGL,NSNGL,FSTCEL,LSTCEL
      REAL VMRIP(MAXGAS),CELLEN,CELDOP,CELTEM,CELPRE,LOPRES,LOTEMP,
     1HIPRES,HITEMP
      INTEGER NPHAS,MAXPAS
      PARAMETER (MAXPAS=100)
      REAL PHASE(MAXPAS),PHPRES(MAXPAS),PHTEMP(MAXPAS)
      CHARACTER*100 IPFILE
      INTEGER I,J,K,IDIP,ISOIP
C--------------------------------------------------------------
      CELLEN=-1.
      CELDOP=0.
C     keyword loop
10    READ(2,1,END=8)TEXT
1     FORMAT(A)
      CALL REMSP(TEXT)
      CALL UPCASE(TEXT)
      IF(TEXT(1:1).EQ.' ')THEN
        RETURN
       ELSE IF(TEXT(1:5).EQ.'GASES')THEN
        READ(TEXT(6:),*)INCGAS
        WRITE(*,53)INCGAS
53      FORMAT(' including,',I3,' cell gases')
        DO 51 I=1,INCGAS
        READ(2,*)IDIP,ISOIP,VMRIP(I)
        CALL ADDGAS(IDIP,ISOIP,LOCID(I))
51      CONTINUE
       ELSE IF(TEXT(1:3).EQ.'DOP')THEN
        READ(TEXT(4:),*)CELDOP
       ELSE IF(TEXT(1:6).EQ.'LENGTH')THEN
        READ(TEXT(7:),*)CELLEN
       ELSE IF(TEXT(1:4).EQ.'SNGL')THEN
        IF(CELLEN.LT.0.)THEN
          CALL WTEXT('no cell length defined')
          STOP
          END IF
        READ(TEXT(5:),*)NSNGL
        FSTCEL=NPATH+1
        LSTCEL=NPATH+NSNGL
        DO 52 ISNGL=1,NSNGL
        READ(2,*)CELPRE,CELTEM
        NPATH=NPATH+1
        ERRLIM(NPATH)=ERRDEF
        IMOD(NPATH)=5
        NLAYER=NLAYER+1
        PRESS(NLAYER)=CELPRE
        TEMP(NLAYER)=CELTEM
        DO 100 I=1,INCGAS
        PP(NLAYER,LOCID(I))=VMRIP(I)*PRESS(NLAYER)
        AMOUNT(NLAYER,LOCID(I))=7.3398E21*VMRIP(I)*
     1  PRESS(NLAYER)*CELLEN/TEMP(NLAYER)
100     CONTINUE
        NLAYIN(NPATH)=1
        LAYINC(1,NPATH)=NLAYER
        DOP(NLAYER)=CELDOP
        SCALE(1,NPATH)=1.
52      CONTINUE
C       now defining the cell calculation
        NCALC=NCALC+1
        ITYPE(NCALC)=128
        NINTP(NCALC)=2
        NREALP(NCALC)=0
        ICALD(1,NCALC)=FSTCEL
        ICALD(2,NCALC)=LSTCEL
C
       ELSE IF(TEXT(1:8).EQ.'PMR TWOP')THEN
        READ(2,*)LOPRES,LOTEMP
        READ(2,*)HIPRES,HITEMP
        NPATH=NPATH+1
        ERRLIM(NPATH)=ERRDEF
        IMOD(NPATH)=6
        NLAYER=NLAYER+2
        PRESS(NLAYER-1)=LOPRES
        TEMP(NLAYER-1)=LOTEMP
        PRESS(NLAYER)=HIPRES
        TEMP(NLAYER)=HITEMP
        DO 107 I=1,INCGAS
        PP(NLAYER-1,LOCID(I))=VMRIP(I)*PRESS(NLAYER-1)
        AMOUNT(NLAYER-1,LOCID(I))=7.3398E21*VMRIP(I)*
     1  PRESS(NLAYER-1)*CELLEN/TEMP(NLAYER-1)
        PP(NLAYER,LOCID(I))=VMRIP(I)*PRESS(NLAYER)
        AMOUNT(NLAYER,LOCID(I))=7.3398E21*VMRIP(I)*
     1  PRESS(NLAYER)*CELLEN/TEMP(NLAYER)
107     CONTINUE
        NLAYIN(NPATH)=2
        LAYINC(1,NPATH)=NLAYER-1
        LAYINC(2,NPATH)=NLAYER
        DOP(NLAYER)=CELDOP
        DOP(NLAYER-1)=CELDOP
        SCALE(1,NPATH)=1.
        SCALE(2,NPATH)=1.
C       now creating the wideband path
        NPATH=NPATH+1
        ERRLIM(NPATH)=ERRDEF
        IMOD(NPATH)=12
        NLAYIN(NPATH)=2
        LAYINC(1,NPATH)=NLAYER-1
        LAYINC(2,NPATH)=NLAYER
        SCALE(1,NPATH)=1.
        SCALE(2,NPATH)=1.
C       now defining the cell calculation
        NCALC=NCALC+1
        ITYPE(NCALC)=129
        NINTP(NCALC)=2
        NREALP(NCALC)=0
        ICALD(1,NCALC)=NPATH-1
        ICALD(2,NCALC)=NPATH
C
       ELSE IF(TEXT(1:8).EQ.'PMR FILE')THEN
C       first reading in the cycle details
        READ(2,1)IPFILE
        CALL REMSP(IPFILE)
        CALL LOCASE(IPFILE)
        OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
        NPHAS=1
102     READ(1,*,END=103)PHPRES(NPHAS),PHTEMP(NPHAS)
        NPHAS=NPHAS+1
        GOTO 102
103     CONTINUE
        CLOSE(UNIT=1)
        NPHAS=NPHAS-1
        FSTCEL=NPATH+1
        LSTCEL=NPATH+NPHAS
        DO 108 K=1,NPHAS
        PHASE(K)=FLOAT(K-1)*360./FLOAT(NPHAS)
        NPATH=NPATH+1
        ERRLIM(NPATH)=ERRDEF
        IMOD(NPATH)=5
        NLAYER=NLAYER+1
        PRESS(NLAYER)=PHPRES(K)
        TEMP(NLAYER)=PHTEMP(K)
        DO 109 I=1,INCGAS
        PP(NLAYER,LOCID(I))=VMRIP(I)*PRESS(NLAYER)
        AMOUNT(NLAYER,LOCID(I))=7.3398E21*VMRIP(I)*
     1  PRESS(NLAYER)*CELLEN/TEMP(NLAYER)
109     CONTINUE
        NLAYIN(NPATH)=1
        LAYINC(1,NPATH)=NLAYER
        DOP(NLAYER)=CELDOP
        SCALE(1,NPATH)=1.
108     CONTINUE
C       now defining the cell calculation
        NCALC=NCALC+1
        ITYPE(NCALC)=130
        NINTP(NCALC)=2
        NREALP(NCALC)=NPHAS
        ICALD(1,NCALC)=FSTCEL
        ICALD(2,NCALC)=LSTCEL
        DO 106 J=1,NPHAS
        RCALD(J,NCALC)=PHASE(J)
106     CONTINUE
       ELSE IF(TEXT(1:3).EQ.'SCR')THEN
        READ(2,*)CELPRE,CELTEM
        NPATH=NPATH+1
        ERRLIM(NPATH)=ERRDEF
        IMOD(NPATH)=13
        NLAYER=NLAYER+1
        PRESS(NLAYER)=CELPRE
        TEMP(NLAYER)=CELTEM
        DO 207 I=1,INCGAS
        PP(NLAYER,LOCID(I))=VMRIP(I)*PRESS(NLAYER)
        AMOUNT(NLAYER,LOCID(I))=7.3398E21*VMRIP(I)*
     1  PRESS(NLAYER)*CELLEN/TEMP(NLAYER)
207     CONTINUE
        NLAYIN(NPATH)=1
        LAYINC(1,NPATH)=NLAYER
        DOP(NLAYER)=CELDOP
        SCALE(1,NPATH)=1.
C       now creating the wideband path
        NPATH=NPATH+1
        ERRLIM(NPATH)=ERRDEF
        IMOD(NPATH)=14
        NLAYIN(NPATH)=1
        LAYINC(1,NPATH)=NLAYER
        SCALE(1,NPATH)=1.
C       now defining the cell calculation
        NCALC=NCALC+1
        ITYPE(NCALC)=131
        NINTP(NCALC)=2
        NREALP(NCALC)=0
        ICALD(1,NCALC)=NPATH-1
        ICALD(2,NCALC)=NPATH
       END IF
      GOTO 10
C
8     RETURN
      END
