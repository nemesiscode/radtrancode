      SUBROUTINE READDUSTOD(OPFILE1,NCONT1,XOD)
C     $Id:
C-----------------------------------------------------------------------------
C_TITLE:  RDLBLD: reads in arrays for use by GENLBL or GENLBL output
C         processors
C_ARGS:
C
C_KEYS:   
C
C_DESCR:  reads in arrays as output by WRLBLD
C
C_FILES:  UNIT 1: driver file, already open when called
C
C_CALLS:
C
C_BUGS:
C
C_HIST:  6feb87 SBC ORIGINAL VERSION
C        20sep91 SBC modified to use new line data routines
C        18nov92 SBC modified to use EMTEMP and BASET
C          1dec92 SBC new format for layer data etc. see WRLDLD
C----------------------------------------------------------------------------
C     note dbcom defines the linedata base variables. it is not normally stored
C     in the same directory as the rest of the code
      INCLUDE '../includes/arrdef.f' 
      INCLUDE '../includes/dbcom.f' 
C---------------------------------------------------------------------------
      INCLUDE '../includes/pathcom.f'
C-----------------------------------------------------------------------------
      INTEGER I,J,NCONT1
      CHARACTER*100 DXD1,DXD2,OPFILE1
      CHARACTER*100 KEYFILTMP
      REAL XOD(MAXCON)
      CALL FILE(OPFILE1,OPFILE1,'drv')
      OPEN(1,FILE=OPFILE1,STATUS='OLD')
      READ(1,*)
      READ(1,*)VMIN,DELV,NPOINT,FWHM
      READ(1,*)WING,VREL
      READ(1,546)KEYFIL
546   FORMAT(A)
334   FORMAT(1X,A30)
      READ(1,*)ICONV,FLAGH2P,NCONT,FLAGC
      IF(NCONT.GT.0)READ(1,334)XSCFIL
      READ(1,*)NLAYER,NPATH,NGAS
      DO 504 I=1,NGAS
      READ(1,*)IDGAS(I)
      READ(1,*)ISOGAS(I),IPROC(I)
504   CONTINUE
C     inputting details of each layer
      INLTE=0
      IEXTRA=0
      READ(1,546)DXD1
      IF(DXD1(1:4).EQ.'NLTE')THEN
       DXD2(1:30)=DXD1(5:34)
       READ(DXD2,*)INLTE
       IEXTRA=IEXTRA+1
      END IF
      DO 508 I=1,IEXTRA
       READ(1,*)
508   CONTINUE
      READ(1,*)
      READ(1,*)
      READ(1,*)
      DO 505 I=1,NLAYER
      READ(1,*)J,BASEH(I),DELH(I),BASEP(I),BASET(I),TOTAM(I),
     1PRESS(I),TEMP(I),DOP(I)
      READ(1,*)(AMOUNT(I,J),PP(I,J),J=1,NGAS)
C     ********************************************
      IF(NCONT.GT.0)READ(1,*)(CONT(J,I),J=1,NCONT)
      IF(FLAGH2P.EQ.1)READ(1,*)HFP(I)
      IF(FLAGC.EQ.1)READ(1,*)HFC(I),(IFC(J,I),J=1,NCONT)
C     ********************************************
505   CONTINUE
      DO 520 I=1,NPATH
      READ(1,*)NLAYIN(I),IMOD(I),ERRLIM(I)
      DO 522 J=1,NLAYIN(I)
      READ(1,*)K,LAYINC(J,I),EMTEMP(J,I),SCALE(J,I)
522   CONTINUE
520   CONTINUE
      READ(1,*)NFILT
      DO 526 I=1,NFILT
      READ(1,*)FILTER(I),VFILT(I)
526   CONTINUE
      READ(1,540)OPFILE
540   FORMAT(1X,A)
      READ(1,*)NCALC
      DO 541 I=1,NCALC
      READ(1,*)ITYPE(I),NINTP(I),NREALP(I),NCHP(I)
      DO 542 J=1,NINTP(I)
      READ(1,*)ICALD(J,I)
542   CONTINUE
      DO 543 J=1,NREALP(I)
      READ(1,*)RCALD(J,I)
543   CONTINUE
      DO 544 J=1,NCHP(I)
      READ(1,545)CCALD(J,I)
545   FORMAT(1X,1A30)
544   CONTINUE
541   CONTINUE

      
      DO J=1,MAXCON
       XOD(J)=0.0
      ENDDO

      DO 506 I=1,NLAYER
       DO 667 J=1,NCONT
        XOD(J)=XOD(J)+CONT(J,I)
667    ENDDO
506   CONTINUE

      NCONT1=NCONT

      RETURN
      END
