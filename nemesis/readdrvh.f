      SUBROUTINE READDRVH(RUNNAME,HEIGHT)
C     $Id: readdrvh.f,v 1.2 2007/07/10 09:59:47 irwin Exp $
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
      INCLUDE '../radtran/includes/dbcom.f' 
C---------------------------------------------------------------------------
      INCLUDE '../radtran/includes/arrdef.f' 
      INCLUDE '../radtran/includes/pathcom.f'
C-----------------------------------------------------------------------------
      INTEGER I,J
      CHARACTER*100 DXD1,DXD2,RUNNAME
      CHARACTER*56 KEYFILTMP
      REAL HEIGHT(100)
      CALL FILE(RUNNAME,RUNNAME,'drv')

      OPEN(1,FILE=RUNNAME,STATUS='OLD')
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
 
      CLOSE(1)

      DO 506 I=1,NLAYER
       HEIGHT(I)=BASEH(I)
506   CONTINUE
  
      RETURN

      END
