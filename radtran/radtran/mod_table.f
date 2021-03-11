      PROGRAM MOD_TABLE
C     $Id: cut_table.f,v 1.1 2011-06-17 15:25:43 irwin Exp $
C***********************************************************************
C_TITL:	CUT_TABLE.f
C
C_DESC: Allows user to changes gas ID of a k-table and also DELV and FWHM
C       NB. For a channel k-table do not change DELV to be > 0.
C       NB. For a continuous k-table do not change DELV to be <=0
C
C_ARGS:	See the definitions below.
C
C_FILE:	UNIT=LUN	output file
C	unit=lun0
C
C_CALL:	file
C	remsp
C	pindex	
C
C_HIST:	10feb95	PGJI	ORIGINAL VERSION
C	13jan11 PGJI	Modified to cut data
C***************************** VARIABLES *******************************

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
      INCLUDE '../includes/dbcom.f'
C NOTE: dbcom defines the linedata base variables. it is not normally
C stored in the same directory as the rest of the code
      INCLUDE '../includes/pathcom.f'

      INTEGER LUN,LUN0,LUN1      
      PARAMETER (LUN=2,LUN0=30,LUN1=31)
C MAXOUT the maximum number of output points

      INTEGER PINDEX,CP,CT,LOOP,I1,I2

      INTEGER IREC,IREC0
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK),VCEN(MAXBIN),X1,X2
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),KVAL

      CHARACTER*100 KTAFIL,OPFILE1,OPFILE2,OUTFIL
      CHARACTER*1 ANS
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet
C******************************** CODE *********************************
      idiag=1

      CALL PROMPT('Enter input filename : ')
      READ(5,23)OPFILE1
23    FORMAT(A)
      CALL FILE(OPFILE1,KTAFIL,'kta')


C     Assume 4-byte words per record
      NW=1

      IRECL = NW*ISYS()
      WRITE(*,*)'irecl = ',irecl

      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL)


      READ(LUN0,REC=9)IDGAS(1)
      READ(LUN0,REC=10)ISOGAS(1)


      print*,'Current gas ID,ISO = ',IDGAS(1),ISOGAS(1)
      print*,'Enter new values : '
      read*, I1,I2
      IDGAS(1)=I1
      ISOGAS(1)=I2

      WRITE(LUN0,REC=9)IDGAS(1)
      WRITE(LUN0,REC=10)ISOGAS(1)

      READ(LUN0,REC=4)DELV
      READ(LUN0,REC=5)FWHM

      print*,'Current gas DELV, FWHM = ',DELV,FWHM
      print*,'Enter new values : '
      read*, DELV,FWHM

      WRITE(LUN0,REC=4)DELV
      WRITE(LUN0,REC=5)FWHM


      CLOSE(LUN0)

      END
