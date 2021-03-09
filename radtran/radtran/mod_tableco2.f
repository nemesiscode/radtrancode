      PROGRAM MOD_TABLECO2
C     $Id: read_table.f,v 1.12 2007-06-28 15:29:12 irwin Exp $
C***********************************************************************
C_TITL:	MOD_TABLECO2
C
C_DESC:	Modifies CO2 table
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
C***************************** VARIABLES *******************************

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
      INCLUDE '../includes/dbcom.f'
C NOTE: dbcom defines the linedata base variables. it is not normally
C stored in the same directory as the rest of the code
      INCLUDE '../includes/pathcom.f'

      INTEGER LUN,LUN0
      PARAMETER (LUN=2,LUN0=30)

      INTEGER PINDEX,CP,CT,LOOP,IPFORM

      INTEGER IREC,IREC0,CT2
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK),VCEN(MAXBIN)
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),TABLE(MAXK,MAXK,MAXG)
      REAL TEMP2(MAXK,MAXK),TN(MAXK),TX,X1,X2
      REAL KOUT,SMIN
      CHARACTER*200 KTAFIL,OPFILE1
      CHARACTER*1 ANS
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C******************************** CODE *********************************
      idiag=1

      CALL PROMPT('Enter k-table filename to modify: ')
      READ(5,23)OPFILE1
23    FORMAT(A)
      CALL FILE(OPFILE1,KTAFIL,'kta')

      WRITE(*,*)'Enter number of 4-byte words per record.'
      CALL PROMPT('old files = 2, new files = 1 : ')
      READ*,NW
      IF(NW.LT.1.OR.NW.GT.2)THEN
        WRITE(*,*)'NW must be 1 or 2 only! -- stopping program.'
        STOP
      ENDIF

      IRECL = NW*ISYS()
      WRITE(*,*)'irecl = ',irecl

      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL)
      READ(LUN0,REC=1)IREC0
      READ(LUN0,REC=2)NPOINT
      READ(LUN0,REC=3)VMIN
      READ(LUN0,REC=4)DELV
      READ(LUN0,REC=5)FWHM
      READ(LUN0,REC=6)NP
      READ(LUN0,REC=7)NT
      READ(LUN0,REC=8)NG
      READ(LUN0,REC=9)IDGAS(1)
      READ(LUN0,REC=10)ISOGAS(1)
      IREC = 11
      WRITE(*,*)' '
      WRITE(*,*)'NPOINT, IREC0 = ',NPOINT,IREC0

      IREC = IREC0

      Print*,'Enter limiting strength of k-table : '
      READ*,SMIN

      IREC=IREC0

      DO 10 IV=1,NPOINT
        DO 20 J=1,NP
          DO 30 K=1,NT
            DO 40 LOOP=1,NG
              READ(LUN0,REC=IREC)K_G(LOOP)
              IF(K_G(LOOP).LT.SMIN)K_G(LOOP)=SMIN
              WRITE(LUN0,REC=IREC)K_G(LOOP)
              IREC = IREC + 1
40          CONTINUE
30        CONTINUE
20      CONTINUE
10    CONTINUE


      CLOSE(LUN0)
 
      END
C***********************************************************************
C***********************************************************************

