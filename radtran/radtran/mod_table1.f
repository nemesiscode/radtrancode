      PROGRAM MOD_TABLE1
C     $Id: cut_table.f,v 1.1 2011-06-17 15:25:43 irwin Exp $
C***********************************************************************
C_TITL:	CUT_TABLE.f
C
C_DESC: Changes k-table from channel to continuous
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

      INTEGER LUN,LUN0,LUN1,IX
      PARAMETER (LUN=2,LUN0=30,LUN1=31)
C MAXOUT the maximum number of output points

      INTEGER PINDEX,CP,CT,LOOP,I1,I2

      INTEGER IREC,IREC0,IREC1,IREC01
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK),VCEN(MAXBIN),X1,X2
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),KVAL

      CHARACTER*100 KTAFIL,OPFILE1,OPFILE2,OUTFIL
      CHARACTER*1 ANS
C******************************** CODE *********************************

      CALL PROMPT('Enter input filename : ')
      READ(5,23)OPFILE1
23    FORMAT(A)
      CALL FILE(OPFILE1,KTAFIL,'kta')

      CALL PROMPT('Enter output filename : ')
      READ(5,23)OPFILE2
      CALL FILE(OPFILE2,OUTFIL,'kta')


C     Assume 4-byte words per record
      NW=1

      IRECL = NW*ISYS()
      WRITE(*,*)'irecl = ',irecl

      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL)

      OPEN(UNIT=LUN1,FILE=OUTFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
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

      IF(DELV.GT.0)THEN
       IX=1
      ELSE
       IX=0
      ENDIF
      IF(IX.EQ.0)THEN
       PRINT*,'DELV = ',DELV
       PRINT*,'Enter new DELV'
       READ*,DELV
       PRINT*,'Enter new FWHM'
       READ*,FWHM
       IREC01=11 + 2*NG + 2 + NP + NT + 2
      ELSE
       DELV = -1.
       FWHM = 0.
       PRINT*,'Enter channel centres'
       READ*,(VCEN(I),I=1,NPOINT)
       IREC01=11 + 2*NG + 2 + NP + NT + 2 + NPOINT
      ENDIF

      print*,'A'
      WRITE(LUN1,REC=1)IREC01
      WRITE(LUN1,REC=2)NPOINT
      WRITE(LUN1,REC=3)VMIN
      WRITE(LUN1,REC=4)DELV
      WRITE(LUN1,REC=5)FWHM
      WRITE(LUN1,REC=6)NP
      WRITE(LUN1,REC=7)NT
      WRITE(LUN1,REC=8)NG
      WRITE(LUN1,REC=9)IDGAS(1)
      WRITE(LUN1,REC=10)ISOGAS(1)

      IREC = 11
      DO 299 J=1,NG
        READ(LUN0,REC=IREC)G_ORD(J)
        WRITE(LUN1,REC=IREC)G_ORD(J)
        IREC = IREC + 1
299   CONTINUE
      DO 399 J=1,NG
        READ(LUN0,REC=IREC)DEL_G(J)
        WRITE(LUN1,REC=IREC)DEL_G(J)
        IREC=IREC+1
399   CONTINUE
      IREC = 11 + 2*NG + 2
      DO 301 J=1,NP
        READ(LUN0,REC=IREC)PRESS1(J)
        WRITE(LUN1,REC=IREC)PRESS1(J)
        IREC = IREC + 1
301   CONTINUE
      DO 302 J=1,NT
        READ(LUN0,REC=IREC)TEMP1(J)
        WRITE(LUN1,REC=IREC)TEMP1(J)
        IREC = IREC + 1
302   CONTINUE
      IF(DELV.LT.0)THEN
       DO 303 J=1,NPOINT
        WRITE(LUN1,REC=IREC)VCEN(J)
        IREC = IREC + 1
303    CONTINUE
      ENDIF


      IREC = IREC0
      IREC1 = IREC01
      DO 10 I=1,NPOINT
       DO 20 J=1,NP
        DO 30 K=1,NT
          DO 40 LOOP=1,NG
            READ(LUN0,REC=IREC)KVAL
            WRITE(LUN1,REC=IREC1)KVAL
            IREC = IREC + 1
            IREC1 = IREC1 + 1
40        CONTINUE
30      CONTINUE
20     CONTINUE
10    CONTINUE

      CLOSE(LUN0)
      CLOSE(LUN1)

      END
