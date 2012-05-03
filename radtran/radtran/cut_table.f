      PROGRAM CUT_TABLE
C     $Id: cut_table.f,v 1.1 2011-06-17 15:25:43 irwin Exp $
C***********************************************************************
C_TITL:	CUT_TABLE.f
C
C_DESC:	extracts specified wavelength range from a look-up table
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

      INTEGER LUN,LUN0,LUN1,MPOINT      
      PARAMETER (LUN=2,LUN0=30,LUN1=31)
C MAXOUT the maximum number of output points

      INTEGER PINDEX,CP,CT,LOOP

      INTEGER IREC,IREC0,NPOINT1,IREC1
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK),VCEN(MAXBIN),X1,X2
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),KVAL,VMIN1

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


27    PRINT*,'Current range = ',VMIN,VMIN+(NPOINT-1)*DELV,DELV,NPOINT

      IF(DELV.LT.0)THEN
       Print*,'DELV is less than zero. Reset?'
       READ(5,23)ANS
       IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')THEN
        CALL PROMPT('Enter new value : ')
        READ*,DELV
       ENDIF
      ENDIF

      CALL PROMPT('Enter output wavelength/wavenumber range : ')
      READ*,X1,X2

      IF(X1.LT.VMIN)X1=VMIN
      VMIN1 = VMIN + DELV*INT((X1-VMIN)/DELV)
      print*,X1,X2,VMIN1,DELV
      NPOINT1 = 1+INT((X2-VMIN1)/DELV)
      print*,NPOINT1
      PRINT*,'Extract range : ',VMIN1, VMIN1+(NPOINT1-1)*DELV,DELV,
     &  NPOINT1

      CALL PROMPT('OK (Y/N)? : ')
      READ(5,23)ANS

      IF(ANS.NE.'Y'.AND.ANS.NE.'y')GOTO 27

      WRITE(LUN1,REC=1)IREC0
      WRITE(LUN1,REC=2)NPOINT1
      WRITE(LUN1,REC=3)VMIN1
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



      N1 = 1 + NINT((VMIN1 - VMIN)/DELV)
      IREC = IREC0 + NP*NT*NG*(N1 - 1)
      IREC1 = IREC0

      WRITE(*,*)'IREC, IREC1 = ',IREC,IREC1

      DO 10 I=1,NPOINT1
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
