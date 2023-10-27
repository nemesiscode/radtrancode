      PROGRAM CONCAT_LBL_TABLE
C     $Id:
C***********************************************************************
C_TITL:	CONCAT_LBL_TABLE.f
C
C_DESC:	Concatenates two LBL-look-up tables
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
C       24oct23 PGJI    Updated for LBL look-up tables.
C***************************** VARIABLES *******************************

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
      INCLUDE '../includes/dbcom.f'
C NOTE: dbcom defines the linedata base variables. it is not normally
C stored in the same directory as the rest of the code
      INCLUDE '../includes/pathcom.f'

      INTEGER LUN,LUN1,LUN2,LUN0
      PARAMETER (LUN=2,LUN0=30)
      PARAMETER (LUN1=31,LUN2=32)
C MAXOUT the maximum number of output points

      INTEGER PINDEX,CP,CT,LOOP

      INTEGER IREC,IREC0,I1,I2,J1,J2
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK)
      REAL P2,T2,PRESS2(MAXK),TEMP2(MAXK)
      REAL PTEMP1(MAXK,MAXK),PTEMP2(MAXK,MAXK)

      CHARACTER*100 KTAFIL1,KTAFIL2,OPFILE1,OPFILE2,OPFILE3,OUTFIL
      CHARACTER*1 ANS
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C******************************** CODE *********************************

      idiag=1

      CALL PROMPT('Enter input filename (1) : ')
      READ(5,23)OPFILE1
23    FORMAT(A)
      CALL FILE(OPFILE1,KTAFIL1,'lta')

      CALL PROMPT('Enter input filename (2) : ')
      READ(5,23)OPFILE2
      CALL FILE(OPFILE2,KTAFIL2,'lta')

      CALL PROMPT('Enter output filename : ')
      READ(5,23)OPFILE3
      CALL FILE(OPFILE3,OUTFIL,'kta')


      IRECL1 = ISYS()
      IRECL2 = ISYS()
      IRECL = ISYS()
      WRITE(*,*)'irecl = ',irecl

      OPEN(UNIT=LUN1,FILE=KTAFIL1,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL1)
      OPEN(UNIT=LUN2,FILE=KTAFIL2,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL2)
      OPEN(UNIT=LUN0,FILE=OUTFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=IRECL)
      READ(LUN1,REC=1)IREC01
      READ(LUN1,REC=2)NPOINT1
      READ(LUN1,REC=3)VMIN1
      READ(LUN1,REC=4)DELV1
      READ(LUN1,REC=5)NP1
      READ(LUN1,REC=6)NT1
      READ(LUN1,REC=7)IDGAS1
      READ(LUN1,REC=8)ISOGAS1

      READ(LUN2,REC=1)IREC02
      READ(LUN2,REC=2)NPOINT2
      READ(LUN2,REC=3)VMIN2
      READ(LUN2,REC=4)DELV2
      READ(LUN2,REC=5)NP2
      READ(LUN2,REC=6)NT2
      READ(LUN2,REC=7)IDGAS2
      READ(LUN2,REC=8)ISOGAS2

      PRINT*,'VMIN1,DELV1,NPOINT1,VMAX : ',VMIN1,DELV1,NPOINT1,
     &   VMIN1 + (NPOINT1-1)*DELV1
      PRINT*,'VMIN2,DELV2,NPOINT2,VMAX : ',VMIN2,DELV2,NPOINT2,
     &   VMIN2 + (NPOINT2-1)*DELV2


      IF(DELV1.NE.DELV2)THEN
       print*,'DELV1, DELV2 = ',DELV1,DELV2
       print*,'Two tables need to be equally spaced'
       print*,'Aborting'
       STOP
      ENDIF

      IF(NP1.NE.NP2.OR.NT1.NE.NT2)THEN
       PRINT*,'Number of T/P points not same'
       PRINT*,'NP1, NP2 = ',NP1,NP2
       PRINT*,'NT1, NT2 = ',NT1,NT2
       STOP
      ENDIF

      IF(IDGAS1.NE.IDGAS2)THEN
       PRINT*,'IDGAS1 <> IDGAS2',IDGAS1,IDGAS2
       STOP
      ENDIF

      IF(ISOGAS1.NE.ISOGAS2)THEN
       PRINT*,'ISOGAS1 <> ISOGAS2',ISOGAS1,ISOGAS2
       STOP
      ENDIF


      IF(VMIN1.GT.VMIN2)THEN
       print*,'Input files in wrong order'
       print*,'VMIN1,VMIN2',VMIN1,VMIN2
       STOP
      ENDIF

      IREC=9

      DO 311 J=1,NP1
        READ(LUN1,REC=IREC)PRESS1(J)
        READ(LUN2,REC=IREC)PRESS2(J)
        IF(PRESS1(J).NE.PRESS2(J))THEN
         PRINT*,'PRESS1(J) <> PRESS2(J)',PRESS1(J),PRESS2(J)
         PRINT*,'Continue (Y/N)?'
         READ(5,23)ANS
         IF(ANS.NE.'Y'.AND.ANS.NE.'y')THEN
          STOP
         ENDIF
        ENDIF 
        IREC = IREC + 1
311   CONTINUE

      IF(NT1.GT.0)THEN
       DO 312 J=1,NT1
        READ(LUN1,REC=IREC)TEMP1(J)
        READ(LUN2,REC=IREC)TEMP2(J)
        IF(TEMP1(J).NE.TEMP2(J))THEN
         PRINT*,'TEMP1(J) <> TEMP2(J)',TEMP1(J),TEMP2(J)
         PRINT*,'Continue (Y/N)?'
         READ(5,23)ANS
         IF(ANS.NE.'Y'.AND.ANS.NE.'y')THEN
          STOP
         ENDIF
        ENDIF 
        IREC = IREC + 1
312    CONTINUE
      ELSE
       DO I=1,NP1
        DO J=1,ABS(NT1)
         READ(LUN1,REC=IREC)PTEMP1(I,J)
         READ(LUN2,REC=IREC)PTEMP2(I,J)
C         print*,I,J,PTEMP1(I,J),PTEMP2(I,J)
         IF(PTEMP1(I,J).NE.PTEMP2(I,J))THEN
          PRINT*,'PTEMP1(I,J) <> PTEMP2(I,J)',PTEMP1(I,J),PTEMP2(I,J)
          PRINT*,'Continue (Y/N)?'
          READ(5,23)ANS
          IF(ANS.NE.'Y'.AND.ANS.NE.'y')THEN
           STOP
          ENDIF
         ENDIF 
         IREC = IREC + 1
        ENDDO
       ENDDO
      ENDIF


      PRINT*,'Table1'
      PRINT*,'Index, Wavelength or Wavenumber'
      DO I=1,NPOINT1
        PRINT*,I,VMIN1+(I-1)*DELV1
      ENDDO

      PRINT*,'Enter starting and finishing wavelength indices of that'
      PRINT*,'part of Table 1 you want to extract and concatenate,'
      CALL PROMPT('i.e., enter indesired I1,I2 : ')
      READ*,I1,I2

      PRINT*,'Table2'
      PRINT*,'Index, Wavelength or Wavenumber'
      DO I=1,NPOINT2
        PRINT*,I,VMIN2+(I-1)*DELV2
      ENDDO

      PRINT*,'Enter starting and finishing wavelength indices of that'
      PRINT*,'part of Table 2 you want to extract and concatenate,'
      CALL PROMPT('i.e., enter desired J1,J2 : ')
      READ*,J1,J2

      NPOINT=1+I2-I1 + 1+J2-J1
      VMIN = VMIN1+(I1-1)*DELV1

      IF(NT1.LT.0)THEN
       IREC0=10 + NP1 + NP1*ABS(NT1) + 2
      ELSE
       IREC0=10 + NP1 + NT1 + 2
      ENDIF

      WRITE(*,*)'VMIN, DELV = ',VMIN,DELV1
      WRITE(*,*)'NP, NT = ',NP1,NT1
      WRITE(*,*)'Gas ID, ISO = ',IDGAS1,ISOGAS1
      WRITE(*,*)' '


      WRITE(LUN0,REC=1)IREC0
      WRITE(LUN0,REC=2)NPOINT
      WRITE(LUN0,REC=3)VMIN
      WRITE(LUN0,REC=4)DELV1
      WRITE(LUN0,REC=5)NP1
      WRITE(LUN0,REC=6)NT1
      WRITE(LUN0,REC=7)IDGAS1
      WRITE(LUN0,REC=8)ISOGAS1

      IREC=9

      WRITE(*,*)'Pressures : '
      DO 301 J=1,NP1
        print*,J,press1(J)
        WRITE(LUN0,REC=IREC)PRESS1(J)
        IREC = IREC + 1
301   CONTINUE

      WRITE(*,*)'Temperatures : '
      IF(NT1.GT.0)THEN
       DO 302 J=1,NT1
         print*,J,temp1(J)
         WRITE(LUN0,REC=IREC)TEMP1(J)
         IREC = IREC + 1
302    CONTINUE
      ELSE
       DO I=1,NP1
        DO J=1,ABS(NT1)
         WRITE(LUN0,REC=IREC)PTEMP1(I,J)
         IREC = IREC + 1
        ENDDO
       ENDDO
      ENDIF

      WRITE(*,*)'Wavelengths/wavenumbers'
      DO 403 J=I1,I2
        print*,I,VMIN1+(J-1)*DELV1
403   CONTINUE
      DO 503 J=J1,J2
        print*,I,VMIN2+(J-1)*DELV2
503   CONTINUE           

      N1=ABS(NT)

      IREC=IREC0

      DO 30 K=1,N1
       DO 20 J=1,NP1
        DO 297 I=I1,I2
           IREC1=IREC01 + (I-1)*NP1*N1+(J-1)*N1+K-1
           READ(LUN1,REC=IREC1)XK
           WRITE(LUN0,REC=IREC)XK
           IREC = IREC + 1
297     CONTINUE
        DO 298 I=J1,J2
           IREC2=IREC02 + (I-1)*NP1*N1+(J-1)*N1+K-1
           READ(LUN2,REC=IREC2)XK
           WRITE(LUN0,REC=IREC)XK
           IREC = IREC + 1
298     CONTINUE
20     CONTINUE
30    CONTINUE

      CLOSE(LUN0)
      CLOSE(LUN1)
      CLOSE(LUN2)

      END
