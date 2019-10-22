      PROGRAM CONCAT_TABLE
C     $Id:
C***********************************************************************
C_TITL:	CONCAT_TABLE.f
C
C_DESC:	Concatenates two k-look-up table
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

      INTEGER LUN,LUN1,LUN2,LUN0
      PARAMETER (LUN=2,LUN0=30)
      PARAMETER (LUN1=31,LUN2=32)
C MAXOUT the maximum number of output points

      INTEGER PINDEX,CP,CT,LOOP

      INTEGER IREC,IREC0,I1,I2,J1,J2
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK),VCEN1(MAXBIN)
      REAL P2,T2,PRESS2(MAXK),TEMP2(MAXK),VCEN2(MAXBIN)
      REAL VCEN(MAXBIN),VEXPECT
      REAL G_ORD1(MAXG),K_G1(MAXG),DEL_G1(MAXG)
      REAL G_ORD2(MAXG),K_G2(MAXG),DEL_G2(MAXG)
      REAL TABLE(MAXK,MAXK,MAXG)

      CHARACTER*100 KTAFIL1,KTAFIL2,OPFILE1,OPFILE2,OPFILE3,OUTFIL
      CHARACTER*1 ANS

C******************************** CODE *********************************

      CALL PROMPT('Enter input filename (1) : ')
      READ(5,23)OPFILE1
23    FORMAT(A)
      CALL FILE(OPFILE1,KTAFIL1,'kta')
      PRINT*,'Enter number of 4-byte words per record.'
      CALL PROMPT('old files = 2, new files = 1 : ')
      READ*,NW1

      CALL PROMPT('Enter input filename (2) : ')
      READ(5,23)OPFILE2
      CALL FILE(OPFILE2,KTAFIL2,'kta')
      PRINT*,'Enter number of 4-byte words per record.'
      CALL PROMPT('old files = 2, new files = 1 : ')
      READ*,NW2

      CALL PROMPT('Enter output filename : ')
      READ(5,23)OPFILE3
      CALL FILE(OPFILE3,OUTFIL,'kta')
      PRINT*,'Enter number of 4-byte words per record.'
      CALL PROMPT('old files = 2, new files = 1 : ')
      READ*,NW


      IRECL1 = NW1*ISYS()
      IRECL2 = NW2*ISYS()
      IRECL = NW*ISYS()
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
      READ(LUN1,REC=5)FWHM1
      READ(LUN1,REC=6)NP1
      READ(LUN1,REC=7)NT1
      READ(LUN1,REC=8)NG1
      READ(LUN1,REC=9)IDGAS1
      READ(LUN1,REC=10)ISOGAS1

      READ(LUN2,REC=1)IREC02
      READ(LUN2,REC=2)NPOINT2
      READ(LUN2,REC=3)VMIN2
      READ(LUN2,REC=4)DELV2
      READ(LUN2,REC=5)FWHM2
      READ(LUN2,REC=6)NP2
      READ(LUN2,REC=7)NT2
      READ(LUN2,REC=8)NG2
      READ(LUN2,REC=9)IDGAS2
      READ(LUN2,REC=10)ISOGAS2

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
      IF(FWHM1.NE.FWHM2)THEN
       print*,'FWHM1, FWHM2 = ',FWHM1,FWHM2
       print*,'Two tables should have same FWHM'
       print*,'Aborting'
       STOP
      ENDIF
      IF(NP1.NE.NP2.OR.NT1.NE.NT2)THEN
       PRINT*,'Number of T/P points not same'
       PRINT*,'NP1, NP2 = ',NP1,NP2
       PRINT*,'NT1, NT2 = ',NT1,NT2
       STOP
      ENDIF

      IF(NG1.NE.NG2)THEN
       PRINT*,'NG1 <> NG2',NG1,NG2
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

      IREC=11
      DO 296 J=1,NG1
        READ(LUN1,REC=IREC)G_ORD1(J)
        READ(LUN2,REC=IREC)G_ORD2(J)
        IREC=IREC+1
        IF(G_ORD1(J).NE.G_ORD2(J)) THEN
         PRINT*,'J, G_ORD1 <> G_ORD2',J, G_ORD1(J),G_ORD2(J)
         PRINT*,'Continue (Y/N)?'
         READ(5,23)ANS
         IF(ANS.NE.'Y'.AND.ANS.NE.'y')THEN
          STOP
         ENDIF
        ENDIF
296   CONTINUE

      DO 398 J=1,NG1
        READ(LUN1,REC=IREC)DEL_G1(J)
        READ(LUN2,REC=IREC)DEL_G2(J)
        IREC=IREC+1
        IF(DEL_G1(J).NE.DEL_G2(J)) THEN
         PRINT*,'J,DEL_G1 <> DEL_G2',J,DEL_G1(J),DEL_G2(J)
         PRINT*,'Continue (Y/N)?'
         READ(5,23)ANS
         IF(ANS.NE.'Y'.AND.ANS.NE.'y')THEN
          STOP
         ENDIF
        ENDIF
398   CONTINUE

      IREC = 11 + 2*NG1 + 2

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
312   CONTINUE

      IREC1=IREC
      IREC2=IREC
      PRINT*,'Table1'
      PRINT*,'Index, Wavelength or Wavenumber'
      IF(DELV1.GT.0)THEN
       DO I=1,NPOINT1
        PRINT*,I,VMIN1+(I-1)*DELV1
       ENDDO
      ELSE
       DO I=1,NPOINT1
        READ(LUN1,REC=IREC1)VCEN1(I)
        PRINT*,I,VCEN1(I)
        IREC1=IREC1+1
       ENDDO
      ENDIF
      PRINT*,'Enter starting and finishing wavelength indices of that'
      PRINT*,'part of Table 1 you want to extract and concatenate,'
      CALL PROMPT('i.e., enter indesired I1,I2 : ')
      READ*,I1,I2

      PRINT*,'Table2'
      PRINT*,'Index, Wavelength or Wavenumber'
      IF(DELV2.GT.0)THEN
       DO I=1,NPOINT2
        PRINT*,I,VMIN2+(I-1)*DELV2
       ENDDO
      ELSE
       DO I=1,NPOINT2
        READ(LUN2,REC=IREC2)VCEN2(I)
        PRINT*,I,VCEN2(I)
        IREC2=IREC2+1
       ENDDO
      ENDIF

      PRINT*,'Enter starting and finishing wavelength indices of that'
      PRINT*,'part of Table 2 you want to extract and concatenate,'
      CALL PROMPT('i.e., enter desired J1,J2 : ')
      READ*,J1,J2


      NPOINT=1+I2-I1 + 1+J2-J1
      IF(DELV1.GT.0)THEN
       VMIN = VMIN1+(I1-1)*DELV1
      ELSE
       VMIN = VCEN1(I1)
      ENDIF

      IREC0=11 + 2*NG1 + 2 + NP1 + NT1 + 2
      IF(DELV1.LT.0)THEN
        IREC0=IREC0+NPOINT
      ENDIF
      IREC=11

      WRITE(*,*)'VMIN = ',VMIN
      WRITE(*,*)'DELV, FWHM = ',DELV1,FWHM1
      WRITE(*,*)'NP, NT, NG = ',NP1,NT1,NG1
      WRITE(*,*)'Gas ID, ISO = ',IDGAS1,ISOGAS1
      WRITE(*,*)' '


      WRITE(LUN0,REC=1)IREC0
      WRITE(LUN0,REC=2)NPOINT
      WRITE(LUN0,REC=3)VMIN
      WRITE(LUN0,REC=4)DELV1
      WRITE(LUN0,REC=5)FWHM1
      WRITE(LUN0,REC=6)NP1
      WRITE(LUN0,REC=7)NT1
      WRITE(LUN0,REC=8)NG1
      WRITE(LUN0,REC=9)IDGAS1
      WRITE(LUN0,REC=10)ISOGAS1

      WRITE(*,*)'G - ordinates'
      DO 299 J=1,NG1
        PRINT*,J,G_ORD1(J)
        WRITE(LUN0,REC=IREC)G_ORD1(J)
        IREC = IREC + 1
299   CONTINUE


      WRITE(*,*)'G - ordinates, weights'
      DO 399 J=1,NG1
        print*,J,DEL_G1(J)
        WRITE(LUN0,REC=IREC)DEL_G1(J)
        IREC=IREC+1
399   CONTINUE

      IREC = 11 + 2*NG1 + 2
      WRITE(*,*)'Pressures : '
      DO 301 J=1,NP1
        print*,J,press1(J)
        WRITE(LUN0,REC=IREC)PRESS1(J)
        IREC = IREC + 1
301   CONTINUE

      WRITE(*,*)'Temperatures : '
      DO 302 J=1,NT1
        print*,J,temp1(J)
        WRITE(LUN0,REC=IREC)TEMP1(J)
        IREC = IREC + 1
302   CONTINUE

      I=1
      WRITE(*,*)'Wavelengths/wavenumbers'
      IF(DELV1.LT.0.0)THEN
       DO 303 J=I1,I2
        WRITE(LUN0,REC=IREC)VCEN1(J)
        print*,I,VCEN1(J)
        IREC=IREC+1
        I=I+1
303    CONTINUE
       DO 304 J=J1,J2
        WRITE(LUN0,REC=IREC)VCEN2(J)
        print*,I,VCEN2(J)
        IREC=IREC+1
        I=I+1
304    CONTINUE
      ELSE
       DO 403 J=I1,I2
        print*,I,VMIN1+(J-1)*DELV1
        I=I+1
403    CONTINUE
       DO 503 J=J1,J2
        print*,I,VMIN2+(J-1)*DELV2
        I=I+1
503    CONTINUE           
      ENDIF

      IREC = IREC0
      IREC1 = IREC01+(I1-1)*NP1*NT1*NG1

      DO 297 I=I1,I2
       DO 20 J=1,NP1
         DO 30 K=1,NT1
           DO 40 LOOP=1,NG1
             READ(LUN1,REC=IREC1)TABLE(J,K,LOOP)
             IREC1 = IREC1 + 1
             IF(I.GE.I1.AND.I.LE.I2)THEN
              WRITE(LUN0,REC=IREC)TABLE(J,K,LOOP)
              IREC = IREC + 1
             ENDIF
40         CONTINUE
30       CONTINUE
20     CONTINUE
297   CONTINUE

      IREC2 = IREC02+(J1-1)*NP2*NT2*NG2

      DO 298 I=J1,J2
       DO 21 J=1,NP2
         DO 31 K=1,NT2
           DO 41 LOOP=1,NG2
             READ(LUN2,REC=IREC2)TABLE(J,K,LOOP)
             IREC2 = IREC2 + 1
             IF(I.GE.J1.AND.I.LE.J2)THEN
              WRITE(LUN0,REC=IREC)TABLE(J,K,LOOP)
              IREC = IREC + 1
             ENDIF
41         CONTINUE
31       CONTINUE
21     CONTINUE
298   CONTINUE

      CLOSE(LUN0)
      CLOSE(LUN1)
      CLOSE(LUN2)

      END
