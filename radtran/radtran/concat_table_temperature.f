      PROGRAM CONCAT_TABLE_TEMPERATURE
C     $Id:
C***********************************************************************
C_TITL:	CONCAT_TABLE.f
C
C_DESC:	Concatenates two k-look-up tables, which have the same wavelength
C       and pressure grids, but different temperature grids.
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
C	 7nov12 JB	Revised from Concat_table to deal with k-tables that
C			cover same temperature and pressure ranges, but 
C			different temperatures.
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
      INTEGER IREC,IREC0,I1,I2,J1,J2, NT,KT1,KT2,LT1,LT2
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK),VCEN1(MAXBIN)
      REAL P2,T2,PRESS2(MAXK),TEMP2(MAXK),VCEN2(MAXBIN)
      REAL VCEN(MAXBIN),VEXPECT
      REAL G_ORD1(MAXG),K_G1(MAXG),DEL_G1(MAXG)
      REAL G_ORD2(MAXG),K_G2(MAXG),DEL_G2(MAXG)
      REAL TABLE(MAXK,2*MAXK,MAXG)
      REAL TABLE1(MAXK,MAXK,MAXG), TABLE2(MAXK,MAXK,MAXG)

      CHARACTER*100 KTAFIL1,KTAFIL2,OPFILE1,OPFILE2,OPFILE3,OUTFIL
      CHARACTER*1 ANS

C******************************** CODE *********************************

      CALL PROMPT('Enter input filename (1, cooler T) : ')
      READ(5,23)OPFILE1
23    FORMAT(A)
      CALL FILE(OPFILE1,KTAFIL1,'kta')
      PRINT*,'Enter number of 4-byte words per record.'
      CALL PROMPT('old files = 2, new files = 1 : ')
      READ*,NW1

      CALL PROMPT('Enter input filename (2, hotter T) : ')
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
C      IF(NP1.NE.NP2.OR.NT1.NE.NT2)THEN
C       PRINT*,'Number of T/P points not same'
C       PRINT*,'NP1, NP2 = ',NP1,NP2
C       PRINT*,'NT1, NT2 = ',NT1,NT2
C       STOP
C      ENDIF

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


C      IF(VMIN1.GT.VMIN2)THEN
C       print*,'Input files in wrong order'
C       print*,'VMIN1,VMIN2',VMIN1,VMIN2
C       STOP
C      ENDIF

      PRINT*,'Table1'
      DO I=1,NPOINT1
       PRINT*,I,VMIN1+(I-1)*DELV1
      ENDDO
      CALL PROMPT('Enter desired I1,I2 : ')
      READ*,I1,I2

      PRINT*,'Table2'
      DO I=1,NPOINT2
       PRINT*,I,VMIN2+(I-1)*DELV2
      ENDDO
      CALL PROMPT('Enter desired J1,J2 : ')
      READ*,J1,J2
      
      IF((J2-J1).NE.(I2-I1))THEN
       PRINT*,'NPOINTS1 <> NPOINTS2'
       PRINT*,'Tables must have same wavelength range'
       STOP
      ENDIF

      
      IF((VMIN1+(I1-1)*DELV1.ne.VMIN2+(J1-1)*DELV1))then
	  PRINT*,'Selected wavelength ranges are different!'
          print*,'I1,J1 = ',I1,J1
          print*,VMIN1+(I1-1)*DELV1,VMIN2+(J1-1)*DELV1
          STOP
      ENDIF
      IF((VMIN1+(I2-1)*DELV1.ne.VMIN2+(J2-1)*DELV1))then
	  PRINT*,'Selected wavelength ranges are different!'
          print*,'I2,J2 = ',I2,J2
          print*,VMIN1+(I2-1)*DELV1,VMIN2+(J2-1)*DELV1
          STOP
      ENDIF

      NPOINT=1+I2-I1
      VMIN = VMIN1+(I1-1)*DELV1

      IREC=11
      
      PRINT*,'Select temperature range for table1 (indices):'
      CALL PROMPT('Enter desired KT1,KT2 : ')
      READ*,KT1,KT2
      
      PRINT*,'Select temperature range for table2 (indices):'
      CALL PROMPT('Enter desired LT1,LT2 : ')
      READ*,LT1,LT2
      
      NT=(KT2-KT1)+1+(LT2-LT1)+1
      PRINT*,'NT=',NT
      PRINT*, 'NT1=', NT1
      PRINT*, 'NT2=', NT2
      

      WRITE(*,*)'VMIN = ',VMIN
      WRITE(*,*)'DELV, FWHM = ',DELV1,FWHM1
      WRITE(*,*)'NP, NT, NG = ',NP1,NT,NG1
      WRITE(*,*)'Gas ID, ISO = ',IDGAS1,ISOGAS1
      WRITE(*,*)' '
      IREC0=11 + 2*NG1 + 2 + NP1 + NT + 2
      IF(DELV1.LT.0)THEN
        IREC0=IREC0+NPOINT
      ENDIF
 
      WRITE(LUN0,REC=1)IREC0
      WRITE(LUN0,REC=2)NPOINT
      WRITE(LUN0,REC=3)VMIN
      WRITE(LUN0,REC=4)DELV1
      WRITE(LUN0,REC=5)FWHM1
      WRITE(LUN0,REC=6)NP1
      WRITE(LUN0,REC=7)NT
      WRITE(LUN0,REC=8)NG1
      WRITE(LUN0,REC=9)IDGAS1
      WRITE(LUN0,REC=10)ISOGAS1

      DO 299 J=1,NG1
        READ(LUN1,REC=IREC)G_ORD1(J)
        READ(LUN2,REC=IREC)G_ORD2(J)
        WRITE(LUN0,REC=IREC)G_ORD1(J)
        print*,J,G_ORD1(J)
        IREC = IREC + 1
        IF(G_ORD1(J).NE.G_ORD2(J)) THEN
         PRINT*,'G_ORD1 <> G_ORD2',G_ORD1(J),G_ORD2(J)
         STOP
        ENDIF
299   CONTINUE


      WRITE(*,*)'G - ordinates, weights'
      DO 399 J=1,NG1
        READ(LUN1,REC=IREC)DEL_G1(J)
        READ(LUN2,REC=IREC)DEL_G2(J)
        print*,J,DEL_G1(J)
        WRITE(LUN0,REC=IREC)DEL_G1(J)
        IREC=IREC+1
        IF(DEL_G1(J).NE.DEL_G2(J)) THEN
         PRINT*,'DEL_G1 <> DEL_G2',DEL_G1(J),DEL_G2(J)
         STOP
        ENDIF
399   CONTINUE


      IREC = 11 + 2*NG1 + 2
      WRITE(*,*)'Pressures : '
      DO 301 J=1,NP1
        READ(LUN1,REC=IREC)PRESS1(J)
        READ(LUN2,REC=IREC)PRESS2(J)
        print*,J,press1(J)
        WRITE(LUN0,REC=IREC)PRESS1(J)
        IF(PRESS1(J).NE.PRESS2(J))THEN
         PRINT*,'PRESS1(J) <> PRESS2(J)',PRESS1(J),PRESS2(J)
        ENDIF 
        IREC = IREC + 1
301   CONTINUE


      WRITE(*,*)'Temperatures : '
      DO 302 J=1,NT1
        READ(LUN1,REC=IREC)TEMP1(J)
        IF(J.ge.KT1.AND.J.le.KT2)THEN
        print*,J,temp1(J)
        print*,'IREC=',IREC
        WRITE(LUN0,REC=IREC)TEMP1(J)
        ENDIF
        IREC = IREC + 1
        IREC01=IREC
302   CONTINUE

      IREC = IREC-NT1
      WRITE(*,*)'Temperatures : '
      DO 303 J=1,NT2
        READ(LUN2,REC=IREC)TEMP2(J)
        IF(J.ge.LT1.AND.J.le.LT2) THEN
        print*,J,temp2(J)
        print*,'IREC=',IREC
        WRITE(LUN0,REC=IREC+(KT2-KT1)+1-LT1+1)TEMP2(J)
        ENDIF
        IREC = IREC + 1 
303   CONTINUE

      IF (TEMP1(1).gt.TEMP2(1)) then
  	print*,'Files entered in wrong order!'
  	stop
      ENDIF


      IREC = IREC0
      IREC01= 11 + 2*NG1 + 2 + NP1 + NT1 + 2
      IREC02= 11 + 2*NG1 + 2 + NP1 + NT2 + 2
      IREC1 = IREC01+(I1-1)*(NP1*(NG1)*NT1)
      IREC2 = IREC02+(J1-1)*(NP1*(NG1)*NT2)
    

      DO 297 I=1,NPOINT1
       DO 20 K=1,NT1
         DO 30 J=1,NP1
           DO 40 LOOP=1,NG1
             READ(LUN1,REC=IREC1)TABLE1(J,K,LOOP) 
       IF((K.ge.KT1.AND.K.le.KT2).and.(I.GE.I1.AND.I.LE.I2))THEN
             WRITE(LUN0,REC=IREC)TABLE1(J,K,LOOP)
             IREC=IREC+1
             ENDIF
             IREC1 = IREC1 + 1
40         CONTINUE
30       CONTINUE
20     CONTINUE
       DO 22 K=1, NT2
        DO 32 J=1,NP1
         DO 42 LOOP=1,NG1	
           READ(LUN2,REC=IREC2)TABLE2(J,K,LOOP)
        IF((K.ge.LT1.AND.K.le.LT2).and.(I.GE.I1.AND.I.LE.I2))THEN
            WRITE(LUN0,REC=IREC)TABLE2(J,K,LOOP)
            IREC=IREC+1
           ENDIF
           IREC2 = IREC2 + 1
42       CONTINUE
32      CONTINUE
22     CONTINUE
297   CONTINUE

      CLOSE(LUN0)
      CLOSE(LUN1)
      CLOSE(LUN2)

      END
