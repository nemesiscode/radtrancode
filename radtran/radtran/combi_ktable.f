      PROGRAM COMBI_KTABLE
C     $Id:
C***********************************************************************
C_TITL:	COMBI_KTABLE.f
C
C_DESC:	Concatenates four k-look-up tables to match same wavelength steps
C       as a fifth.
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
C_HIST:	08oct23	PGJI	ORIGINAL VERSION
C***************************** VARIABLES *******************************

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
      INCLUDE '../includes/dbcom.f'
C NOTE: dbcom defines the linedata base variables. it is not normally
C stored in the same directory as the rest of the code
      INCLUDE '../includes/pathcom.f'

      INTEGER LUN,LUNI(4),LUNREF,LUN0

      INTEGER PINDEX,CP,CT,LOOP

      INTEGER IREC,IREC0,IDUM
      REAL P(4),T(4),PRESSREF(MAXK),TEMPREF(MAXK),VCENREF(MAXBIN)
      REAL PTEMPREF(MAXK,MAXK),VCENI(4,MAXBIN)
      REAL VCEN(MAXBIN),VEXPECT
      REAL G_ORDREF(MAXG),K_G(MAXG),DEL_GREF(MAXG)
      REAL TABLE(MAXK,MAXK,MAXG)
      REAL VMINI(4),DELVI(4),FWHMI(4)

      INTEGER NPOINTI(4),NP(4),NT(4),NG(4),IDGASI(4),ISOGASI(4)
      INTEGER NPOINTREF,NPREF,NTREF,NGREF,IDGASREF,ISOGASREF
      INTEGER IDFIL(MAXBIN),IDOFF(MAXBIN)

      CHARACTER*100 KTAFIL1,KTAFIL2,OPFILE1,OPFILE2,OUTFIL
      CHARACTER*100 KTAFIL3,KTAFIL4,OPFILE3,OPFILE4
      CHARACTER*100 OPFILEREF,KTAFILREF,OPFILEOUT
      CHARACTER*1 ANS
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C******************************** CODE *********************************

      idiag=1
23    FORMAT(A)
      NW=1

      LUNI(1)=10
      LUNI(2)=11
      LUNI(3)=12
      LUNI(4)=13
      LUNREF=14
      LUN0=15

      CALL PROMPT('input k-table to extract wavelengths: ')
      READ(5,23)OPFILEREF
      CALL FILE(OPFILEREF,KTAFILREF,'kta')

      CALL PROMPT('Enter input filename (1) : ')
      READ(5,23)OPFILE1
      CALL FILE(OPFILE1,KTAFIL1,'kta')
      CALL PROMPT('Enter input filename (2) : ')
      READ(5,23)OPFILE2
      CALL FILE(OPFILE2,KTAFIL2,'kta')
      CALL PROMPT('Enter input filename (3) : ')
      READ(5,23)OPFILE3
      CALL FILE(OPFILE3,KTAFIL3,'kta')
      CALL PROMPT('Enter input filename (4) : ')
      READ(5,23)OPFILE4
      CALL FILE(OPFILE4,KTAFIL4,'kta')

      CALL PROMPT('Enter output filename : ')
      READ(5,23)OPFILEOUT
      CALL FILE(OPFILEOUT,OUTFIL,'kta')


      IRECL = NW*ISYS()

      OPEN(UNIT=LUNREF,FILE=KTAFILREF,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL)
      OPEN(UNIT=LUNI(1),FILE=KTAFIL1,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL)
      OPEN(UNIT=LUNI(2),FILE=KTAFIL2,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL)
      OPEN(UNIT=LUNI(3),FILE=KTAFIL3,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL)
      OPEN(UNIT=LUNI(4),FILE=KTAFIL4,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL)
      OPEN(UNIT=LUN0,FILE=OUTFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=IRECL)

      DO IFILE=1,4
       READ(LUNI(IFILE),REC=1)IDUM
       READ(LUNI(IFILE),REC=2)NPOINTI(IFILE)
       READ(LUNI(IFILE),REC=3)VMINI(IFILE)
       READ(LUNI(IFILE),REC=4)DELVI(IFILE)
       READ(LUNI(IFILE),REC=5)FWHMI(IFILE)
       READ(LUNI(IFILE),REC=6)NP(IFILE)
       READ(LUNI(IFILE),REC=7)NT(IFILE)
       READ(LUNI(IFILE),REC=8)NG(IFILE)
       READ(LUNI(IFILE),REC=9)IDGASI(IFILE)
       READ(LUNI(IFILE),REC=10)ISOGASI(IFILE)
      ENDDO

      READ(LUNREF,REC=1)IRECREF
      READ(LUNREF,REC=2)NPOINTREF
      READ(LUNREF,REC=3)VMINREF
      READ(LUNREF,REC=4)DELVREF
      READ(LUNREF,REC=5)FWHMREF
      READ(LUNREF,REC=6)NPREF
      READ(LUNREF,REC=7)NTREF
      READ(LUNREF,REC=8)NGREF
      READ(LUNREF,REC=9)IDGASREF
      READ(LUNREF,REC=10)ISOGASREF

      IREC=11
      DO 296 J=1,NGREF
        READ(LUNREF,REC=IREC)G_ORDREF(J)
        IREC=IREC+1
296   CONTINUE

      DO 398 J=1,NGREF
        READ(LUNREF,REC=IREC)DEL_GREF(J)
        IREC=IREC+1
398   CONTINUE

      IREC = 11 + 2*NGREF + 2

      DO 311 J=1,NP1
        READ(LUNREF,REC=IREC)PRESSREF(J)
        IREC = IREC + 1
311   CONTINUE

      IF(NTREF.GT.0)THEN
       DO 312 J=1,NT1
        READ(LUNREF,REC=IREC)TEMPREF(J)
        IREC = IREC + 1
312    CONTINUE
      ELSE
       DO I=1,NPREF
        DO J=1,ABS(NTREF)
         READ(LUNREF,REC=IREC)PTEMPREF(I,J)
         IREC = IREC + 1
        ENDDO
       ENDDO
      ENDIF

      IF(NTREF.LT.0)THEN
       IREC0=11 + 2*NGREF + 2 + NPREF + NPREF*ABS(NTREF) + 2
      ELSE
       IREC0=11 + 2*NGREF + 2 + NPREF + NTREF + 2
      ENDIF

      IREC=IREC0
      DO I=1,NPOINTREF
        READ(LUNREF,REC=IREC)VCENREF(I)
        IREC=IREC+1
      ENDDO

      DO 313 IFILE=1,4
       IREC=IREC0
       DO I=1,NPOINTI(IFILE)
        READ(LUNI(IFILE),REC=IREC)VCENI(IFILE,I)
        DO J=1,NPOINTREF
         IF(VCENI(IFILE,I).EQ.VCENREF(J))THEN
          IDFIL(J)=IFILE
          IDOFF(J)=I
         ENDIF
        ENDDO
        IREC=IREC+1
       ENDDO
313   CONTINUE


      DO 314 I=1,NPOINTREF
       PRINT*,I,VCENREF(I),IDFIL(I),IDOFF(I),VCENI(IDFIL(I),
     1   IDOFF(I))     
314   CONTINUE

      stop

      WRITE(LUN0,REC=1)IREC0
      WRITE(LUN0,REC=2)NPOINTREF
      WRITE(LUN0,REC=3)VMINREF
      WRITE(LUN0,REC=4)DELV1REF
      WRITE(LUN0,REC=5)FWHM1REF
      WRITE(LUN0,REC=6)NPREF
      WRITE(LUN0,REC=7)NTREF
      WRITE(LUN0,REC=8)NGREF
      WRITE(LUN0,REC=9)IDGASREF
      WRITE(LUN0,REC=10)ISOGASREF

      WRITE(*,*)'G - ordinates'
      DO 299 J=1,NGREF
        PRINT*,J,G_ORDREF(J)
        WRITE(LUN0,REC=IREC)G_ORDREF(J)
        IREC = IREC + 1
299   CONTINUE


      WRITE(*,*)'G - ordinates, weights'
      DO 399 J=1,NGREF
        print*,J,DEL_GREF(J)
        WRITE(LUN0,REC=IREC)DEL_GREF(J)
        IREC=IREC+1
399   CONTINUE

      IREC = 11 + 2*NGREF + 2
      WRITE(*,*)'Pressures : '
      DO 301 J=1,NPREF
        print*,J,pressref(J)
        WRITE(LUN0,REC=IREC)PRESSREF(J)
        IREC = IREC + 1
301   CONTINUE

      WRITE(*,*)'Temperatures : '
      IF(NTREF.GT.0)THEN
       DO 302 J=1,NTREF
         print*,J,tempREF(J)
         WRITE(LUN0,REC=IREC)TEMPREF(J)
         IREC = IREC + 1
302    CONTINUE
      ELSE
       DO I=1,NPREF
        DO J=1,ABS(NTREF)
         WRITE(LUN0,REC=IREC)PTEMPREF(I,J)
         IREC = IREC + 1
        ENDDO
       ENDDO
      ENDIF

      I=1
      WRITE(*,*)'Wavelengths/wavenumbers'
      IF(DELV1.LE.0.0)THEN
       DO 303 J=1,NPOINTREF
        WRITE(LUN0,REC=IREC)VCENREF(J)
        print*,I,VCENREF(J)
        IREC=IREC+1
        I=I+1
303    CONTINUE
      ELSE
       print*,'Abort 1'
       stop
      ENDIF

      IREC = IREC0+NPOINTREF

      DO 297 I=1,NPOINTREF
       IFILE=IDFIL(I)
       I1=IDOFF(I)

       IREC1=IREC0+NPOINTI(IFILE)+(I1-1)*NPREF*ABS(NTREF)*NGREF
       
       DO 20 J=1,NPREF
         DO 30 K=1,ABS(NTREF)
           DO 40 LOOP=1,NGREF
             READ(LUNI(IFILE),REC=IREC1)TABLE(J,K,LOOP)
             IREC1 = IREC1 + 1
             WRITE(LUN0,REC=IREC)TABLE(J,K,LOOP)
             IREC = IREC + 1
40         CONTINUE
30       CONTINUE
20     CONTINUE
297   CONTINUE

      CLOSE(LUNREF)
      CLOSE(LUN0)

      DO I=1,4
       CLOSE(LUNI(I))
      ENDDO

      END
