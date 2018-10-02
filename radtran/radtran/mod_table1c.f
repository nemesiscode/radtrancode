      PROGRAM MOD_TABLE1C
C     $Id:
C***********************************************************************
C_TITL:	CUT_TABLE.f
C
C_DESC: Changes k-table from channel to continuous. Update of mod_table1b, which
C       is able to deal with k-table format where calculation tables are not
C       fixed, but vary with temperature.
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

      implicit none
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
      INCLUDE '../includes/dbcom.f'
C NOTE: dbcom defines the linedata base variables. it is not normally
C stored in the same directory as the rest of the code
      INCLUDE '../includes/pathcom.f'

      INTEGER LUN,LUN0,LUN1,IX,NW,IRECL,ISYS
      INTEGER NP,NT,J,I,N1,N2,NG
      REAL VV1,VV2,X
      PARAMETER (LUN=2,LUN0=30,LUN1=31)
C MAXOUT the maximum number of output points

      INTEGER PINDEX,CP,CT,LOOP,I1,I2

      INTEGER IREC,IREC0,IREC1,IREC01
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK),VCEN(MAXBIN),X1,X2
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),KVAL,TEMP2(MAXK,MAXK)

      CHARACTER*100 KTAFIL,OPFILE1,OPFILE2,OUTFIL
      CHARACTER*1 ANS
C******************************** CODE *********************************

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
      IF(DELV.GT.0.0)THEN
       WRITE(*,*)'VMIN, VMAX = ',VMIN,dble(VMIN + (NPOINT - 1)*DELV)
      ELSE
       WRITE(*,*)'VMIN = ',VMIN
      ENDIF
      WRITE(*,*)'DELV, FWHM = ',dble(DELV),dble(FWHM)
      WRITE(*,*)'NP, NT, NG = ',NP,NT,NG
      WRITE(*,*)'Gas ID, ISO = ',IDGAS(1),ISOGAS(1)
      WRITE(*,*)' '
      DO 299 J=1,NG
        READ(LUN0,REC=IREC)G_ORD(J)
        IREC = IREC + 1
299   CONTINUE
      WRITE(*,*)'G - ordinates, weights'
      DO 399 J=1,NG
        READ(LUN0,REC=IREC)DEL_G(J)
        WRITE(*,*)g_ord(j),del_g(j)
        IREC=IREC+1
399   CONTINUE
      IREC = 11 + 2*NG + 2
      WRITE(*,*)'Pressures : '
      DO 301 J=1,NP
        READ(LUN0,REC=IREC)PRESS1(J)
        WRITE(*,*)press1(j)
        PRESS1(J) = LOG(PRESS1(J))
        IREC = IREC + 1
301   CONTINUE
      WRITE(*,*)'Temperatures : '
      IF(NT.GT.0)THEN
       DO 302 J=1,NT
         READ(LUN0,REC=IREC)TEMP1(J)
         WRITE(*,*)temp1(j)
         IREC = IREC + 1
302    CONTINUE
      ELSE
       DO 307 I=1,NP
        DO 308 J=1,ABS(NT)
         READ(LUN0,REC=IREC)TEMP2(I,J)
C         WRITE(*,*)temp2(i,j)
         IREC = IREC + 1
308     CONTINUE
        WRITE(*,*)EXP(PRESS1(I)),(TEMP2(I,J),J=1,ABS(NT))
307    CONTINUE
      ENDIF

C     Read in central wavelengths if non-uniform grid
      IF(DELV.LE.0.0)THEN
       PRINT*,'Channel centres'
       print*,'Press a key to continue'
       READ(5,23)ANS
       DO 303 J=1,NPOINT
        READ(LUN0,REC=IREC)VCEN(J)
        PRINT*,J,VCEN(J)
        IREC=IREC+1
303    CONTINUE
      ENDIF

      IREC = IREC0

      CALL PROMPT('Enter wavelength range to set to zero : ')
      
      IF(DELV.GT.0)THEN
        CALL PROMPT('Enter wavenumber range to set tp zero : ')
        READ*,VV1,VV2
        N1 = 1 + NINT((VV1 - VMIN)/DELV)
        N2 = 1 + NINT((VV2 - VMIN)/DELV)
        WRITE(*,*)'Bin = ',N1,' Wavenumber = ',(VMIN + (N1-1)*DELV)
        WRITE(*,*)'Bin = ',N2,' Wavenumber = ',(VMIN + (N2-1)*DELV)
      ELSE
        CALL PROMPT('Enter channel range : ')
        READ*,N1,N2
      ENDIF

      NBIN=1+N2-N1
      X=0.
      IREC = IREC0 + NP*ABS(NT)*NG*(N1 - 1)
      DO 10 J=1,NBIN*NP*ABS(NT)*NG
            WRITE(LUN0,REC=IREC)X
            IREC = IREC + 1
10    CONTINUE

      CLOSE(LUN0)

      END
