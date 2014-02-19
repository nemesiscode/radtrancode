      PROGRAM READ_TABLE
C     $Id: read_table.f,v 1.12 2007-06-28 15:29:12 irwin Exp $
C***********************************************************************
C_TITL:	READ_TABLE.f
C
C_DESC:	Reads an absorption coefficient look-up table
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

      INTEGER PINDEX,CP,CT,LOOP

      INTEGER IREC,IREC0
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK),VCEN(MAXBIN)
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),TABLE(MAXK,MAXK,MAXG)

      CHARACTER*100 KTAFIL,OPFILE1

C******************************** CODE *********************************

      CALL PROMPT('Enter input filename : ')
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
      DO 302 J=1,NT
        READ(LUN0,REC=IREC)TEMP1(J)
        WRITE(*,*)temp1(j)
        IREC = IREC + 1
302   CONTINUE
C     Read in central wavelengths if non-uniform grid
      IF(DELV.LT.0.0)THEN
       PRINT*,'Channel centres'
       DO 303 J=1,NPOINT
        READ(LUN0,REC=IREC)VCEN(J)  
        PRINT*,J,VCEN(J)  
        IREC=IREC+1
303    CONTINUE
      ENDIF

      IREC = IREC0

      IF(DELV.GT.0)THEN
       CALL PROMPT('Enter wavenumber [cm^-1] : ')
       READ*,VV
       N1 = 1 + NINT((VV - VMIN)/DELV)
       WRITE(*,*)'Bin = ',N1,' Wavenumber = ',(VMIN + (N1-1)*DELV)
      ELSE
       CALL PROMPT('Enter channel number : ')
       READ*,N1
      ENDIF

      IREC = IREC0 + NP*NT*NG*(N1 - 1)
      WRITE(*,*)'IREC = ',IREC
      DO 20 J=1,NP
        DO 30 K=1,NT
          DO 40 LOOP=1,NG
            READ(LUN0,REC=IREC)TABLE(J,K,LOOP)
            IREC = IREC + 1
40        CONTINUE
30      CONTINUE
20    CONTINUE



      WRITE(*,*)'Enter Pressure [atm] and Temperature [K]'
      READ*,P1,T1
      PMAX = PRESS1(NP)
      PMIN = PRESS1(1)
      P1 = LOG(P1)
      IF(P1.LT.PMIN)THEN
        WRITE(*,*)'**WARNING** P < PMIN ==> P, PMIN = ',P,PMIN
        WRITE(*,*)' '
        WRITE(*,*)'Setting P equal to PMIN.'
        P1 = PMIN
      ENDIF
      IF(P1.GT.PMAX)THEN
        WRITE(*,*)'**WARNING** P > PMAX ==> P, PMAX = ',P,PMAX
        WRITE(*,*)' '
        WRITE(*,*)'Setting P equal to PMAX.'
        P1 = PMAX
      ENDIF
      TMAX = TEMP1(NT)
      TMIN = TEMP1(1)
      IF(T1.LT.TMIN)THEN
        WRITE(*,*)'**WARNING** T < TMIN ==> T, TMIN = ',T,TMIN
        WRITE(*,*)' '
        WRITE(*,*)'Setting T equal to TMIN.'
        T1 = TMIN
      ENDIF
      IF(T1.GT.TMAX)THEN
        WRITE(*,*)'**WARNING** T > TMAX ==> T, TMAX = ',T,TMAX
        WRITE(*,*)' '
        WRITE(*,*)'Setting T equal to TMAX.'
        T1 = TMAX
      END IF

      CP = PINDEX(PRESS1,NP,P1)
      CT = PINDEX(TEMP1,NT,T1)

      WRITE(*,*)'Pressure index (CP), temperature index (CT) = ',CP,CT
      WRITE(*,*)' '
      IF(PRESS1(CP+1).NE.PRESS1(CP))THEN
       T = (P1 - PRESS1(CP))/(PRESS1(CP+1) - PRESS1(CP))
      ELSE
       T=0.
      ENDIF
      IF(TEMP1(CT+1).NE.TEMP1(CT))THEN
       U = (T1 - TEMP1(CT))/(TEMP1(CT+1) - TEMP1(CT))
      ELSE
       U=0.
      ENDIF
      WRITE(*,*)'(P - PRESS1(CP))/(PRESS1(CP+1) - PRESS1(CP)) = ',T
      WRITE(*,*)'(T - TEMP1(CT))/(TEMP1(CT+1) - TEMP1(CT)) = ',U
      WRITE(*,*)' '
      WRITE(*,*)'    G_ORD    |      K_G     |   K_G*26850'

      IF(TABLE(CP,CT,NG).GT.0.0)THEN
        DO 80 LOOP=1,NG
C          print*,EXP(PRESS1(CP)),26850.0*TABLE(CP,CT,LOOP),
C     1     EXP(PRESS1(CP+1)),26850.0*TABLE(CP+1,CT,LOOP)
          Y1 = LOG(TABLE(CP,CT,LOOP))
          Y2 = LOG(TABLE(CP+1,CT,LOOP))
          Y3 = LOG(TABLE(CP+1,CT+1,LOOP))
          Y4 = LOG(TABLE(CP,CT+1,LOOP))
          K_G(LOOP) = EXP((1.0 - T)*(1.0 - U)*Y1 + T*(1.0 - U)*Y2 +
     1    T*U*Y3 + (1.0 - T)*U*Y4)
          WRITE(*,*)G_ORD(LOOP),K_G(LOOP),K_G(LOOP)*26850.0
80      CONTINUE
      ELSE
        DO 81 LOOP=1,NG
          K_G(LOOP) = 0.0
          WRITE(*,*)G_ORD(LOOP),K_G(LOOP),K_G(LOOP)*26850.0
81      CONTINUE
      ENDIF

      OPEN(12,FILE='read_table.dat',form='unformatted', 
     1   status='unknown')
       WRITE(12)G_ORD
       WRITE(12)TABLE
      CLOSE(12)
        
      END
C***********************************************************************
C***********************************************************************

