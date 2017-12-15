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

      INTEGER PINDEX,CP,CT,LOOP,IPFORM

      INTEGER IREC,IREC0,CT2
      REAL P1,T1,PRESS1(MAXK),TEMP1(MAXK),VCEN(MAXBIN)
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),TABLE(MAXK,MAXK,MAXG)
      REAL TEMP2(MAXK,MAXK),TN(MAXK),TX,X1,X2
      REAL KOUT
      CHARACTER*200 KTAFIL,OPFILE1
      CHARACTER*1 ANS

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

304   PRINT*,'Output data for a particular wavelength (1)'
      PRINT*,'or output spectrum for particular P,T (2)'
      CALL PROMPT('?:')
      READ*,IPFORM
      IF(IPFORM.LT.1.OR.IPFORM.GT.2)GOTO 304  

      IF(IPFORM.EQ.1)THEN
       IF(DELV.GT.0)THEN
        CALL PROMPT('Enter wavenumber [cm^-1] : ')
        READ*,VV
        N1 = 1 + NINT((VV - VMIN)/DELV)
        WRITE(*,*)'Bin = ',N1,' Wavenumber = ',(VMIN + (N1-1)*DELV)
       ELSE
        CALL PROMPT('Enter channel number : ')
        READ*,N1
       ENDIF

       IREC = IREC0 + NP*ABS(NT)*NG*(N1 - 1)
       WRITE(*,*)'IREC = ',IREC
       DO 20 J=1,NP
        DO 30 K=1,ABS(NT)
          DO 40 LOOP=1,NG
            READ(LUN0,REC=IREC)TABLE(J,K,LOOP)
C            print*,j,k,loop,irec,table(j,k,loop)
            IREC = IREC + 1
40        CONTINUE
30      CONTINUE
20     CONTINUE



       WRITE(*,*)'Enter Pressure [atm] and Temperature [K]'
       READ*,P1,T1
       PMAX = PRESS1(NP)
       PMIN = PRESS1(1)
       P1 = LOG(P1)
       IF(P1.LT.PMIN)THEN
        WRITE(*,*)'**WARNING** P < PMIN ==> P, PMIN = ',P1,PMIN
        WRITE(*,*)' '
        WRITE(*,*)'Setting P equal to PMIN.'
        P1 = PMIN
       ENDIF
       IF(P1.GT.PMAX)THEN
        WRITE(*,*)'**WARNING** P > PMAX ==> P, PMAX = ',P1,PMAX
        WRITE(*,*)' '
        WRITE(*,*)'Setting P equal to PMAX.'
        P1 = PMAX
       ENDIF
       CP = PINDEX(PRESS1,NP,P1)
       IF(CP.LT.1)CP=1
       IF(CP.GE.NP)CP=NP-1
       V = (P1 - PRESS1(CP))/(PRESS1(CP+1) - PRESS1(CP))

       WRITE(*,*)'Pressure index CP,NP = ',CP,NP
       WRITE(*,*)'P,P(CP),P(CP+1)=',P1,PRESS1(CP),PRESS1(CP+1)
       WRITE(*,*)'(P - PRESS1(CP))/(PRESS1(CP+1) - PRESS1(CP)) = ',
     &   V

       IF(NT.GT.0)THEN
         TMAX = TEMP1(NT)
         TMIN = TEMP1(1)
         IF(T1.LT.TMIN)THEN
          WRITE(*,*)'**WARNING** T < TMIN ==> T, TMIN = ',T1,TMIN
          WRITE(*,*)' '
          WRITE(*,*)'Setting T equal to TMIN.'
          T1 = TMIN
         ENDIF
         IF(T1.GT.TMAX)THEN
          WRITE(*,*)'**WARNING** T > TMAX ==> T, TMAX = ',T1,TMAX
          WRITE(*,*)' '
          WRITE(*,*)'Setting T equal to TMAX.'
          T1 = TMAX
         END IF

         CT = PINDEX(TEMP1,NT,T1)
         IF(CT.LT.1)CT=1
         IF(CT.GE.NT)CT=NT-1

         WRITE(*,*)'temperature index CT,NT = ',CT,NT
         U = (T1 - TEMP1(CT))/(TEMP1(CT+1) - TEMP1(CT))
         WRITE(*,*)'T,T(CT),T(CT+1)=',T1,TEMP1(CT),TEMP1(CT+1)
         WRITE(*,*)'(T1 - TEMP1(CT))/(TEMP1(CT+1) - TEMP1(CT)) = ',
     &    U

       ELSE

         DO I=1,ABS(NT)
          TN(I)=TEMP2(CP,I)
         ENDDO
         TX=T1
         IF(TX.LT.TN(1))TX=TN(1)
         IF(TX.GT.TN(ABS(NT)))TX=TN(ABS(NT))
         CT = PINDEX(TN,ABS(NT),TX)
         IF(CT.LT.1)CT=1
         IF(CT.GE.ABS(NT))CT=ABS(NT)-1
         U = (TX - TN(CT))/(TN(CT+1) - TN(CT))

         WRITE(*,*)'temperature index CT,NT = ',CT,ABS(NT)
         WRITE(*,*)'T,T(CP,CT),T(CP,CT+1)=',TX,TEMP2(CP,CT),
     &		TEMP2(CP,CT+1)
         WRITE(*,*)'(TX-TEMP2(CP,CT))/(TEMP2(CP,CT+1)-
     &TEMP2(CP,CT))=',U

         TX=T1
         DO I=1,ABS(NT)
          TN(I)=TEMP2(CP+1,I)
         ENDDO
         IF(TX.LT.TN(1))TX=TN(1)
         IF(TX.GT.TN(ABS(NT)))TX=TN(ABS(NT))
         CT2 = PINDEX(TN,ABS(NT),TX)
         IF(CT2.LT.1)CT2=1
         IF(CT2.GE.ABS(NT))CT2=ABS(NT)-1
         U2 = (TX - TN(CT2))/(TN(CT2+1) - TN(CT2))

         WRITE(*,*)'temperature index CT2,NT = ',CT2,ABS(NT)
         WRITE(*,*)'T,T(CP+1,CT),T(CP+1,CT+1)=',TX,TEMP2(CP+1,CT2),
     &		TEMP2(CP+1,CT2+1)
         WRITE(*,*)'(TX-TEMP2(CP+1,CT2))/(TEMP2(CP+1,CT2+1)-
     &TEMP2(CP+1,CT2))=',U2

       ENDIF
       WRITE(*,*)' '
       WRITE(*,*)'    G_ORD    |      K_G     |   K_G*26850'

       IF(TABLE(CP,CT,NG).GT.0.0)THEN
         DO 80 LOOP=1,NG
          IF(NT.GT.0)THEN
            Y1 = LOG(TABLE(CP,CT,LOOP))
            Y2 = LOG(TABLE(CP+1,CT,LOOP))
            Y3 = LOG(TABLE(CP+1,CT+1,LOOP))
            Y4 = LOG(TABLE(CP,CT+1,LOOP))
            K_G(LOOP) = EXP((1.0 - V)*(1.0 - U)*Y1 + V*(1.0 - U)*Y2 +
     1    V*U*Y3 + (1.0 - V)*U*Y4)
          ELSE
            Y1 = LOG(TABLE(CP,CT,LOOP))
            Y2 = LOG(TABLE(CP+1,CT2,LOOP))
            Y3 = LOG(TABLE(CP+1,CT2+1,LOOP))
            Y4 = LOG(TABLE(CP,CT+1,LOOP))
            X1=(1.0-U)*Y1 + U*Y4
            X2=(1.0-U2)*Y2 + U2*Y3
            X = (1.0-V)*X1 + V*X2
            K_G(LOOP)=EXP(X)
          ENDIF
          WRITE(*,*)G_ORD(LOOP),K_G(LOOP),K_G(LOOP)*26850.0
80      CONTINUE
       ELSE
        DO 81 LOOP=1,NG
          K_G(LOOP) = 0.0
          WRITE(*,*)G_ORD(LOOP),K_G(LOOP),K_G(LOOP)*26850.0
81      CONTINUE
       ENDIF

       OPEN(12,FILE='read_table.dat',status='unknown')
        WRITE(12,*)NG,NP,ABS(NT)
        WRITE(12,*)G_ORD
        WRITE(12,*)PRESS1
        IF(NT.GT.0)THEN
         WRITE(12,*)TEMP1
        ELSE
         WRITE(12,*)TEMP2
        ENDIF
        WRITE(12,*)TABLE
       CLOSE(12)
       
      ELSE
  
c      write a simple opacity spectrum that can be plotted        
       OPEN(16,FILE='read_table_kspec.dat',status='unknown')
c      tabulated output         
       OPEN(12,FILE='read_table_wv.dat',status='unknown')
       
        WRITE(12,*)NPOINT
        IF(DELV.GT.0.0)THEN
         DO I = 1,NPOINT
          VCEN(I)=VMIN+(I-1)*DELV
         ENDDO
        ENDIF
        WRITE(12,*)(VCEN(I),I=1,NPOINT)
        WRITE(12,*)NG
        WRITE(12,*)(G_ORD(I),I=1,NG)

        WRITE(*,*)'Enter Pressure [atm] and Temperature [K]'
        READ*,P1,T1
        PMAX = PRESS1(NP)
        PMIN = PRESS1(1)
        P1 = LOG(P1)
        IF(P1.LT.PMIN)THEN
         WRITE(*,*)'**WARNING** P < PMIN ==> P, PMIN = ',P1,PMIN
         WRITE(*,*)' '
         WRITE(*,*)'Setting P equal to PMIN.'
         P1 = PMIN
        ENDIF
        IF(P1.GT.PMAX)THEN
         WRITE(*,*)'**WARNING** P > PMAX ==> P, PMAX = ',P1,PMAX
         WRITE(*,*)' '
         WRITE(*,*)'Setting P equal to PMAX.'
         P1 = PMAX
        ENDIF
        CP = PINDEX(PRESS1,NP,P1)
        
        IF(NT.GT.0)THEN
         TMAX = TEMP1(NT)
         TMIN = TEMP1(1)
        ELSE
         TMAX = TEMP2(CP,ABS(NT))
         TMIN = TEMP2(CP,1)
        ENDIF
        IF(T1.LT.TMIN)THEN
          WRITE(*,*)'**WARNING** T < TMIN ==> T, TMIN = ',T1,TMIN
          WRITE(*,*)' '
          WRITE(*,*)'Setting T equal to TMIN.'
         T1 = TMIN
        ENDIF
        IF(T1.GT.TMAX)THEN
         WRITE(*,*)'**WARNING** T > TMAX ==> T, TMAX = ',T1,TMAX
         WRITE(*,*)' '
         WRITE(*,*)'Setting T equal to TMAX.'
         T1 = TMAX
        END IF

        IF(NT.GT.0)THEN
         CT = PINDEX(TEMP1,NT,T1)
        ELSE
         DO I=1,ABS(NT)
          TN(I)=TEMP2(CP,I)
         ENDDO
         CT = PINDEX(TN,ABS(NT),T1)
        ENDIF
 
        WRITE(*,*)'Pressure index (CP), temperature index (CT) = ',CP,CT
        WRITE(*,*)'Tabulated log(pressure),temperature : '

        IF(NT.GT.0)THEN
         WRITE(*,*)PRESS1(CP),TEMP1(CT)
        ELSE
         WRITE(*,*)PRESS1(CP),TEMP2(CP,CT)
        ENDIF

        DO I=1,NPOINT

         IREC = IREC0 + NP*ABS(NT)*NG*(I - 1)
         DO 21 J=1,NP
          DO 31 K=1,ABS(NT)
           DO 41 LOOP=1,NG
            READ(LUN0,REC=IREC)TABLE(J,K,LOOP)
            IREC = IREC + 1
41         CONTINUE
31        CONTINUE
21       CONTINUE

         WRITE(12,*)(TABLE(CP,CT,LOOP),LOOP=1,NG)

c        integrate k over all ng points
	   kout = 0.
	   do LOOP=1,ng
	      kout = kout + TABLE(CP,CT,LOOP)*DEL_G(LOOP)
	   enddo
         write(16,*) vcen(i),kout

        ENDDO

       CLOSE(12)
       CLOSE(16)

      ENDIF

 
      END
C***********************************************************************
C***********************************************************************

