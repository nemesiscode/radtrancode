      SUBROUTINE CALCXMAP(IPFILE,NV,XMAP,HSWITCH,GNAME)
C     ****************************************************************
C     Subroutine to create matrix to map between user-defined gradient
C     variables and internal variables relating to the .prf profile.
C
C     Pat Irwin	1/7/03	Original version
C     Pat Irwin 29/2/12	Updated for Radtrans2.0
C
C     ****************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS),RHUM
      REAL H1(MAXPRO),P1(MAXPRO),T1(MAXPRO),VMR1(MAXPRO,MAXGAS)
      REAL XV(MAXPRO),MOLWT,LATITUDE,TEMP
      INTEGER NPRO,NVMR,ID(MAXGAS),ISO(MAXGAS),NPRO1,J,IFORM
      INTEGER IPLANET,I,N,K,NN,NCONT,NPARAM,NV,IV
      INTEGER IDV(MAXV),ICOL,COLV(MAXV),IPARAM,ILEN,IADD
      REAL ALTV(MAXV),XMAP(MAXV,MAXGAS+2+MAXCON,MAXPRO)
      REAL HMIN,HMAX,F
C     H is the height in kilometres above some NOMINAL zero.
C     P is the pressure in atmospheres (not bar).
C     T is the temperature in Kelvin.
C     VMR holds the NVMR volume mixing ratio profiles for each of the gases.
C     There are NPRO points in each profile.
C     ID and ISO hold the local identifier and isotope identifier for
C     each gas. Note that this program does not check that you only include
C     each gas once or that the identifiers are valid.
C
C----------------------------------------------------------------------------
      CHARACTER*40 GNAME(MAXV),TNAME
      LOGICAL DSWITCH,HSWITCH
      CHARACTER*100 IPFILE,MODFIL,DMODFIL,HMODFIL,BUFFER,TEXT
C----------------------------------------------------------------------------

      DO IV=1,MAXV
       DO IPARAM=1,MAXGAS+2+MAXCON
        DO I=1,MAXPRO
         XMAP(IV,IPARAM,I)=0.0
        ENDDO
       ENDDO
      ENDDO

1     FORMAT(A)

      DSWITCH=.FALSE.
      HSWITCH=.FALSE.
      CALL FILE(IPFILE,IPFILE,'pat')
      OPEN(UNIT=2,FILE=IPFILE,STATUS='OLD')
2     READ(2,1,END=9)TEXT
      CALL REMSP(TEXT)
      CALL UPCASE(TEXT)
      IF(TEXT(1:5).EQ.'MODEL')THEN
       MODFIL=TEXT(6:)
       CALL LOCASE(MODFIL)
      ELSE IF(TEXT(1:10).EQ.'DUST MODEL')THEN
       DMODFIL=TEXT(11:)
       CALL LOCASE(DMODFIL)
       DSWITCH=.TRUE.
      ELSE IF(TEXT(1:13).EQ.'FPARAH2 MODEL')THEN
       HMODFIL=TEXT(14:)
       CALL LOCASE(HMODFIL)
       HSWITCH=.TRUE.
      ENDIF
      GOTO 2
9     CONTINUE
      CLOSE(2)

      PRINT*,'T/P profile is : ',MODFIL

      IF(DSWITCH)THEN
       PRINT*,'Dust profile is : ',DMODFIL
      ELSE
       PRINT*,'No dust profile defined'
      ENDIF

      IF(HSWITCH)THEN
       PRINT*,'para-H2 profile is : ',HMODFIL
      ELSE
       PRINT*,'No para-H2 profile defined'
      ENDIF


      CALL FILE(MODFIL,MODFIL,'prf')
      OPEN(UNIT=1,FILE=MODFIL,STATUS='OLD')
C     First skip header
54    READ(1,1)BUFFER
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)IFORM
      READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
      DO 20 I=1,NVMR
        READ(1,*)ID(I),ISO(I)
20    CONTINUE
C     reading the first block of profiles
      READ(1,*)
      N=MIN(NVMR,3)
C     N is the maximum VMR which can be read in from the next block
      DO 30 I=1,NPRO
      IF(IFORM.EQ.0)THEN
         READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,N)
      ELSE IF(IFORM.EQ.1)THEN
         READ(1,*)H(I),T(I),(VMR(I,J),J=1,N)
      ELSE
         CALL WTEXT('invalid format')
         STOP
      END IF
30    CONTINUE
C     reading in additional blocks if any
C     N VMR profiles have been read in so far
33    IF(NVMR.GT.N)THEN
        READ(1,*)
C       profiles up to VMR(?,K) to be read from this block
        K=MIN(NVMR,(N+6))
        DO 32 I=1,NPRO
        READ(1,*)(VMR(I,J),J=N+1,K)
32      CONTINUE
        N=K
        GOTO 33
      END IF
      CLOSE(UNIT=1)

C     all processing below assumes that heights are in ascending order
C     so sorting just in case
      DO 12 J=1,NPRO
      DO 12 I=1,NPRO-1
       IF(ABS(H(I)-H(I+1)).LT.0.01)THEN
        WRITE(*,14)
14      FORMAT(' identical height values found')
C       STOP
       END IF
       IF(H(I).GT.H(I+1))THEN
        TEMP=H(I+1)
        H(I+1)=H(I)
        H(I)=TEMP
        TEMP=P(I+1)
        P(I+1)=P(I)
        P(I)=TEMP
        TEMP=T(I+1)
        T(I+1)=T(I)
        T(I)=TEMP
        DO 15 K=1,NVMR
         TEMP=VMR(I+1,K)
         VMR(I+1,K)=VMR(I,K)
         VMR(I,K)=TEMP
15      CONTINUE
       END IF
12    CONTINUE

      CALL FILE(DMODFIL,DMODFIL,'prf')
      OPEN(UNIT=3,FILE=DMODFIL,STATUS='OLD')
55     READ(3,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 55
       READ(BUFFER,*)NN, NCONT
      CLOSE(UNIT=3)


      WRITE(101)NPRO
101   FORMAT('There are ',I4,' levels')
      HMIN=H(1)
      HMAX=H(NPRO)
      PRINT*,'Minimum, maximum height = ',HMIN,HMAX
      WRITE(102)NVMR
102   FORMAT('There are ',I3,' gases')
      PRINT*,'Gases are : '
      DO I=1,NVMR
       PRINT*,I,ID(I),ISO(I)
      ENDDO
      IF(DSWITCH)THEN
       WRITE(6,103)NCONT
103    FORMAT('There are ',I3,' particle types')
      ENDIF

      CALL PROMPT('Enter number of gradients required : ')
      READ*,NV

      NPARAM=NVMR+1+NCONT
      IF(HSWITCH)NPARAM=NPARAM+1

      DO 99 IV=1,NV
       PRINT*,'Variable : ',IV 
       PRINT*,'Enter ID (choice of following)'
       WRITE(6,201)NVMR
201    FORMAT('1  to ',I2,' : gas v.m.r.')
       WRITE(6,202)NVMR+1
202    FORMAT(I2,'       : temperature')
       IF(DSWITCH)THEN
        WRITE(6,203)NVMR+2,NVMR+2+NCONT-1
203     FORMAT(I2,' to ',I2,' : aerosol opacity')
       ENDIF
       IF(HSWITCH)THEN
        WRITE(6,204)NPARAM-1,NPARAM
204     FORMAT(I2,' to ',I2,' : para-H2 fraction')
       ENDIF

       READ*,IDV(IV)
205    CALL PROMPT('Column sens. (0) or sens. at defined height(1) : ')
       READ*,ICOL
       IF(ICOL.LT.0.OR.ICOL.GT.1)GOTO 205
       COLV(IV)=ICOL
       ALTV(IV)=-999.9
       IF(ICOL.EQ.1)THEN
231     CALL PROMPT('Enter required altitude : ')
        READ*,ALTV(IV)
        IF(ALTV(IV).LT.HMIN.OR.ALTV(IV).GT.HMAX)THEN
         PRINT*,'Height is out of range. Must be between'
         PRINT*,HMIN,HMAX
         GOTO 231
        ENDIF
       ENDIF


       ILEN=1
       IF(COLV(IV).EQ.0)THEN
        TNAME=' Column '
        ILEN = 9
       ENDIF

       IF(IDV(IV).LE.NVMR)THEN
        WRITE(TEXT,301)ID(IDV(IV)),ISO(IDV(IV))
301     FORMAT(' Gas : ',I3,I3)
        IADD=13
       ELSEIF(IDV(IV).EQ.NVMR+1)THEN
        TEXT = ' Temperature'
        IADD=12
       ELSEIF(IDV(IV).GT.NVMR+1.AND.IDV(IV).LE.NVMR+1+NCONT)THEN
	WRITE(TEXT,302)IDV(IV)-(NVMR+1)
302     FORMAT(' Aerosol : ',I3)
        IADD=14
       ELSEIF(IDV(IV).EQ.NVMR+1+NCONT+1)THEN
        TEXT = ' Para-H2 fraction'
        IADD = 17
       ENDIF

       TNAME(ILEN:) = TEXT
       ILEN=ILEN+IADD

       IF(COLV(IV).EQ.1)THEN
        WRITE(TEXT,303)ALTV(IV)
303     FORMAT(' at ',F8.2,' km')
        TNAME(ILEN+1:) = TEXT
       ENDIF

       GNAME(IV)=TNAME
       WRITE(6,1)GNAME(IV)


99    CONTINUE

      DO IV=1,NV
       IF(COLV(IV).EQ.0)THEN
        DO I=1,NPRO
         XMAP(IV,IDV(IV),I)=1.0
        ENDDO
       ELSE
        J=0
        DO I=1,NPRO-1
         IF(H(I).LE.ALTV(IV).AND.H(I+1).GT.ALTV(IV))THEN
          J=I
          F=(ALTV(IV)-H(I))/(H(I+1)-H(I))
         ENDIF
        ENDDO
        IF(J.LT.1)THEN
          J=NPRO-1
          F=1.0
        ENDIF
        XMAP(IV,IDV(IV),J)=1.0-F
        XMAP(IV,IDV(IV),J+1)=F
       ENDIF      
      ENDDO

      RETURN

      END
      
