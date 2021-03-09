      PROGRAM CONV_SPEC
C     $Id: conv_spec.f,v 1.4 2011-06-17 15:40:25 irwin Exp $
C-----------------------------------------------------------------------------
C_TITLE:  CONV: to average output of GENLBL
C
C_ARGS:   none
C
C_KEYS:   ATMO,SPEC,VMS,PROG
C
C_DESCR:  
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:  15may90  SBC ORIGINAL VERSION derived from outspec
C        17may90  SBC replaced triangular integration with two square ones
C        25aug92  SBC mod to allow multiple paths
C                     included filters
C-----------------------------------------------------------------------------
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C-----------------------------------------------------------------------------
      INTEGER MAXPT
      PARAMETER(MAXPT=150000)
      CHARACTER*100 FILE1,IPFILE
      REAL OUTPUT(MAXPT),ERROR(MAXPT)
      REAL Y(MAXPT)
      INTEGER I,K,L,IREC,IREC0,ISPACE
      LOGICAL BIT
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      idiag=1

      CALL PROMPT('driving file name? :')
      READ(5,503)OPFILE
503   FORMAT(A)
      CALL FILE(OPFILE,IPFILE,'drv') 
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
      CALL RDLBLD
C
      I=2*ISYS()
      OPEN(UNIT=2,FILE=OPFILE,STATUS='OLD',ACCESS='DIRECT',
     1RECL=I)
      READ(2,REC=1)IREC0
      READ(2,REC=2)NPOINT
      READ(2,REC=3)VMIN
      READ(2,REC=4)DELV
      READ(2,REC=5)NPATH
      IF(NPOINT.GT.MAXPT)THEN
        PRINT*,'Conv_Spec. MAXPT not big enough',NPOINT,MAXPT
        STOP
      END IF
C
      PRINT*,'Enter both the required FWHM, and the'
      CALL PROMPT('Wavespace [0=wavenumber, 1=microns]:')
      READ(*,*)WID,ISPACE
      PRINT*,'Enter ICONV [1=square, 2=triangle,'
      CALL PROMPT('3=gauss or 4=sinc2 ]:')
      READ(*,*)ICONV
      CALL PROMPT('Enter FWHM of existing spectrum :')
      READ(*,*)FWHM


C     reading in the data
      IREC=IREC0
      DO 200 I=1,NPATH*NPOINT
      READ(2,REC=IREC)OUTPUT(I),ERROR(I)
      IREC=IREC+1
200   CONTINUE
C


      DO 140 K=1,NPATH
      print*,ITYPE(K)
      IF(ITYPE(K).EQ.256.OR.
     &  (ITYPE(K).LT.128.AND.BIT(1,ITYPE(K))) ) THEN
       ISPEC = 1
       Print*,'Path ',K,'Radiance'
      ELSE
       ISPEC = 0
       PRINT*,'Path ',K,'Transmission/absorption'
      ENDIF
      DO 110 I=1,NPOINT
      L=K+(I-1)*NPATH
      Y(I)=OUTPUT(L)
110   CONTINUE

      CALL NEWCONV(NPOINT,Y,ICONV,ISPACE,WID,VMIN,DELV,FWHM,
     1ISPEC)
      DO 120 I=1,NPOINT
      L=K+(I-1)*NPATH
      OUTPUT(L)=Y(I)
120   CONTINUE

      DO 131 I=1,NPOINT
      L=K+(I-1)*NPATH
      Y(I)=ERROR(L)*ERROR(L)
131   CONTINUE
      CALL NEWCONV(NPOINT,Y,ICONV,ISPACE,WID,VMIN,DELV,FWHM,
     1ISPEC)
      DO 135 I=1,NPOINT
      L=K+(I-1)*NPATH
      ERROR(L)=SQRT(Y(I))
135   CONTINUE


140   CONTINUE


C
C     output the header
      CALL FILE(OPFILE,FILE1,'ave')
      CLOSE(UNIT=2)
      I=2*ISYS()
      OPEN(UNIT=2,FILE=FILE1,STATUS='UNKNOWN',ACCESS='DIRECT',
     1RECL=I)
      WRITE(2,REC=1)IREC0
      WRITE(2,REC=2)NPOINT
      WRITE(2,REC=3)VMIN
      WRITE(2,REC=4)DELV
      WRITE(2,REC=5)NPATH
C     and the data
      IREC=IREC0
      DO 201 I=1,NPATH*NPOINT
      WRITE(2,REC=IREC)OUTPUT(I),ERROR(I)
      IREC=IREC+1
201   CONTINUE
C
      STOP
      END
      
