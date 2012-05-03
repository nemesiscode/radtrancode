      PROGRAM LISPEC
C     $Id: li_spec.f,v 1.3 2011-06-17 15:40:27 irwin Exp $
C-----------------------------------------------------------------------------
C_TITLE:  OUTSPEC: to list output of GENLBL calculations to file or terminal
C
C_ARGS:   none
C
C_KEYS:   ATMO,SPEC,VMS,PROG
C
C_DESCR:  lists any one of the spectra in a file produced by MAINLBL
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:   5feb87  SBC ORIGINAL VERSION
C        15may90  SBC modified order of parameters in output array
C-----------------------------------------------------------------------------
      INTEGER MAXPT
      PARAMETER (MAXPT=150001)
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INTEGER I,J,K,L,IREC,IREC0,NOUT,LUN
      LOGICAL ABSORB
      CHARACTER Q
      CHARACTER*100 NAME,IPFILE,FILE2
      REAL OUTPUT(MAXPT),ERROR(MAXPT),V
      CALL PROMPT('driving file name? :')
      READ(*,503)OPFILE
503   FORMAT(A)
      CALL FILE(OPFILE,IPFILE,'drv') 
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
      CALL RDLBLD
C
C     check if want to look at averaged file
      CALL PROMPT('averaged?')
      READ(*,201)Q
201   FORMAT(1A1)
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')THEN
        CALL FILE(OPFILE,FILE2,'ave')
       ELSE
        FILE2=OPFILE
       END IF
      I=2*ISYS()
      print*,'Record Length = ',I
      OPEN(UNIT=2,FILE=FILE2,STATUS='OLD',ACCESS='DIRECT',RECL=I) 
      READ(2,REC=1)IREC0
      READ(2,REC=2)NPOINT
      READ(2,REC=3)VMIN
      READ(2,REC=4)DELV
      READ(2,REC=5)NPATH

      PRINT*,'header: IREC0,NPOINT,VMIN,DELV,NPATH',
     1  IREC0,NPOINT,VMIN,DELV,NPATH

      IF(NPOINT.GT.MAXPT)THEN
        WRITE(*,301)
301     FORMAT(' not enough space')
        STOP
        END IF
C
C-----------------------------------------------------
      NOUT=NPATH*NPOINT
      IREC=IREC0
      CALL PROMPT('output to file? [Y/N]')
      READ(*,4)Q
4     FORMAT(1A1)
      LUN=6
      IF(Q.EQ.'Y'.OR.Q.EQ.'y')THEN
        CALL PROMPT('filename?')
        READ(*,2)NAME
2       FORMAT(A)
        OPEN(UNIT=3,FILE=NAME,STATUS='NEW')
        LUN=3
        END IF
      CALL PROMPT('output 1-calc ?')
      READ(*,151)Q
151   FORMAT(1A1)
      ABSORB=.FALSE.
      IF(Q.EQ.'Y'.OR.Q.EQ.'y')ABSORB=.TRUE.
      DO 200 I=1,NOUT
      READ(2,REC=IREC)OUTPUT(I),ERROR(I)
      IF(ABSORB)OUTPUT(I)=1.-OUTPUT(I)
      IREC=IREC+1
200   CONTINUE
      DO 210 I=1,NPOINT
      K=(I-1)*NPATH
      V=VMIN+(I-1)*DELV
      L=MIN((K+5),(K+NPATH))
      WRITE(LUN,211)V,(OUTPUT(J),J=K+1,L)
      IF((NPATH).GT.5)WRITE(LUN,212)(OUTPUT(J),J=K+6,K+NPATH)
211   FORMAT(1X,F13.6,' ',5E13.6)
212   FORMAT(12X,5E13.6)
210   CONTINUE
      STOP
      END

