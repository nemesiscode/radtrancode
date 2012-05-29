C           PROGRAM AGROUND
C     $Id: aground.f,v 1.4 2011-06-17 15:40:25 irwin Exp $
C-----------------------------------------------------------------------------
C_TITLE:  AGROUND  
C
C_ARGS:   none
C
C_KEYS:   ATMO,SPEC,VMS,PROG
C
C_CALLS:
C
C_BUGS:
C
C_DESC:	Adds ground emission times transmission to atmospheric emission.
C	Works with original Radtran .out files rather than intermediate 
C       ASCII files.
C
C	Pat Irwin	13/1/00	Original
C	Pat Irwin	26/4/12	Updated comments
C
C-----------------------------------------------------------------------------
C     note that include is not F77 and is included only during development
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C-----------------------------------------------------------------------------
      INTEGER MAXPT
      PARAMETER(MAXPT=150000)
      INTEGER IREC1,IREC2,IREC0,I,NPLOT,NINTER,NCHX,NCHY
      REAL VMAX,YMIN,YMAX,XMIN,XMAX,DX,DY
      REAL Y(MAXPT),BRIGHT,PLOG10,ERROR(MAXPT)
      double precision X(MAXPT),ddelv,dvmin
      character*20 tmp
      REAL YEM(MAXPT),YTR(MAXPT),GEMISS,TGROUND
      CHARACTER*100 XLABEL,YLABEL,FILE2
      CHARACTER Q
      CHARACTER*100 IPFILE,XFORM,YFORM
      INTEGER ISPACE
      LOGICAL FIRST,LOG,ASCOUT,EMISS,MICRON
      LOGICAL BRTEMP

      BRTEMP=.FALSE.
      EMISS=.TRUE.
      
      CALL PROMPT('driving file name?')
      READ(*,503)OPFILE
503   FORMAT(1A30)
      CALL FILE(OPFILE,IPFILE,'drv')
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
      CALL RDLBLD
C
C     check if want to look at averaged file
      CALL PROMPT('averaged?')
      READ(*,201)Q
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')THEN
        CALL FILE(OPFILE,FILE2,'ave')
       ELSE
        FILE2=OPFILE
      END IF
      I=2*ISYS()
      OPEN(UNIT=2,FILE=FILE2,STATUS='OLD',ACCESS='DIRECT',RECL=I)
      READ(2,REC=1)IREC0
      READ(2,REC=2)NPOINT
      READ(2,REC=3)VMIN
      READ(2,REC=4)DELV
      READ(2,REC=5)NPATH
      IF(NPOINT.GT.MAXPT)THEN
        WRITE(*,*)'NPOINT',NPOINT,'MAXPT',MAXPT
        WRITE(*,301)
301     FORMAT('Aground: NPOINT > MAXPT. Abort')
        STOP
      END IF

C-----------------------------------------------------
c	some fudge to fix the wavnumber scale for very small spacings! 	
C-----------------------------------------------------
      open(66,file='.aground.tmp',status='unknown')
      write(66,*) vmin,delv
      close(66)
	open(66,file='.aground.tmp',status='old')
      read(66,*) dvmin,ddelv
      close(66)
	print*,dvmin,ddelv
C-----------------------------------------------------

C
C-----------------------------------------------------
10    CONTINUE

      CALL DISCAL

      CALL PROMPT('Enter ID of atm emission calc : ')
      READ*,IEMCALC

      CALL PROMPT('Enter ID of atm transmission calc : ')
      READ*,ITRCALC

      IF(IEMCALC.LT.1.OR.ITRCALC.LT.1)THEN
       PRINT*,'ID out of range. Enter again'
       GOTO 10
      ENDIF

      IF(IEMCALC.GT.NCALC.OR.ITRCALC.GT.NCALC)THEN
       PRINT*,'ID out of range. Enter again'
       GOTO 10
      ENDIF

      CALL PROMPT('Enter ground temperature and emissivity : ')

      READ*,TGROUND,GEMISS

      VMAX=VMIN+(NPOINT-1)*DELV
      CALL PROMPT('plot brightness temperature? : ')
      READ(*,201)Q
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')BRTEMP=.TRUE.
      CALL PROMPT('Plot against wavenumber (1) or microns (2) : ')
      READ*,ISPACE
      IF(ISPACE.EQ.2)THEN
       MICRON=.TRUE.
      ELSE
       MICRON=.FALSE.
      END IF

      DO 110 I=1,NPOINT
      IREC1=IREC0+(I-1)*NPATH+ITRCALC-1
      IREC2=IREC0+(I-1)*NPATH+IEMCALC-1

      READ(2,REC=IREC1)YTR(I),ERROR(I)
      READ(2,REC=IREC2)YEM(I),ERROR(I)

      
cc      X(I)=VMIN+(I-1)*DELV
      X(I)=DVMIN+dble(I-1)*DDELV

      Y(I) = YEM(I) + GEMISS*PLANCK(real(X(I)),TGROUND)*YTR(I)

      IF(BRTEMP)THEN
        IF(Y(I).GT.0)THEN
         Y(I)=BRIGHT(real(X(I)),Y(I))
        ELSE
         Y(I)=0.
        END IF
      ELSE
       IF(EMISS.AND.MICRON)Y(I)=Y(I)*real(X(I)*X(I))/1E4
      END IF
      IF(MICRON)X(I)=1E4/X(I)
110   CONTINUE

      CALL PROMPT('log plot? [Y/N]')
      READ(*,201)Q
      CALL UPCASE(Q)
      IF(Q.EQ.'Y')THEN
        LOG=.TRUE.
        YMAX=-50.0
        YMIN=50.0
       ELSE
        LOG=.FALSE.
        YMAX=Y(1)
        YMIN=YMAX
      END IF
      XMAX=SNGL(X(1))
      XMIN=XMAX
      DO 304 I=1,NPOINT
      IF(LOG)Y(I)=PLOG10(Y(I))
      IF(Y(I).NE.-999.9)THEN
       YMAX=MAX(YMAX,Y(I))
       YMIN=MIN(YMIN,Y(I))
      ENDIF
      XMAX=MAX(XMAX,SNGL(X(I)))
      XMIN=MIN(XMIN,SNGL(X(I)))
304   CONTINUE
      IF(MICRON)THEN
       WRITE(*,184)XMIN,XMAX
      ELSE
       WRITE(*,204)XMIN,XMAX
      END IF
184   FORMAT(' wavelength ( x ) limits are: ',2F10.3)
204   FORMAT(' wavenumber ( x ) limits are: ',2F10.3)
      WRITE(*,305)YMIN,YMAX
305   FORMAT(' y limits are: ',2E12.5)
      IF(YMIN.EQ.YMAX)THEN
        Print*,'Ylimits are the same - no data'
        CALL PROMPT('input xmin,xmax,ymin,ymax')
        READ(*,*)XMIN,XMAX,YMIN,YMAX
      ELSE
        CALL PROMPT('force plotting axes? [y/n]')
        READ(*,201)Q
201     FORMAT(1A1)
        CALL UPCASE(Q)
        IF(Q.EQ.'Y')THEN
         CALL PROMPT('input xmin,xmax,ymin,ymax')
         READ(*,*)XMIN,XMAX,YMIN,YMAX
        ELSE
C        XMAX=VMAX
C        XMIN=VMIN
         YMAX=YMAX + (YMAX-YMIN)*0.05
        END IF
      END IF

C     *** start of debug lines ***
C      PRINT*,'ABSORB = ', ABSORB
C     *** end of debug lines ***

      IF(MICRON)THEN
       XLABEL='Wavelength (microns)'
      ELSE
       XLABEL='wavenumbers cm-1'
      END IF
      IF(EMISS)THEN
       IF(BRTEMP)THEN
        YLABEL='Brightness Temperature (K)'
       IF(LOG)YLABEL='Log10(Brightness Temperature (K))'
       ELSE IF(MICRON)THEN
        YLABEL='Radiance W cm-2 sr-1 micron-1'
        IF(LOG)YLABEL='Log10(Radiance W cm-2 sr-1 micron-1)'
       ELSE
        YLABEL='Radiance W cm-2 sr-1 cm'
        IF(LOG)YLABEL='Log10(Radiance W cm-2 sr-1 cm)'
       END IF
      END IF
      NPLOT=NPOINT
C------------------------------------------------------------------

C     ------------------------------------------------------------------------
C
C     plotting section

      ASCOUT=.TRUE.
      IF(ASCOUT)OPEN(12,FILE='plspec.idl',STATUS='UNKNOWN')
      FIRST=.TRUE.
      IF(ASCOUT)THEN
       WRITE(12,*)NPLOT
       WRITE(12,*)'XLIM,YLIM : '
       WRITE(12,*)XMIN,XMAX,YMIN,YMAX
       WRITE(12,*)'XLABEL, YLABEL : '
       WRITE(12,322)XLABEL
       WRITE(12,322)YLABEL
322    FORMAT(A)
323    FORMAT(1X,F16.6,' ',E16.6)
       DO 998 I=1,NPLOT
                IF(X(I).LT.XMIN)GOTO 998
                IF(X(I).GT.XMAX)GOTO 998
                WRITE(12,323)X(I),Y(I)
998    CONTINUE
      END IF
      IF(ASCOUT)CLOSE(12)
      STOP
      END
C------------------------------------------------------------------------------


      REAL FUNCTION PLOG10(Y)
      REAL Y

      IF(Y.GT.0)THEN
       PLOG10=LOG10(Y)
      ELSE
       PRINT*,'PLOG10. Trying to take log of ',Y
       PRINT*,'Setting to -999.9'
       PLOG10=-999.9
      ENDIF
      RETURN
      END
