      PROGRAM PLOTLI
C     $Id: pl_lines.f,v 1.5 2011-06-17 15:53:02 irwin Exp $
C---------------------------------------------------------------------------
C
C_TITLE: PLOT_LINES: plots out line details for a particular spectral region
C
C_ARGS:  None.
C
C_KEYS:  PROG, VMS, CAL .
C
C_DESCR: Prompts for line data key, gas identifiers, wavelength range and
C        plotting details. Plots out line strengths versus wavenumber.
C
C_FILES: unit DBLUN line data and gas files
C
C_CALLS:  LINES         reads line data
C         RDKEY         reads in details from the key file
C         PROMPT        user prompt utility
C         REMSP         removes leading spaces from text string
C         RDGAS         reads in gas information
C         RDISO         reads in isotope information
C         FNDWAV        searches database for wavenumber entry
C         RDLINE        performs internal read on line data record
C         ASKYN         prompts user for yes/no answer
C         GOBLIN graphics routines
C
C_BUGS:
C
C_HIST:
C        originally was ISO.FOR used for mars
C        plots out line locations for particular gases and isotopes
C        modified 15jun87 to use routine LINES as used by GENLBL so that
C        current line data always available
C        18sep89 SBC added standard header, renamed from LINEPLT
C
C_END:
C
C--------------------------------------------------------------------------
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
      INTEGER LINLIM
      PARAMETER (LINLIM=10000000)
      INTEGER I,L,ISOGAS,IDGAS,ILAST,NLIN,NGAS,IFIRST,MAXLIN
      REAL YMIN,YMAX,DX,DY
      DOUBLE PRECISION VMIN,VMAX,DELV,DVMIN,DVLIN,VLIN(LINLIM)
      CHARACTER*100 BUFFER
      INTEGER ISCALE,IXFORM,IYFORM,IPTF
      CHARACTER*12 XFORM,YFORM
      REAL SCALE,PARTF,TCORS1,TCORS2,TS1,TS2,TSTIM,DPEXP
      REAL ALIN(LINLIM),ELIN(LINLIM),
     1 SBLIN(LINLIM),TDW(LINLIM),TDWS(LINLIM),PSHIFT(LINLIM),
     2 DOUBV(LINLIM)
      DOUBLE PRECISION SLIN(LINLIM),LNABSCO
      CHARACTER*15 LLQ(LINLIM)
      INTEGER IDLIN(LINLIM)
C     variables relating to gases
      LOGICAL TDYPLT,LINEAR,POINTS,ASKYN
C     TDYPLT if true then only plots strongest line in any wavenumber
C     region smaller than 1/2000th of the entire plot
C     This is to avoid  multiple lines working their way through the paper
C     on pen plotters
C     LINEAR is true for a linear scale in strength
C     POINTS is true plots only a point for each line
C     POINTS true and TDYPLT true isn't sensible
      DATA TDYPLT/.FALSE./
      DATA POINTS/.FALSE./
      DATA LINEAR/.FALSE./
C
      CALL PROMPT('data base key?')
      READ(*,102)KEYFIL
102   FORMAT(A)
      CALL RDKEY(2)
      CALL RDGAS
      CALL RDISO
      CALL PROMPT('gas id, isotope code?')
      READ(*,*)IDGAS,ISOGAS
      
      IPTF=0
      CALL PROMPT('Enter calculation temperature (296 default) : ')
      READ(*,*)TEMP
      
      TCORS1 = PARTF(IDGAS,ISOGAS,TEMP,IPTF)
      TCORS2 = 1.439*(TEMP - 296.)/(296.*TEMP)

      LINEAR=ASKYN('linear plot? [Y/N]')
      IF(POINTS)TDYPLT=.FALSE.
      CALL PROMPT('vmin,vmax?')
      READ(*,*)VMIN,VMAX
      DELV=VMAX-VMIN
C     DVMIN is the minimum spacing between plotted lines to avoid overwriting
C     and very long plots
      DVMIN=DELV/2000
      MAXLIN=LINLIM
      IFIRST=1
      NGAS=1
      PRINT*,'MAXLIN=',MAXLIN
      CALL LINESS(VMIN,DELV,MAXLIN,VLIN,SLIN,ALIN,ELIN,IDLIN,SBLIN,
     1PSHIFT,DOUBV,TDW,TDWS,LLQ,NLIN,IFIRST,ILAST,NGAS,IDGAS,ISOGAS)
      WRITE(*,110)NLIN
110   FORMAT(' found ',I8,' lines')
      IF(NLIN.EQ.0)STOP

      DO 207 I=1,NLIN

       TS1 = (1.0 - DPEXP(-1.439*VLIN(I)/TEMP))
       TS2 = (1.0 - DPEXP(-1.439*VLIN(I)/296.0))
       TSTIM=1.0
       IF(TS2.NE.0.)TSTIM=TS1/TS2
C       print*,SLIN(I)
       LNABSCO=LOG(SLIN(I))+LOG(TCORS1)+TCORS2*ELIN(I)+LOG(TSTIM)
       SLIN(I)=DEXP(LNABSCO)
C       print*,LNSABSCO,LOG(TCORS1),TCORS2*ELIN(I),LOG(TSTIM)
C       print*,SLIN(I)
C       stop
207   CONTINUE

      IF(LINEAR)THEN
        DO 203 I=1,NLIN
        YMAX=SNGL(MAX(DBLE(YMAX),SLIN(I)))
203     CONTINUE
        ISCALE=INT(1000+LOG10(YMAX))-1000
        SCALE=10.**ISCALE
        YMAX=YMAX*1.1/SCALE
        YMIN=0.
        ISCALE=ISCALE-47
       ELSE
        YMAX=-17.
        YMIN=-30.
        DO 204 I=1,NLIN
         IF(SLIN(I).GT.0)THEN
          YMIN=MIN(YMIN,FLOAT(INT(1000+LOG10(SLIN(I)))-1047))
         END IF
204     CONTINUE
        DY=1.+FLOAT(INT((YMAX-YMIN)/14))
        END IF
C
      call prompt('Wavenumber (1) or wavelength (2) : ')
      read*,iwave

      if (iwave.eq.2)then
       x1 = 1e4/vmax
       x2 = 1e4/vmin
       vmin = x1
       vmax = x2
       do i=1,nlin
        vlin(i)=1e4/vlin(i)
       end do
       print*,vmin,vmax
       dvmin = 0.0
      end if

      OPEN(23,FILE='pl_lines.asc',STATUS='unknown')

      if(iwave.eq.2)then
       WRITE(23,51)
      else
       WRITE(23,50)
      end if
50    FORMAT('wavenumbers (cm-1)')
51    FORMAT('wavelength (microns)')

      IF(ISOGAS.NE.0)THEN
      WRITE(23,44)IDGAS,GASNAM(IDGAS),ISOGAS,
     1DBISO(ISOGAS,IDGAS),NLIN
44      FORMAT('gas:',I2,' ( ',1A6,' ), Isotope:',I5,'  (',
     1      I4,'),  ',I8,' lines')
      ELSE
      WRITE(23,47)IDGAS,GASNAM(IDGAS),ISOGAS,NLIN
47      FORMAT('gas:',I2,' ( ',1A6,' ), Isotope:',I5,'  (ALL) ,  '
     1      ,I8,' lines')
      END IF

      WRITE(23,46)DBNAME
46      FORMAT('data base:',A)

      IF(ISOGAS.NE.0)THEN
        WRITE(23,45)
45      FORMAT('note: strengths corrected for atmospheric abundance')
      END IF
C
      IF(TDYPLT)THEN
        L=1
        DO 300 I=2,NLIN
        DVLIN=ABS(VLIN(I)-VLIN(L))
        IF(DVLIN.LT.DVMIN)THEN
          SLIN(L)=MAX(SLIN(I),SLIN(L))
          GOTO 300
        ELSE
          L=L+1
          SLIN(L)=SLIN(I)
          VLIN(L)=VLIN(I)
        END IF
300        CONTINUE
      ELSE
        L=NLIN
      END IF

      WRITE(23,*)L
      DO 200 I=1,L
      IF(.NOT.LINEAR)THEN
        IF(SLIN(I).GT.0.)THEN
          SLIN(I)=LOG10(SLIN(I))-47
        ELSE
          GOTO 200
        END IF
       ELSE
         SLIN(I)=SLIN(I)*1e-20
         SLIN(I)=SLIN(I)*1e-27
C        SLIN(I)=SLIN(I)/DBLE(SCALE)
      END IF
      WRITE(23,*)VLIN(I),SNGL(SLIN(I))
200   CONTINUE

      CLOSE(23)

      STOP
      END
