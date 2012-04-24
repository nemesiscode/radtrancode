      SUBROUTINE LOADLINEBUFFER(NXTREC,VMIN,DELV,MAXLIN,NLINR,VLIN,
     1 SLIN,ALIN,ELIN,IDLIN,SBLIN,PSHIFT,DOUBV,TDW,TDWS,LLQ,NLIN,
     2 FSTLIN,LSTLIN,NGAS,IDGAS,ISOGAS,IFULL,V2)
C     $Id: 
C***********************************************************************
C_TITL:	LOADLINEBUFFER.f
C
C_DESC:	Reads in line data bins for program LBL. Reads in line data for
C	wavenumber interval vmin to vmin+delv, ignores other than the NGAS
C	gases with identifiers listed in IDGAS.
C	Only returns parameters used by LBL. i.e. not suitable for
C	non-lte calculations. Routine prompts for line data key if non
C	already set.
C
C_ARGS:	Input variables:
C	NXTREC		INT	Next record to read
C	VMIN		REAL	Lowest wavenumber of bin.
C	DELV		REAL	Size of bin in wavenumbers.
C	MAXLIN		INT	Maximum number of lines which can be held.
C	NLINR		INT	Maximum number of lines to read here.
C	VLIN(MAXLIN)	REAL	Line position [cm-1].
C	SLIN(MAXLIN)	REAL	Line strength [cm-1/molecule/cm2 at 296K]
C	ALIN(MAXLIN)	REAL	Lorentz width [cm-1].
C	ELIN(MAXLIN)	REAL	Lower state energy [cm-1].
C	IDLIN(MAXLIN)	INT	Gas identifier POSITION IN IDGAS ARRAY!!
C	SBLIN(MAXLIN)	REAL	Self broadening coefficient air-broadened
C				- self-broadened width.
C	TDW(MAXLIN)	REAL	Temperature coefficient of air-broadened
C				width.
C	TDWS(MAXLIN)	REAL	Temperature coefficient of self-broadened
C				width.
C	LLQ(MAXLIN)	CHARA*15 Lower state local quanta.
C	NLIN		INT	On exit holds the number of lines read in.
C	FSTLIN		INT	1st line position to load.
C	LSTLIN		INT	On exit holds last line position
C				(=FSTLIN-1 if none).
C	NGAS		INT	Number of gases of interest
C	IDGAS(NGAS)	INT	Local identifiers of the NGAS gases
C	ISOGAS(NGAS)	INT	Isotopic id. zero if all.
C
C_OUT:  Output variables:
C       IFULL		INT	Set to 1 if buffer full, 0 otherwise
C	V2		REAL	Wavenumber of last line read in.
C	NXTREC		INT	Next record to be read in.
C
C_FILE:	unit=1 (DBLUN) line data and gas files
C
C_CALL:	RDKEY         reads in details from the key file
C	REMSP         removes leading spaces from text string
C	RDGAS         reads in gas information
C	RDISO         reads in isotope information
C	FNDWAV        searches database for wavenumber entry
C	RDLINE        performs internal read on line data record
C
C_HIST:	20jun86	SBC	ORIGINAL VERSION
C	10jun87	SBC	Modified to include new HITRAN data base. Uses
C			code largely copied from HITRAN.FOR but not the
C			original subroutines which don't really fit the
C			program structure here.
C	17nov87	SBC	Removed RELABU from parameter list - now stored
C			in FORMAT.FOR Also ISOAFG defined in FORMAT now
C	18sep89	SBC	Switched to relative pointers to conform to D
C			Edwards format.
C	22feb91	SBC	Altered to use new ascii data base format.
C***************************** VARIABLES *******************************

      IMPLICIT NONE

C ../includes/dbcom.f stores the line database variables (e.g. RELABU).
      INCLUDE '../includes/dbcom.f' 

      INTEGER I,irec,LINE,MAXREC,MINREC,NXTREC,J
      INTEGER MAXLIN,FSTLIN,LSTLIN,NLIN,IFULL,NLINR
      INTEGER IDLIN(MAXLIN),NGAS,IDGAS(NGAS),ISOGAS(NGAS)

      REAL VMIN,DELV,VLIM,V2
      REAL VLIN(MAXLIN),SLIN(MAXLIN),ALIN(MAXLIN),ELIN(MAXLIN)
      REAL SBLIN(MAXLIN),TDW(MAXLIN),TDWS(MAXLIN),PSHIFT(MAXLIN)
      REAL DOUBV(MAXLIN)
      CHARACTER*15 LLQ(MAXLIN)
      CHARACTER*256 BUFFER
      LOGICAL WARNTD

      REAL TINW(MAXDGAS)
      DATA (TINW(J),J=1,42)/0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     $     0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     $     0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     $     0.5,0.5/


C******************************** CODE *********************************

      WARNTD=.TRUE.



C     Assume line database file is already open.

     

      VLIM= VMIN + DELV

      LINE= FSTLIN
      NLIN = 0
      IFULL = 0
      IREC=NXTREC-1

111   IREC=IREC+1

      READ(DBLUN,110,REC=IREC)BUFFER(1:DBRECL)
110   FORMAT(A)
C     WRITE(6,110)BUFFER(1:DBRECL)
C     print*,'Line ID : ',LNID,LOCID(LNID)
      CALL RDLINE(BUFFER)
      IF(LNWAVE.LT.VLIM)THEN
         DO 114 I=1,NGAS
C NOTE: that a line can be included more than once
          IF(IDGAS(I).EQ.LOCID(LNID))THEN
            IF(ISOGAS(I).EQ.0)THEN
              SLIN(LINE)= LNSTR
            ELSE
              IF(LNISO.NE.DBISO(ISOGAS(I),LOCID(LNID)))GOTO 114
              SLIN(LINE)= LNSTR/RELABU(ISOGAS(I),IDGAS(I))
C              WRITE(*,*)' dbiso = ',DBISO(ISOGAS(I),LOCID(LNID))
C              WRITE(*,*)' relabu = ',RELABU(ISOGAS(I),IDGAS(I))
            ENDIF
            IDLIN(LINE)= I
            VLIN(LINE)= LNWAVE
            IF(LNWIDA.GT.0)THEN
C Set the air-broadened halfwidth equal to the Lorentz halfwidth
              ALIN(LINE)= LNWIDA
            ELSE
          WRITE(*,*)' LOADLINEBUFFER.f :: Warning: Lorentz halfwidth=0.'
              WRITE(*,*)' Setting to 0.075.'
              ALIN(LINE)= 0.075
            ENDIF
            PSHIFT(LINE)= LNPSH
            DOUBV(LINE)= LDOUBV
            TDW(LINE)= LNTDEP
            TDWS(LINE)= LNTDEPS
            IF(DBRECL.EQ.160)THEN
             LLQ(LINE)=LNLLQ04
            ELSE
             LLQ(LINE)(1:9)= LNLLQ
             LLQ(LINE)(10:15)='      '
            ENDIF
            IF(LNWIDS.GT.1.E-20)THEN
C NOTE: SBLIN is the correction to air broadening so that zero is valid
              SBLIN(LINE)= LNWIDA - LNWIDS
            ELSE
              SBLIN(LINE)= 0.0
            ENDIF
            ELIN(LINE)= LNLSE
            LINE= LINE + 1
            NLIN=NLIN+1
C            WRITE(*,*)LNWAVE,LNSTR,LNWIDA,LNLSE
            IF(LINE.GT.MAXLIN)THEN	!Loop round
              LINE=1
            ENDIF
            IF(NLIN.GT.NLINR)THEN
              WRITE(*,*)' LOADLINEBUFFER.f ::  NLIN > NLINR.'
              WRITE(*,*)' '
              WRITE(*,*)' Line buffer is full'
              WRITE(*,*)' '
              WRITE(*,*)' MAXLIN,NLINR = ',MAXLIN,NLINR
              WRITE(*,*)' Wavelength of last line read: ',VLIN(LINE)
              IFULL=1
              GOTO 112
            ENDIF
          ENDIF
114      CONTINUE
         GOTO 111
      ENDIF
112   V2=VLIN(LINE)
      LSTLIN= LINE - 1
      NXTREC=IREC

      DO 51 J=1,NLIN
        IF(TDW(J).LT.1.E-30)THEN
          TDW(J)= TINW(IDLIN(J))
          IF(WARNTD)THEN
            WRITE(*,*)' LOADLINEBUFFER :: Warning: setting default line'
            WRITE(*,*)' widths.'
            WARNTD= .FALSE.
          ENDIF
        ENDIF

        IF(TDWS(J).LT.1.E-30)THEN
          TDWS(J)= TINW(IDLIN(J))
          IF(WARNTD)THEN
            WRITE(*,*)' LOADLINEBUFFER :: Warning: setting default line'
            WRITE(*,*)' widths.'
            WARNTD= .FALSE.
          ENDIF
        ENDIF
51    CONTINUE

      RETURN

      END
