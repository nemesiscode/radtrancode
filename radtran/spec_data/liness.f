      SUBROUTINE LINESS(VMIN,DELV,MAXLIN,VLIN,SLIN,ALIN,ELIN,IDLIN,
     1 SBLIN,PSHIFT,DOUBV,TDW,TDWS,LLQ,NLIN,FSTLIN,LSTLIN,
     2 NGAS,IDGAS,ISOGAS)
C     $Id: liness.f,v 1.6 2011-06-17 14:58:50 irwin Exp $
C***********************************************************************
C_TITL:	LINESS.f
C
C_DESC:	Reads in line data bins for program LBL. Reads in line data for
C	wavenumber interval vmin to vmin+delv, ignores other than the NGAS
C	gases with identifiers listed in IDGAS.
C	Only returns parameters used by LBL. i.e. not suitable for
C	non-lte calculations. Routine prompts for line data key if non
C	already set.
C
C_ARGS:	Input variables:
C	VMIN		REAL	Lowest wavenumber of bin.
C	DELV		REAL	Size of bin in wavenumbers.
C	MAXLIN		INT	Maximum number of lines which can be held.
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

      INTEGER I,irec,LINE,MAXREC,MINREC
      INTEGER MAXLIN,FSTLIN,LSTLIN,NLIN
      INTEGER IDLIN(MAXLIN),NGAS,IDGAS(NGAS),ISOGAS(NGAS)

      REAL VMIN,DELV,VLIM
      REAL VLIN(MAXLIN),ALIN(MAXLIN),ELIN(MAXLIN)
      REAL SBLIN(MAXLIN),TDW(MAXLIN),TDWS(MAXLIN),PSHIFT(MAXLIN)
      REAL DOUBV(MAXLIN)
      DOUBLE PRECISION SLIN(MAXLIN)
      CHARACTER*15 LLQ(MAXLIN)
      CHARACTER*256 BUFFER

C******************************** CODE *********************************

C      WRITE(*,*)' LINESS.f :: vmin, delv = ',vmin,delv
C      WRITE(*,*)' LINESS.f :: maxlin, fstlin= ',maxlin,fstlin
C      WRITE(*,*)' LINESS.f :: ngas = ',ngas
C      DO i=1,ngas
C        WRITE(*,*)' LINESS.f :: idgas, isogas = ',idgas,isogas
C      ENDDO
C      WRITE(*,*)' LINESS.f :: dblun = ',DBLUN
C      WRITE(*,*)' LINESS.f :: keyfil = ',KEYFIL

C First open database .key file
      IF(KEYFIL(1:3).EQ.'   ')THEN
        WRITE(*,*)' LINESS.f :: Reading .key file again.'
        CALL RDKEY(1)
        CALL RDGAS
        CALL RDISO
      ENDIF
      OPEN(UNIT=DBLUN,FILE=DBFILE,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='DIRECT',RECL=DBRECL)

      CALL FNDWAV(VMIN)
      MINREC= DBREC
      VLIM= VMIN + DELV
      CALL FNDWAV(VLIM)
      MAXREC= DBREC - 1
      WRITE(*,*)' LINESS.f :: MINREC,MAXREC = ',MINREC,MAXREC

      LINE= FSTLIN
      DO 111 IREC=MINREC,MAXREC
        READ(DBLUN,110,REC=IREC)BUFFER(1:DBRECL)
110     FORMAT(A)
C        WRITE(6,110)BUFFER(1:DBRECL)
C        print*,'Line ID : ',LNID,LOCID(LNID)
        CALL RDLINE(BUFFER)

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
              WRITE(*,*)' LINESS.f :: Warning: Lorentz halfwidth=0.'
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
C            WRITE(*,*)LNWAVE,LNSTR,LNWIDA,LNLSE
            IF(LINE.GT.MAXLIN)THEN	!If lines run past end stop
              WRITE(*,*)' LINESS.f :: Error: LINE > MAXLIN.'
              WRITE(*,*)' Stopping program.'
              WRITE(*,*)' '
              WRITE(*,*)' Either reduce the spectra range or increase'
              WRITE(*,*)' MAXLIN.'
              WRITE(*,*)' '
              WRITE(*,*)' LINE, MAXLIN = ',LINE,MAXLIN
              WRITE(*,*)' Wavelength of last line read: ',VLIN(MAXLIN)
              STOP
            ENDIF
          ENDIF
114     CONTINUE
111   CONTINUE
      NLIN= LINE - FSTLIN
      LSTLIN= LINE - 1

C      print*,fstlin,lstlin
      CLOSE(UNIT=DBLUN)

      RETURN

      END
