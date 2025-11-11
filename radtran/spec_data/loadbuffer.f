      SUBROUTINE LOADBUFFER(VMIN,VMAX,FSTREC,MAXLIN,MAXBIN,IB,NGAS,
     1 IDGAS,ISOGAS,VBOT,WING,NLINR,VLIN,SLIN,ALIN,ELIN,IDLIN,SBLIN,
     2 PSHIFT,DOUBV,TDW,TDWS,LLQ,NXTREC,FSTLIN,LSTLIN,LASTBIN)

C     $Id:
C***********************************************************************
C_TITL:	LOADBUFFER.f
C
C_DESC:	Reads in line data buffer. Reads in a maximum of MAXLIN lines
C	within wavenumber interval vmin to vmax, ignoring lines for gases
C       other than the NGAS gases with identifiers listed in IDGAS.
C
C_ARGS:	Input variables:
C	VMIN		REAL*8	Lowest wavenumber of region of interest
C	VMAX		REAL*8	Highest wavenumber of region of interest
C	FSTREC		INT     First record in line database to query 
C	MAXLIN		INT	Max number of lines to read in
C	MAXBIN		INT 	Max number of continuum bins
C	IB		INT	Buffer to fill (1 or 2)
C	NGAS		INT	Number of gases of interest
C	IDGAS(NGAS)	INT	Local identifiers of the NGAS gases
C	ISOGAS(NGAS)	INT	Isotopic ID (zero if all)
C	VBOT		REAL	Bottom wavenumber of continuum bins
C	WING		REAL	Width of continuum bins.
C
C_ARGS: Output variables:
C	NLINR		INT	Number of lines actually read in
C	VLIN(2,MAXLIN)	REAL*8	Line position [cm-1].
C	SLIN(2,MAXLIN)	REAL	Line strength [cm-1/molecule/cm2 at 296K]
C	ALIN(2,MAXLIN)	REAL	Lorentz width [cm-1].
C	ELIN(2,MAXLIN)	REAL	Lower state energy [cm-1].
C	IDLIN(2,MAXLIN)	INT	Gas identifier POSITION IN IDGAS ARRAY!!
C	SBLIN(2,MAXLIN)	REAL	Self broadening coefficient air-broadened
C				- self-broadened width.
C	TDW(2,MAXLIN)	REAL	Temperature coefficient of air-broadened
C				width.
C	TDWS(2,MAXLIN)	REAL	Temperature coefficient of self-broadened
C				width.
C	LLQ(2,MAXLIN)	CHARA*15 Lower state local quanta.
C	NXTREC		INTEGER	Next record in database to query
C	FSTLIN(2,MAXBIN) INT	First line read into in each continuum bin
C	LSTLIN(2,MAXBIN) INT	Last line read into in each continuum bin
C	LASTBIN(2)        INT	Bin position of last line in buffer
C
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
C_HIST:	15apr11	PGJI	ORIGINAL VERSION
C***************************** VARIABLES *******************************

      IMPLICIT NONE

C ../includes/dbcom.f stores the line database variables (e.g. RELABU).
      INCLUDE '../includes/dbcom.f' 

      INTEGER I,IREC,LINE,FSTREC,NXTREC,J,IB
      INTEGER MAXLIN,NLINR,MAXBIN,CURBIN
      INTEGER IDLIN(2,MAXLIN),NGAS,IDGAS(NGAS),ISOGAS(NGAS)
      INTEGER FSTLIN(2,MAXBIN),LSTLIN(2,MAXBIN),NBINX
      INTEGER LASTBIN(2),FIRSTBIN(2),NBINY
      DOUBLE PRECISION VMIN,VMAX
      DOUBLE PRECISION VLIN(2,MAXLIN)
      DOUBLE PRECISION SLIN(2,MAXLIN)
      REAL ALIN(2,MAXLIN),ELIN(2,MAXLIN),FH2
      REAL SBLIN(2,MAXLIN),TDW(2,MAXLIN)
      REAL TDWS(2,MAXLIN),PSHIFT(2,MAXLIN)
      REAL DOUBV(2,MAXLIN),VBOT,WING
      CHARACTER*15 LLQ(2,MAXLIN)
      CHARACTER*256 BUFFER
      LOGICAL WARNTD
      PARAMETER (FH2=0.865)
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      REAL TINW(MAXDGAS)
      DATA (TINW(J),J=1,42)/0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     $     0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     $     0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
     $     0.5,0.5/


C******************************** CODE *********************************

      WARNTD=.TRUE.


C     Assume line database file is already open.

      if(idiag.gt.0)print*,'loadbuffer: FSTREC=',FSTREC

C      print*,'VMIN,VMAX = ',VMIN,VMAX
C      print*,'MAXLIN,MAXBIN,IB = ',MAXLIN,MAXBIN,IB
C      print*,'NGAS',NGAS
C      do i=1,NGAS
C       print*,idgas(i),isogas(i)
C      enddo
C      print*,'VBOT,WING = ',VBOT,WING
   
      IREC=FSTREC-1
      LINE=0

      if(idiag.gt.0)print*,'IREC,LINE = ',IREC,LINE
    
      DO I=1,MAXBIN
       FSTLIN(IB,I)=-1
       LSTLIN(IB,I)=-1
      ENDDO
      FIRSTBIN(IB)=-1
      LASTBIN(IB)=-1
      DO I=1,MAXLIN
       VLIN(IB,I)=0.
       SLIN(IB,I)=0.
       ALIN(IB,I)=0.
       ELIN(IB,I)=0.
       IDLIN(IB,I)=0
       SBLIN(IB,I)=0.
       PSHIFT(IB,I)=0.
       DOUBV(IB,I)=0.
       TDW(IB,I)=0.
       TDWS(IB,I)=0.
       LLQ(IB,I)=' '
      ENDDO
111   IREC=IREC+1

C      print*,'IREC = ',irec
      IF(IREC.GT.DBSIZ)GOTO 102

      READ(DBLUN,110,REC=IREC)BUFFER(1:DBRECL)
110   FORMAT(A)

      CALL RDLINE(BUFFER)
   
C      WRITE(*,110)BUFFER(1:DBRECL)

      IF(LNWAVE.LT.VMAX)THEN
         DO 114 I=1,NGAS
C NOTE: that a line can be included more than once
          IF(IDGAS(I).EQ.LOCID(LNID))THEN
            
            LINE= LINE + 1


            IF(ISOGAS(I).EQ.0)THEN
              SLIN(IB,LINE)= LNSTR
            ELSE
              IF(LNISO.NE.DBISO(ISOGAS(I),LOCID(LNID)))THEN
C               Reset counter - this is a different isotope
                LINE=LINE-1
                GOTO 114
              ENDIF
              SLIN(IB,LINE)= LNSTR/RELABU(ISOGAS(I),IDGAS(I))
C              WRITE(*,*)' dbiso = ',DBISO(ISOGAS(I),LOCID(LNID))
C              WRITE(*,*)' relabu = ',RELABU(ISOGAS(I),IDGAS(I))
            ENDIF


            IDLIN(IB,LINE)= I

            VLIN(IB,LINE)= LNWAVE
            ELIN(IB,LINE)= LNLSE
            PSHIFT(IB,LINE)= LNPSH
            DOUBV(IB,LINE)= LDOUBV



            CURBIN = 1+INT((LNWAVE-VBOT)/WING)
            LSTLIN(IB,CURBIN)=LINE
            IF(FSTLIN(IB,CURBIN).EQ.-1)THEN
             FSTLIN(IB,CURBIN)=LINE
            ENDIF

            IF(LNWIDA.GT.0)THEN
C Set the air-broadened halfwidth equal to the Lorentz halfwidth
              ALIN(IB,LINE)= LNWIDA
            ELSE
              WRITE(*,*)' LOADBUFFER.f :: Warning: Lorentz halfwidth=0.'
              WRITE(*,*)' Setting to 0.075.'
              ALIN(IB,LINE)= 0.075
            ENDIF

            IF(LNWIDS.GT.1.E-20)THEN
C NOTE: SBLIN is the correction to air broadening so that zero is valid
              SBLIN(IB,LINE)= LNWIDA - LNWIDS
            ELSE
              SBLIN(IB,LINE)= 0.0
            ENDIF

            TDW(IB,LINE)= LNTDEP
            IF(LNTDEP.LT.1.E-30)THEN
             TDW(IB,LINE)= TINW(LOCID(LNID))
             IF(WARNTD)THEN
              WRITE(*,*)' LOADBUFFER :: Warning: setting default line'
              WRITE(*,*)' widths.'
              WARNTD= .FALSE.
             ENDIF
            ENDIF

            TDWS(IB,LINE)= LNTDEPS
            IF(LNTDEPS.LT.1.E-30)THEN
             TDWS(IB,LINE)= TINW(LOCID(LNID))
             IF(WARNTD)THEN
              WRITE(*,*)' LOADBUFFER :: Warning: setting default line'
              WRITE(*,*)' widths.'
              WARNTD= .FALSE.
             ENDIF
            ENDIF


            IF(LNWIDA1.GT.1.E-20)THEN
C            Composite 'Oxford'-ExoMOL format, which lists H2, He and Self-broadening.
C	     FH2 set as a parameter above. Edit as required. 
             ALIN(IB,LINE) = FH2*LNWIDA + (1.0-FH2)*LNWIDA1
             SBLIN(IB,LINE) = ALIN(IB,LINE) - LNWIDS            
             TDW(IB,LINE) = FH2*LNTDEP + (1-FH2)*LNTDEP1
            ENDIF

            IF(DBRECL.EQ.160)THEN
             LLQ(IB,LINE)=LNLLQ04
            ELSE
             LLQ(IB,LINE)(1:9)= LNLLQ
             LLQ(IB,LINE)(10:15)='      '
            ENDIF
            IF(LINE.EQ.MAXLIN)THEN
              if(idiag.gt.0)print*,'Reached end of buffer'
              NLINR=MAXLIN
              NXTREC=IREC+1
              GOTO 101
            ENDIF

          ENDIF
114      CONTINUE
         GOTO 111

      ENDIF

102   CONTINUE

      NLINR=LINE
      NXTREC=IREC

101   CONTINUE
      NBINX = 1+ INT((VMAX-VBOT)/WING) 
      DO I=1,NBINX
       IF(FSTLIN(IB,I).GT.0)THEN
        LASTBIN(IB)=I
        IF(FIRSTBIN(IB).LT.0)THEN
         FIRSTBIN(IB)=I
        ENDIF
       ENDIF
      ENDDO
      if(idiag.gt.0)then
       print*,'IB,NLINR,NXTREC,FIRSTBIN,LASTBIN',IB,NLINR,NXTREC,
     1 FIRSTBIN(IB),LASTBIN(IB)
      endif
      NBINY = 1+LASTBIN(IB)-FIRSTBIN(IB)
      IF(FIRSTBIN(IB).GT.0)THEN
       if(idiag.gt.0)print*,'Bins covered = ',NBINY
       IF(NBINY.LT.3)THEN
        print*,'Error in loadbuffer.f -'
        print*,'Line buffers should cover more than 2 wing bins.'
        print*,'You need to increase MAXLIN in lincomc.f and'
        print*,'recompile.'
        stop
       ENDIF
      ENDIF

      if(idiag.gt.0)print*,'loadbuffer end'

      RETURN

      END
