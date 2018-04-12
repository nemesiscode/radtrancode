      PROGRAM LISTLI
C     $Id: li_lines.f,v 1.5 2011-06-17 14:57:59 irwin Exp $
C---------------------------------------------------------------------------
C
C_TITLE: LIST_LINES: list spectral lines from binary data base
C
C_ARGS:  None.
C
C_KEYS:  PROG, VMS, CAL .
C
C_DESCR: Prompts for line data key, gases and wavenumber range. Lists main
C        line data parameters for these gases.
C
C_FILES: unit DBLUN  line data and gas files
C
C_CALLS:  LINES         reads line data
C         RDKEY         reads in details from the key file
C         PROMPT        user prompt utility
C         REMSP         removes leading spaces from text string
C         RDGAS         reads in gas information
C         RDISO         reads in isotope information
C         FNDWAV        searches database for wavenumber entry
C         RDLINE        performs internal read on line data record
C        			
C_BUGS:  
C
C_HIST:  18sep89 SBC ORIGINAL VERSION copied from LIST.FOR
C        22feb91 SBC modified to use new data base format
C        17jun92 SBC changed strength scaling to e-47 to allow more range
C         8dec92 SBC switched to using GETLIN instead of LINES
C	  5mar97 PGJI switched back to LINES
C
C_END:
C
C--------------------------------------------------------------------------
      INCLUDE '../includes/dbcom.f' 
C--------------------------------------------------------------
      INTEGER LINLIM
      PARAMETER (LINLIM=200000)
      REAL VLIN(LINLIM),ALIN(LINLIM),ELIN(LINLIM),
     1 SBLIN(LINLIM),TDW(LINLIM),TDWS(LINLIM),PSHIFT(LINLIM),
     2 DOUBV(LINLIM)
      DOUBLE PRECISION SLIN(LINLIM)
      CHARACTER*15 LLQ(LINLIM)
      INTEGER IDLIN(LINLIM),IDGAS(MAXDGAS),ISOGAS(MAXDGAS)
      INTEGER NGAS,IGAS,MAXLIN,IFIRST,NLIN,ILAST,I,DISEXP
      INTEGER NEXTRA
      CHARACTER*1 ANS

      DOUBLE PRECISION VMIN,VMAX,DELV
      CHARACTER*12 MANTIS
C
      CALL PROMPT('data base key?')
      READ(*,102)KEYFIL
      CALL REMSP(KEYFIL)
102   FORMAT(A)
      CALL RDKEY(2)
      CALL RDGAS
      CALL RDISO
      CALL PROMPT('how many gases?')
      READ(*,*)NGAS
      DO 108 IGAS=1,NGAS
      CALL PROMPT('gas id, isotope code?')
      READ(*,*)IDGAS(IGAS),ISOGAS(IGAS)
      IF(ISOGAS(IGAS).EQ.0)THEN
       PRINT*,IDGAS(IGAS),ISOGAS(IGAS),
     &		'Output uncorrected line-strengths'
      ELSE
       PRINT*,IDGAS(IGAS),ISOGAS(IGAS),
     &		'Output line-strengths/terrestrial relative abund.'
      ENDIF
108   CONTINUE

    
      CALL PROMPT('vmin,vmax?')
      READ(*,*)VMIN,VMAX
      MAXLIN=LINLIM
      NEXTRA=LINLIM
      IFIRST=1
      DELV=VMAX-VMIN

      CALL LINESS(VMIN,DELV,MAXLIN,VLIN,SLIN,ALIN,ELIN,IDLIN,SBLIN,
     1PSHIFT,DOUBV,TDW,TDWS,LLQ,NLIN,IFIRST,ILAST,NGAS,IDGAS,ISOGAS)

      WRITE(*,110)NLIN
110   FORMAT(' found ',I8,' lines')

      CALL PROMPT('Print to screen(1) or to file (2) : ')
      READ*,IPRINT

      IF(IPRINT.EQ.1)THEN
       IW=6
      ELSE
       OPEN(12,file='li_lines.out',status='unknown')
       IW=12
       WRITE(IW,110)NLIN
       WRITE(IW,*)'NGAS = ',NGAS
       DO I=1,NGAS
        WRITE(IW,*)IDGAS(I),ISOGAS(I)
       ENDDO
      ENDIF
      WRITE(IW,202)
      DO 200 I=1,NLIN

        IF(MOD(I,20).EQ.1)THEN
	 IF(IW.EQ.6.AND.I.NE.1)THEN
          READ(5,101)ANS
          WRITE(IW,202)
         ENDIF
	ENDIF
C       now need to rescale strengths by 1.E-47.
C       can't do this trivially since 1.e-47 is not a valid number and
C       the result would be useless for weak lines anyway
        WRITE(MANTIS,203)SLIN(I)
203     FORMAT(E12.5)
        DISEXP=INT((1000+DLOG10(SLIN(I))))-1046
      WRITE(IW,201)VLIN(I),MANTIS(1:8),DISEXP,ALIN(I),ALIN(I)-SBLIN(I),
     1 ELIN(I),TDW(I),TDWS(I),IDGAS(IDLIN(I)),ISOGAS(IDLIN(I))
201     FORMAT(1X,F10.3,1X,1A8,'E',I4,F7.4,F7.4,F8.1,F5.2,F5.2,I4,I4)

200   CONTINUE      

      IF(IPRINT.EQ.2)THEN
       CLOSE(12)
      ENDIF

101   FORMAT(A)

202   FORMAT(1X,5X,'V',4X,6X,'S',6X,3X,'A',3X,3X,'SB',3X,3X,'E',3X,
     1' TDW ',1X,'TDWS',1X,'ID',2X,'ISO')

      STOP
      END
