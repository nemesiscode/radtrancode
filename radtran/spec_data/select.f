      PROGRAM SELECT
C     $Id: select.f,v 1.4 2011-06-17 14:52:18 irwin Exp $
************************************************************************
************************************************************************
C_TITL:	SELECT.f
C
C_DESC:	Copies a subset of a line data base to a new data base. Copies
C	linedata bases (eg HITRAN or GEISA) from an existing data base
C	into a sequential ascii file which can be turned into a new data
C	base using makedb. All lines or only selected lines between given
C	wavenumber limits are copies. Selection can be performed by gas 
C	id, gas and isotope or by strength limit. The strength limit for a
C	particular gas or isotope is set so that the sum of strengths of
C	ommitted lines is less than n% of the total sum of strengths. This
C	criteria is calculated by finding the number of lines in each 
C	decade of strength.
C	NOTE: the strength limit can be set independently for each isotope
C	or can be set for all of the isotopes to be included. In the
C	latter case terrestrial isotopic abundances are assumed, i.e the
C	strengths are used uncorrected. The logic behind this is that if
C	the isotopic abundance may vary you want to treat them
C	independently anyway.
C
C_ARGS:	See the definitions below.
C
C_FILE: unit=2	Input database file.
C	unit=3	Output file.
C
C_CALL:	RDGAS		Reads in the gas file information.
C	RDISO		Reads in the isotope file information.
C	RDLINE		Reads line data from buffer.
C	ASKYN		Prompts user for yes/no answer.
C	WTEXT		Writes test to screen.
C	PROMPT		Prompts user for input.
C	RDKEY		Reads in details from the key file.
C	FNDWAV		Searches data base for wavenumber entry.
C	READ_DEL	Reads in Delaye water broadening data.
C	READ_YAM	Reads in Yamamoto CO2 broadening data.
C
C_HIST:	12sep91	SBC	Original version
C	26apr92	SBC	Removed list directed internal read - not strict
C			F77.
C	10feb92	PGJI	Added routine to read in delaye data and change
C			call to EDLINE.
C	26feb92	PGJI	Added routine to output Schofield PMR2 compatable
C			data.
C       29aug07 CAN     Found a bug: the final write statement was outside
C                       the 'IF' statement, meaning that the line strength
C                       cut-off stuff was not working - just writing out
C                       all lines regardless.
C
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INCLUDE '../includes/dbcom.f' 
C ../includes/dbcom.f stores the linedata base variables.

      INTEGER MINSTR,MAXSTR
      REAL LIMSTR
      PARAMETER (MINSTR=-37,MAXSTR=38,LIMSTR=10.**MINSTR)

      INTEGER NLINES(MINSTR:MAXSTR,MAXISO,MAXDGAS)
      INTEGER TOTLIN(MINSTR:MAXSTR,MAXDGAS)
      INTEGER I,J,K,ID,ISO,IBIN,LINE,FSTLIN,LSTLIN,NBIN
      INTEGER FIRST(2),LAST(2),NPAR,IERROR,READI

      REAL VMIN,VMAX,BINSIZ,LIMIT(MAXISO,MAXDGAS),PERCEN,VLOW,VHIGH
C VMIN: Wavenumber [cm-1] minimum.
C VMAX: Wavenumber [cm-1] maximum.
C BINSIZ: Size [cm-1] of bins for limit selection.
      REAL TOTSTR,SUMSTR
      REAL ENERGY(190),ACO2(190),NCO2(190),AH2O(190),NH2O(190)
      REAL AN2(190),NN2(190),AO2(190),NO2(190)
      REAL YACO2(200),YNCO2(200),YAN2(200),YNN2(200),ECO2(200),EN2(200)

      CHARACTER*6 QIDENT(190)
      CHARACTER*100 OPNAME
      CHARACTER*100 TEXT
      CHARACTER*256 BUFFER

      LOGICAL ASKYN,INCGAS(MAXISO,MAXDGAS),ALLISO(MAXISO,MAXDGAS),INC
C INCGAS: Flag to show if a gas is included.
C ALLISO: Flag to show if using all isotopes for limit calculation.

C******************************** CODE *********************************

C Open database ...
      CALL PROMPT('name of data base key')
      READ(*,1)KEYFIL
      CALL REMSP(KEYFIL)
1     FORMAT(A)
      CALL RDKEY(2)
      CALL RDGAS
      CALL RDISO
      CALL READ_DEL(QIDENT,ENERGY,ACO2,NCO2,AH2O,NH2O,AN2,NN2,AO2,NO2)
      CALL READ_YAM(YACO2,YNCO2,YAN2,YNN2,ECO2,EN2)

      DO I=1,MAXDGAS
        DO J=1,MAXISO
          INCGAS(J,I) = .FALSE.
          ALLISO(J,I) = .FALSE.
        ENDDO
      ENDDO

C Read in spectral parameters ...
2     CALL PROMPT('Enter wavenumber limits')
      READ(*,*)VMIN,VMAX
      IF(VMAX.LE.VMIN)GOTO 2

      WRITE(*,*)'Lines whose summed strength is less than n% of'
      WRITE(*,*)'the total sum strength in a bin may be omitted'
22    CALL PROMPT('Enter bin width and cut-off limit [%]')
      READ(*,*)BINSIZ,PERCEN
      PERCEN = PERCEN*0.01
      IF(BINSIZ.LT.0.01)THEN
        WRITE(*,*)'Bin size is too small'
        GOTO 22
      ENDIF

      CALL PROMPT('Enter new file name')
      READ(*,21)OPNAME
21    FORMAT(A)
      OPEN(UNIT=3,FILE=OPNAME,STATUS='UNKNOWN')
      WRITE(BUFFER,114)
114   FORMAT(' # data records written by routine SELECT')
      WRITE(3,111)BUFFER(1:DBRECL)
      WRITE(BUFFER,201)
201   FORMAT(' # original data base files:')
      WRITE(3,111)BUFFER(1:DBRECL)
      WRITE(BUFFER,202)DBFILE
202   FORMAT(' #  ',A)
      WRITE(3,111)BUFFER(1:DBRECL)
      WRITE(BUFFER,202)KEYFIL
      WRITE(3,111)BUFFER(1:DBRECL)
      CALL EDSET

C Select the gases ...
10    CONTINUE

      CALL WTEXT('--------------------------------------------------')
      CALL WTEXT('enter gas details')
      CALL WTEXT('enter gas id and isotope [0=all isotopes]')
      CALL WTEXT('ALL includes all gases, GO starts execution,')
      CALL WTEXT(' ? gives summary')
      CALL PROMPT('ID ISO ?')
      READ(*,11)TEXT
11    FORMAT(A)
      CALL UPCASE(TEXT)
      CALL REMSP(TEXT)
      IF(TEXT(1:2).EQ.'GO')GOTO 20
      IF(TEXT(1:3).EQ.'ALL')THEN
        INC = ASKYN('all isotopes for strength cut-off limit?')
        DO 13 I=1,MAXDGAS
          IF(LOCID(I).NE.0)THEN
            DO 12 J=1,DBNISO(LOCID(I))
              INCGAS(J,LOCID(I)) = .TRUE.
              ALLISO(J,LOCID(I)) = INC
12          CONTINUE
          ENDIF
13      CONTINUE
        GOTO 10
      ENDIF

      IF(TEXT(1:1).EQ.'?')GOTO 15
C NOTE: separate command line processor used below rather than internal
C read because list directed internal read is considered not F77 by some
C compilers (notably Prospero Fortran for MS-DOS)
      NPAR = 2
      CALL CLP(TEXT,NPAR,FIRST,LAST)
      ID = READI(TEXT(FIRST(1):LAST(1)),IERROR)
      IF(IERROR.NE.0)GOTO 15
      ISO = READI(TEXT(FIRST(2):LAST(2)),IERROR)
      IF(IERROR.NE.0)GOTO 15
      WRITE(*,*)'ID ISO = ',ID,ISO
      IF(ISO.EQ.0)THEN
        DO 14 J=1,DBNISO(ID)
          INCGAS(J,ID) = .TRUE.
          ALLISO(J,ID) = .TRUE.
14      CONTINUE
      ELSE
        IF(INCGAS(ISO,ID))THEN
          INCGAS(ISO,ID) = .FALSE.
        ELSE
          INCGAS(ISO,ID) = .TRUE.
          ALLISO(ISO,ID) = ASKYN('all isotopes for limit?')
        ENDIF
      ENDIF
      GOTO 10
15    CONTINUE
      CALL WTEXT(' ')
      CALL WTEXT('  ID name     isotope, included?, use all for limit?')
      K = 0
      DO 16 I=1,MAXDGAS
        DO 18 J=1,DBNISO(I)
          IF(INCGAS(J,I))GOTO 19
18      CONTINUE
        GOTO 16

19      CONTINUE

        WRITE(TEXT,17)I,GASNAM(I),(J,INCGAS(J,I),ALLISO(J,I),
     1  J=1,DBNISO(I))
        K = K + 1
17      FORMAT(I3,1X,1A8,13(I3,L1,L1))
16    CONTINUE
      IF(K.EQ.0)CALL WTEXT('no gases included')
      CALL WTEXT(TEXT)
      GOTO 10
20    CONTINUE

      OPEN(DBLUN,FILE=DBFILE,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='DIRECT',RECL=DBRECL)

      WRITE(BUFFER,203)
203   FORMAT(' # original data base header if any:')
      WRITE(3,111)BUFFER(1:DBRECL)
      LINE = 1
204   READ(DBLUN,111,REC=LINE)BUFFER(1:DBRECL)
      IF(BUFFER(2:2).EQ.'#'.OR.BUFFER(1:1).EQ.'#')THEN
        WRITE(3,111)BUFFER(1:DBRECL)
        LINE = LINE + 1
        GOTO 204
      ENDIF
      WRITE(BUFFER,205)
205   FORMAT(' # selection criteria:')
      WRITE(3,111)BUFFER(1:DBRECL)
      WRITE(BUFFER,206)VMIN,VMAX,BINSIZ
206   FORMAT(' # wavenumber range =',F10.3,' -',F10.3,' with',F8.3,
     1 ' cm-1 bins')
      WRITE(3,111)BUFFER(1:DBRECL)

      DO 207 I=1,MAXDGAS
         IF(DBNISO(I).LT.1)GOTO 207
         WRITE(BUFFER,208)GASNAM(I),(INCGAS(J,I),J=1,DBNISO(I))
 208     FORMAT(' #',1A8,' INCGAS:',20(1X,L1))
         WRITE(3,111)BUFFER(1:DBRECL)
         WRITE(BUFFER,209)GASNAM(I),(ALLISO(J,I),J=1,DBNISO(I))
 209     FORMAT(' #',1A8,' ALLISO:',20(1X,L1))
         WRITE(3,111)BUFFER(1:DBRECL)
 207  CONTINUE
      
      
C     For each bin ...
      NBIN = INT((VMAX-VMIN)/BINSIZ) + 1
      
C     First setting LSTLIN to one record before the first in the first bin
C     (FNDWAV returns first record in DBREC) each bin then starts off at one
C     record after LSTLIN and then sets it to the last record in the bin.
      CALL FNDWAV(VMIN)
      LSTLIN = DBREC - 1
      DO 100 IBIN=1,NBIN
         VLOW = VMIN + FLOAT(IBIN-1)*BINSIZ
         VHIGH = VMIN + FLOAT(IBIN)*BINSIZ
         WRITE(*,121)IBIN,VLOW,VHIGH
 121     FORMAT(I5,F12.3,'-',F12.3)
C     Find wavenumber region in data base
         FSTLIN = LSTLIN + 1
         CALL FNDWAV(VHIGH)
         LSTLIN = DBREC - 1

         IF(PERCEN.GT.0.0)THEN
C     Load NLINES array
C     NLINES (strength decade, iso, ngas) is the number of lines in 
C     each strength decade
            DO 120 I=1,MAXDGAS
               DO 122 J=1,DBNISO(I)
                  DO 123 K=MINSTR,MAXSTR
                     NLINES(K,J,I) = 0
 123              CONTINUE
 122           CONTINUE
 120        CONTINUE
            
            DO 110 LINE=FSTLIN,LSTLIN
               READ(DBLUN,111,REC=LINE)BUFFER(1:DBRECL)
 111           FORMAT(A)
               CALL RDLINE(BUFFER)
               IF(LNSTR.LT.LIMSTR)THEN
                  WRITE(TEXT,115)LNWAVE,LNSTR,LNID,LNISO
 115              FORMAT('WARNING - strength too low for storage',
     1                 F12.6,E12.5,2I3)
                  CALL WTEXT(TEXT)
                  GOTO 110
               ENDIF

               I = INT(ALOG10(LNSTR))
               IF(I.GT.MAXSTR)THEN
                  WRITE(*,113)LNWAVE,LNSTR,LNID,LNISO
 113              FORMAT('WARNING - strength too high for storage',
     1                 F12.6,E12.5,2I3)
                  STOP
               ENDIF
               
C     see if this line is the right isotope
               K = -1
               DO 112 J=1,DBNISO(LOCID(LNID))
                  IF(LNISO.EQ.DBISO(J,LOCID(LNID)))K=J
 112           CONTINUE
               
C     if yes, then add 1 to count
               IF(K.NE.-1)THEN
c                  if (LNID.EQ.6) print*, LNWAVE,LNSTR,LNID,LNISO        
                  NLINES(I,K,LOCID(LNID))=NLINES(I,K,LOCID(LNID))+1
               ENDIF
               
 110        CONTINUE
            
C count the total lines for all isotopes for this strength decade, gas id
            DO 130 I=1,MAXDGAS
               DO 132 J=MINSTR,MAXSTR
                  TOTLIN(J,I) = 0
                  DO 131 K=1,DBNISO(I)
                     IF(INCGAS(K,I))TOTLIN(J,I) = 
     C                    TOTLIN(J,I) + NLINES(J,K,I)
 131              CONTINUE
 132           CONTINUE
 130        CONTINUE
            
C Compute limit for each isotope ...
        DO 140 I=1,MAXDGAS
          DO 146 J=1,DBNISO(I)
            IF(INCGAS(J,I))THEN
              IF(ALLISO(J,I))THEN
                TOTSTR = 0.
                DO 144 K=MINSTR,MAXSTR
                  TOTSTR = TOTSTR + FLOAT(TOTLIN(K,I))*10.**K
144             CONTINUE
                LIMIT(J,I) = 10.**MINSTR
                SUMSTR = 0.
                DO 145 K=MINSTR,MAXSTR
                  SUMSTR = SUMSTR + FLOAT(TOTLIN(K,I))*10.**K
                  IF(SUMSTR.GT.PERCEN*TOTSTR)GOTO 143
                  LIMIT(J,I) = 10.**K
145             CONTINUE
              ELSE
                TOTSTR = 0.
                DO 141 K=MINSTR,MAXSTR
                  TOTSTR = TOTSTR + FLOAT(NLINES(K,J,I))*10.**K
c                  print*, FLOAT(NLINES(K,J,I))*10.**K, totstr
141             CONTINUE
                LIMIT(J,I) = 10.**MINSTR
                SUMSTR = 0.
                DO 142 K=MINSTR,MAXSTR
                  SUMSTR = SUMSTR + FLOAT(NLINES(K,J,I))*10.**K

c                  print*, FLOAT(NLINES(K,J,I))*10.**K, sumstr,
c     C                 percen*totstr, totstr, limit(j,i)
                  IF(SUMSTR.GT.PERCEN*TOTSTR)GOTO 143
                  LIMIT(J,I) = 10.**K
142             CONTINUE
              ENDIF

143           CONTINUE

            ENDIF
146       CONTINUE
140     CONTINUE

C Print out computed limits for each gas ...
      WRITE(*,*)'Limiting strengths for each gas:'
      DO 134 I=1,MAXDGAS
        WRITE(*,*)I,DBNISO(I),(LIMIT(K,I),K=1,DBNISO(I))
134   CONTINUE

      ENDIF

C Copy required lines ...
      DO 150 LINE=FSTLIN,LSTLIN
        READ(DBLUN,111,REC=LINE)BUFFER(1:DBRECL)
        CALL RDLINE(BUFFER)
        K = -1
        DO 152 J=1,DBNISO(LOCID(LNID))
          IF(LNISO.EQ.DBISO(J,LOCID(LNID)))K = J
152     CONTINUE
        IF(K.NE.-1.) THEN
           IF(INCGAS(K,LOCID(LNID)))THEN
              IF(PERCEN.GT.0.0)THEN
                 IF (LNSTR.GE.LIMIT(K,LOCID(LNID)))THEN
c                    print*, 'comp here ', lnstr, limit(k,locid(lnid))
                    CALL EDLINE(BUFFER,QIDENT,ACO2,NCO2,AH2O,NH2O,YACO2,
     1                   YNCO2,YAN2,YNN2)
                    WRITE(3,111)BUFFER(1:DBRECL)
                 ENDIF
              ELSE
                 CALL EDLINE(BUFFER,QIDENT,ACO2,NCO2,AH2O,NH2O,YACO2,
     1                YNCO2,YAN2,YNN2)
                 WRITE(3,111)BUFFER(1:DBRECL)
              ENDIF
           ENDIF
        ENDIF
150   CONTINUE

      IF(VHIGH.GE.VMAX)GOTO 101

100   CONTINUE
101   CONTINUE

      STOP

      END
************************************************************************
************************************************************************
