      PROGRAM SUMMARY
C     $Id: summary.f,v 1.2 2011-06-17 14:49:53 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  SUMMARY: outputs a summery of numbers of lines and strengths
C         for all species in given intervals
C
C_KEYS:   PROG,LINEDATA
C
C_DESCR:  outputs numbers of lines and maximum and minimum strengths for each
C         isotope of each species in each user defined bin over a given
C         wavenumber interval
C
C_FILES :
C         unit 1 - misc
C         unit 2 - data base
C
C_CALLS:  RDKEY     reads key file
C         FNDWAV    searches data base for wavenumber entry
C         RDLINE    uses internal read to translate line data record
C         PROMPT    prompts user for input
C         REMSP     removes leading spaces from text string
C         WTEXT     writes text to screen
C
C_BUGS:
C
C_HIST:
C         18jun92 SBC Original version
C
C_END:
C--------------------------------------------------------------
      INCLUDE '../includes/dbcom.f' 
C--------------------------------------------------------------
      CHARACTER*256 BUFFER
      REAL VMIN,VMAX,VBIN,BINSIZ
      DOUBLE PRECISION SMAX(MAXISO,MAXDGAS),SMIN(MAXISO,MAXDGAS)
      INTEGER NLINES(MAXISO,MAXDGAS),IBIN,NBIN,I,J,IMIN,IMAX,TOTLIN,OVTOT
      CALL PROMPT('data base key?')
      READ(*,102)KEYFIL
102   FORMAT(A)
      I=2
      CALL RDKEY(I)
      CALL RDGAS
      CALL RDISO

      OPEN(UNIT=DBLUN,FILE=DBFILE,STATUS='OLD',RECL=DBRECL,
     1ACCESS='DIRECT',FORM='FORMATTED')
      CALL PROMPT('wavenumber range?')
      READ(*,*)VMIN,VMAX
      CALL PROMPT('bin size?')
      READ(*,*)BINSIZ
      NBIN=1+INT((VMAX-VMIN)/BINSIZ)
      IF(NBIN.GT.1.AND.(VMIN+FLOAT(NBIN-1)*BINSIZ-VMAX).GT.-0.1)
     1NBIN=NBIN-1
      CALL FNDWAV(VMIN)
      OPEN(UNIT=3,FILE='summary.out',STATUS='UNKNOWN')
      WRITE(3,220)DBNAME
220   FORMAT(' data base: ',A)
      WRITE(3,221)DBFILE
221   FORMAT(' from DB file: ',A)
      WRITE(3,222)KEYFIL
222   FORMAT(' as defined in key: ',A)
      DO 10 IBIN=1,NBIN
      DO 30 I=1,MAXDGAS
      DO 30 J=1,MAXISO
      NLINES(J,I)=0
      SMAX(J,I)=0.
      SMIN(J,I)=1.E32
30    CONTINUE
      VBIN=VMIN+FLOAT(IBIN-1)*BINSIZ
      VMAX=VMIN+FLOAT(IBIN)*BINSIZ
      WRITE(*,202)IBIN,VBIN,VMAX
20    READ(DBLUN,12,REC=DBREC)BUFFER(1:DBRECL)
12    FORMAT(A)
      DBREC=DBREC+1
      CALL RDLINE(BUFFER)
      DO 213 I=1,DBNISO(LOCID(LNID))
      IF(LNISO.EQ.DBISO(I,LOCID(LNID)))GOTO 214
213   CONTINUE
      WRITE(*,215)LNID,LOCID(LNID),LNISO
215   FORMAT(' %unable to identify gas',I3,
     1' (Hitran ',I3,')  isotope',I4)
      GOTO 20
214   NLINES(I,LOCID(LNID))=NLINES(I,LOCID(LNID))+1
      SMAX(I,LOCID(LNID))=MAX(SMAX(I,LOCID(LNID)),LNSTR)
      SMIN(I,LOCID(LNID))=MIN(SMIN(I,LOCID(LNID)),LNSTR)
      IF(LNWAVE.LT.VMAX)GOTO 20
      WRITE(3,201)
201   FORMAT('--------------------------------------------------------')
      WRITE(3,202)IBIN,VBIN,VMAX
202   FORMAT(' BIN:',I3,3X,F10.3,' -',F10.3,' cm-1')
      OVTOT=0
      DO 204 I=1,MAXDGAS
      TOTLIN=0
      DO 209 J=1,MAXISO
      IF(NLINES(J,I).LT.1)GOTO 209
      IF(TOTLIN.EQ.0)THEN
        WRITE(3,203)
203     FORMAT('         ',' id',' iso','    smin     ',
     1  '    smax     ',' nlines')
        END IF
      TOTLIN=TOTLIN+NLINES(J,I)
      WRITE(BUFFER,207)SMIN(J,I),SMAX(J,I)
207   FORMAT(2E13.5)
      IF(SMIN(J,I).GT.0.)THEN
        IMIN=INT(1000+DLOG10(SMIN(J,I)))-1046
        WRITE(BUFFER(11:13),208)IMIN
       ELSE
        WRITE(BUFFER(11:13),212)
212     FORMAT('***')
        END IF
      IF(SMAX(J,I).GT.0.)THEN
        IMAX=INT(1000+DLOG10(SMAX(J,I)))-1046
        WRITE(BUFFER(24:26),208)IMAX
       ELSE
        WRITE(BUFFER(11:13),212)
        END IF
208   FORMAT(I3.2)
      WRITE(3,205)GASNAM(I),I,J,BUFFER(1:26),NLINES(J,I)
205   FORMAT(1X,1A8,I3,I4,1A26,I7)
209   CONTINUE
      IF(TOTLIN.GT.0)WRITE(3,210)TOTLIN,GASNAM(I)
210   FORMAT(' total of',I8,1X,A,' lines')
      OVTOT=OVTOT+TOTLIN
204   CONTINUE
      WRITE(3,211)OVTOT
211   FORMAT(' overall total of',I9,' lines in bin')
10    CONTINUE
      WRITE(3,201)
      WRITE(*,206)
206   FORMAT(' output written to summary.out')
      STOP
      END
