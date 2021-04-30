      PROGRAM CP_LINES_SEQ_MODWIDJ
C     $Id:
C--------------------------------------------------------------
C_TITLE:  COPY: copies a subset of sequential-access line data base to a new data base
C
C_KEYS:   PROG,LINEDATA
C
C_DESCR:  copies line data bases (eg HITRAN or GEISA) from an existing data
C         base into a sequential ascii file which can be turned into a new data
C         base using makedb. Unlike SELECT.FOR it copies ALL lines
C         between given wavenumber limits rather than selecting them by
C         gas, isotope and strength.
C
C_FILES :
C         unit 2 - data base files
C         unit 3 - output file
C
C_BUGS:   
C
C_HIST:
C         26apr92 SBC Original version
C	  10apr18 PGJI Updated for sequential access
C_END:
C--------------------------------------------------------------
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
C
      DOUBLE PRECISION VMIN,VMAX,STRMIN
      REAL FH2,JWID,WCOEFF(5),JWIDS
      INTEGER JL,I
      CHARACTER*256 BUFFER,TBUFF
      CHARACTER*100 OPNAME,AFILE
C
 
      PRINT*,'Enter name of file containing self-broadened widths'
      CALL PROMPT('as function of J : ')
      READ(5,1)AFILE
      OPEN(12,FILE=AFILE,STATUS='OLD')
       READ(12,*)(WCOEFF(I),I=1,5)
       READ(12,*)JMAX,KMAX
      CLOSE(12)

C     open existing sequential access database
      CALL PROMPT('Enter name of sequential-access data base:')
      READ(*,1)KEYFIL
1     FORMAT(A)
      OPEN(DBLUN,FILE=KEYFIL,STATUS='OLD',ACTION='READ')

C     read in parameters
2     CALL PROMPT('wavenumber limits to extract?')
      READ(*,*)VMIN,VMAX
      IF(VMAX.LE.VMIN)GOTO 2

      CALL PROMPT('Enter name of output file:')
      READ(*,1)OPNAME
      OPEN(UNIT=3,FILE=OPNAME,STATUS='UNKNOWN')

      CALL PROMPT('Enter minimum line strength : ')
      READ*,STRMIN

      STRMIN=STRMIN*1e20
      STRMIN=STRMIN*1e27

      CALL PROMPT('Enter record length : ')
      READ*,DBRECL

C     For some weird reason the first release of NH3 ExoMOL had the
C     H2-broadening coefficients in the foreign-broadened field and
C     He-broadening in the self-broadened field
C      CALL PROMPT('Enter assumed H2 vmr : ')
C      READ*,FH2


      PRINT*,'Finding first line in database with v >= ',VMIN
100   READ(DBLUN,1,END=999)BUFFER
      CALL RDLINE(BUFFER)
      IF(LNWAVE.LT.VMIN)GOTO 100
      
      PRINT*,'Starting transfer'
      IF(LNSTR.GE.STRMIN)THEN
C        WRITE(6,1)BUFFER(1:DBRECL)
C        READ(LNLLQ04(1:3),*)JL
C        PRINT*,JL
C        SWID=JWID(WCOEFF,JMAX,JL)
C        PRINT*,SWID
C        WRITE(TBUFF,302)SWID
C        BUFFER(41:45)=TBUFF(1:5)
C        WRITE(6,1)BUFFER(1:DBRECL)
C        STOP
C        WRITE(3,1)BUFFER(1:DBRECL)

        WRITE(6,1)BUFFER(1:DBRECL)
        print*,LNUGQI04
        print*,LNLGQI04
        print*,LNULQ04
        print*,LNLLQ04
        READ(LNLLQ04(1:3),*)JL
        READ(LNLLQ04(4:6),*)KL
        PRINT*,JL,KL
        SWID=JWIDS(WCOEFF,JMAX,KMAX,JL,KL)
        PRINT*,SWID
        WRITE(TBUFF,302)SWID
        BUFFER(41:45)=TBUFF(1:5)
        WRITE(6,1)BUFFER(1:DBRECL)
        STOP
        WRITE(3,1)BUFFER(1:DBRECL)
      ENDIF
110   READ(DBLUN,1,END=999)BUFFER
      CALL RDLINE(BUFFER)
      IF(LNWAVE.GT.VMAX)GOTO 998
      
      IF(LNSTR.GE.STRMIN)THEN
C        READ(LNLLQ(1:3),*)JL
C        SWID=JWID(WCOEFF,JMAX,JL)
C        WRITE(TBUFF,302)SWID
C        BUFFER(41:45)=TBUFF(1:5)
C        WRITE(3,1)BUFFER(1:DBRECL)

        READ(LNLLQ(1:3),*)JL
        READ(LNLLQ(4:6),*)KL
        SWID=JWIDS(WCOEFF,JMAX,KMAX,JL,KL)
        WRITE(TBUFF,302)SWID
        BUFFER(41:45)=TBUFF(1:5)
        WRITE(3,1)BUFFER(1:DBRECL)
      ENDIF
      GOTO 110

301   FORMAT(F5.4)
302   FORMAT(F5.3)

998   PRINT*,'Transfer complete'
      CLOSE(3)
      CLOSE(DBLUN)
      STOP  

999   PRINT*,'Reached end of database'
      CLOSE(3)
      CLOSE(DBLUN)

      END

