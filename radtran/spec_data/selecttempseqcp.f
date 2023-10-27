      PROGRAM SELECTTEMPSEQCP
C     $Id: selecttemp.f,v 1.2 2011-06-17 14:51:54 irwin Exp $
************************************************************************
************************************************************************
C_TITL:	SELECTTEMP.f
C
C_DESC:	Copies a subset of a sequential line data base to a new data base. 
C       Copies linedata bases (eg HITRAN or GEISA) from an existing data base
C	into a sequential ascii file which can be turned into a new data
C	base using makedb. All lines or only selected lines between given
C	wavenumber limits are copies. Selection can be performed by gas 
C	id, gas and isotope or by strength limit at a particular
C       temperature. The strength limit for a particular gas or isotope 
C       at that temperature is set so that the sum of strengths of
C	ommitted lines is less than n% of the total sum of strengths. This
C	criteria is calculated by finding the number of lines in each 
C	decade of strength.
C
C       A correction is then made to the remaining lines.
C
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
C	09jun10 PGJI	Revised version to select lines relevant for
C                       a given temperature.
C
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INCLUDE '../includes/dbcom.f' 
C ../includes/dbcom.f stores the linedata base variables.

      INTEGER MINSTR,MAXSTR,NSAV,IDF,ISOF
      DOUBLE PRECISION LIMSTR,CC,LNSTR1,XS,STRMIN,STRMINX
      PARAMETER (MINSTR=-323,MAXSTR=308,CC=10.)

      INTEGER NLINES(MINSTR:MAXSTR,MAXISO,MAXDGAS)
      INTEGER TOTLIN(MINSTR:MAXSTR,MAXDGAS),MAXBIN
      PARAMETER(MAXBIN=10000)
      DOUBLE PRECISION STRLOSE,STRKEEP,STROUT,WEIGHT
      DOUBLE PRECISION TSTROUT(MAXBIN),TXLSE(MAXBIN)
      DOUBLE PRECISION TXWIDA(MAXBIN),TXWIDS(MAXBIN)
      DOUBLE PRECISION TXTDEP(MAXBIN)



      REAL XWIDA,XWIDS,XLSE,XTDEP
      INTEGER JBIN
      INTEGER I,J,K,ID,ISO,IBIN,LINE,FSTLIN,LSTLIN,NBIN,IPTF
      INTEGER FIRST(2),LAST(2),NPAR,IERROR,READI,NLIN,NKEEP,NLOSE

      DOUBLE PRECISION VMIN,VMAX,BINSIZ,VLOW,VHIGH,VMEAN
      DOUBLE PRECISION LIMIT(MAXISO,MAXDGAS),TOTSTR,SUMSTR
C VMIN: Wavenumber [cm-1] minimum.
C VMAX: Wavenumber [cm-1] maximum.
C BINSIZ: Size [cm-1] of bins for limit selection.
      REAL TCALC,TCORS1,TCORS2,TSTIM
      DOUBLE PRECISION LNABSCO
      REAL ENERGY(190),ACO2(190),NCO2(190),AH2O(190),NH2O(190)
      REAL AN2(190),NN2(190),AO2(190),NO2(190)
      REAL YACO2(200),YNCO2(200),YAN2(200),YNN2(200),ECO2(200),EN2(200)
      REAL TEMP,TS1,TS2,PARTF
      CHARACTER*6 QIDENT(190)
      CHARACTER*100 OPNAME,OPSTR
      CHARACTER*100 TEXT
      CHARACTER*256 BUFFER,BUFFERSAV

      LOGICAL ASKYN,INCGAS(MAXISO,MAXDGAS),ALLISO(MAXISO,MAXDGAS),INC
C INCGAS: Flag to show if a gas is included.
C ALLISO: Flag to show if using all isotopes for limit calculation.

C******************************** CODE *********************************
      LIMSTR=CC**MINSTR

C Open database ...
      CALL PROMPT('Enter name of sequential access line data: ')
      READ(*,1)DBFILE
1     FORMAT(A)
      CALL PROMPT('Enter name of gas file : ')
      READ(*,1)GASFIL
      CALL PROMPT('Enter name of isotopes file : ')
      READ(*,1)ISOFIL
      CALL PROMPT('Enter DBRECL : ')
      READ*,DBRECL

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

      WRITE(*,*)'Lines whose strength is less than cut-off will'
      WRITE(*,*)'be assigned to continuum'
22    CALL PROMPT('Enter bin width and line strength cut-off')
      READ(*,*)BINSIZ,STRMINX
      IF(BINSIZ.LT.0.01)THEN
        WRITE(*,*)'Bin size is too small'
        GOTO 22
      ENDIF
      STRMIN=STRMINX*1E20
      STRMIN=STRMIN*1E27

      CALL PROMPT('Enter calculation temperature : ')
      READ(*,*)TEMP

      CALL PROMPT('Enter new file  name')
      READ(*,21)OPNAME
21    FORMAT(A)

      OPEN(UNIT=3,FILE=OPNAME,STATUS='UNKNOWN')
      CALL FILE(OPNAME,OPSTR,'lco')
      OPEN(UNIT=4,FILE=OPSTR,STATUS='UNKNOWN')


      WRITE(BUFFER,114)
114   FORMAT(' # data records written by routine SELECTTEMPSEQCP')
      WRITE(3,1)BUFFER(1:DBRECL)
      WRITE(BUFFER,201)
201   FORMAT(' # original data base file:')
      WRITE(3,1)BUFFER(1:DBRECL)
      WRITE(BUFFER,202)DBFILE
202   FORMAT(' #  ',A)
      WRITE(3,1)BUFFER(1:DBRECL)

      WRITE(4,*)'VMIN, VMAX = ',VMIN,VMAX
      WRITE(4,*)'BINSIZ,STRMIN = ',BINSIZ,STRMINX

      CALL PROMPT('Enter IPTF (Partition function flag) : ')
      READ*,IPTF

      WRITE(4,*)'Calc. Temp, IPTF = ',TEMP,IPTF

      CALL PROMPT('Enter ID,ISO to put in header : ')
      READ(5,*)IDF,ISOF
      WRITE(4,*)'ID, ISO = ',IDF,ISOF

      OPEN(DBLUN,FILE=DBFILE,STATUS='OLD')

      WRITE(BUFFER,203)
203   FORMAT(' # original data base header if any:')
      WRITE(3,1)BUFFER(1:DBRECL)
204   READ(DBLUN,1,END=999)BUFFER(1:DBRECL)
      IF(BUFFER(2:2).EQ.'#'.OR.BUFFER(1:1).EQ.'#')THEN
        WRITE(3,1)BUFFER(1:DBRECL)
        GOTO 204
      ENDIF
      WRITE(BUFFER,205)
205   FORMAT(' # selection criteria:')
      WRITE(3,1)BUFFER(1:DBRECL)
      WRITE(BUFFER,206)VMIN,VMAX,BINSIZ
206   FORMAT(' # wavenumber range =',F10.3,' -',F10.3,' with',F8.3,
     1 ' cm-1 bins')
      WRITE(3,1)BUFFER(1:DBRECL)

      DO 207 I=1,MAXDGAS
         IF(DBNISO(I).LT.1)GOTO 207
         WRITE(BUFFER,208)GASNAM(I),(INCGAS(J,I),J=1,DBNISO(I))
208      FORMAT(' #',1A8,' INCGAS:',20(1X,L1))
         WRITE(3,1)BUFFER(1:DBRECL)
         WRITE(BUFFER,209)GASNAM(I),(ALLISO(J,I),J=1,DBNISO(I))
209      FORMAT(' #',1A8,' ALLISO:',20(1X,L1))
         WRITE(3,1)BUFFER(1:DBRECL)
207   CONTINUE
      
      
      TCORS2=1.439*(TEMP-296.)/(296.*TEMP)
C     For each bin ...
      NBIN = INT((VMAX-VMIN)/BINSIZ)
      DO IBIN=1,NBIN
       TSTROUT(IBIN)=0.
       TXLSE(IBIN)=0.
       TXWIDA(IBIN)=0.
       TXWIDS(IBIN)=0.
       TXTDEP(IBIN)=0.
      ENDDO

      JBIN=-1

      WRITE(4,*)'NBIN = ',NBIN

C     Find point in database where wavenumber is equal to VMIN
      PRINT*,'Finding first line in database with v >= ',VMIN
151   READ(DBLUN,1,END=999)BUFFER(1:DBRECL)
      CALL RDLINE(BUFFER)
      IF(LNWAVE.LT.VMIN)GOTO 151

      PRINT*,'Starting transfer'
      print*,LNSTR,STRMIN

      TS1 = 1.0-EXP(-1.439*SNGL(LNWAVE)/TEMP)
      TS2 = 1.0-EXP(-1.439*SNGL(LNWAVE)/296.0)
      TSTIM=1.0
      IF(TS2.NE.0) TSTIM = TS1/TS2
      K = -1
      DO J=1,DBNISO(LOCID(LNID))
        IF(LNISO.EQ.DBISO(J,LOCID(LNID)))K=J
      ENDDO                
      TCORS1=PARTF(LOCID(LNID),K,TEMP,IPTF)
      LNABSCO = LOG(LNSTR)+LOG(TCORS1)+
     & (TCORS2*LNLSE)+LOG(TSTIM)
      LNSTR1=DEXP(LNABSCO)

      IF(LNSTR1.GE.STRMIN)THEN
        WRITE(3,1)BUFFER(1:DBRECL)
      ELSE
        IBIN=1+INT((LNWAVE-VMIN)/BINSIZ)

        WEIGHT=LNSTR1*1E-20
        TSTROUT(IBIN)=TSTROUT(IBIN)+WEIGHT
        TXWIDA(IBIN)=TXWIDA(IBIN)+LNWIDA*WEIGHT
        TXWIDS(IBIN)=TXWIDS(IBIN)+LNWIDS*WEIGHT
        TXLSE(IBIN)=TXLSE(IBIN)+LNLSE*WEIGHT
        TXTDEP(IBIN)=TXTDEP(IBIN)+LNTDEP*WEIGHT

      ENDIF

110   READ(DBLUN,1,END=999)BUFFER
      CALL RDLINE(BUFFER)
      IF(LNWAVE.GT.VMAX)GOTO 998

      TS1 = 1.0-EXP(-1.439*SNGL(LNWAVE)/TEMP)
      TS2 = 1.0-EXP(-1.439*SNGL(LNWAVE)/296.0)
      TSTIM=1.0
      IF(TS2.NE.0) TSTIM = TS1/TS2
      K = -1
C      print*,LNID,LNISO
      DO J=1,DBNISO(LOCID(LNID))
         IF(LNISO.EQ.DBISO(J,LOCID(LNID)))K=J
      ENDDO
C      print*,RELABU(K,LOCID(LNID))
      IF(ISOF.EQ.0.OR.(LNID.EQ.IDF.AND.LNISO.EQ.ISOF))THEN
       TCORS1=PARTF(LOCID(LNID),K,TEMP,IPTF)
       LNABSCO = LOG(LNSTR)+LOG(TCORS1)+
     &  (TCORS2*LNLSE)+LOG(TSTIM)
       LNSTR1=DEXP(LNABSCO)
       IF(ISOF.GT.0)LNSTR1=LNSTR1/RELABU(K,LOCID(LNID))
       IF(LNSTR1.GE.STRMIN)THEN
         WRITE(3,1)BUFFER(1:DBRECL)
       ELSE
         IBIN=1+INT((LNWAVE-VMIN)/BINSIZ)
         IF(IBIN.NE.JBIN)THEN
          print*,'Processing bin ',IBIN,' out of ',NBIN
          JBIN=IBIN
         ENDIF
         WEIGHT=LNSTR1*1E-20
         TSTROUT(IBIN)=TSTROUT(IBIN)+WEIGHT
         TXWIDA(IBIN)=TXWIDA(IBIN)+LNWIDA*WEIGHT
         TXWIDS(IBIN)=TXWIDS(IBIN)+LNWIDS*WEIGHT
         TXLSE(IBIN)=TXLSE(IBIN)+LNLSE*WEIGHT
         TXTDEP(IBIN)=TXTDEP(IBIN)+LNTDEP*WEIGHT

       ENDIF
      ENDIF
      GOTO 110
998   PRINT*,'Transfer complete'
999   CONTINUE
      CLOSE(3)
      CLOSE(DBLUN)

      DO IBIN=1,NBIN
        IF(TSTROUT(IBIN).GT.0)THEN
         TXWIDA(IBIN)=TXWIDA(IBIN)/TSTROUT(IBIN)
         TXWIDS(IBIN)=TXWIDS(IBIN)/TSTROUT(IBIN)
         TXLSE(IBIN)=TXLSE(IBIN)/TSTROUT(IBIN)
         TXTDEP(IBIN)=TXTDEP(IBIN)/TSTROUT(IBIN)
         TSTROUT(IBIN)=TSTROUT(IBIN)*1E-27
        ELSE
         TXWIDA(IBIN)=0.1
         TXWIDS(IBIN)=0.1
         TXLSE(IBIN)=888.0
         TXTDEP(IBIN)=0.5
        ENDIF
      ENDDO


      DO 100 IBIN=1,NBIN
        VLOW = VMIN + FLOAT(IBIN-1)*BINSIZ
        VHIGH = VMIN + FLOAT(IBIN)*BINSIZ
        VMEAN=0.5*(VLOW+VHIGH)
        STROUT = TSTROUT(IBIN)
        XWIDA = TXWIDA(IBIN)
        XWIDS = TXWIDS(IBIN)
        XLSE = TXLSE(IBIN)
        XTDEP = TXTDEP(IBIN)

        WRITE(4,*)VMEAN,STROUT,XLSE,XWIDA,XWIDS,XTDEP

100   CONTINUE

      CLOSE(4)


      STOP

      END
************************************************************************
************************************************************************
