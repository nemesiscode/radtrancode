      PROGRAM SELECTTEMPSEQTOPBUFFLOOP
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

      INTEGER MINSTR,MAXSTR,ISWAP,IDF,ISOF
      DOUBLE PRECISION LIMSTR,CC,LNSTR1,XS
      PARAMETER (MINSTR=-323,MAXSTR=308,CC=10.)

      INTEGER NLINES(MINSTR:MAXSTR,MAXISO,MAXDGAS)
      INTEGER TOTLIN(MINSTR:MAXSTR,MAXDGAS),ILUN(20),JLUN(20)
      DOUBLE PRECISION STRLOSE,STRKEEP,STROUT,WEIGHT
      REAL XWIDA,XWIDS,XLSE,XTDEP
      INTEGER NSAV,MAXLIN,NTEMP,ITEMP
      PARAMETER(MAXLIN = 2000000)
      INTEGER IDX(MAXLIN),IFLAG(MAXLIN),LE,JE,I1,I2,NHEAD
      INTEGER I,J,K,ID,ISO,IBIN,LINE,FSTLIN,LSTLIN,NBIN,IPTF
      INTEGER FIRST(2),LAST(2),NPAR,IERROR,READI,NLIN,NKEEP,NLOSE
      DOUBLE PRECISION SSL(MAXLIN),SSLTMP(MAXLIN)

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
      REAL TEMP,TS1,TS2,PARTF,TMIN,TMAX,TEMP1(20)
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

      WRITE(*,*)'Lines whose summed strength is less than n% of'
      WRITE(*,*)'the total sum strength in a bin may be omitted'
22    CALL PROMPT('Enter bin width and number of lines to extract')
      READ(*,*)BINSIZ,NSAV
      IF(BINSIZ.LT.0.01)THEN
        WRITE(*,*)'Bin size is too small'
        GOTO 22
      ENDIF

      CALL PROMPT('Enter number of temperatures : ')
      READ(*,*)NTEMP

      CALL PROMPT('Enter Tmin, Tmax : ')
      READ(*,*)TMIN,TMAX

      DO I=1,NTEMP
       TEMP1(I)=TMIN+(I-1)*(TMAX-TMIN)/FLOAT(NTEMP-1)
      ENDDO
      IF(NTEMP.EQ.1)TEMP1(1)=TMIN

      CALL PROMPT('Enter IPTF (Partition function flag) : ')
      READ*,IPTF


      CALL PROMPT('Enter new file root name')
      READ(*,21)OPNAME
21    FORMAT(A)

      DO I=1,LEN(OPNAME)
       IF(OPNAME(I:I).NE.' ')J=I
      ENDDO

      OPNAME(J+1:J+2)='00'

      LE=J+2

      print*,LE,LE

      OPEN(12,FILE='TEMP.DAT',STATUS='UNKNOWN')
      NHEAD=0
      OPEN(DBLUN,FILE=DBFILE,STATUS='OLD')
204   READ(DBLUN,1,END=999)BUFFER(1:DBRECL)
      IF(BUFFER(2:2).EQ.'#'.OR.BUFFER(1:1).EQ.'#')THEN
        WRITE(12,1)BUFFER(1:DBRECL)
        NHEAD=NHEAD+1
        GOTO 204
      ENDIF
      CLOSE(12)

      NBIN = INT((VMAX-VMIN)/BINSIZ)

      CALL PROMPT('Enter ID,ISO to put in header : ')
      READ(5,*)IDF,ISOF

      DO 1001 ITEMP=1,NTEMP
       ILUN(ITEMP)=50+ITEMP-1
       JLUN(ITEMP)=80+ITEMP-1
       TEMP=TEMP1(ITEMP)
       print*,ITEMP,TEMP
       I1=INT(ITEMP/10)
       I2=ITEMP-10*I1
       PRINT*,I1,I2
       JE=LE-1
       OPNAME(JE:JE)=CHAR(I1+48)
       OPNAME(LE:LE)=CHAR(I2+48)

       CALL FILE(OPNAME,OPSTR,'par')
       OPEN(UNIT=ILUN(ITEMP),FILE=OPSTR,STATUS='UNKNOWN')
       CALL FILE(OPNAME,OPSTR,'lco')
       OPEN(UNIT=JLUN(ITEMP),FILE=OPSTR,STATUS='UNKNOWN')

       WRITE(BUFFER,114)
114    FORMAT(' # data records written by SELECTTEMPSEQ*')
       WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)
       WRITE(BUFFER,201)
201    FORMAT(' # original data base file:')
       WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)
       WRITE(BUFFER,202)DBFILE
202    FORMAT(' #  ',A)
       WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)

       WRITE(JLUN(ITEMP),*)'VMIN, VMAX = ',VMIN,VMAX
       WRITE(JLUN(ITEMP),*)'BINSIZ,NSAV = ',BINSIZ,NSAV
       WRITE(JLUN(ITEMP),*)'Calc. Temp, IPTF = ',TEMP,IPTF
       WRITE(JLUN(ITEMP),*)'ID, ISO = ',IDF,ISOF
       WRITE(JLUN(ITEMP),*)'NBIN = ',NBIN
       WRITE(BUFFER,203)
203    FORMAT(' # original data base header if any:')
       WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)

       OPEN(12,FILE='TEMP.DAT',STATUS='OLD')
       DO I=1,NHEAD
        READ(12,1)BUFFER
        WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)
       ENDDO
       CLOSE(12)

       WRITE(BUFFER,205)
205    FORMAT(' # selection criteria:')
       WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)
       WRITE(BUFFER,206)VMIN,VMAX,BINSIZ
206    FORMAT(' # wavenumber range =',F10.3,' -',F10.3,' with',F8.3,
     1 ' cm-1 bins')
       WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)

       DO 207 I=1,MAXDGAS
         IF(DBNISO(I).LT.1)GOTO 207
         WRITE(BUFFER,208)GASNAM(I),(INCGAS(J,I),J=1,DBNISO(I))
208      FORMAT(' #',1A8,' INCGAS:',20(1X,L1))
         WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)
         WRITE(BUFFER,209)GASNAM(I),(ALLISO(J,I),J=1,DBNISO(I))
209      FORMAT(' #',1A8,' ALLISO:',20(1X,L1))
         WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)
207    CONTINUE

1001  CONTINUE


C     Find point in database where wavenumber is equal to VMIN
      PRINT*,'Finding first line in database with v >= ',VMIN
151   READ(DBLUN,1,END=999)BUFFER(1:DBRECL)
      CALL RDLINE(BUFFER)
      IF(LNWAVE.LT.VMIN)GOTO 151

      PRINT*,'NBIN = ',NBIN

      DO 100 IBIN=1,NBIN

        print*,'ibin = ',ibin
        VLOW = VMIN + FLOAT(IBIN-1)*BINSIZ
        VHIGH = VMIN + FLOAT(IBIN)*BINSIZ
        VMEAN=0.5*(VLOW+VHIGH)
        WRITE(*,121)IBIN,VLOW,VHIGH
121     FORMAT(I5,F12.3,'  -',F12.3)

        NLIN=1
C       Write line into memory
117     CONTINUE
C        print*,nlin,lnwave,lnstr
        CALL ADD_LINE_BUFF(BUFFER,NLIN)
        READ(DBLUN,1,END=999)BUFFER(1:DBRECL)
        CALL RDLINE(BUFFER)
        IF(LNWAVE.LT.VHIGH)THEN
         NLIN=NLIN+1
         GOTO 117
        ENDIF
        BUFFERSAV=BUFFER

        PRINT*,'NLIN,NSAV = ',NLIN,NSAV


        DO 1000 ITEMP=1,NTEMP     

         TEMP=TEMP1(ITEMP)
         TCORS2=1.439*(TEMP-296.)/(296.*TEMP)

         DO 303 LINE=1,NLIN
          CALL READ_LINE_BUFF(BUFFER,LINE)
          CALL RDLINE(BUFFER)

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
     &         (TCORS2*LNLSE)+LOG(TSTIM)
          LNSTR1=DEXP(LNABSCO)
          SSL(LINE)=LNSTR1
          IDX(LINE)=LINE
          IFLAG(LINE)=0
C         print*,nlin,ssl(nlin),idx(nlin),iflag(nlin)
303      CONTINUE

         IF(NLIN.GT.NSAV)THEN
C          Now sort the lines into order of strength
C          print*,'SORT, NLIN = ',NLIN

C          If we've already sorted the lines for a previous temperature
C          start from that sorted list as it'll be faster to reshuffle for the 
C          new temperature from there.           
           IF(ITEMP.GT.1)THEN
            DO I=1,NLIN
             SSLTMP(I)=SSL(IDX(I))
            ENDDO
            DO I=1,NLIN
             SSL(I)=SSLTMP(I)
            ENDDO
           ENDIF

101        CONTINUE
           ISWAP=0
           DO I=1,NLIN-1
            IF(SSL(I).LT.SSL(I+1))THEN
             XS=SSL(I)
             SSL(I)=SSL(I+1)
             SSL(I+1)=XS
             J=IDX(I)
             IDX(I)=IDX(I+1)
             IDX(I+1)=J
             ISWAP=1
            ENDIF
           ENDDO
           IF(ISWAP.GT.0)GOTO 101

C          Mark lines to extract
           DO I=1,NSAV
            IFLAG(IDX(I))=1
C            print*,I,SSL(I),IDX(I)
           ENDDO
         ENDIF

C        Initialise continuum
         STROUT=0.0
         XWIDA=0.
         XWIDS=0.
         XLSE=0.
         XTDEP=0.

         DO 160 LINE=1,NLIN
           CALL READ_LINE_BUFF(BUFFER,LINE)
      
           IF(NLIN.GT.NSAV.AND.IFLAG(LINE).EQ.0)THEN

            CALL RDLINE(BUFFER)

            TS1 = 1.0-EXP(-1.439*SNGL(LNWAVE)/TEMP)
            TS2 = 1.0-EXP(-1.439*SNGL(LNWAVE)/296.0)
            TSTIM=1.0
            IF(TS2.NE.0) TSTIM = TS1/TS2
            K = -1
            DO J=1,DBNISO(LOCID(LNID))
             IF(LNISO.EQ.DBISO(J,LOCID(LNID)))K=J
            ENDDO                
            TCORS1=PARTF(LOCID(LNID),K,TEMP,IPTF)
            LNABSCO = LOG(LNSTR)+LOG(TCORS1)+(TCORS2*LNLSE)+LOG(TSTIM)
            LNSTR1=DEXP(DBLE(LNABSCO))

            WEIGHT=LNSTR1*1E-20
            STROUT=STROUT+WEIGHT
C            print*,weight,strout
            XWIDA=XWIDA+LNWIDA*WEIGHT
            XWIDS=XWIDS+LNWIDS*WEIGHT
            XLSE=XLSE+LNLSE*WEIGHT
            XTDEP=XTDEP+LNTDEP*WEIGHT
           ELSE
            WRITE(ILUN(ITEMP),1)BUFFER(1:DBRECL)
           ENDIF 
160      CONTINUE

         IF(STROUT.GT.0)THEN
           XWIDA=XWIDA/STROUT
           XWIDS=XWIDS/STROUT
           XLSE=XLSE/STROUT
           XTDEP=XTDEP/STROUT
           STROUT=STROUT*1E-27
         ELSE
           XWIDA=0.1
           XWIDS=0.1
           XLSE=888.0
           XTDEP=0.5
         ENDIF
         WRITE(JLUN(ITEMP),*)VMEAN,STROUT,XLSE,XWIDA,
     1     XWIDS,XTDEP

1000    CONTINUE

        IF(VHIGH.GT.VMAX)GOTO 102

        BUFFER=BUFFERSAV
        CALL RDLINE(BUFFER)

100   CONTINUE
102   CONTINUE
999   CONTINUE

      CLOSE(DBLUN)
      DO ITEMP=1,NTEMP
       CLOSE(ILUN(ITEMP))
       CLOSE(JLUN(ITEMP))
      ENDDO

      END
************************************************************************
************************************************************************
