      REAL FUNCTION PARTF(ID,ISO,TEMP,IPTF)
C     $Id: partf.f,v 1.3 2011-06-17 15:53:02 irwin Exp $
C-----------------------------------------------------------------------------
C_TITLE:  PARTF: computes total partition function
C
C_ARGS:   ID:INTEGER     local gas identifier
C         ISO:INTEGER    local isotope identifier
C         TEMP:REAL      temperature in Kelvin
C	  IPTF:INTEGER   Flag to indicate alternative processing
C
C_KEYS:   FUNCT,ATMOS,LBL
C
C_DESCR:  uses four term polynomial fits to compute ration of total partition
C         functions for use by line by line programs.
C         Representation is as described in Gimache et al 1990 and distributed
C         with HITRAN 91
C
C_DESCR1: The partition function is used to scale the line strengths to 
C	  different temperatures. It has two components, rotational and
C	  vibrational. Thus the line strength ratio is:
C
C	  S(T)      Z_v(T0)*Z_r(T0)
C	  ---    =  --------------  x other terms.
C         S(T0)      Z_v(T)*Z_r(T)
C
C         The rotational partition function is relatively easy to
C         calculate and the above may be simplified to:
C
C                       n
C	  S(T)      (T0)   Z_v(T0)
C	  ---    =  (--) x ------  x other terms.
C         S(T0)     ( T)   Z_v(T)
C
C         where n = 1.0 for linear molecules and 1.5 for non-linear
C         molecules
C
C         To fold in the vibrational part of the partition function
C         we use a 4th order polynomial scheme:
C
C	  S(T)      A + B*T0 + C*T0^2 + D*T0^3
C	  ---    =  --------------------------  x other terms.
C         S(T0)     A + B*T  + C*T^2  + D*T^3
C         
C         where A,B,C,D are coefficients which are read from the
C	  'gas_info.dat' file.
C
C         for linear molecules with no great vibrational contrib to the 
C	  partition function we set:
C
C		A=0
C		B=1
C		C=0
C		D=0
C	  
C         for linear molecules with no great vibrational contrib to the 
C         partition function we set:
C
C		A= -85.1670
C		B = 6.644465 
C		C = 0.0453985
C		D = -2.93911e-5
C
C         which is a reasonable fit to T**1.5 from 30  to 350K
C
C	  If IPARTF=1 then the code scans look-up tables of the partition
C         functions of different gases to see if any match and if so uses 
C         these data instead.
C 
C_CALLS:  
C
C_BUGS:   not all gases have valid fits
C
C_HIST:  23oct92  SBC  ORIGINAL VERSION with values copied from GENLN2 QTIPS
C        30JUN10  JML  Updated CH4 coefficients for three temperature ranges
C	 27OCT14  PGJI Modified IPTF=1 to indicate the need to read in tables
C   		       of tabulated partition functions.
C
C-----------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE '../includes/dbcom.f'
      INCLUDE '../includes/partfcom.f'
C-----------------------------------------------------------------------------
      INTEGER ID,ISO,JEXTRA,IPTF,I,J,NTAB
      REAL TEMP,XTEMP(MTABEXTRA),YTEMP(MTABEXTRA)
      CHARACTER*100 PARTFIL
C
      REAL T0,T02,T03,P1,WT
      PARAMETER (T0=296.,T02=T0*T0,T03=T0*T02)
C
      REAL PARTT0,PTMP

      IF(IPTF.EQ.0) THEN
       IF(ISO.GT.0)THEN
        PARTT0=DBQTA(ISO,ID)+DBQTB(ISO,ID)*T0
     1   +DBQTC(ISO,ID)*T02+DBQTD(ISO,ID)*T03
C
        PTMP=PARTT0/(DBQTA(ISO,ID) + TEMP*( DBQTB(ISO,ID) +
     1   TEMP*( DBQTC(ISO,ID) + TEMP*DBQTD(ISO,ID) ) ) )

       ELSE
        PARTT0=0.
        PTMP=0.
        DO 255 J=1,DBNISO(ID)
C         print*,RELABU(J,ID),(DBQTA(J,ID)+DBQTB(J,ID)*T0
C     1   +DBQTC(J,ID)*T02+DBQTD(J,ID)*T03),(DBQTA(J,ID) + 
C     3   TEMP*( DBQTB(J,ID) +
C     1   TEMP*( DBQTC(J,ID) + TEMP*DBQTD(J,ID) ) ) ),
C     3   (DBQTA(J,ID)+DBQTB(J,ID)*T0
C     1   +DBQTC(J,ID)*T02+DBQTD(J,ID)*T03)/(DBQTA(J,ID) +
C     3   TEMP*( DBQTB(J,ID) +
C     1   TEMP*( DBQTC(J,ID) + TEMP*DBQTD(J,ID) ) ) )

         PARTT0=PARTT0+RELABU(J,ID)*(DBQTA(J,ID)+DBQTB(J,ID)*T0
     1   +DBQTC(J,ID)*T02+DBQTD(J,ID)*T03)
         PTMP=PTMP+RELABU(J,ID)*(DBQTA(J,ID) + TEMP*( DBQTB(J,ID) +
     1   TEMP*( DBQTC(J,ID) + TEMP*DBQTD(J,ID) ) ) )

255     CONTINUE
        PTMP=PARTT0/PTMP
       ENDIF

      ELSE

       IF(IREADEXTRA.NE.-1)THEN
        PRINT*,'Extra partition functions database has not been'
        PRINT*,'initialised. Initialise here'
        CALL INIT_PARTF
       ENDIF

       IF(ISO.GT.0)THEN
        JEXTRA=-1
        DO 10 I=1,NPARTEXTRA
         IF(IDEXTRA(I).EQ.ID.AND.
     &		(ISOEXTRA(I).EQ.ISO))THEN
          JEXTRA=I
          NTAB=NTABEXTRA(I)
          DO 20 J=1,NTAB
           XTEMP(J)=TEMPEXTRA(I,J)
           YTEMP(J)=PARTFEXTRA(I,J)
20        CONTINUE
          CALL VERINT(XTEMP,YTEMP,NTAB,PTMP,TEMP)
         ENDIF

10      CONTINUE

        IF(JEXTRA.LT.0)THEN
         PRINT*,'Problem in PARTF - gas not found in partfextra.dat'
         PRINT*,'ID,ISO = ',ID,ISO
         PRINT*,'Stopping'
         STOP
        ENDIF
       ELSE
        JEXTRA=-2
        PTMP=0.
        WT=0.
C       We'll keep a tab on the total relative abundances in case the
C       number of isotopes in tabulated partition function file is less
C       than the total number of isotopes for the gas in the database.
        DO 15 I=1,NPARTEXTRA
         IF(IDEXTRA(I).EQ.ID)THEN
          NTAB=NTABEXTRA(I)
          DO 25 J=1,NTAB
           XTEMP(J)=TEMPEXTRA(I,J)
           YTEMP(J)=PARTFEXTRA(I,J)
25        CONTINUE
          CALL VERINT(XTEMP,YTEMP,NTAB,P1,TEMP)
C          print*,ID,ISOEXTRA(I),P1,RELABU(ISOEXTRA(I),ID),WT
          PTMP=PTMP+RELABU(ISOEXTRA(I),ID)*P1
          WT=WT+RELABU(ISOEXTRA(I),ID)
         ENDIF
15      CONTINUE

        IF(WT.LE.0)THEN
         PRINT*,'Problem in PARTF - no isotopes found in partfextra.dat'
         PRINT*,'for gas requested'
         PRINT*,'ID = ',ID
         PRINT*,'Stopping'
         STOP
        ENDIF

        PTMP=PTMP/WT

       ENDIF 
      ENDIF


      IF(PTMP.LT.0.0)THEN
       PRINT*,'PROBLEM IN PARTF - RESULT IS NEGATIVE'
       PRINT*,'ID,ISO,IPTF = ',ID,ISO,IPTF
       PRINT*,'TEMP = ',TEMP
       PRINT*,'DBQTA,DBQTB,DBQTC,DBQTD',DBQTA(ISO,ID),DBQTB(ISO,ID),
     &DBQTC(ISO,ID),DBQTD(ISO,ID)
       PRINT*,'PARTT0 = ',PARTT0
       PRINT*,'PARTF = ',PTMP
       PRINT*,'IPTF = ',IPTF
       PRINT*,'JEXTRA = ',JEXTRA
       PRINT*,'Stopping'
       STOP
      ENDIF

      PARTF=PTMP

      RETURN
      END
