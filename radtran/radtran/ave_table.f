      PROGRAM AVE_TABLE
C     $Id: ave_table.f,v 1.12 2011-06-17 15:40:25 irwin Exp $
C***********************************************************************
C_TITL:	AVE_TABLE.f
C
C_DESC: Reads in a k-table and smooths it to desired resolution. For best
C	results, the smoothed ktable should be at least Nyquist-sampled
C	(e.g. input DELV = 2.0, input FWHM = 5.0).
C
C_ARGS:	See the definitions below.
C
C_FILE:	unit=30		(=LUN0), Input <>.kta file.
C	unit=31		(=LUN1), Output <>.kta file.
C
C_CALL:	
C	prompt		Prompts the user for input.
C	file		Force a file extension.
C
C_HIST:	19/4/01	PGJI	ORIGINAL VERSION.
C	???????	PGJI	Major modification after 'features' arose with CH4
C			at 1300cm-1.
C	6nov03	PGJI	
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
C ../includes/bincom.f stores the line bin variables (including NLINES,
C FSTLIN, LSTLIN) and band parameters.
      INCLUDE '../includes/dbcom.f'
C ../includes/dbcom.f stores the line database variables.
      INCLUDE '../includes/pathcom.f'
C ../includes/parcom.f stores the parameter values such as MAXLAY,
C MAXPATH, LIMGAS, etc.

      INTEGER LUN0,LUN1,MBIN,NC,IC0,JC,I1
      PARAMETER (LUN0=30,LUN1=31,MBIN=40)
C MAXOUT: the maximum number of output points.

      INTEGER J,J1,N
      REAL COEF(4),XA(4),YA(4),XMIN,XMAX,X1,X2
      INTEGER I,IP,IT,LOOP,IBIN
C I, IP, IT, LOOP: Incremators.
C IBIN: (=INT(0.5*NBIN))
      INTEGER IREC,IREC0,IREC2,IRECL,ISYS
C IREC: record incremator corresponding to input variable TABLE
C IREC2: record incremator corresponding to output TABLE
C IRECL: record length of input file; number of 4-byte words per record --
C =2 for old files and =1 for the newer files.
      INTEGER NP,NT,NG
C NP: Number of pressures.
C NT: Number of temperatures.
C NG: Number of ordinates in k-distribution.

      REAL TOT_TIME
C TOT_TIME: The total system time it took to complete program execution.
      real rate
      integer c1,c2,cr,TIME1,TIME2,cm
C TIME: Temporary variable returned by GETTIME containing the system time.
C TIME1: System time at the beginning of program execution.
C TIME2: System time at the end of program execution.

      INTEGER NPOINT0,IDGAS1,ISOGAS1,IC
C IDGAS1: The local gas identifier.
C ISOGAS1: The local gas-isotopic identifier; if zero all isotopes of the
C gas are included.
      REAL V,V1,XK,F
      REAL DELV0,FWHM0,XOFF,XREF,SUMC
      REAL VMIN0,VMAX,VMAX0,V_START,V_END,FRAC

      REAL TEMP1(MAXK),PRESS1(MAXK)
C TEMP1: Temperature [Kelvin].
C PRESS1: Pressure [atm].
      REAL G_ORD(MAXG),DEL_G(MAXG)
C G_ORD: Gauss-Legendre ordinates for calculating the k-distribution.
C DEL_G: Gauss-Legendre weights for integration.
      REAL TABLE(MAXK,MAXK,MAXG),TABLES(MBIN,MAXK,MAXK,MAXG)
      REAL TMP,GW(MBIN*20),GK(MBIN*20),GSW(MBIN*20)

      CHARACTER*100 KTAFIL,OPFILE1
      CHARACTER*1 ANS

      LOGICAL RESORT
C******************************** CODE *********************************

C Obtain the system time for use in determining the total time used by the
C program since execution.
      CALL system_clock(count_rate=cr)
      CALL system_clock(count_max=cm)
      rate = REAL(cr)

      CALL system_clock(TIME1)

      IRECL = ISYS()

      CALL PROMPT('Enter input filename : ')
      READ(5,23)OPFILE1
23    FORMAT(A)
      CALL FILE(OPFILE1,KTAFIL,'kta')

      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='OLD',ACCESS='DIRECT',
     1 RECL=IRECL)
      READ(LUN0,REC=1)IREC0
      READ(LUN0,REC=2)NPOINT
      READ(LUN0,REC=3)VMIN
      READ(LUN0,REC=4)DELV
      READ(LUN0,REC=5)FWHM
      READ(LUN0,REC=6)NP
      READ(LUN0,REC=7)NT
      READ(LUN0,REC=8)NG
      READ(LUN0,REC=9)IDGAS1
      READ(LUN0,REC=10)ISOGAS1
      IREC = 11
      VMAX = VMIN + (NPOINT - 1)*DELV
      WRITE(*,*)' '
      WRITE(*,*)'Input file ... '
      WRITE(*,*)'  NPOINT = ',NPOINT
      WRITE(*,*)'  VMIN, VMAX = ',VMIN,VMAX
      WRITE(*,*)'  DELV, FWHM = ',DELV,FWHM
      WRITE(*,*)'  NP, NT, NG = ',NP,NT,NG
      WRITE(*,*)'  Gas ID, ISO = ',IDGAS1,ISOGAS1
      WRITE(*,*)' '
      DO 299 J=1,NG
        READ(LUN0,REC=IREC)G_ORD(J)
        IREC = IREC + 1
299   CONTINUE
      WRITE(*,*)'  G-ordinates, weights : '
      DO 399 J=1,NG
        READ(LUN0,REC=IREC)DEL_G(J)
        WRITE(*,*)' ',G_ORD(j),DEL_G(j)
        IREC = IREC + 1
399   CONTINUE
      IREC = 11 + 2*NG + 2
      WRITE(*,*)'  Pressures : '
      DO 301 J=1,NP
        READ(LUN0,REC=IREC)PRESS1(J)
        WRITE(*,*)' ',PRESS1(j)
        PRESS1(J) = LOG(PRESS1(J))
        IREC = IREC + 1
301   CONTINUE
      WRITE(*,*)'  Temperatures : '
      DO 302 J=1,NT
        READ(LUN0,REC=IREC)TEMP1(J)
        WRITE(*,*)' ',TEMP1(j)
        IREC = IREC + 1
302   CONTINUE
      IREC = IREC0

      IF(DELV.NE.FWHM.AND.(DELV*2.0).NE.FWHM)THEN
        WRITE(*,*)'AVE_TABLE.f :: **ERROR** Table to be averaged must'
        WRITE(*,*)'have either FWHM=DELV, or FWHM=2*DELV.'
        WRITE(*,*)' '
        WRITE(*,*)'Stopping program.'
        STOP
      ENDIF

      CALL PROMPT('Enter output (smoothed) filename : ')
      READ(5,23)OPFILE1

      CALL PROMPT('Enter desired vmin, vmax, delv, fwhm(square) : ')
      READ*,VMIN0,VMAX0,DELV0,FWHM0

C     Snap FWHM0 and DELV0 to multiples of existing grid 
      FWHM0 = FWHM*INT(FWHM0/FWHM + 0.5)
      DELV0 = DELV*INT(DELV0/DELV + 0.5)
      NBIN = INT(FWHM0/FWHM)

      WRITE(*,*)'Snapped FHWM, DELV and NBIN',FWHM0,DELV0,NBIN
      IF(NBIN.GT.MBIN)THEN
        WRITE(*,*)'AVE_TABLE.f :: **ERROR** NBIN > MBIN'
        WRITE(*,*)'Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'NBIN, MBIN = ',NBIN,MBIN
        STOP
      ENDIF

      IBIN = INT(0.5*NBIN)
      IF(NBIN.EQ.2*IBIN)THEN                  ! if even number of bins
        WRITE(*,*)'Even number of averaging bins/output bin'
        XOFF = 0.5*FWHM	 ! Central output wavenumber is slightly offset
      ELSE
        print*,'Odd number of averaging bins/output bin'
        XOFF = 0.0
      ENDIF                                   

      WRITE(*,*)' '
      WRITE(*,*)'NBIN, IBIN, XOFF = ',NBIN,IBIN,XOFF

C     Add offsets to wavenumber range
      VMIN0 = VMIN0+XOFF 
      VMAX0 = VMAX0+XOFF 

C     Check range to make sure we don't run off either end of table
      V_START = (VMIN0+XOFF) - IBIN*DELV    ! central nu of lowest bin
      V_END = (VMAX0-XOFF)+IBIN*DELV	   ! central nu of highest bin

      IF(V_START.LT.VMIN.OR.V_END.GT.VMAX)THEN
        WRITE(*,*)'AVE_TABLE.F :: **ERROR** VMIN and/or VMAX outside'
        WRITE(*,*)'allowed limits'
        WRITE(*,*)' '
        WRITE(*,*)'Allowed ...'
        IF(V_START.LT.VMIN)THEN
          VMIN0 = VMIN0 + (VMIN-V_START) 
        ENDIF
        IF(V_END.GT.VMAX)THEN
          VMAX0 = VMAX0 - (V_END-VMAX) 
        ENDIF
        WRITE(*,*)'VMIN0 : ',VMIN0
        WRITE(*,*)'VMAX0 : ',VMAX0
      ENDIF

      NPOINT0 = 1 + INT((VMAX0 - VMIN0)/DELV0)

      WRITE(*,*)' '
      WRITE(*,*)'Output file ...'
      WRITE(*,*)'  NPOINT = ',NPOINT0
      WRITE(*,*)'  VMIN, VMAX = ',VMIN0,VMIN0 + (NPOINT0 - 1)*DELV0
      WRITE(*,*)'  DELV, FWHM = ',DELV0,FWHM0

      CALL FILE(OPFILE1,KTAFIL,'kta')
      OPEN(UNIT=LUN1,FILE=KTAFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=IRECL)
      WRITE(LUN1,REC=1)IREC0
      WRITE(LUN1,REC=2)NPOINT0
      WRITE(LUN1,REC=3)VMIN0
      WRITE(LUN1,REC=4)DELV0
      WRITE(LUN1,REC=5)FWHM0
      WRITE(LUN1,REC=6)NP
      WRITE(LUN1,REC=7)NT
      WRITE(LUN1,REC=8)NG
      WRITE(LUN1,REC=9)IDGAS1
      WRITE(LUN1,REC=10)ISOGAS1
      IREC = 11
      DO 199 J=1,NG
        WRITE(LUN1,REC=IREC)G_ORD(J)
        IREC = IREC + 1
199   CONTINUE
      DO 200 J=1,NG
        WRITE(LUN1,REC=IREC)DEL_G(J)
        IREC = IREC + 1
200   CONTINUE
      IREC = 11 + 2*NG + 2
      DO 201 J=1,NP
        WRITE(LUN1,REC=IREC)EXP(PRESS1(J))
        IREC = IREC + 1
201   CONTINUE
      DO 202 J=1,NT
        WRITE(LUN1,REC=IREC)TEMP1(J)
        IREC = IREC + 1
202   CONTINUE


      DO 1000 I=1,NPOINT0
        WRITE(*,*)'AVE_TABLE.f :: CALCULATING : ',I,' OF ',NPOINT0

        V = VMIN0+(I-1)*DELV0
        IREC2 = IREC0 + NP*NT*NG*(I - 1)

        DO IP=1,NP
          DO IT=1,NT
           DO LOOP=1,NG
            TABLE(IP,IT,LOOP) = 0.0
           ENDDO
          ENDDO
        ENDDO

        DO 101 J=1,NBIN
          V1 = (V+XOFF) + (J-1-IBIN)*FWHM   ! central nu of next bin
          J1 = 1 + INT((V1 - VMIN)/DELV)
C          PRINT*,V,I,J,V1,J1

C Read in the table for the wavenumber region concerned.
          IREC = IREC0 + NP*NT*NG*(J1 - 1)

          DO 901 IP=1,NP
            DO 902 IT=1,NT
              DO 903 LOOP=1,NG
                READ(LUN0,REC=IREC)TABLES(J,IP,IT,LOOP)
                IREC = IREC + 1
903           CONTINUE
902         CONTINUE
901       CONTINUE
101     CONTINUE

        DO 801 IP=1,NP
          DO 802 IT=1,NT

            NC = 0
            DO 803 J=1,NBIN 
              DO 804 LOOP=1,NG
                NC = NC + 1
                GW(NC) = DEL_G(LOOP)/(1.0*NBIN)
                GK(NC) = TABLES(J,IP,IT,LOOP)
C                print*,nc,gw(nc),gk(nc)
804           CONTINUE  
803         CONTINUE

C Now sort ...
777         RESORT = .FALSE.
            DO 805 IC=1,NC-1
              IF(GK(IC).GT.GK(IC+1))THEN
                TMP = GK(IC)
                GK(IC) = GK(IC+1)
                GK(IC+1) = TMP
                TMP = GW(IC)
                GW(IC) = GW(IC+1)
                GW(IC+1) = TMP
                RESORT = .TRUE.
              ENDIF
805         CONTINUE
            IF(RESORT) GOTO 777

            GSW(1) = GW(1)
            DO 806 IC=2,NC
              GSW(IC) = GSW(IC-1) + GW(IC)
806         CONTINUE

C            do ic=1,nc
C             print*,gsw(ic),gw(ic),gk(ic)
C            enddo

            LOOP = 1
            XREF = DEL_G(LOOP)
            SUMC = 0.0
            DO 807 IC=1,NC
              IF(GSW(IC).GE.XREF)THEN
                IF(IC.GT.1)THEN
                 FRAC = (GSW(IC)-XREF)/(GSW(IC)-GSW(IC-1))
                ELSE
                 FRAC = (GSW(IC)-XREF)/GSW(IC)
                ENDIF
                SUMC = SUMC + (1.0-FRAC)*GK(IC)*GW(IC)
                TABLE(IP,IT,LOOP) = SUMC/DEL_G(LOOP)
                IF(LOOP.LT.NG) THEN
                  SUMC = FRAC*GK(IC)*GW(IC)
                  LOOP = LOOP + 1
                  XREF = XREF + DEL_G(LOOP)
                ENDIF
              ELSE
                SUMC = SUMC + GK(IC)*GW(IC)
              ENDIF
807         CONTINUE
            IF(TABLE(IP,IT,NG).EQ.0.0)THEN
              TABLE(IP,IT,LOOP) = SUMC/DEL_G(NG)     ! Catch last point
cc              WRITE(*,*)'AVE_TABLE.f :: NG, IP, IT, V = ',NG,IP,IT,V 
            ENDIF  

            DO 808 LOOP=1,NG
              WRITE(LUN1,REC=IREC2)TABLE(IP,IT,LOOP)
C              print*,IP,IT,LOOP,TABLE(IP,IT,LOOP)
              IREC2 = IREC2 + 1
808         CONTINUE
C            read(5,23)ans

802       CONTINUE
801     CONTINUE

1000  CONTINUE

      CLOSE(LUN0)
      CLOSE(LUN1) 

      WRITE(*,*)' AVE_TABLE.F :: calculation complete.'
      CALL system_clock(TIME2)
      TOT_TIME = (TIME2 - TIME1)/rate
      WRITE(*,244)TOT_TIME
244   FORMAT(/' Elapsed time (sec) = ',F8.1)

      END
C***********************************************************************
C***********************************************************************



      SUBROUTINE FILE(NAME1,NAME2,EXT)
C     $Id: ave_table.f,v 1.12 2011-06-17 15:40:25 irwin Exp $
C***********************************************************************
C_DESC:	Forces correct VMS style file extension for a filename
C	(i.e. assumes a <4 character extension after a dot separator). 
C
C_HIST:	1/1/90    Original Version:       SBC
C	3/10/94   Updated Header          PGJI
C
C***********************************************************************
      INTEGER I,L,LE
      CHARACTER*(*) NAME1,NAME2,EXT

C First copy the file name to the second array and remove any leading 
C spaces
      NAME2 = NAME1
      CALL REMSP(NAME2)
C Now stepping back from the end of the name through the last four
C characters
      LE = LEN(NAME2)
      IF(LE.LT.1)RETURN
      DO 40 L=LE,1,-1
      IF(NAME2(L:L).NE.' ')THEN
        LE = L
        GOTO 50
        ENDIF
40    CONTINUE
      RETURN
50    CONTINUE
      L = LE - 3
      IF(L.LT.1)L = 1
      DO 10 I=LE,L,-1
C If a directory tree delimeter is found just add the extension to the end 
C of the name
      IF(NAME2(I:I).EQ.'/'.OR.NAME2(I:I).EQ.']')GOTO 30
C If a dot is found overwrite any existsing extension
      IF(NAME2(I:I).EQ.'.')THEN
        NAME2(I+1:I+3) = EXT
        GOTO 20
        ENDIF
10    CONTINUE
C If no dot is found in last four characters add the extension to the end 
C of the filename
30    NAME2(LE+1:LE+1) = '.'
      NAME2(LE+2:LE+4) = EXT

20    RETURN

      END 
C***********************************************************************
C***********************************************************************



      SUBROUTINE REMSP(TEXT)
C     $Id: ave_table.f,v 1.12 2011-06-17 15:40:25 irwin Exp $
C***********************************************************************
C_DESC:	Removes leading spaces from text string
C
C_HIST:	1/1/90	SBC	Original Version
C	3/10/94	PGJI
C***********************************************************************
      CHARACTER TEXT*(*)
      INTEGER I,J,K,L

      J = LEN(TEXT)
      DO 100 I=1,J
        IF(TEXT(I:I).NE.' ')GOTO 10
100   CONTINUE
      RETURN
10    CONTINUE
      IF(I.EQ.1)RETURN
      DO 20 K=I,J
        L = K - I + 1
        TEXT(L:L) = TEXT(K:K)
20    CONTINUE
      DO 30 K=J-I+2,J
        TEXT(K:K) = ' '
30    CONTINUE
      END


      INTEGER FUNCTION PINDEX(X,NX,X1)
      INTEGER NX
      REAL X(NX),X1
      IF(X1.GT.X(NX))THEN
       PRINT*,'Warning in INDEX. X1>XMAX'
       X1 = X(NX)
      END IF
      IF(X1.LT.X(1))THEN
       PRINT*,'Warning in INDEX. X1>XMIN'
       X1 = X(1)
      ENDIF
 
      J = 0

      DO 10 I=1,NX-1
       IF(X1.GE.X(I))THEN
        J = I
       ELSE
        GOTO 20
       ENDIF
10    CONTINUE
20    CONTINUE

      PINDEX = J

      END
C***********************************************************************
C***********************************************************************

