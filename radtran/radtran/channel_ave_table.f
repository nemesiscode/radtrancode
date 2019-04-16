      PROGRAM CHANNEL_AVE_TABLE
C***********************************************************************
C_TITL:	CHANNEL_AVE_TABLE.f
C
C_DESC: Reads in a k-table and converts it to a channel-averaged table.
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
C	6nov03		PGJI
C      20mar19		PGJI converted to do channel-averaging	
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
      PARAMETER (LUN0=30,LUN1=31,MBIN=1000)
C MAXOUT: the maximum number of output points.

      INTEGER J,J1,J2,N
      REAL COEF(4),XA(4),YA(4),XMIN,XMAX,X1,X2
      INTEGER I,IP,IT,LOOP,IBIN
C I, IP, IT, LOOP: Incremators.
C IBIN: (=INT(0.5*NBIN))
      INTEGER IREC,IREC0,IREC2,IRECL,ISYS
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
      REAL V,V1,F,V2,SUMK,KX,SUM1
      REAL DELV0,FWHM0,XOFF,XREF,SUMC
      REAL VMIN0,VMAX,VMAX0,V_START,V_END,FRAC

      REAL TEMP1(MAXK),PRESS1(MAXK)
C TEMP1: Temperature [Kelvin].
C PRESS1: Pressure [atm].
      REAL G_ORD(MAXG),DEL_G(MAXG)
C G_ORD: Gauss-Legendre ordinates for calculating the k-distribution.
C DEL_G: Gauss-Legendre weights for integration.
      REAL TMP,GW(MBIN*20),GK(MBIN*20)
      REAL XK,KG(MAXG),YY
      REAL VCEN(100),VFIL(100,1000),FIL(100,1000)
      INTEGER NSUB(100),NCONV,K,IREC1
      REAL VMIN1,DELV1,XF(1000),YF(1000)

      CHARACTER*100 KTAFIL,OPFILE1,FILFILE
      CHARACTER*1 ANS

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

      CALL PROMPT('Enter output (smoothed) filename : ')
      READ(5,23)OPFILE1

      CALL PROMPT('Enter name of filter file : ')
      READ(5,23)FILFILE

      OPEN(12,FILE=FILFILE,STATUS='old')
       READ(12,*)NCONV
       DO 442 K=1,NCONV
          read(12,*)vcen(K)
          read(12,*)nsub(K)
          do j=1,nsub(k)
           read(12,*)vfil(k,j),fil(k,j)
C           print*,vfil(k,j),fil(k,j)
          enddo
442    CONTINUE
      CLOSE(12)


      IREC1=IREC0+NCONV
      VMIN1=VCEN(1)
      VMAX=VCEN(NCONV)
      DELV1=-1.0
      FWHM=0.0
      CALL FILE(OPFILE1,KTAFIL,'kta')
      OPEN(UNIT=LUN1,FILE=KTAFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=IRECL)
      WRITE(LUN1,REC=1)IREC1
      WRITE(LUN1,REC=2)NCONV
      WRITE(LUN1,REC=3)VMIN1
      WRITE(LUN1,REC=4)DELV1
      WRITE(LUN1,REC=5)FWHM
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
      DO 303 J=1,NCONV
        WRITE(LUN1,REC=IREC)VCEN(J)
        IREC=IREC+1
303   CONTINUE


      DO 1000 I=1,NCONV

        v1 = vfil(i,1)
        v2 = vfil(i,nsub(i))
        j1 = 1 + int((v1-vmin)/delv)
        j2 = 1 + int((v2-vmin)/delv)
        
        print*,i,nsub(i)
        print*,v1,v2
        print*,j1,j2,1+j2-j1

        do j=1,nsub(i)
         xf(j)=vfil(i,j)
         yf(j)=fil(i,j)
         sumk=sumk+yf(j)
         print*,j,xf(j),yf(j)
        enddo

        WRITE(*,*)'CHANNEL_AVE_TABLE.f :: CALCULATING : ',
     &    I,' OF ',NCONV


        DO 901 IP=1,NP
         DO 902 IT=1,NT

          NC=0
          SUMK=0.
          DO 900 J=J1,J2

            V = VMIN + (J-1)*DELV
            CALL VERINT(XF,YF,NSUB(I),YY,V)
        
            IREC=IREC0 + NP*NT*NG*(J-1) + NT*NG*(IP-1) + NG*(IT-1)

            DO 903 LOOP=1,NG
                READ(LUN0,REC=IREC)KX
                IREC = IREC + 1
                NC=NC+1
                GK(NC)=KX
                GW(NC)=DEL_G(LOOP)*YY
                SUMK=SUMK+GW(NC)
903         CONTINUE
            
900       CONTINUE

          IF (NC.GT.MBIN*20)THEN
           PRINT*,'Error in Channel_ave_table - NC too big'
           PRINT*,NC,MBIN*20,MBIN
           STOP
          ENDIF

          SUM1=0.
          DO LOOP=1,NC
           GW(LOOP)=GW(LOOP)/SUMK
           SUM1=SUM1+GW(LOOP)
          ENDDO


C         Now sort k-coefficients
          CALL RANK(DEL_G,NG,GK,GW,NC,KG)


          IREC2=IREC1 + NP*NT*NG*(I-1) + NT*NG*(IP-1) + NG*(IT-1)

          DO 808 LOOP=1,NG
             WRITE(LUN1,REC=IREC2)KG(LOOP)
             IREC2 = IREC2 + 1
808       CONTINUE


902      CONTINUE
901     CONTINUE

1000  CONTINUE

      CLOSE(LUN0)
      CLOSE(LUN1) 

      WRITE(*,*)'CHANNEL_AVE_TABLE.F :: calculation complete.'
      CALL system_clock(TIME2)
      TOT_TIME = (TIME2 - TIME1)/rate
      WRITE(*,244)TOT_TIME
244   FORMAT(/' Elapsed time (sec) = ',F8.1)

      END
C***********************************************************************
C***********************************************************************

