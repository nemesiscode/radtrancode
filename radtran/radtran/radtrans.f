      PROGRAM RADTRANS
C     $Id: radtrans.f,v 1.8 2011-06-17 15:40:27 irwin Exp $
C***********************************************************************
C_TITL:	RADTRANS.f
C
C_AUTH:	P. G. J. Irwin and S. B. Calcutt,
C	University of Oxford, Department of Physics,
C	Atmospheric, Oceanic and Planetary Physics,
C	Clarendon Laboratory,
C	Parks Rd,
C	Oxford OX1 3PU
C	United Kingdom
C
C_DESC:	General purpose radiative transfer program for general planetary
C	atmospheres. Combines lbl, band and c-k techniques for simulating
C	remote sensing infrared measurements of atmospheres.
C
C	Identical to Radtran EXCEPT that for lbl calculations, the line
C	binning scheme of Radtran is dropped and instead all the lines are
C	read in at the start of the program. This means that smaller
C	intervals sometimes have to be used but the program is simpler and
C	easier to understand.
C
C_ARGS: See the definitions below.
C
C_FILE:	UNIT=17		The log file <opfile>.log .
C	UNIT=1		The driver file <opfile>.drv .
C	UNIT=2		The output file <opfile>.
C
C_CALL:	gettime		Obtains the system time for use in determining
C			the total time used by the program since
C			execution.
C	prompt		Prompts the user for input.
C	file		Forces a file extension.
C	rdlbld		Read the LBL driver file.
C	rdkey		Reads in linedata .key file.
C	rdkey_corrk	Reads in linedata .key file for correlated-k
C			calculations.
C	rdkey_band	Reads in linedata .key file for band calculations.
C	read_bands	Reads in a band parameter file.
C	read_ktable	Reads an absorption coefficient look-up table.
C	rdgas		Reads in gas details for line database.
C	rdiso		Reads in isotope details for line database.
C	read_xsec	Reads in dust/aerosol cross-section files.
C	genrads		Computes multiple path spectra at "infinite" or
C			finite resolution for several gases.
C
C_HIST:	10apr90	SBC	Modified parameters in GEN_LBL call
C	1oct94	PGJI	Major modifications to combine different models.
C	11oct96	PGJI	Revised for most recent c-k and scattering models.
C***************************** VARIABLES *******************************

      IMPLICIT NONE

C NOTE: dbcom defines the linedata base variables. it is not normally 
C stored in the same directory as the rest of the code
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
      INCLUDE '../includes/dbcom.f'
      INCLUDE '../includes/pathcom.f'

      INTEGER LUN,ULOG
      PARAMETER (LUN=2,ULOG=17)

      INTEGER IRECK(10),JJ
      INTEGER ICHECK
      INTEGER I,J,IREC,IREC0,ISYS,LUN0,INORMAL,IRAY
      INTEGER NPA,NTA,NGA,NP,NT,NG,NOUT,IPTF,IMIE,IMIE1

      REAL OUTPUT(MAXOUT),ERROR(MAXOUT),MAXDV
      REAL PRESS1(MAXK),TEMP1(MAXK),PRESSA(MAXK),TEMPA(MAXK)
      REAL GABSA(MAXG),DEL_GA(MAXG),GABS(MAXG),DEL_G(MAXG)

      REAL TOT_TIME
C TOT_TIME: The total system time it took to complete program execution.
      DOUBLE PRECISION TIME,TIME1,TIME2
C TIME: Temporary variable returned by GETTIME containing the system time.
C TIME1: System time at the beginning of program execution.
C TIME2: System time at the end of program execution.

      REAL TMIN,TMAX,PMIN,PMAX
C TMIN: Temperature [Kelbin] minimum.
C TMAX: Temperature [Kelvin] maximum.
C PMIN: Pressure [atm] minimum.
C PMAX: Pressure [atm] maximum.

      CHARACTER*100 DRVFIL,KTAFIL,LOGFIL,CIAFIL,FLAFIL

C DRVFIL: Driver file.
C KTAFIL: Absorption-coefficient file.
C LOGFIL: Log file.
      CHARACTER*100 CORNAME(10)
C CORNAME: Appears to be precurrsor to the .kls file. Contains the names
C of the .kta files for use during calculation.

      INCLUDE '../includes/ciacom.f'
      INCLUDE '../includes/gascom.f'
      INCLUDE '../includes/planrad.f'

      CHARACTER*100 ANAME
      REAL DNU
      INTEGER IPARA      

      LOGICAL BIT

C     Setup HG common block (if required)
      REAL vwave
      REAL xwave(maxsec),xf(maxcon,maxsec),xg1(maxcon,maxsec),
     1 xg2(maxcon,maxsec)
      REAL twave,tico,frac,tnco
      COMMON /hgphas/ xwave,xf,xg1,xg2,tnco,twave,frac,tico
      COMMON /IMIESCAT/ IMIE1
C******************************** CODE *********************************

      xwave(1)=-1.
      jradf=-1
      jloggf=-1


C Obtain the system time for use in determining the total time used by the
C program since execution.

      CALL GETTIME(TIME)
      TIME1= TIME

      CALL PROMPT('driving file name?')
      READ(*,503)OPFILE
503   FORMAT(A)
      CALL FILE(OPFILE,DRVFIL,'drv')
      CALL FILE(OPFILE,LOGFIL,'log')
      CALL FILE(OPFILE,CIAFIL,'cia')

      OPEN(12,FILE=CIAFIL,STATUS='OLD')
       READ(12,503)ANAME
       READ(12,*) DNU
       READ(12,*) IPARA
      CLOSE(12)
      IREAD1=1
      IREAD2=1
      IF(IPARA.EQ.0)THEN
       ANAME1=ANAME
       DNU1=DNU
      ELSE
       ANAME2=ANAME
       DNU2=DNU
       IPARA2=IPARA
      ENDIF

C     Set whether we need to read in a para-H2 profile
      FLAGH2P=0
      if(IPARA.NE.0)FLAGH2P=1

      OPEN(UNIT=ULOG,FILE=LOGFIL,STATUS='UNKNOWN')

      WRITE(ULOG,503)'************************************************'
      WRITE(ULOG,503)'*                                              *'
      WRITE(ULOG,503)'*                   RADTRANS                   *'
      WRITE(ULOG,503)'*                                              *'
      WRITE(ULOG,503)'*   General purpose radiative transfer model   *'
      WRITE(ULOG,503)'*                                              *'
      WRITE(ULOG,503)'*               Planetary Group,               *'
      WRITE(ULOG,503)'*  Atmospheric, Oceanic and Planetary Physics, *'
      WRITE(ULOG,503)'*                   Oxford                     *'
      WRITE(ULOG,503)'*        P.G.J.Irwin and S.B.Calcutt           *'
      WRITE(ULOG,503)'*                                              *'
      WRITE(ULOG,503)'*           Original      11/4/95              *'
      WRITE(ULOG,503)'*           Revised       20/6/96 (and since)  *'
      WRITE(ULOG,503)'*           This version  26/4/12              *'
      WRITE(ULOG,503)'*                                              *'
      WRITE(ULOG,503)'************************************************'
      WRITE(ULOG,505)DRVFIL
505   FORMAT(1X,'Driving filename : ',A)

      OPEN(UNIT=1,FILE=DRVFIL,STATUS='OLD')

C Read in arrays from LUN= 1 (i.e. DRVFIL) such as VMIN, DELV, NPOINT,
C FWHM, WING, VREL, and KEYFIL.
      CALL RDLBLD


      CALL FILE(OPFILE,FLAFIL,'fla')

      CALL READFLAGS(FLAFIL,INORMAL,IRAY,IH2O,ICH4,IO3,INH3,
     1 IPTF,IMIE)
      IMIE1=IMIE

      MAXDV=VREL

      WRITE(ULOG,*)'Spectral Model Defined by : ' 
      WRITE(ULOG,506)VMIN,DELV,NPOINT,FWHM
506   FORMAT(1X,'Vmin = ',F7.2,', DelV = ',F7.3,', NPOINT = ',I6,
     1', FWHM = ',F7.3)
      WRITE(ULOG,507)ICONV,WING,VREL
507   FORMAT(1X,'ICONV = ',I2,', PAR1 = ',F5.1,', PAR2 = ',F5.1)
      WRITE(ULOG,*)' '

      WRITE(ULOG,15)NCALC
15    FORMAT(1X,'There are ',I3,' calculations : ')
      DO 14 I=1,NCALC
        WRITE(ULOG,*)I,ITYPE(I)
        IF(ICONV.GE.10.AND.ICONV.LT.20)THEN
          IF(.NOT.BIT(4,ITYPE(I)))THEN
            IF(.NOT.BIT(3,ITYPE(I)))THEN
              WRITE(*,*)'RADTRANS.f :: Error: Must have Curtis-'
              WRITE(*,*)'Godson paths for Band Models.'
              WRITE(*,*)' '
              WRITE(*,*)'Stopping program.'
              STOP
            ENDIF
          ENDIF
          IF(ITYPE(I).EQ.160.OR.ITYPE(I).EQ.161)THEN
            WRITE(*,*)'RADTRANS.f :: Error: Can_t have combined'
            WRITE(*,*)'cell/atmosphere paths for Band Models.'
            WRITE(*,*)' '
            WRITE(*,*)'Stopping program.'
            STOP
          ENDIF
        ENDIF
14    CONTINUE
      IF(ICONV.LT.12.OR.ICONV.EQ.20)THEN
        CALL RDKEY(LUN)
      ELSE IF(ICONV.EQ.24)THEN
        CALL RDKEY_CORRK(LUN,NGAS,IDGAS,ISOGAS,CORNAME)
      ELSE
        CALL RDKEY_BAND(LUN)
        CALL READ_BANDS
      ENDIF

      IF(ICONV.GE.20.AND.ICONV.LT.24)THEN
        CALL FILE(OPFILE,DGFILE,'gin')
        WRITE(ULOG,508)DGFILE
508     FORMAT(1X,'G-interval file : ',A)
      ENDIF

      CALL FILE(OPFILE,RADFILE,'sca')

      TMIN=TEMP(1)
      TMAX=TEMP(1)
      PMIN=PRESS(1)
      PMAX=PRESS(1)

      DO I=2,NLAYER
        TMIN= MIN(TMIN,TEMP(I))
        TMAX= MAX(TMAX,TEMP(I))
        PMIN= MIN(PMIN,PRESS(I))
        PMAX= MIN(PMAX,PRESS(I))
      ENDDO

      IF(ICONV.EQ.24.AND.NGAS.GT.0)THEN
        LUN0= 40
        KTAFIL= CORNAME(1)
        CALL FILE(KTAFIL,KTAFIL,'kta')
        CALL READ_KTABLE(KTAFIL,LUN0,NPOINT,VMIN,DELV,FWHM,
     1  IDGAS(1),ISOGAS(1),PRESS1,TEMP1,NP,NT,GABS,DEL_G,NG,IREC0)
        IRECK(1)= IREC0

	IF (PMIN.LT.PRESS1(1)) THEN
	  WRITE(*,*)'*** WARNING*** PMIN_PRF < PMIN_KTABLE'
	  WRITE(*,*)PMIN,EXP(PRESS1(1))
	ENDIF
	IF (PMAX.GT.PRESS1(NP)) THEN
	  WRITE(*,*)'*** WARNING *** PMAX_PRF > PMAX_KTABLE: '
	  WRITE(*,*)PMAX,EXP(PRESS1(NP))
	ENDIF
	IF (TMIN.LT.TEMP1(1)) THEN
	  WRITE(*,*)'*** WARNING *** TMIN_PRF < TMIN_KTABLE: '
	  WRITE(*,*)TMIN,TEMP1(1)
	ENDIF
	IF (TMAX.GT.TEMP1(NT)) THEN
	  WRITE(*,*)'*** WARNING *** TMAX_PRF > TMAX_KTABLE: '
	  WRITE(*,*)TMAX,TEMP1(NT)
	ENDIF


        IF(NGAS.GT.1)THEN
          DO 302 I=2,NGAS
            LUN0= 40 + I - 1
            KTAFIL= CORNAME(I)
            CALL FILE(KTAFIL,KTAFIL,'kta')     
            WRITE(*,*)'RADTRANS.f :: KTA file: ',KTAFIL
            CALL READ_KTABLE(KTAFIL,LUN0,NPOINT,VMIN,DELV,FWHM,
     1      IDGAS(I),ISOGAS(I),PRESSA,TEMPA,NPA,NTA,GABSA,
     2      DEL_GA,NGA,IREC0)
            IF(NPA.NE.NP)THEN
              WRITE(*,*)'RADTRANS.f :: Error: NPA <> NP'
	      WRITE(*,*)'Stopping program.'
	      WRITE(*,*)' '
	      WRITE(*,*)'NPA, NP: ',NPA,NP
	      STOP
            ENDIF
            DO J=1,NP
	      IF(PRESSA(J).NE.PRESS1(J))THEN
                WRITE(*,*)'RADTRANS.f :: Error: PRESSA(J) <> PRESS1(J)'
	        WRITE(*,*)'Stopping program.'
	        WRITE(*,*)' '
	        WRITE(*,*)'PRESSA(J), PRESS1(J) = ',PRESSA(J),PRESS1(J)
                STOP
              ENDIF
            ENDDO
            IF(NTA.NE.NT)THEN
              WRITE(*,*)'RADTRANS.f :: Error: NTA <> NT'
              WRITE(*,*)'Stopping program.'
              WRITE(*,*)' '
              WRITE(*,*)'NTA, NT = ',NTA,NT
	      STOP
            ENDIF
            DO J=1,NT
              IF(TEMPA(J).NE.TEMP1(J))THEN
                WRITE(*,*)'RADTRANS.f :: Error: TEMPA(J) <> TEMP1(J)'
	        WRITE(*,*)'Stopping program.'
	        WRITE(*,*)' '
	        WRITE(*,*)'TEMPA(J), TEMP1(J) = ',TEMPA(J),TEMP1(J)
                STOP
              ENDIF
            ENDDO
            IF(NGA.NE.NG)THEN
              WRITE(*,*)'RADTRANS.f :: Error: NGA <> NG'
              WRITE(*,*)'Stopping program.'
              WRITE(*,*)' '
              WRITE(*,*)'NGA, NG = ',NGA,NG
              STOP
            ENDIF
            DO J=1,NG
              IF(GABSA(J).NE.GABS(J))THEN
                WRITE(*,*)'RADTRANS.f :: Error: GABSA(J) <> GABS(J)'
                WRITE(*,*)'Stopping program.'
                WRITE(*,*)' '
                WRITE(*,*)'GABSA(J), GABS(J) = ',GABSA(J),GABS(J)
                STOP
              ENDIF
	      IF(DEL_GA(J).NE.DEL_G(J))THEN
                WRITE(*,*)'RADTRANS.f :: Error: DEL_GA(J) <> DEL_G(J)'
                WRITE(*,*)'Stopping program.'
                WRITE(*,*)' '
                WRITE(*,*)'DEL_GA(J), DEL_G(J) = ',DEL_GA(J),DEL_G(J)
                STOP
              ENDIF
            ENDDO
            IRECK(I)= IREC0
302       CONTINUE
        ENDIF
      ENDIF

      IF(NGAS.GT.0)THEN
        CALL RDGAS
        IF(ICONV.LT.12.OR.ICONV.EQ.20)CALL RDISO
      ENDIF

      NOUT= NPATH*NPOINT

      IF(NOUT.GT.MAXOUT)THEN
        WRITE(*,*)'RADTRANS.f :: Error: NOUT > MAXOUT -- the output'
        WRITE(*,*)'arrays are too small; need to reset and recompile.'
        WRITE(*,*)'Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)'NOUT, MAXOUT = ',NOUT,MAXOUT
        WRITE(*,*)'(NOUT=NPATH*NPOINT) NPATH, NPOINT = ',NPATH,NPOINT
        STOP
      ENDIF

      IF(NCONT.GT.0)THEN
        CALL READ_XSEC(XSCFIL)
        WRITE(ULOG,509)XSCFIL
509     FORMAT(1X,'Scattering X-section file : ',A)
      ENDIF



      WRITE(*,*)' CALLING genrads'
      CALL GENRADS(NLAYER,NPATH,AMOUNT,PP,PRESS,TEMP,DELH,
     1 NLAYIN,LAYINC,SCALE,EMTEMP,NGAS,IDGAS,ISOGAS,INORMAL,
     2 IRAY,IPTF,IPROC,VMIN,DELV,NPOINT,FWHM,ICONV,WING,VREL,MAXDV,
     3 IMOD,INLTE,NFILT,FILTER,VFILT,
     4 NCONT,CONT,NSEC,VSEC,XSEC,
     5 PRESS1,TEMP1,NP,NT,IRECK,GABS,DEL_G,NG,
     6 DOP,ERRLIM,FLAGH2P,HFP,OUTPUT,ERROR)
      WRITE(*,*)' genrads COMPLETE'

201   I= 2*ISYS()
      WRITE(*,*)'RADTRANS.f :: Writing output.'
      WRITE(*,*)'RECL = ',I
      OPEN(UNIT=2,FILE=OPFILE,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 ERR=202,RECL=I)
      IREC0=10   
      WRITE(LUN,REC=1)IREC0
      WRITE(LUN,REC=2)NPOINT
      WRITE(LUN,REC=3)VMIN
      WRITE(LUN,REC=4)DELV
      WRITE(LUN,REC=5)NPATH
      IREC=IREC0
C NOTE: NPOINT point output spectra are ouput in order of path
      DO 200 I=1,NOUT
        WRITE(LUN,REC=IREC)OUTPUT(I),ERROR(I)
        IREC= IREC + 1
200   CONTINUE

      WRITE(*,*)' RADTRANS.f :: calculation complete.'
      WRITE(ULOG,*)' RADTRANS.f :: calculation complete.'
      CALL GETTIME(TIME)
      TIME2= TIME
      TOT_TIME=SNGL(TIME2-TIME1)
      WRITE(*,244)TOT_TIME
      WRITE(ULOG,244)TOT_TIME
244   FORMAT(' Elapsed time (s) = ',F8.1)
      CLOSE(ULOG)
      CLOSE(LUN)
      STOP

202   CONTINUE
      WRITE(*,203)OPFILE
203   FORMAT(' unable to open file ',1A30)
      CALL PROMPT('new name? :')
      READ(*,204)OPFILE
204   FORMAT(1A30)

      GOTO 201

      END
C***********************************************************************
C***********************************************************************

