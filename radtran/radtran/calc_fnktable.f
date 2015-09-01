      PROGRAM CALC_FNKTABLE
C     $Id: calc_fnktable.f,v 1.4 2011-10-26 09:44:08 irwin Exp $
C***********************************************************************
C_TITL:	CALC_FNKTABLE.f
C
C_DESC:	Calculates the absorption coefficient look-up table for a gas of
C	the user's choice.
C
C_ARGS:	See the definitions below.
C
C_FILE:	unit=LUN0	output file "<ktafil>.kta"
C	unit=LUN1	calc_ntable_temp.dat
C
C_CALL:	FILE		Forces correct VMS style file extension for a
C			filename. i.e. assumes a <4 character extension
C			after a dot separator.
C	RDKEY		Reads in line data key file.
C	RDGAS		Reads in details of gases in data base.
C	RDISO		Reads in isotope details for line data base.
C	RDKEY_BAND	Reads in line data key file for band data.
C	READ_BAND	Routine to read in a band parameter file.
C	LBL_KCONT	Subroutine to calculate the line continuum in the
C			range VMIN - VMAX in bins of width WING.
C	LBL_FKNEW	Calculates the cumulative K-Distribution for a
C			spectral interval for a mixture of gases at a
C			number of pressures and temperatures. This is done
C			by first generating the lbl absorption coefficient
C			spectra and then analysing these according to the 
C			equation given by Lacis and Oinas (1991):
C				f(k) = (1/(V2-V1))*SUM(ABS(dV/DK))
C			f(k) is then summed to give the cumulative k
C			distribution.
C	ESUMSET
C	CALC_ESUM5	Fits and exponential sum series to an input
C			transmission curve using the Levenburg-Marquedt
C			method non-lnear least squares method described in
C			Numerical recipes.
C	REMSP		Removes leading spaces from text string.
C
C_HIST: 30.3.00	PGJI	ORIGINAL VERSION
C	4/4/06	NT added in some pre-tabluation of variables to speed
C			up the code in lbl_fknew.
C			Changed code to calculate del_g and g_ord at run time
C			instead of using data statements. these are then
C			passed to calc_fkdist via calc_fknew
C		** NB I have not changed ESUMSET to have Del_g passed to it **
C		** i left that as it was in calc_fntable **
C	26/4/12	PGJI	Updated for new radtrans/nemesis
C
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
C ../includes/bincom.f stores the line bin variables (including NLINES,
C FSTLIN, LSTLIN) and band parameters.
      INCLUDE '../includes/dbcom.f'
C ../includes/dbcom.f stores the line database variables.
      INCLUDE '../includes/lincom.f'
C ../includes/lincom.f stores the linedata variables (including MAXLIN,
C VLIN, SLIN, ALIN, ELIN, SBLIN, TDW, TDWS and that lot).
      INCLUDE '../includes/pathcom.f'
C ../includes/parcom.f stores the parameter values such as MAXLAY,

      INTEGER I,II,IV,J,K,LI,LJ,IWAVE,JV1,JV2,JV3
      INTEGER IREC,IREC0,IRECL,ISYS,ISPEC
      INTEGER IDGAS1,ISOGAS1,IPROC1
C IDGAS1: The local gas identifier.
C ISOGAS1: The local gas-isotopic identifier; if zero all isotopes of the
C gas are included.
C IPROC1: Line wing identifier.

      INTEGER NP,NT,LOOP
C NP: Number of pressures.
C NT: Number of temperatures.
      INTEGER LUN,LUN0,LUN1
      PARAMETER (LUN=2,LUN0=30,LUN1=31)
      REAL QROT
      INTEGER MDATA,NG,NGMAX,MFIL,NFIL
      PARAMETER (MDATA=20,QROT=1.5,NGMAX=21,MFIL=1000)
C NG: Number of ordinates in k-distribution.

      REAL TOT_TIME,TFIL(MFIL),VFIL(MFIL)
C TOT_TIME: The total system time it took to complete program execution.
      DOUBLE PRECISION TIME,TIME1,TIME2
C TIME: Temporary variable returned by GETTIME containing the system time.
C TIME1: System time at the beginning of program execution.
C TIME2: System time at the end of program execution.

      REAL VSTART,VEND,DELVSF,V1
C VSTART: Beginning of spectral range wavenumber [cm-1].
C VEND: End of spectral range wavenumber [cm-1].
C DELVSF: DELV scale factor. Used to avoid non-integer DO loops.

      REAL VMAX,PMIN,PMAX,TMIN,TMAX,VMIN1,VMAX1
C VMAX: Wavenumber [cm-1] maximum (VMIN is already declared elsewhere --
C likely in some COMMON block in one of those ../include files).
C PMIN: Pressure [atm] minimum.
C PMAX: Pressure [atm] maximum.
C TMIN: Temperature [Kelvin] minimum.
C TMAX: Temperature [Kelbin] maximum.

      REAL P1,TE1,DT,DP,XP,TEMP1(MAXK),PRESS1(MAXK)
C P1: Pressure [atm] at level J.
C TE1: Temperature [Kelvin] at level K.
C PRESS1: Pressure [atm].
C TEMP1: Temperature [Kelvin].

      REAL FRAC,MAXDV
C FRAC: Required fraction (0=air broadened,1=self).
C MAXDV: Line wing cut-off parameter [cm-1]: The maximum line width away
C within which to consider the contribution of the line wings.

      REAL U,XE(MDATA),YE(MDATA),SIGE(MDATA)
      REAL SDES,SUM1,ALAMDA1
      REAL KNU0,DELAD,Y0,EL,SFB,CB1,CB2

      REAL G_ORD(ngmax),K_G(ngmax),DEL_G(ngmax),ERRK(ngmax)
C G_ORD: Gauss-Legendre ordinates for calculating the k-distribution.
C K_G: Calculated k-distribution.
C DEL_G: Gauss-Legendre weights for integration.
C **** all these are now defined by zgauleg.f ***

      CHARACTER*100 KTAFIL

C******************************** CODE *********************************

C Obtain the system time for use in determining the total time used by the
C program since execution.
      CALL GETTIME(TIME)
      TIME1 = TIME
C     Set NFIL to -1 to stop calc_fkdist_wavec trying to average over a
C     filter function

      NFIL=-1

      CALL PROMPT('Enter NG : ')
      READ*,NG

c  ** calc g_ord and del_g **
      call zgauleg(g_ord,del_g,ng,ngmax)



      CALL PROMPT('Use Wavelengths(0) or Wavenumbers(1): ')
      READ*,IWAVE

      WRITE(*,*)'Enter minumum : '
      READ*,VMIN

      WRITE(*,*)'Enter FWHM, DELV and NPOINT : '
      READ*,FWHM,DELV,NPOINT
      VMAX = VMIN + (NPOINT - 1)*DELV
      WRITE(*,*)' VMIN --> VMAX by DELV: ',VMIN,VMAX,DELV

      WRITE(*,*)'Enter gas ID,ISO,IPROC : '
      READ*,IDGAS1,ISOGAS1,IPROC1

      WRITE(*,*)'Linedata (1) or banddata (2) : '
      READ*,ISPEC
      IF(ISPEC.EQ.1)THEN
        WRITE(*,*)'Enter WING,VREL : '
        READ*,WING,VREL
      ENDIF

      WRITE(*,*)'Enter number of pressure points ( <= 20 ) : '
      READ*,NP
      WRITE(*,*)'Enter log(pmin), log(pmax) : '
      READ*,PMIN,PMAX
      DP = (PMAX - PMIN)/(NP - 1)
      DO 5 J=1,NP
        XP = PMIN + (J - 1)*DP
        PRESS1(J) = EXP(XP)
        WRITE(*,*)J,PRESS1(J)
5     CONTINUE

      WRITE(*,*)'Enter number of temperature points ( <= 20 ) : '
      READ*,NT
      WRITE(*,*)'Enter Tmin, Tmax : '
      READ*,TMIN,TMAX
      DT = (TMAX - TMIN)/(NT - 1)
      DO 6 J=1,NT
        TEMP1(J) = TMIN + (J - 1)*DT
        WRITE(*,*)J,TEMP1(J)
6     CONTINUE

      WRITE(*,*)'---------- Hydrogen CIA Absorption details ----------'
      WRITE(*,*)'Enter fractional abundance of absorber'
      WRITE(*,*)'0.0 will set the line width to be completely foreign'
      WRITE(*,*)'broadened. 1.0 will set the line width to be'
      WRITE(*,*)'completely self-broadened : '
      READ*,FRAC

      WRITE(*,*)'Enter line wing cut-off (cm^-1) : '
      READ*,MAXDV

      WRITE(*,*)'Enter database name : '
      READ(5,23)OPFILE
23    FORMAT(A)
      CALL FILE(OPFILE,KEYFIL,'key')
      IF(ISPEC.EQ.1)THEN
        CALL RDKEY(LUN)
        CALL RDGAS
        CALL RDISO
      ELSE
        CALL RDKEY_BAND(LUN)
        CALL READ_BANDS
        print*,'READBAND OK'
        CALL RDGAS
      ENDIF

      WRITE(*,*)'Enter output filename : '
      READ(5,23)OPFILE
      CALL FILE(OPFILE,KTAFIL,'kta')
      IRECL = ISYS()
      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=IRECL)
      IREC0 = 11 + 2*NG + 2 + NP + NT + 2
cc      WRITE(*,*)' CALC_FNKTABLE.f :: IREC0 = ',irec0
      WRITE(LUN0,REC=1)IREC0
      WRITE(LUN0,REC=2)NPOINT
      WRITE(LUN0,REC=3)VMIN
      WRITE(LUN0,REC=4)DELV
      WRITE(LUN0,REC=5)FWHM
      WRITE(LUN0,REC=6)NP
      WRITE(LUN0,REC=7)NT
      WRITE(LUN0,REC=8)NG
      WRITE(LUN0,REC=9)IDGAS1
      WRITE(LUN0,REC=10)ISOGAS1
      IREC = 11
      DO 299 J=1,NG
        WRITE(LUN0,REC=IREC)G_ORD(J)
        IREC = IREC + 1
299   CONTINUE
      DO 399 J=1,NG
        WRITE(LUN0,REC=IREC)DEL_G(J)
        IREC = IREC + 1
399   CONTINUE
      IREC = 11 + 2*NG + 2
      DO 301 J=1,NP
        WRITE(LUN0,REC=IREC)PRESS1(J)
        IREC = IREC + 1
301   CONTINUE
      DO 302 J=1,NT
        WRITE(LUN0,REC=IREC)TEMP1(J)
        IREC = IREC + 1
302   CONTINUE

      IF(ISPEC.EQ.1)THEN
C Calculate continuum absorptions

        IF(IWAVE.EQ.0)THEN
          V1 = VMIN
          VMIN1 = 1E4/(VMAX+0.5*FWHM)
          VMAX1 = 1E4/(V1-0.5*FWHM)
        ELSE
          VMIN1=VMIN-0.5*FWHM
          VMAX1=VMAX+0.5*FWHM
        ENDIF


C        print*,vmin1,vmax1
        DO 105 J=1,NP
          DO 102 K=1,NT
            WRITE(*,*)'Calculating Continua: IP, IT = ',J,K
C            print*,press1(j),temp1(k)
            CALL LBL_KCONT(VMIN1,VMAX1,WING,VREL,PRESS1(J),TEMP1(K),
     1      IDGAS1,ISOGAS1,FRAC,IPROC1,J,K,MAXDV,IPTF,IEXO)
102       CONTINUE
105     CONTINUE
      ENDIF

      IREC = IREC0
      i = 1
      ii = 1
      print*,irec,irec0,i,ii

C Calc_ntable_temp.dat is a temporary file whose purpose is to assist
C the user to get back on track should the program crash mid-calc.
C Calc_ntable_temp.dat is discarded at the end of the calculation.
      OPEN(UNIT=LUN1,FILE='calc_ntable_temp.dat',STATUS='UNKNOWN')

C The use of a non-integer DO loop variable or expression is obsolescent
C in Fortran 90 and deleted in Fortran 95. As DELV is often LT 1.0, a
C scale factor is used below to allow IV to be declared as an integer and
C to avoid a non-integer DO loop.
      DELVSF = 1/DELV
      print*,delv,delvsf,ispec
      print*,vmin*delvsf,vmax*delvsf,delv*delvsf
      print*,int(vmin*delvsf+0.1),int(vmax*delvsf+0.1),
     1   int(delv*delvsf+0.1)
      jv1=int(vmin*delvsf+0.1)
      jv2=int(vmax*delvsf+0.1)
      jv3=int(delv*delvsf+0.1)
      DO 10 IV=jv1,jv2,jv3
        I=1+IV-INT(VMIN*DELVSF)
        print*,iv,i,ii
        IF(ISPEC.NE.1)THEN
          KNU0 = BANDPAR(I,1,1)
          DELAD = BANDPAR(I,1,2)
          Y0 = BANDPAR(I,1,3)
          EL = BANDPAR(I,1,4)
          SFB = BANDPAR(I,1,5)
          CB1 = BANDPAR(I,1,6) * 1.e-4
          CB2= BANDPAR(I,1,7)
        ENDIF

C During the complilation of some of the larger k-tables,
C calc_ntable_temp.dat can get extremely large in size. To combat this,
C the IF-THEN loop below was inserted so that the file is regularly
C discarded if it rearches a certain size (ii= 100 should equate to about
C 4 Megs in size).
        IF (ii.GE.100) THEN
          CLOSE(UNIT=LUN1,STATUS='DELETE')
          OPEN(UNIT=LUN1,FILE='calc_ntable_temp.dat',STATUS='UNKNOWN')
          ii = 0
        ENDIF

      WRITE(*,*)'CALC_FNKTABLE.f :: Current Wavenumber = ',IV/DELVSF
      WRITE(LUN1,*)'CALC_FNKTABLE.f :: Current Wavenumber = ',IV/DELVSF
        DO 20 J=1,NP
          P1 = PRESS1(J)
          DO 30 K=1,NT
            TE1 = TEMP1(K)
cc            WRITE(*,*)'Pressure, temperature: ',P1,TE1

            IF(ISPEC.EQ.1)THEN
              VSTART = (IV/DELVSF) - 0.5*FWHM
              VEND = VSTART + FWHM
              IF(IWAVE.EQ.0)THEN
               V1=VSTART
               VSTART=1E4/VEND
               VEND=1E4/V1
              ENDIF
              CALL LBL_FKNEW(IWAVE,VSTART,VEND,P1,TE1,IDGAS1,ISOGAS1,
     1        IPROC1,J,K,FRAC,MAXDV,IPTF,NPOINT)

              DELV=(VEND-VSTART)/FLOAT(NPOINT)

              CALL CALC_FKDIST_WAVEC(IWAVE,VSTART,DELV,NPOINT,
     1  NFIL,VFIL,TFIL,G_ORD,DEL_G,K_G,NGMAX,NG)

            ELSE
              print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              print*,'! esumset needs del_g setting by hand !'
              print*,'! update routine or do it manually    !'
              print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              STOP  
              CALL ESUMSET(KNU0,DELAD,Y0,EL,SFB,cb1,cb2,QROT,P1,TE1,
     1        FRAC,XE,YE,SIGE,MDATA,G_ORD,K_G,NG)
              ALAMDA1 = 2000.
              CALL CALC_ESUM5(MDATA,XE,YE,SIGE,NG,DEL_G,K_G,ERRK,
     1        ALAMDA1)
              sdes = 0.0
              DO 55 li=1,MDATA
                U = XE(LI)
                sum1 = 0.0
                DO 303 LJ=1,NG
                  sum1 = sum1 + EXP(-k_g(LJ)*U)*del_g(LJ)
303             CONTINUE
                sdes = sdes + (sum1 - ye(li))**2
55            CONTINUE
              sdes = SQRT(sdes/MDATA)
cc              WRITE(*,*)'Esum standard deviation [%] = ',sdes*100.
            ENDIF
            WRITE(LUN1,*)(K_G(LOOP),LOOP=1,NG)
            DO 40 LOOP=1,NG
              WRITE(LUN0,REC=IREC)K_G(LOOP)
              IREC = IREC + 1
40          CONTINUE
30        CONTINUE
20      CONTINUE
        i = i + 1
        ii = ii + 1
10    CONTINUE


C-----------------------------------------------------------------------
C
C	Close files and shut down the program.
C
C-----------------------------------------------------------------------
C If the code succeeds in completion, delete "calc_ntable_temp.dat" (LUN1)
C since its main usefulness is when the code crashes prior to completion.
      CLOSE(UNIT=LUN1,STATUS='DELETE')

      WRITE(*,*)' CALC_FNKTABLE.f :: calculation complete.'
      CALL GETTIME(TIME)
      TIME2 = TIME
      TOT_TIME = SNGL(TIME2 - TIME1)
      WRITE(*,244)TOT_TIME
244   FORMAT(/' Elapsed time (s) = ',F8.1)

      END
************************************************************************
************************************************************************
