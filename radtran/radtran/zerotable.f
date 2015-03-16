      PROGRAM ZEROTABLE
C     $Id:
C---------------------------------------------------------------------------
C_TITLE:  ZEROTABLE
C
C_ARGS:
C
C_KEYS:
C
C_DESCR:  generates dummy k-table with coefficients set to zero.
C
C_FILES:  UNIT LUN: output file
C
C---------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
      INCLUDE '../includes/dbcom.f'

      INTEGER I,II,IV,J,K,LI,LJ
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

      INTEGER MDATA,QROT,NG
      PARAMETER (MDATA=20,QROT=1.5,NG=10)
C NG: Number of ordinates in k-distribution.

      REAL TOT_TIME
C TOT_TIME: The total system time it took to complete program execution.
      DOUBLE PRECISION TIME,TIME1,TIME2
C TIME: Temporary variable returned by GETTIME containing the system time.
C TIME1: System time at the beginning of program execution.
C TIME2: System time at the end of program execution.

      REAL VSTART,VEND,DELVSF
C VSTART: Beginning of spectral range wavenumber [cm-1].
C VEND: End of spectral range wavenumber [cm-1].
C DELVSF: DELV scale factor. Used to avoid non-integer DO loops.

      REAL VMAX,PMIN,PMAX,TMIN,TMAX
C VMAX: Wavenumber [cm-1] maximum (VMIN is already declared elsewhere --
C likely in some COMMON block in one of those ../include files).
C PMIN: Pressure [atm] minimum.
C PMAX: Pressure [atm] maximum.
C TMIN: Temperature [Kelvin] minimum.
C TMAX: Temperature [Kelbin] maximum.

      REAL P1,TE1,DT,DP,XP,TEMP1(20),PRESS1(20)
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

      REAL G_ORD(21),K_G(21),DEL_G(21),ERRK(21)
C G_ORD: Gauss-Legendre ordinates for calculating the k-distribution.
C K_G: Calculated k-distribution.
C DEL_G: Gauss-Legendre weights for integration.

      DATA G_ORD/0.013047, 0.067468, 0.160295, 0.283302, 0.425563,
     1           0.574437, 0.716698, 0.839705, 0.932532, 0.986953,
     2           0., 0., 0., 0., 0.,
     3           0., 0., 0., 0., 0.,
     4           0. /

      DATA DEL_G/0.033336, 0.074726, 0.109543, 0.134633, 0.147762,
     1           0.147762, 0.134633, 0.109543, 0.074726, 0.033336,
     2           0., 0., 0., 0., 0.,
     3           0., 0., 0., 0., 0.,
     4           0. /
cc      NG=20
cc      DATA G_ORD/0.003435, 0.018014, 0.043883, 0.080441, 0.126834,
cc     1           0.181973, 0.244566, 0.313147, 0.386107, 0.461736,
cc     2           0.538263, 0.613892, 0.686853, 0.755433, 0.818026,
cc     3           0.873166, 0.919558, 0.956117, 0.981986, 0.996564,
cc     4           0. /

cc      DATA DEL_G/0.008807, 0.020301, 0.031336, 0.041638, 0.050965,
cc     1           0.059097, 0.065844, 0.071048, 0.074586, 0.076377,
cc     2           0.076377, 0.074586, 0.071048, 0.065844, 0.059097,
cc     3           0.050965, 0.041638, 0.031336, 0.020301, 0.008807,
cc     4           0. /

      CHARACTER*100 KTAFIL
      CHARACTER*100 OPFILE
      REAL VMIN,FWHM,DELV
      INTEGER NPOINT

C******************************** CODE *********************************

C Obtain the system time for use in determining the total time used by the
C program since execution.
      CALL GETTIME(TIME)
      TIME1 = TIME

      WRITE(*,*)'Enter wavenumber minumum : '
      READ*,VMIN

      WRITE(*,*)'Enter FWHM, DELV and NPOINT : '
      READ*,FWHM,DELV,NPOINT
      VMAX = VMIN + (NPOINT - 1)*DELV
      WRITE(*,*)' VMIN --> VMAX by DELV: ',VMIN,VMAX,DELV

      WRITE(*,*)'Enter gas ID,ISO,IPROC : '
      READ*,IDGAS1,ISOGAS1,IPROC1

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

      WRITE(*,*)'Enter fractional abundance of absorber'
      WRITE(*,*)'0.0 will set the line width to be completely foreign'
      WRITE(*,*)'broadened. 1.0 will set the line width to be'
      WRITE(*,*)'completely self-broadened : '
      READ*,FRAC


      WRITE(*,*)'Enter output filename : '
      READ(5,23)OPFILE
23    FORMAT(A)
      CALL FILE(OPFILE,KTAFIL,'kta')
      IRECL = ISYS()
      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=IRECL)
      IREC0 = 11 + 2*NG + 2 + NP*NT + 2
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


      IREC = IREC0
 
      DO I=1,NG
       K_G(LOOP)=0.0
      ENDDO

      DO 10 IV=1,NPOINT
        DO 20 J=1,NP
          DO 30 K=1,NT
            DO 40 LOOP=1,NG
              WRITE(LUN0,REC=IREC)K_G(LOOP)
              IREC = IREC + 1
40          CONTINUE
30        CONTINUE
20      CONTINUE
10    CONTINUE


      CLOSE(LUN0)

C-----------------------------------------------------------------------
C
C	Close files and shut down the program.
C
C-----------------------------------------------------------------------

      END
************************************************************************
************************************************************************
