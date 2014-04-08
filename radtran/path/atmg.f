      SUBROUTINE ATMG(TEXT)
C     $Id: atmg.f,v 1.7 2007-06-28 15:08:58 irwin Exp $
C***********************************************************************
C_TITL:	ATMG
C
C_DESC:	Reads atmospheric calculation for path.f and calculates
C	layers.
C
C_ARGS:	See definitions below.
C
C_FILE:	unit 1 - atmospheric model, dust and cell input files
C	unit 2 - the path file [.pat]
C
C_CALL:	remsp	Removes leading spaces from a character string.
C	upcase	Capitolises text string.
C
C_HIST:	26feb93	SBC	Original version
C***********************************************************************

      CHARACTER TEXT*(*)

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C ../includes/pathcom.f holds the variables used by the software when
C calculating atmospheric paths (e.g. NPATH and IMOD).
      INCLUDE '../includes/laycom.f'
C ../includes/laycom.f holds variables used only by the path software
C parameters are passed between routines mostly using common blocks
C because of the extensive use of large arrays. NOTE: laycom uses
C parameters defined in pathcom.
      INCLUDE '../includes/laygrad.f'
C ../includes/laygrad.f holds the variables for use in gradient
C calculations.

C Miscellaneous variables used in code ...
      REAL DTR
      PARAMETER (DTR=3.1415926/180.)
C DTR: conversion factor for degrees to radians
      INTEGER NUSE,LOCLAY(MAXINC),USELAY(MAXINC),NCG,FSTCG,LSTCG,IPATH
      REAL SF(MAXINC),SIN2A,COSA,Z0,EMITT(MAXLAY)
C NUSE: number of atmospheric layers (calculated by layer) to be
C used in this calculation.
C LOCLAY: the layer numbers (relative to the first, i.e. 1 to NLAY) to
C use in the order to be used (first is farthest from observer).
C USELAY: actual layer numbers in main arrays.
C NCG: number of Curtis-Godson paths needed.
C FSTCG: layer number of the first curtis-godson path defined.
C LSTCG: layer number of the last curtis-godson path defined.
C IPATH: number of paths needed for this atmospheric calculation.
C SF: scale factor to apply to each layer for this path.
C SIN2A: square of the sine of the angle from the nadir.
C COSA: cosine of angle from the nadir.
C Z0: distance of the start of the path from the centre of the planet.
      INTEGER I,J,K,L,M,NPATH1
      REAL S0,S1,STMP,XIFC(MAXCON,MAXPRO)
      LOGICAL DEF,SURFACE

C Path file format is atm followed by optional keywords. keyword are:
C         limb             followed by bottom layer to use
C         nadir            followed by angle (degrees) from nadir and bottom
C                          layer to use (normally 1)
C         (no)wf           weighting function
C         (no)cg           curtis godson
C         (no)therm        thermal emission
C	  (no)hemisphere   hemispherical integration
C         (no)scatter      Full scattering calculation
C	  (no)single	   Single scattering calculation
C         (no)absorb       calculate absorption not transmission
C         (no)binbb        use planck function at bin centre in genlbl
C         (no)broad        calculate emission outside of genlbl

C********************************* CODE ********************************

C Obviously you must have read in a model before you can calculate layers
      IF(.NOT.MODEL)THEN
        WRITE(*,*)' ATMG.F :: No model defined for atmosphere calc.'
        WRITE(*,*)' ATMG.F :: Stopping program.'
        STOP
      ENDIF

C Initialising atmosphere flags: NOTE: that the default is a limb path
C from bottom layer
      WF=.FALSE.
      CG=.FALSE.
      THERM=.TRUE.
      HEMISP=.FALSE.
      SCATTER=.FALSE.
      SINGLE=.FALSE.
      BINBB=.TRUE.
      ABSORB=.FALSE.
      LIMB=.TRUE.
      ANGLE=90.
      BOTLAY=1

C Atmospheric path keyword loop
4     READ(2,1,END=5)TEXT
1     FORMAT(A)
      CALL REMSP(TEXT)
      CALL UPCASE(TEXT)
C Checking for negated keywords
      IF(TEXT(1:2).EQ.'NO')THEN
        DEF=.FALSE.
        TEXT(1:2)='  '
        CALL REMSP(TEXT)
      ELSE
        DEF=.TRUE.
      ENDIF
      IF(TEXT(1:1).EQ.' ')THEN
        GOTO 5
      ELSE IF(TEXT(1:2).EQ.'WF')THEN
        WF=DEF
      ELSE IF(TEXT(1:2).EQ.'CG')THEN
        CG=DEF
      ELSE IF(TEXT(1:5).EQ.'THERM')THEN
        THERM=DEF
      ELSE IF(TEXT(1:10).EQ.'HEMISPHERE')THEN
        HEMISP=DEF
      ELSE IF(TEXT(1:7).EQ.'SCATTER')THEN
        SCATTER=DEF
      ELSE IF(TEXT(1:6).EQ.'SINGLE')THEN
        SINGLE=DEF
      ELSE IF(TEXT(1:5).EQ.'BROAD')THEN
        BROAD=DEF
      ELSE IF(TEXT(1:6).EQ.'ABSORB')THEN
        ABSORB=DEF
      ELSE IF(TEXT(1:5).EQ.'BINBB')THEN
        BINBB=DEF
      ELSE IF(TEXT(1:5).EQ.'NADIR')THEN
C "NONADIR" doesn't make sense so DEF must be true
        IF(DEF)THEN
          LIMB=.FALSE.
          READ(TEXT(6:),*)ANGLE,BOTLAY
          IF(ANGLE.GT.90.)THEN
            ANGLE=180-ANGLE
            SURFACE=.TRUE.
          ELSE
            SURFACE=.FALSE.
          ENDIF
        ELSE
          WRITE(*,*)' ATMG.f :: "nonadir" keyword no allowed.'
          WRITE(*,*)' ATMG.f :: Stopping program.'
          STOP
        ENDIF
      ELSE IF(TEXT(1:4).EQ.'LIMB')THEN
C "NOLIMB" doesn't make sense so DEF must be true
        IF(DEF)THEN
          LIMB=.TRUE.
          ANGLE=90.
          READ(TEXT(5:),*)BOTLAY
        ELSE
          WRITE(*,*)' ATMG.f :: "nolimb" keyword no allowed.'
          WRITE(*,*)' ATMG.f :: Stopping program.'
          STOP
        ENDIF
      ELSE
        WRITE(*,*)' ATMG.f :: Unrecognised atmospheric'
        WRITE(*,*)' Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)' Unrecognised keyword: ',TEXT
        STOP
      ENDIF
      GOTO 4
5     CONTINUE

C--------------------------------------------------------------
C
C	Checking that keywords make sense
C
C--------------------------------------------------------------
      IF (THERM.AND.WF)THEN
        WRITE(*,*)' ATMG.f :: Can_t handle a thermal weighting-'
        WRITE(*,*)' function calculation. Stopping program.'
        STOP
      ENDIF

      IF (CG)THEN
        WRITE(*,*)' ATMG.f :: Can_t yet handle Curtis-Godson'
        WRITE(*,*)' paths. Stopping program.'
        STOP
      ENDIF


      IF(THERM.AND.ABSORB)THEN
        WRITE(*,*)' ATMG.f :: Can_t use absorption for thermal'
        WRITE(*,*)' calculations -- resetting.'
        ABSORB=.FALSE.
      ENDIF
      IF(.NOT.THERM)THEN
cc        IF(BROAD)WRITE(*,*)' ATMG.f :: BROAD ignored for non-'
cc        IF(BROAD)WRITE(*,*)' thermal calculation.'
cc        IF(BINBB)WRITE(*,*)' ATMG.f :: BINBB ignored for non-'
cc        IF(BINBB)WRITE(*,*)' thermal calculation.'
        BROAD=.FALSE.
        BINBB=.FALSE.
      ENDIF
      IF(BROAD.AND.BINBB)THEN
        WRITE(*,*)' ATMG.f :: Can_t use BROAD and BINBB -- setting'
        WRITE(*,*)' to BINBB.'
        BROAD=.FALSE.
      ENDIF
      IF(SINGLE.AND.SCATTER)THEN
        CALL WTEXT('can"t have both SINGLE and SCATTER - set SCATTER')
        SINGLE=.FALSE.
      ENDIF
      IF((SCATTER.OR.SINGLE).AND.THERM)THEN
        WRITE(*,*)' ATMG.f :: THERM not required as scattering'
        WRITE(*,*)' includes emission.'
        THERM=.FALSE.
      ENDIF
      IF(HEMISP.AND.(.NOT.THERM))THEN
        WRITE(*,*)' ATMG.f :: HEMISP assumes THERM.'
        WRITE(*,*)' Setting THERM=.TRUE.'
        THERM=.TRUE.
      ENDIF
      IF((SCATTER.OR.SINGLE).AND.CG)THEN
        WRITE(*,*)' ATMG.f :: CG and SCATTER not allowed.'
        WRITE(*,*)' Setting CG=.FALSE.'
        CG=.FALSE.
      ENDIF
      IF(SCATTER.OR.SINGLE)THEN
        IF(LIMB)THEN
          WRITE(*,*)' ATMG.f :: SCATTER and LIMB not catered for.'
          WRITE(*,*)' Stopping program.'
          STOP
        ELSE IF(ANGLE.NE.0)THEN
          WRITE(*,*)' ATMG.f :: NADIR must be zero for scattering'
          WRITE(*,*)' calculations. Setting nadir angle to zero.'
          ANGLE=0.
        ENDIF
      ENDIF
      IF(HEMISP)THEN 
        IF(LIMB)THEN
          WRITE(*,*)' ATMG.f :: Can_t do HEMISP and LIMB.'
          WRITE(*,*)' Stopping program.'
          STOP
        ELSE IF(ANGLE.NE.0)THEN
          WRITE(*,*)' ATMG.f :: NADIR must be zero for hemisphere'
          WRITE(*,*)' calculation. Setting nadir angle to be zero.'
          ANGLE=0.
        ENDIF
      ENDIF
      HT=BASEH(BOTLAY)
      SIN2A=SIN(DTR*ANGLE)**2
      COSA=COS(DTR*ANGLE)
      Z0=RADIUS+HT

C--------------------------------------------------------------
C
C	Now calculate the atmospheric paths
C
C	NOTE: atmospheric layers are numbered from 1 to NLAY with NLAY at
C	the top so you have to add FSTLAY-1 when refering to layers in the
C	main arrays. The layers in the actual calculation are numbered 1
C	to NUSE with NUSE farthest from the observer (since this is more
C	convenient for emission integrals).
C
C--------------------------------------------------------------
      IF(LIMB)THEN
        NUSE=2*(NLAY-BOTLAY+1)
        IF(NUSE.GT.MAXINC)THEN
          WRITE(*,*)' ATMG.f :: Error: NUSE > MAXINC'
          WRITE(*,*)' Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' NUSE, MAXINC: ',NUSE,MAXINC
          STOP
        ENDIF
        DO 111 I=1,NUSE
          IF(I.LE.NUSE/2)THEN
            LOCLAY(I)=NLAY-I+1
          ELSE
            LOCLAY(I)=BOTLAY+(I-NUSE/2-1)
          ENDIF
          USELAY(I)=LOCLAY(I)+FSTLAY-1
          EMITT(I)=TEMP(USELAY(I))
111     CONTINUE
      ELSE IF(SURFACE)THEN
        NUSE=NLAY-BOTLAY+1
        IF(NUSE.GT.MAXINC)THEN
          WRITE(*,*)' ATMG.f :: Error: NUSE > MAXINC'
          WRITE(*,*)' Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' NUSE, MAXINC: ',NUSE,MAXINC
          STOP
        ENDIF
        DO 110 I=1,NUSE
          LOCLAY(I)=I
          USELAY(I)=LOCLAY(I)+FSTLAY-1
          EMITT(I)=TEMP(USELAY(I))
110    CONTINUE
      ELSE
        NUSE=NLAY-BOTLAY+1
        IF(NUSE.GT.MAXINC)THEN
          WRITE(*,*)' ATMG.f :: Error: NUSE > MAXINC',NUSE,MAXINC
          WRITE(*,*)' Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' NUSE, MAXINC: ',NUSE,MAXINC
          STOP
        ENDIF
        DO 112 I=1,NUSE
          LOCLAY(I)=NLAY-I+1
          USELAY(I)=LOCLAY(I)+FSTLAY-1
          EMITT(I)=TEMP(USELAY(I))
112     CONTINUE
      ENDIF

C--------------------------------------------------------------
C
C	Computing the scale factors for the layers
C
C--------------------------------------------------------------
      DO 113 I=1,NUSE
        STMP = (RADIUS+BASEH(USELAY(I)))**2 - SIN2A*Z0**2
C Sometimes, there are rounding errors here that cause the program to try
C to take the square root of a _very_ small negative number. This quietly
C fixes that and hopefully doesn't break anything else ...
        IF (STMP .LT. 0) STMP = 0
        S0 = SQRT(STMP) - Z0 * COSA
        IF(LOCLAY(I).EQ.NLAY)THEN
          S1=SQRT((RADIUS+H(NPRO))**2 - SIN2A*Z0**2) - Z0 * COSA
          SF(I)=(S1-S0)/(H(NPRO)-BASEH(USELAY(I)))
        ELSE
          S1=SQRT( (RADIUS+BASEH(1+USELAY(I)))**2 - SIN2A*Z0**2 )
     1    -Z0*COSA
           SF(I)=(S1-S0)/(BASEH(1+USELAY(I))-BASEH(USELAY(I)))
         ENDIF
113   CONTINUE


C--------------------------------------------------------------
C
C	First calculate any curtis-godson paths needed
C
C--------------------------------------------------------------
      IF(CG)THEN
C Elsewhere in the Radtran suite layers are numbered from the bottom of
C the atmosphere up (ie bottom layer is layer number 1), but here the
C layers are numbered from the top down (ie top layer is layer number 1)
C Therefore a warning should be printed to alert the user.
        WRITE(*,*)' NOTE: The layers in this program are numbered from'
        WRITE(*,*)' the top down (ie the layer closest to the'
        WRITE(*,*)' spacecraft is layer 1). This is contrary to the'
        WRITE(*,*)' convention used elsewhere in the Radtran suite,'
        WRITE(*,*)' where the deepest layer is known as layer 1.'

        NCG=1
C If calculating a weighting function or thermal emission, need layers
C from each point in the atmosphere
        IF(WF.OR.THERM)NCG=NUSE
        FSTCG=NLAYER+1
        LSTCG=NLAYER+NCG
        DO 310 M=1,NCG
          K=M + NUSE - NCG
C K is the layer furthest from the observer for THIS CG path. Note that
C the last layer calculated is the longest, i.e. the furthest from the
C observer, so if NCG=1, this is the path calculated (i.e. if NCG=1 then
C K=NUSE, otherwise K runs from 1 to NUSE).

C Initialising the variables for this CG layer
          NLAYER=NLAYER + 1
          IF(NLAYER.GT.MAXLAY)THEN
            WRITE(*,*)' ATMG.f :: Error: NLAYER > MAXLAY'
            WRITE(*,*)' Stopping program: either reduce the number of'
            WRITE(*,*)' layers, or increase MAXLAY.'
            WRITE(*,*)' '
            WRITE(*,*)' NLAYER, MAXLAY: ',NLAYER,MAXLAY
            STOP
          ENDIF

          DELH(NLAYER)=0.
          TEMP(NLAYER)=0.
          PRESS(NLAYER)=0.
          HFP(NLAYER)=0.
          HFC(NLAYER)=0.
          TOTAM(NLAYER)=0.
          DO 28 I=1,NGAS
            AMOUNT(NLAYER,I)=0.
            PP(NLAYER,I)=0.
28        CONTINUE
          DO 328 I=1,NCONT
            XIFC(I,NLAYER)=0.0
            CONT(I,NLAYER)=0.
328       CONTINUE
          IF(LIMB)THEN
            EMITT(NLAYER)=TEMP(USELAY(K))
            IF(K.GT.(NUSE/2))THEN
              IF(LOCLAY(K).LT.NLAY)THEN
                BASET(NLAYER)=BASET(1+USELAY(K))  
                BASEH(NLAYER)=BASEH(1+USELAY(K))  
                BASEP(NLAYER)=BASEP(1+USELAY(K))  
              ELSE
                BASET(NLAYER)=T(NPRO)
                BASEH(NLAYER)=H(NPRO)
                BASEP(NLAYER)=P(NPRO)
              ENDIF
            ELSE
              BASET(NLAYER)=BASET(USELAY(K))
              BASEH(NLAYER)=BASEH(USELAY(K))
              BASEP(NLAYER)=BASEP(USELAY(K))
            ENDIF
          ELSE
            EMITT(NLAYER)=TEMP(USELAY(K))
            BASET(NLAYER)=BASET(USELAY(K))
            BASEH(NLAYER)=BASEH(USELAY(K))
            BASEP(NLAYER)=BASEP(USELAY(K))
          ENDIF
          WRITE(*,*)' ATMG.f :: Creating CG layer ',NLAYER
          WRITE(*,*)' by adding atmospheric layers 1 to ',K
          DO 201 J=1,K
C i.e. including layers from the observer (layer 1) to layer K
            L=USELAY(J)
            TOTAM(NLAYER)=TOTAM(NLAYER)+TOTAM(L)*SF(J)
            PRESS(NLAYER)=PRESS(NLAYER)+PRESS(L)*TOTAM(L)*SF(J)
            HFP(NLAYER)=HFP(NLAYER)+HFP(L)*TOTAM(L)*SF(J)
            HFC(NLAYER)=HFC(NLAYER)+HFC(L)*TOTAM(L)*SF(J)
            TEMP(NLAYER)=TEMP(NLAYER)+TEMP(L)*TOTAM(L)*SF(J)
            DELH(NLAYER)=DELH(NLAYER)+DELH(L)
            DO 202 I=1,NVMR
              AMOUNT(NLAYER,ATMGAS(I))=
     1        AMOUNT(NLAYER,ATMGAS(I))+AMOUNT(L,ATMGAS(I))*SF(J)
              PP(NLAYER,ATMGAS(I))= PP(NLAYER,ATMGAS(I)) + 
     1        PP(L,ATMGAS(I))*AMOUNT(L,ATMGAS(I))*SF(J)
202         CONTINUE
            DO 302 I=1,NCONT
              XIFC(I,NLAYER)=XIFC(I,NLAYER)+IFC(I,L)*TOTAM(L)*SF(J)
              CONT(I,NLAYER)=CONT(I,NLAYER)+CONT(I,L)*SF(J)
302         CONTINUE
201       CONTINUE
          PRESS(NLAYER)=PRESS(NLAYER)/TOTAM(NLAYER)
          HFP(NLAYER)=HFP(NLAYER)/TOTAM(NLAYER)
          HFC(NLAYER)=HFC(NLAYER)/TOTAM(NLAYER)
          TEMP(NLAYER)=TEMP(NLAYER)/TOTAM(NLAYER)
          DELH(NLAYER)=0.0	! Correction to stop CIACON error
          DO I=1,NCONT
           IFC(I,NLAYER)=INT(XIFC(I,NLAYER)/TOTAM(NLAYER)+0.5)
          ENDDO

          DO 203 I=1,NVMR
            IF(AMOUNT(NLAYER,ATMGAS(I)).GT.0.0)THEN
              PP(NLAYER,ATMGAS(I))=
     1        PP(NLAYER,ATMGAS(I))/AMOUNT(NLAYER,ATMGAS(I))
            ELSE
              WRITE(*,*)' ATMG.f :: Warning: ...'
              WRITE(*,*)' ATMG.f :: AMOUNT(',NLAYER,ATMGAS(I),')=0.0'
            ENDIF
203       CONTINUE
310     CONTINUE
      ENDIF

C--------------------------------------------------------------
C
C	Now add the paths
C
C--------------------------------------------------------------
      NPATH1=NPATH + 1
      IPATH=1
C Need multiple paths if calculating a weighting function or if performing
C a thermal integration outside genlbl
      IF(WF)IPATH=NUSE
      IF(THERM.AND.BROAD)IPATH=NUSE
      DO 29 J=1,IPATH
cc        WRITE(*,*)' ATMG.f :: j, ipath: ',J,IPATH
        NPATH=NPATH+1
        IF(NPATH.GT.MAXPAT)THEN
          WRITE(*,*)' ATMG.f :: Error: NPATH > MAXPAT'
          WRITE(*,*)' Stopping program: reduce the number of paths'
          WRITE(*,*)' and/or recompile.'
          WRITE(*,*)' '
          WRITE(*,*)' NPATH, MAXPAT: ',NPATH,MAXPAT
          STOP
        ENDIF
        ERRLIM(NPATH)=ERRDEF
C Setting the correct type for the genlbl path
        IMOD(NPATH)=0
        IF(THERM)THEN
          IF(.NOT.BROAD)THEN
            IF(BINBB)THEN
              IF(HEMISP)THEN
                IMOD(NPATH)=18
              ELSE
                IMOD(NPATH)=3
              ENDIF
            ELSE
              IF(HEMISP)THEN
                IMOD(NPATH)=17
              ELSE
                IMOD(NPATH)=2
              ENDIF
            ENDIF
          ENDIF
        ELSE
          IF(ABSORB)IMOD(NPATH)=1
        ENDIF
        IF(SCATTER)IMOD(NPATH)=15
        IF(SINGLE)IMOD(NPATH)=16
        IF(CG)THEN
          IF(IMOD(NPATH).GE.2.AND.IMOD(NPATH).LE.3)
     1    IMOD(NPATH)=IMOD(NPATH)+8
          NLAYIN(NPATH)=J+NCG-IPATH
C NOTE: you always want to include CG layers from 1 to NLAYIN above so
C that if IPATH=1 and NCG=1 you just include CG layers 1. If IPATH=1 and
C NCG=NUSE you include CG layers 1 to NCG. If IPATH=NCG=NUSE you include
C layers 1-1, 1-2,... 1-NCG. IPATH=NUSE and NCG=1 is not allowed.
          DO 311 I=1,NLAYIN(NPATH)
            LAYINC(I,NPATH)=FSTCG+I-1
            SCALE(I,NPATH)=1.
            EMTEMP(I,NPATH)=EMITT(LAYINC(I,NPATH))
311       CONTINUE
          IF(WF)THEN
            NLAYIN(NPATH)=1
            LAYINC(1,NPATH)=FSTCG+J-1
            SCALE(1,NPATH)=1.
            EMTEMP(1,NPATH)=EMITT(LAYINC(1,NPATH))
          ENDIF
        ELSE
          NLAYIN(NPATH)=J+NUSE-IPATH
C NLAYIN is chosen so that if IPATH=1, use layers 1 to NUSE but if
C IPATH=NUSE then include paths 1 to J. i.e. 1 to 1, 1 to 2 ... up to 1 to
C NUSE.
          DO 22 I=1,NLAYIN(NPATH)
            LAYINC(I,NPATH)=USELAY(I)
            EMTEMP(I,NPATH)=EMITT(I)
            SCALE(I,NPATH)=SF(I)
22        CONTINUE
        ENDIF
29    CONTINUE

C--------------------------------------------------------------
C
C	Having calculated the atmospheric paths, now output the
C	calculation.
C
C--------------------------------------------------------------
      NCALC=NCALC+1
      IF(NCALC.GT.MAXCAL)THEN
        WRITE(*,*)' ATMG.f :: Error: NCALC > MAXCAL',NCALC,MAXCAL
        WRITE(*,*)' Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)' NCALC, MAXCAL: ',NCALC,MAXCAL
        STOP
      ENDIF
      IF(LIMB)THEN
        ITYPE(NCALC)=64
      ELSE
        ITYPE(NCALC)=0
      ENDIF
      IF(ABSORB)ITYPE(NCALC)=ITYPE(NCALC)+1
      IF(THERM)ITYPE(NCALC)=ITYPE(NCALC)+2
      IF(WF)ITYPE(NCALC)=ITYPE(NCALC)+4
      IF(CG)ITYPE(NCALC)=ITYPE(NCALC)+8
      IF(SCATTER.OR.SINGLE)ITYPE(NCALC)=256
      NINTP(NCALC)=2
      ICALD(1,NCALC)=NPATH1
      ICALD(2,NCALC)=NPATH
      NREALP(NCALC)=2
      RCALD(1,NCALC)=ANGLE
      RCALD(2,NCALC)=HT

      RETURN

      END
