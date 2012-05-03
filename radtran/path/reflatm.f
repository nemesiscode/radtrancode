      SUBROUTINE REFLATM(TEXT)
C     $Id: reflatm.f,v 1.2 2003-10-31 10:58:53 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  REFLATM: reads atmospheric calculation for path.f and calculates layers
C
C_KEYS:   RADTRAN,SUBR
C
C_DESCR:  
C
C_ARGS:   
C
C_FILES : unit 1 - atmospheric model, dust and cell input files
C         unit 2 - the path file [.pat]
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:   26feb93 SBC Original version
C
C_END:
C--------------------------------------------------------------
      CHARACTER TEXT*(*)
C--------------------------------------------------------------
C     Variables to hold calculated layers and the details of each paths and
C     calculation requested.
C     pathcom holds the variables used bu the genlbl software too
C     laycom holds variables used only by the path software
C     parameters are passed between routines mostly using common blocks
C     because of the extensive use of large arrays
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/laycom.f'
C     note that laycom uses parameters defined in pathcom
C--------------------------------------------------------------
C     miscellaneous variables used in code
      REAL DTR
      PARAMETER (DTR=3.1415926/180.)
C     DTR is conversion factor for degrees to radians
      INTEGER NUSE,LOCLAY(MAXINC),USELAY(MAXINC),NCG,FSTCG,LSTCG,IPATH
      REAL SF(MAXINC),SIN2A,COSA,Z0
      REAL SIN2A1,SIN2A2,COSA1,COSA2
C     NUSE is the number of atmospheric layers (calculated by layer)
C     to be used in this calculation
C     LOCLAY are the layer numbers (relative to the first, i.e. 1 to NLAY) to
C     use in the order to be used (first is farthest from observer)
C     USELAY are the corresponding actual layer numbers in main arrays
C     NCG is the number of Curtis-Godson paths needed
C     FSTCG is the layer number of the first curtis-godson path defined
C     LSTCG the number of the last
C     IPATH is the number of paths needed for this atmospheric calculation
C     SF is the scale factor to apply to each layer for this path
C     SIN2A = square of the sine of the angle from the nadir
C     COSA = the cosine
C     Z0 = the distance of the start of the path from the centre of the planet
      LOGICAL DEF,SURFACE
      INTEGER I,J,K,L,M,NPATH1
      REAL S0,S1
C--------------------------------------------------------------
C      path file format is atm followed by optional keywords. keyword are:
C         limb             followed by bottom layer to use
C         nadir            followed by angle (degrees) from nadir and bottom
C                          layer to use (normally 1)
C         (no)wf           weighting function
C         (no)cg           curtis godson
C         (no)therm        thermal emission
C         (no)absorb       calculate absorption not transmission
C         (no)binbb        use planck function at bin centre in genlbl
C         (no)broad        calculate emission outside of genlbl
C
      print*,'Reflecting layer atmospheric path'
      IF(.NOT.MODEL)THEN
C         obviously you must have read in a model before you can calculate
C         layers
        CALL WTEXT('no model defined for atmosphere calc.')
        STOP
      END IF
C
C     initialising atmosphere flags
      CG=.FALSE.
      ABSORB=.FALSE.

C     atmospheric path keyword loop
4     READ(2,1,END=5)TEXT
1     FORMAT(A)
      CALL REMSP(TEXT)
      CALL UPCASE(TEXT)
C     checking for negated keywords
      IF(TEXT(1:2).EQ.'NO')THEN
        DEF=.FALSE.
        TEXT(1:2)='  '
        CALL REMSP(TEXT)
       ELSE
        DEF=.TRUE.
      END IF
      IF(TEXT(1:1).EQ.' ')THEN
        GOTO 5
      ELSE IF(TEXT(1:6).EQ.'ANGLES')THEN
       READ(TEXT(7:),*)ANGLE1,ANGLE2,BOTLAY
      ELSE IF(TEXT(1:2).EQ.'CG')THEN
        CG=DEF
      ELSE IF(TEXT(1:6).EQ.'ABSORB')THEN
        ABSORB=DEF
      ELSE
        CALL WTEXT('unrecognised atmospheric keyword')
        STOP
      END IF
      GOTO 4
C
5     CONTINUE
C     checking that keywords make sense
C
      HT=BASEH(BOTLAY)
      SIN2A1=SIN(DTR*ANGLE1)**2
      COSA1=COS(DTR*ANGLE1)
      SIN2A2=SIN(DTR*ANGLE2)**2
      COSA2=COS(DTR*ANGLE2)
      Z0=RADIUS+HT
C
C     now calculating atmospheric paths
C     
C     calculating which layers to use
C     note atmospheric layers are numbered from 1 to NLAY with NLAY at the top
C     so have to add FSTLAY-1 when refering to layers in the main arrays.
C     The layers in the actual calculation are numbered 1 to NUSE with NUSE
C     farthest from the observer (since this is more convenient for emission
C     integrals).
      NUSE=2*(NLAY-BOTLAY+1)
      DO 111 I=1,NUSE
        IF(I.LE.NUSE/2)THEN
          LOCLAY(I)=NLAY-I+1
         ELSE
          LOCLAY(I)=BOTLAY+(I-NUSE/2-1)
          END IF
        USELAY(I)=LOCLAY(I)+FSTLAY-1
111     CONTINUE
C
C     computing the scale factors for the layers
      DO 113 I=1,NUSE
      IF(I.LE.NUSE/2)THEN
       SIN2A=SIN2A1
       COSA=COSA1
      ELSE
       SIN2A=SIN2A2
       COSA=COSA2
      END IF
      S0=SQRT( (RADIUS+BASEH(USELAY(I)))**2 - SIN2A*Z0**2 )
     1-Z0*COSA
      IF(LOCLAY(I).EQ.NLAY)THEN
        S1=SQRT((RADIUS+H(NPRO))**2 - SIN2A*Z0**2 ) 
     1  -Z0*COSA
        SF(I)=(S1-S0)/(H(NPRO)-BASEH(USELAY(I)))
       ELSE
        S1=SQRT( (RADIUS+BASEH(1+USELAY(I)))**2 - SIN2A*Z0**2 )
     1  -Z0*COSA
        SF(I)=(S1-S0)/(BASEH(1+USELAY(I))-BASEH(USELAY(I)))
        END IF
113   CONTINUE
C

C     first calculating any curtis-godson paths needed
      IF(CG)THEN
        NCG=1
C       if calculating a weighting function or thermal emission, need layers
C       from each point in the atmosphere
        FSTCG=NLAYER+1
        LSTCG=NLAYER+NCG
        K = NUSE
C       K is the layer furthest from the observer for THIS CG path.
C       note that the last layer calculated is the longest, i.e. the furthest
C       from the observer, so if NCG=1, this is the path calculated.
C       i.e. if NCG=1 then K=NUSE, otherwise K runs from 1 to NUSE
C
C       initialising the variables for this CG layer
        NLAYER=NLAYER+1
        DELH(NLAYER)=0.
        TEMP(NLAYER)=0.
        PRESS(NLAYER)=0.
        TOTAM(NLAYER)=0.
        DO 28 I=1,NGAS
        AMOUNT(NLAYER,I)=0.
        PP(NLAYER,I)=0.
28      CONTINUE
        DO 328 I=1,NCONT
        CONT(I,NLAYER)=0.
328     CONTINUE
C     Elsewhere in the Radtran suite layers are numbered
C     from the bottom of the atmosphere up (ie bottom layer
C     is layer number 1), but here the layers are numbered
C     from the top down (ie top layer is layer number 1)
C     Therefore a warning should be printed to alert the user.
        WRITE(*,2051)
 2051   FORMAT('Note: The layers in this program are numbered from
     & the top down,'/'ie the layer closest to the spacecraft is
     & layer 1. This is'/'contrary to the convention used elsewhere in
     & the Radtran suite, where'/'the deepest layer is known as
     & layer 1')


        BASET(NLAYER)=BASET(BOTLAY)
        BASEP(NLAYER)=BASEP(BOTLAY)
        BASEH(NLAYER)=BASEH(BOTLAY)

        WRITE(*,205)NLAYER,K
205     FORMAT(' creating CG layer:',I3,
     1  ' by adding atmospheric layers 1 to',I3)
        DO 201 J=1,K
C       i.e. including layers from the observer (layer 1) to layer K
        L=USELAY(J)
        TOTAM(NLAYER)=TOTAM(NLAYER)+TOTAM(L)*SF(J)
        PRESS(NLAYER)=PRESS(NLAYER)+PRESS(L)*TOTAM(L)*SF(J)
        TEMP(NLAYER)=TEMP(NLAYER)+TEMP(L)*TOTAM(L)*SF(J)
        DELH(NLAYER)=DELH(NLAYER)+DELH(L)
        DO 202 I=1,NVMR
        AMOUNT(NLAYER,ATMGAS(I))=
     1  AMOUNT(NLAYER,ATMGAS(I))+AMOUNT(L,ATMGAS(I))*SF(J)
        PP(NLAYER,ATMGAS(I))=
     1  PP(NLAYER,ATMGAS(I))+PP(L,ATMGAS(I))*AMOUNT(L,ATMGAS(I))*SF(J)
202     CONTINUE
        DO 302 I=1,NCONT
        CONT(I,NLAYER)=CONT(I,NLAYER)+CONT(I,L)*SF(J)
302     CONTINUE
201     CONTINUE
        PRESS(NLAYER)=PRESS(NLAYER)/TOTAM(NLAYER)
        TEMP(NLAYER)=TEMP(NLAYER)/TOTAM(NLAYER)
        DO 203 I=1,NVMR
         IF(AMOUNT(NLAYER,ATMGAS(I)).GT.0)THEN
           PP(NLAYER,ATMGAS(I))=
     1       PP(NLAYER,ATMGAS(I))/AMOUNT(NLAYER,ATMGAS(I))
         ELSE
           PP(NLAYER,ATMGAS(I))=0.0
         ENDIF
203     CONTINUE
        END IF
C
C     now add the paths
      NPATH1=NPATH+1
      IPATH=1
      DO 29 J=1,IPATH
      WRITE(*,207)J
207   FORMAT(' path',I3)
      NPATH=NPATH+1
      ERRLIM(NPATH)=ERRDEF
C     setting the correct type for the radtran path

      IMOD(NPATH)=0
      IF(ABSORB)IMOD(NPATH)=1

      IF(CG)THEN
        NLAYIN(NPATH)=J+NCG-IPATH
C       note that always want to include CG layers from 1 to NLAYIN above
C       so that if IPATH=1 and NCG=1 you just include CG layers 1
C               if IPATH=1 and NCG=NUSE you include CG layers 1 to NCG
C               if IPATH=NCG=NUSE you include layers 1-1, 1-2,... 1-NCG
C       IPATH=NUSE and NCG=1 is not allowed
        DO 311 I=1,NLAYIN(NPATH)
        LAYINC(I,NPATH)=FSTCG+I-1
        SCALE(I,NPATH)=1.
        EMTEMP(I,NPATH)=BASET(LAYINC(I,NPATH))
311     CONTINUE
       ELSE
        NLAYIN(NPATH)=J+NUSE-IPATH
C       NLAYIN chosen so that if IPATH=1, use layers 1 to NUSE but
C       if IPATH=NUSE then include paths 1 to J. i.e. 1 to 1, 1 to 2...
C       up to 1 to NUSE
        DO 22 I=1,NLAYIN(NPATH)
        LAYINC(I,NPATH)=USELAY(I)
        IF(I.GT.(NUSE/2))THEN
            IF(LOCLAY(I).EQ.NLAY)THEN
              EMTEMP(I,NPATH)=T(NPRO)
             ELSE
              EMTEMP(I,NPATH)=BASET(LAYINC(I,NPATH)+1)
             END IF
        ELSE
            EMTEMP(I,NPATH)=BASET(LAYINC(I,NPATH))
        END IF
        SCALE(I,NPATH)=SF(I)
22      CONTINUE
        END IF

29    CONTINUE
C
C     have calculated the atmospheric paths, now outputing the calculation
      NCALC=NCALC+1
      ITYPE(NCALC)=32
      IF(ABSORB)ITYPE(NCALC)=ITYPE(NCALC)+1
      IF(CG)ITYPE(NCALC)=ITYPE(NCALC)+8
      NINTP(NCALC)=2
      ICALD(1,NCALC)=NPATH1
      ICALD(2,NCALC)=NPATH
      NREALP(NCALC)=3
      RCALD(1,NCALC)=ANGLE1
      RCALD(2,NCALC)=ANGLE2
      RCALD(3,NCALC)=HT
C
      RETURN
      END
