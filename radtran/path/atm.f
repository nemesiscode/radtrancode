
      SUBROUTINE ATM(TEXT)
C     $Id: atm.f,v 1.13 2011-06-17 14:47:04 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  ATM: reads atmospheric calculation for path.f and calculates layers
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
      INCLUDE '../includes/planrad.f'
C     note that laycom uses parameters defined in pathcom
C--------------------------------------------------------------
C     miscellaneous variables used in code
      REAL DTR,F
      PARAMETER (DTR=3.1415926/180.)
C     DTR is conversion factor for degrees to radians
      INTEGER NUSE,LOCLAY(MAXINC),USELAY(MAXINC),NCG,FSTCG,LSTCG,IPATH
      REAL SF(MAXINC),SIN2A,COSA,Z0,EMITT(MAXLAY)
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
      REAL S0,S1,STMP,XIFC(MAXCON,MAXPRO)
C--------------------------------------------------------------
C      path file format is atm followed by optional keywords. keyword are:
C         limb             followed by bottom layer to use
C         nadir            followed by angle (degrees) from nadir and bottom
C                          layer to use (normally 1)
C         (no)wf           weighting function
C         (no)netflux      Net flux calculation
C         (no)cg           curtis godson
C         (no)therm        thermal emission
C	  (no)hemisphere   Integrate emission into hemisphere
C         (no)scatter      Full scattering calculation
C         (no)nearlimb     Near-limb scattering calculation
C         (no)single	   Single scattering calculation
C         (no)absorb       calculate absorption not transmission
C         (no)binbb        use planck function at bin centre in genlbl
C         (no)broad        calculate emission outside of genlbl

C
      IF(.NOT.MODEL)THEN
C         obviously you must have read in a model before you can calculate
C         layers
        CALL WTEXT('no model defined for atmosphere calc.')
        STOP
      END IF
      
C     Update radius if it's being fitted  
      if(jradf.gt.0)radius=radius2
C
C     initialising atmosphere flags
C     note that default is a limb path from bottom layer
      WF=.FALSE.
      CG=.FALSE.
      THERM=.TRUE.
      HEMISP=.FALSE.
      NEARLIMB=.FALSE.
      SCATTER=.FALSE.
      SURFACE=.FALSE.
      SINGLE=.FALSE.
      BINBB=.TRUE.
      ABSORB=.FALSE.
      LIMB=.TRUE.
      NETFLUX=.FALSE.
      ANGLE=90.
      BOTLAY=1
C
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
       ELSE IF(TEXT(1:2).EQ.'WF')THEN
        WF=DEF
       ELSE IF(TEXT(1:7).EQ.'NETFLUX')THEN
        NETFLUX=DEF
       ELSE IF(TEXT(1:2).EQ.'CG')THEN
        CG=DEF
       ELSE IF(TEXT(1:5).EQ.'THERM')THEN
        THERM=DEF
       ELSE IF(TEXT(1:10).EQ.'HEMISPHERE')THEN
        HEMISP=DEF
       ELSE IF(TEXT(1:8).EQ.'NEARLIMB')THEN
        NEARLIMB=DEF
       ELSE IF(TEXT(1:6).EQ.'SINGLE')THEN
        SINGLE=DEF
       ELSE IF(TEXT(1:7).EQ.'SCATTER')THEN
        SCATTER=DEF
       ELSE IF(TEXT(1:5).EQ.'BROAD')THEN
        BROAD=DEF
       ELSE IF(TEXT(1:6).EQ.'ABSORB')THEN
        ABSORB=DEF
       ELSE IF(TEXT(1:5).EQ.'BINBB')THEN
        BINBB=DEF
       ELSE IF(TEXT(1:5).EQ.'NADIR')THEN
C       NONADIR doesn't make sense so DEF must be true
        IF(DEF)THEN
          LIMB=.FALSE.
          IPZEN=0
          READ(TEXT(6:),*,END=101)ANGLE,BOTLAY,IPZEN
          GOTO 102
101       READ(TEXT(6:),*)ANGLE,BOTLAY
102       IF(ANGLE.GT.90.)THEN
           ANGLE=180-ANGLE
           SURFACE = .TRUE.
          ELSE
           SURFACE = .FALSE.
          END IF
         ELSE
          CALL WTEXT('nonadir keyword no allowed')
          STOP
         END IF
       ELSE IF(TEXT(1:4).EQ.'LIMB')THEN
C       NOLIMB doesn't make sense so DEF must be true
        IF(DEF)THEN
          LIMB=.TRUE.
          ANGLE=90.
          READ(TEXT(5:),*)BOTLAY
         ELSE
          CALL WTEXT('nolimb keyword no allowed')
          STOP
          END IF
       ELSE
        CALL WTEXT('unrecognised atmospheric keyword : ')
        PRINT*,TEXT
        STOP
        END IF
      GOTO 4
C
5     CONTINUE
C     checking that keywords make sense
C      IF (THERM.AND.WF)THEN
C         CALL WTEXT('Can"t handle a thermal weighting-function calc.')
C         STOP
C       END IF
      IF(THERM.AND.ABSORB)THEN
        CALL WTEXT('can"t use absorption for thermal calcs - resetting')
        ABSORB=.FALSE.
      END IF
      IF(SINGLE.AND.SCATTER)THEN
        CALL WTEXT('can"t have both SINGLE and SCATTER - set SCATTER')
        SINGLE=.FALSE.
      ENDIF
      IF(.NOT.THERM)THEN
C        IF(BROAD)CALL WTEXT('BROAD ignored for non thermal calculation')
C        IF(BINBB)CALL WTEXT('BINBB ignored for non thermal calculation')
        BROAD=.FALSE.
        BINBB=.FALSE.
        END IF
      IF(BROAD.AND.BINBB)THEN
        CALL WTEXT('can"t use BROAD and BINBB - setting to BINBB')
        BROAD=.FALSE.
      END IF
      IF((SCATTER.OR.SINGLE).AND.THERM)THEN
        CALL WTEXT('THERM not required. Scattering includes emission')
        THERM=.FALSE.
      END IF
      IF(HEMISP.AND.(.NOT.THERM))THEN
        CALL WTEXT('HEMISP assumes THERM')
        THERM=.TRUE.
      ENDIF
      IF((SCATTER.OR.SINGLE).AND.CG)THEN
        CALL WTEXT('CG and SCATTER not allowed. Setting CG=.FALSE.')
        CG=.FALSE.
      END IF
      IF(SCATTER.OR.SINGLE)THEN
        IF (LIMB) THEN
          IF(SINGLE)THEN
           CALL WTEXT('SINGLE and LIMB not catered for.')
           STOP
          ENDIF
        ELSE IF (ANGLE.NE.0) THEN
          CALL WTEXT('atm.f : ')
          CALL WTEXT('NADIR must be zero for scattering calculation.')
          CALL WTEXT('   Setting nadir angle to zero.')
          ANGLE=0.
        ENDIF
      END IF
      IF(HEMISP)THEN
        IF (LIMB) THEN
          CALL WTEXT('Can not do HEMISP and LIMB')
          STOP
        ELSE IF (ANGLE.NE.0) THEN
          CALL WTEXT('NADIR must be zero for hemisphere calculation.')
          CALL WTEXT('   Setting nadir angle to zero.')
          ANGLE=0.
        ENDIF
      END IF
C

      IF(IPZEN.EQ.1)THEN
C      Compute zenith angle of ray at bottom of bottom layer, assuming it
C      has been defined at the 0km level
       Z0=RADIUS+BASEH(BOTLAY)
       ANGLE=(1./DTR)*ASIN(RADIUS*SIN(DTR*ANGLE)/Z0)
      ELSEIF(IPZEN.EQ.2)THEN
C      Compute zenith angle of ray at bottom of bottom layer, assuming it
C      has been defined at the top of the atmosphere
       Z0=RADIUS+BASEH(NLAY)+DELH(NLAY)
C      Calculate tangent altitude of ray at lowest point
       HTAN=Z0*SIN(ANGLE*DTR)-RADIUS
       PRINT*,'Near-limb path does not reach bottom layer'
       PRINT*,'Tangent altitude, radius is : ',HTAN,HTAN+RADIUS
       IF(HTAN.LE.BASEH(BOTLAY))THEN
C      Calculate zenith angle at bottom of lowest layer
        ANGLE=(1./DTR)*ASIN(Z0*SIN(DTR*ANGLE)/(RADIUS+BASEH(BOTLAY)))
       ELSE
C       We need to model this ray as a tangent path.
        LIMB=.TRUE.
        ANGLE=90.
C       Find number of bottom layer. Snap to layer with nearest base height
C       to computed tangent height.
        DO I=1,NLAY
         IF(BASEH(I).LT.HTAN)BOTLAY=I
        ENDDO
        IF(BOTLAY.LT.NLAY)THEN
         F=(HTAN-BASEH(BOTLAY))/(BASEH(BOTLAY+1)-BASEH(BOTLAY))
         IF(F.GT.0.5)BOTLAY=BOTLAY+1
        ENDIF
        print*,'botlay',botlay,baseh(botlay),baseh(botlay+1)
       ENDIF
      ENDIF

      SIN2A=SIN(DTR*ANGLE)**2
      COSA=COS(DTR*ANGLE)
      Z0=RADIUS+BASEH(BOTLAY)

      print*,'ATMG Check: ',LIMB,NADIR,ANGLE,SIN2A,COSA,Z0


C     now calculating atmospheric paths
C     
C     calculating which layers to use
C     note atmospheric layers are numbered from 1 to NLAY with NLAY at the top
C     so have to add FSTLAY-1 when refering to layers in the main arrays.
C     The layers in the actual calculation are numbered 1 to NUSE with NUSE
C     farthest from the observer (since this is more convenient for emission
C     integrals).

      IF(LIMB)THEN
        NUSE=2*(NLAY-BOTLAY+1)
        IF(NUSE.GT.MAXINC)THEN
         PRINT*,'Error in atm.f NUSE > MAXINC',NUSE,MAXINC
         STOP
        ENDIF
        DO 111 I=1,NUSE
         IF(I.LE.NUSE/2)THEN
           LOCLAY(I)=NLAY-I+1
         ELSE
           LOCLAY(I)=BOTLAY+(I-NUSE/2-1)
         END IF
         USELAY(I)=LOCLAY(I)+FSTLAY-1
  	 EMITT(I)=TEMP(USELAY(I))
111     CONTINUE
      ELSEIF(SURFACE)THEN
        NUSE=NLAY-BOTLAY+1
        IF(NUSE.GT.MAXINC)THEN
         PRINT*,'Error in atm.f NUSE > MAXINC',NUSE,MAXINC
         STOP
        ENDIF
        DO 1103 I=1,NUSE
         LOCLAY(I)=I
         USELAY(I)=LOCLAY(I)+FSTLAY-1
	 EMITT(I)=TEMP(USELAY(I))
1103    CONTINUE
      ELSE
        NUSE=NLAY-BOTLAY+1
        IF(NUSE.GT.MAXINC)THEN
         PRINT*,'Error in atm.f NUSE > MAXINC',NUSE,MAXINC
         STOP
        ENDIF
        DO 112 I=1,NUSE
         LOCLAY(I)=NLAY-I+1
         USELAY(I)=LOCLAY(I)+FSTLAY-1
	 EMITT(I)=TEMP(USELAY(I))
112     CONTINUE
      END IF


C     computing the scale factors for the layers
      DO 113 I=1,NUSE
         STMP = (RADIUS+BASEH(USELAY(I)))**2 - SIN2A*Z0**2
C        Sometimes, there are rounding errors here that cause the program
C        to try to take the square root of a _very_ small negative number.
C        This quietly fixes that and hopefully doesn't break anything
C        else....
         IF (STMP .LT. 0) STMP = 0
         S0 = SQRT(STMP) - Z0 * COSA
         IF(LOCLAY(I).EQ.NLAY)THEN
            S1=SQRT((RADIUS+H(NPRO))**2 - SIN2A*Z0**2) - Z0 * COSA
            SF(I)=(S1-S0)/(H(NPRO)-BASEH(USELAY(I)))
         ELSE
            S1=SQRT( (RADIUS+BASEH(1+USELAY(I)))**2 - SIN2A*Z0**2 )
     1           -Z0*COSA
            SF(I)=(S1-S0)/(BASEH(1+USELAY(I))-BASEH(USELAY(I)))
         END IF
 113  CONTINUE
C
C     first calculating any curtis-godson paths needed
C


      IF(CG)THEN
C
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
C
        NCG=1
C       if calculating a weighting function or thermal emission, need layers
C       from each point in the atmosphere
        IF(WF.OR.THERM)NCG=NUSE
        FSTCG=NLAYER+1
        LSTCG=NLAYER+NCG
        DO 310 M=1,NCG
        K=M+NUSE-NCG
C       K is the layer furthest from the observer for THIS CG path.
C       note that the last layer calculated is the longest, i.e. the furthest
C       from the observer, so if NCG=1, this is the path calculated.
C       i.e. if NCG=1 then K=NUSE, otherwise K runs from 1 to NUSE
C
C       initialising the variables for this CG layer
        NLAYER=NLAYER+1
        IF(NLAYER.GT.MAXLAY)THEN
         PRINT*,'Overflow in atm.f  NLAYER > MAXLAY'
         PRINT*,NLAYER,MAXLAY
         PRINT*,'Reduce number of layers, or increase MAXLAY'
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
28      CONTINUE
        DO 328 I=1,NCONT
        XIFC(I,NLAYER)=0.0
        CONT(I,NLAYER)=0.
328     CONTINUE
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
             END IF
           ELSE
             BASET(NLAYER)=BASET(USELAY(K))
             BASEH(NLAYER)=BASEH(USELAY(K))
             BASEP(NLAYER)=BASEP(USELAY(K))
           END IF
         ELSE
          EMITT(NLAYER)=TEMP(USELAY(K))
          BASET(NLAYER)=BASET(USELAY(K))
          BASEH(NLAYER)=BASEH(USELAY(K))
          BASEP(NLAYER)=BASEP(USELAY(K))
          END IF
        WRITE(*,205)NLAYER,K
205     FORMAT(' creating CG layer:',I3,
     1  ' by adding atmospheric layers 1 to',I3)
        DO 201 J=1,K
C       i.e. including layers from the observer (layer 1) to layer K
        L=USELAY(J)
        TOTAM(NLAYER)=TOTAM(NLAYER)+TOTAM(L)*SF(J)
        PRESS(NLAYER)=PRESS(NLAYER)+PRESS(L)*TOTAM(L)*SF(J)
        HFP(NLAYER)=HFP(NLAYER)+HFP(L)*TOTAM(L)*SF(J)
        HFC(NLAYER)=HFC(NLAYER)+HFC(L)*TOTAM(L)*SF(J)
        TEMP(NLAYER)=TEMP(NLAYER)+TEMP(L)*TOTAM(L)*SF(J)
        DELH(NLAYER)=DELH(NLAYER)+DELH(L)
        DO 202 I=1,NVMR
        AMOUNT(NLAYER,ATMGAS(I))=
     1  AMOUNT(NLAYER,ATMGAS(I))+AMOUNT(L,ATMGAS(I))*SF(J)
        PP(NLAYER,ATMGAS(I))=
     1  PP(NLAYER,ATMGAS(I))+PP(L,ATMGAS(I))*AMOUNT(L,ATMGAS(I))*SF(J)
202     CONTINUE
        DO 302 I=1,NCONT
         XIFC(I,NLAYER)=XIFC(I,NLAYER)+IFC(I,L)*TOTAM(L)*SF(J)
        CONT(I,NLAYER)=CONT(I,NLAYER)+CONT(I,L)*SF(J)
302     CONTINUE
201     CONTINUE
        PRESS(NLAYER)=PRESS(NLAYER)/TOTAM(NLAYER)
        HFP(NLAYER)=HFP(NLAYER)/TOTAM(NLAYER)
        HFC(NLAYER)=HFC(NLAYER)/TOTAM(NLAYER)
        TEMP(NLAYER)=TEMP(NLAYER)/TOTAM(NLAYER)
	DELH(NLAYER)=0.0		! Correction to stop CIACON error
        DO I=1,NCONT
         IFC(I,NLAYER)=INT(XIFC(I,NLAYER)/TOTAM(NLAYER)+0.5)
        ENDDO
        DO 203 I=1,NVMR
         IF(AMOUNT(NLAYER,ATMGAS(I)).GT.0.0)THEN
          PP(NLAYER,ATMGAS(I))=
     1    PP(NLAYER,ATMGAS(I))/AMOUNT(NLAYER,ATMGAS(I))
         ELSE
          Print*,'atm.f: WARNING. AMOUNT(',NLAYER,ATMGAS(I),') = 0.0'
         END IF
203     CONTINUE
310     CONTINUE
      END IF
C
C     now add the paths
      NPATH1=NPATH+1
      IPATH=1
C     need multiple paths if calculating a weighting function or if performing
C     a thermal integration outside genlbl
      IF(WF)IPATH=NUSE
      IF(THERM.AND.BROAD)IPATH=NUSE
      IF(NETFLUX)IPATH=NUSE
      DO 29 J=1,IPATH
C      WRITE(*,207)J
207   FORMAT(' path',I3)
      NPATH=NPATH+1
      IF(NPATH.GT.MAXPAT)THEN
       PRINT*,'Error in atm.f'
       PRINT*,'NPATH > MAXPAT',NPATH,MAXPAT
       PRINT*,'Reduce number of paths or recompile'
       STOP
      ENDIF
      ERRLIM(NPATH)=ERRDEF
C     setting the correct type for the genlbl path
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
         END IF
        END IF
      ELSE
       IF(ABSORB)IMOD(NPATH)=1
      ENDIF


      IF(SCATTER)IMOD(NPATH)=15
      IF(NETFLUX)THEN
       IF(SCATTER)THEN
        IMOD(NPATH)=24
       ELSE
        IMOD(NPATH)=21
       ENDIF
      ENDIF
      IF(SCATTER.AND.LIMB)THEN
C       Assumes int. rad. field calc.
        IMOD(NPATH)=23  
      ENDIF
      IF(NEARLIMB)IMOD(NPATH)=23
      IF(SINGLE)IMOD(NPATH)=16
      IF(CG)THEN
        IF(IMOD(NPATH).GE.2.AND.IMOD(NPATH).LE.3)
     1    IMOD(NPATH)=IMOD(NPATH)+8

        IF(IMOD(NPATH).GE.17.AND.IMOD(NPATH).LE.18)
     1    IMOD(NPATH)=IMOD(NPATH)+2

        NLAYIN(NPATH)=J+NCG-IPATH
C       note that always want to include CG layers from 1 to NLAYIN above
C       so that if IPATH=1 and NCG=1 you just include CG layers 1
C               if IPATH=1 and NCG=NUSE you include CG layers 1 to NCG
C               if IPATH=NCG=NUSE you include layers 1-1, 1-2,... 1-NCG
C       IPATH=NUSE and NCG=1 is not allowed
        DO 311 I=1,NLAYIN(NPATH)
        LAYINC(I,NPATH)=FSTCG+I-1
        SCALE(I,NPATH)=1.
        EMTEMP(I,NPATH)=EMITT(LAYINC(I,NPATH))
311     CONTINUE
        IF(WF)THEN
          NLAYIN(NPATH)=1
          LAYINC(1,NPATH)=FSTCG+J-1
          SCALE(1,NPATH)=1.
          EMTEMP(1,NPATH)=EMITT(LAYINC(1,NPATH))
        ENDIF
      ELSE
        NLAYIN(NPATH)=J+NUSE-IPATH
C       NLAYIN chosen so that if IPATH=1, use layers 1 to NUSE but
C       if IPATH=NUSE then include paths 1 to J. i.e. 1 to 1, 1 to 2...
C       up to 1 to NUSE
        DO 22 I=1,NLAYIN(NPATH)
           LAYINC(I,NPATH)=USELAY(I)
  	   EMTEMP(I,NPATH)=EMITT(I)
           SCALE(I,NPATH)=SF(I)
22      CONTINUE
      END IF

29    CONTINUE
C
C     have calculated the atmospheric paths, now outputing the calculation
      NCALC=NCALC+1
      IF(NCALC.GT.MAXCAL)THEN
       PRINT*,'Error in atm.f NCALC > MAXCAL',NCALC,MAXCAL
       STOP
      ENDIF
      IF(LIMB)THEN
       ITYPE(NCALC)=64
      ELSE
       ITYPE(NCALC)=0
      END IF
      IF(ABSORB)ITYPE(NCALC)=ITYPE(NCALC)+1
      IF(THERM)ITYPE(NCALC)=ITYPE(NCALC)+2
      IF(WF)ITYPE(NCALC)=ITYPE(NCALC)+4
      IF(CG)ITYPE(NCALC)=ITYPE(NCALC)+8
      IF(SCATTER.OR.SINGLE)ITYPE(NCALC)=256
      NINTP(NCALC)=3
      ICALD(1,NCALC)=NPATH1
      ICALD(2,NCALC)=NPATH
      ICALD(3,NCALC)=BOTLAY
      NREALP(NCALC)=2
      RCALD(1,NCALC)=ANGLE
      RCALD(2,NCALC)=HT
C
      RETURN
      END
