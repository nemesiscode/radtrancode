      SUBROUTINE CLOUDLAYER(TEXT)
C     $Id: cloudlayer.f,v 1.1.1.1 2000-08-17 09:26:56 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  LAYER: calculates atmospheric layers for path.f
C
C_KEYS:   RADTRAN,SUBR
C
C_DESCR:  
C
C_ARGS:   
C
C_FILES : unit 2 - the path file [.pat]
C
C_CALLS:  
C
C_BUGS:
C
C_HIST:   26feb93 SBC Original version
C
C_END:
C--------------------------------------------------------------
      IMPLICIT NONE 
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
      INCLUDE '../includes/constdef.f'
C     note that laycom uses parameters defined in pathcom
C--------------------------------------------------------------
C     miscellaneous variables used in code
      INTEGER NINT
      REAL SUM,TMP
      PARAMETER (NINT=101)
C     integrations are performed using simpson's rule, W holds the weights
      REAL W(NINT)
      REAL DTR
      PARAMETER (DTR=3.1415926/180.)
C     DTR is conversion factor for degrees to radians
      REAL PBOT,LNPBOT,PNOW,TNOW,S,S0,S1,SMAX,SIN2A,COSA,Z0,DUDS
C     PBOT is the pressure at the bottom of the lowest layer
C     LNPBOT is ln(PBOT)
C     PNOW is the current working pressure
C     TNOW = the temperature at current height
C     S is the distance along the atmospheric path
C     S0 and S1 are the values of S at the bottom and top of the layer
C     SMAX is the distance S along the path to the top of the model atmosphere
C     SIN2A = square of the sine of the angle from the nadir
C     COSA = the cosine
C     Z0 = the distance of the start of the path from the centre of the planet
C     DUDS = DU/DS = the number of molecules per cm2 per km along the path
      INTEGER I,J,K,L,LAYINT
      REAL HEIGHT,VMR1,DELS,CPBOT(4),CHBOT(4)
      REAL LPR(4),COD1(4),COD2(4),CSH1(4),CSH2(4)
      INTEGER NLAYG(3),MFINE,NFINE,NLAYTOP,NLAYBOT,I1,K1
      PARAMETER (MFINE=1000)
      REAL TFINE(MFINE),HFINE(MFINE),PFINE(MFINE),CTEMP(MFINE)
      REAL CDENS(7,MFINE),LNP,XP,DZ,SCALEH,SH
      REAL CLH(7),CSH(7),COD(7),LP(MAXPRO),CNOW(7)
             
C--------------------------------------------------------------
C
C     setting defaults for the layer parameters defined in laycom.f
      LAYTYP=2
      LAYINT=1
      LAYHT=H(1)
      LAYANG=0.
      NLAY=20
      NFINE=1000
2     READ(2,1,END=3)TEXT
1     FORMAT(A)
      CALL UPCASE(TEXT)
      CALL REMSP(TEXT)
      IF(TEXT(1:1).EQ.' ')THEN
        GOTO 3
       ELSE IF(TEXT(1:5).EQ.'NLAYG')THEN
        READ(TEXT(7:),*)(NLAYG(J),J=1,3)
       ELSE IF(TEXT(1:7).EQ.'NLAYTOP')THEN
        READ(TEXT(8:),*)NLAYTOP
       ELSE IF(TEXT(1:7).EQ.'NLAYBOT')THEN
        READ(TEXT(8:),*)NLAYBOT
       ELSE IF(TEXT(1:6).EQ.'CLOUD')THEN
        DO J=1,3
         READ(2,*,END=3)CPBOT(J),COD1(J),CSH1(J),COD2(J),CSH2(J)
        ENDDO
        READ(2,*)CPBOT(4),COD2(4),CSH2(4)  
       ELSE IF(TEXT(1:6).EQ.'LAYANG')THEN
        READ(TEXT(7:),*)LAYANG
       ELSE IF(TEXT(1:5).EQ.'LAYHT')THEN
        READ(TEXT(6:),*)LAYHT
       ELSE IF(TEXT(1:6).EQ.'LAYINT')THEN
        READ(TEXT(7:),*)LAYINT
       ELSE
        CALL WTEXT('invalid layer keyword')
        print*,text
        STOP
        END IF
      GOTO 2
3     CONTINUE
C
C     Simpsons rule weights
      DO 130 I=1,NINT
      W(I)=2.
      IF(I.EQ.2*(I/2))W(I)=4.
130   CONTINUE
      W(1)=1.
      W(NINT)=1.
C

      SIN2A=SIN(DTR*LAYANG)**2
      COSA=COS(DTR*LAYANG)
      Z0=RADIUS+LAYHT

      print*,'npro = ',npro

C
C     computing the bases of each layer
      CALL VERINT(H,P,NPRO,PBOT,LAYHT)
      LNPBOT=ALOG(PBOT)
      SMAX=SQRT((RADIUS+H(NPRO))**2-(SIN2A*(Z0)**2))
     1-COSA*(Z0)
C     Assume LAYTYP=1. splitting by equal log pressure
      print*,'CLOUDLAYER: MAXCON = ',MAXCON
      DO I=1,MAXCON
       DO J=1,MAXLAY
        CONT(I,J)=0
       END DO
      END DO
      I=0
      DO 104 J=1,NLAYBOT
       I=I+1
       LNP = LNPBOT + FLOAT(J-1)*(LOG(CPBOT(1))-LNPBOT)/FLOAT(NLAYBOT)
       XP = EXP(LNP)
       BASEP(I)=XP
       CALL VERINT(P,H,NPRO,BASEH(I),XP)
104   CONTINUE
      LPR(1)=LOG(CPBOT(1))
      LPR(2)=LOG(CPBOT(2))
      IF(CPBOT(3).GT.CPBOT(4))THEN
       LPR(3) = LOG(CPBOT(3))
       LPR(4) = LOG(CPBOT(4))
      ELSE
       LPR(3) = LOG(CPBOT(4))
       LPR(4) = LOG(CPBOT(3))
      ENDIF

      IF(CPBOT(3).EQ.CPBOT(4))NLAYG(3)=0


      DO K=1,4
       CALL VERINT(P,H,NPRO,CHBOT(K),CPBOT(K))
      ENDDO

      DO 1005 K = 1,3
        DO 655 J=1,NLAYG(K)
          I=I+1
          LNP = LPR(K) + FLOAT(J-1)*(LPR(K+1)-LPR(K))/FLOAT(NLAYG(K))
          XP = EXP(LNP)
          BASEP(I)=XP
          CALL VERINT(P,H,NPRO,BASEH(I),XP)
655      CONTINUE
1005    CONTINUE

        DO 668 J=1,NLAYTOP
          I=I+1
          LNP = LPR(4) + FLOAT(J-1)*(LOG(P(NPRO))-LPR(4))/FLOAT(NLAYTOP)
          XP = EXP(LNP)
          BASEP(I)=XP
          CALL VERINT(P,H,NPRO,BASEH(I),XP)
668     CONTINUE
      NLAY=I
      DO 200 I=1,NLAY-1
      DELH(I)=BASEH(I+1)-BASEH(I)
200   CONTINUE
      DELH(NLAY)=H(NPRO)-BASEH(NLAY)

      NCONT=7
      DO K=1,3
       CLH(2*(K-1)+1) = CHBOT(K)
       CLH(2*(K-1)+2) = CHBOT(K)
       CSH(2*(K-1)+1) = CSH1(K)
       CSH(2*(K-1)+2) = CSH2(K)
       COD(2*(K-1)+1) = COD1(K)
       COD(2*(K-1)+2) = COD2(K)
      ENDDO
      CLH(7) = CHBOT(4)
      CSH(7) = CSH2(4)
      COD(7) = COD2(4)

C     SET UP A FINE SCALE GRID
      DO I=1,NPRO
       LP(I) = LOG(P(I))
      ENDDO
      DO I=1,NFINE
       DO K=1,7
        CDENS(K,I)=0.0
       ENDDO
      ENDDO

      DO I=1,NFINE
       HFINE(I) = LAYHT + FLOAT(I-1)*(H(NPRO)-LAYHT)/(FLOAT(NFINE)-1) 
       CALL VERINT(H,LP,NPRO,XP,HFINE(I))
       PFINE(I) = EXP(XP)
       CALL VERINT(H,T,NPRO,TFINE(I),HFINE(I))
      ENDDO

      DO I=1,NFINE
       IF(I.LT.NFINE)THEN
        DZ = HFINE(I+1)-HFINE(I)
        SCALEH = -DZ/LOG(PFINE(I+1)/PFINE(I))
       ENDIF
       DO K=1,7
        IF(HFINE(I).GE.CLH(K))THEN
         IF(I.EQ.1)THEN
          CDENS(K,I)=1.0
         ELSEIF(CDENS(K,I-1).EQ.0.0)THEN
          CDENS(K,I)=1.0
         ELSE
          SH = SCALEH*CSH(K)
          CDENS(K,I)=CDENS(K,I-1)*EXP(-DZ/SH)
         ENDIF
        ENDIF
       ENDDO
      ENDDO


C
C     computing details of each layer
      DO 110 L=1,NLAY
      I=L+NLAYER
C       print*,L,I
C     find the temperature associated with the base of this layer
      CALL VERINT(H,T,NPRO,BASET(I),BASEH(I))
C     and the distance from the begining of the path and the ratio of the
C     distance of the path through to the vertical extent
      S0=SQRT((RADIUS+BASEH(I))**2-SIN2A*(Z0)**2)
     1  -(Z0)*COSA
      IF(L.LT.NLAY)THEN
        S1=SQRT((RADIUS+BASEH(I+1))**2-SIN2A*(Z0)**2)
     1  -(Z0)*COSA
        LAYSF(I)=(S1-S0)/(BASEH(I+1)-BASEH(I))
      ELSE
        S1=SMAX
        LAYSF(I)=(S1-S0)/(H(NPRO)-BASEH(I))
      END IF
C     initialise the variables describing the layer
      PRESS(I)=0.
      TEMP(I)=0.
      TOTAM(I)=0.
      DO 121 J=1,NVMR
      AMOUNT(I,J)=0.
      PP(I,J)=0.
121   CONTINUE




C     now find the layer details depending upon type specified LAYINT
      IF(LAYINT.EQ.0)THEN
C       just using the values at the centre of the layer
        S=0.5*(S0+S1)
        HEIGHT=SQRT((S+(Z0)*COSA)**2+SIN2A*(Z0)**2)
     1  -RADIUS
        CALL VERINT(H,P,NPRO,PRESS(I),HEIGHT)
        CALL VERINT(H,T,NPRO,TEMP(I),HEIGHT)
        DO K=1,NCONT
         IF(HEIGHT.GE.CLH(K)) THEN
          DO I1=1,NFINE
           CTEMP(I1) = CDENS(K,I1)
          ENDDO
          CALL VERINT(HFINE,CTEMP,NFINE,TMP,HEIGHT)
          CONT(K,I)=TMP
         ENDIF
        END DO

C        print*,I,HEIGHT,(CONT(K,I),K=1,7)

        DUDS=MODBOLTZ*PRESS(I)/TEMP(I)
        TOTAM(I)=DUDS*(S1-S0)
        DO 128 J=1,NVMR
        CALL VERINT(H,VMR(1,J),NPRO,VMR1,HEIGHT)
        AMOUNT(I,J)=AMOUNT(I,J)+VMR1*TOTAM(I)
        PP(I,J)=PP(I,J)+VMR1*PRESS(I)
128     CONTINUE


       ELSE IF(LAYINT.EQ.1)THEN
C       computing a Curtis-Godson equivalent path for a gas with constant
C       mixing ratio
        DELS=(S1-S0)/FLOAT(NINT-1)
        DO 120 K=1,NINT
        S=S0+FLOAT(K-1)*DELS
        HEIGHT=SQRT((S+(Z0)*COSA)**2+SIN2A*(Z0)**2)
     1  -RADIUS
        CALL VERINT(H,P,NPRO,PNOW,HEIGHT)
        CALL VERINT(H,T,NPRO,TNOW,HEIGHT)
        DO K1=1,NCONT
         IF(HEIGHT.GT.CLH(K1)) THEN
          DO I1=1,NFINE
           CTEMP(I1) = CDENS(K1,I1)
          ENDDO
          CALL VERINT(HFINE,CTEMP,NFINE,CNOW(K1),HEIGHT)
         ELSE
          CNOW(K1)=0.0
         ENDIF
        END DO

C       calculating the number of molecules per km per cm2
        DUDS=MODBOLTZ*PNOW/TNOW
        TOTAM(I)=TOTAM(I)+DUDS*W(K)
        TEMP(I)=TEMP(I)+TNOW*DUDS*W(K)
        PRESS(I)=PRESS(I)+PNOW*DUDS*W(K)
        DO K1=1,NCONT
         CONT(K1,I) = CONT(K1,I)+W(K)*CNOW(K1)
        ENDDO
        DO 123 J=1,NVMR
        CALL VERINT(H,VMR(1,J),NPRO,VMR1,HEIGHT)
        AMOUNT(I,J)=AMOUNT(I,J)+VMR1*DUDS*W(K)
        PP(I,J)=PP(I,J)+VMR1*PNOW*DUDS*W(K)
123     CONTINUE
120     CONTINUE


        TEMP(I)=TEMP(I)/TOTAM(I)
        PRESS(I)=PRESS(I)/TOTAM(I)
        DO 173 J=1,NVMR
        AMOUNT(I,J)=AMOUNT(I,J)*DELS/3.
        PP(I,J)=PP(I,J)/TOTAM(I)
173     CONTINUE
        TOTAM(I)=TOTAM(I)*DELS/3.
       ELSE
        CALL WTEXT('unidentified layer integration type LAYINT')
        STOP
        END IF


C     now scaling all the amounts back to a vertical path
      TOTAM(I)=TOTAM(I)/LAYSF(I)
      DO 129 J=1,NVMR
      AMOUNT(I,J)=AMOUNT(I,J)/LAYSF(I)
129   CONTINUE

110   CONTINUE
      

      print*,'Normalising cloud opacities'
      DO K=1,7
       SUM=0.0
       DO I=1,NLAY
        SUM=SUM+CONT(K,I)
       ENDDO
       print*,K,COD(K),SUM
       DO I=1,NLAY
        CONT(K,I)=CONT(K,I)*COD(K)/SUM
       ENDDO
      ENDDO


      FSTLAY=NLAYER+1
      NLAYER=NLAYER+NLAY
      LSTLAY=NLAYER
      LAYERS=.TRUE.
C
      RETURN
      END
