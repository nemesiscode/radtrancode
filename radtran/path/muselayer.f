      SUBROUTINE MUSELAYER(TEXT)
C     $Id: muselayer.f,v 1.2 2002-07-24 09:30:07 irwin Exp $
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
      PARAMETER (NINT=101)
C     integrations are performed using simpson's rule, W holds the weights
      REAL W(NINT)
      REAL DTR,TGR
      PARAMETER (DTR=3.1415926/180.)
C     DTR is conversion factor for degrees to radians
      REAL PBOT,PNOW,TNOW,S,S0,S1,SMAX,SIN2A,COSA,Z0,DUDS
C     PBOT is the pressure at the bottom of the lowest layer
C     PNOW is the current working pressure
C     LNPNOW = ln(PNOW)
C     TNOW = the temperature at current height
C     S is the distance along the atmospheric path
C     S0 and S1 are the values of S at the bottom and top of the layer
C     SMAX is the distance S along the path to the top of the model atmosphere
C     SIN2A = square of the sine of the angle from the nadir
C     COSA = the cosine
C     Z0 = the distance of the start of the path from the centre of the planet
C     DUDS = DU/DS = the number of molecules per cm2 per km along the path
      INTEGER I,J,K,L,LAYINT
      REAL HEIGHT,VMR1,DELS,CDENS(10,10),CLBOT(10),CLTOP(10)
      INTEGER NCLAY(10)
      REAL A,B,C
C--------------------------------------------------------------
C
C     setting defaults for the layer parameters defined in laycom.f
      LAYTYP=2
      LAYINT=1
      LAYHT=H(1)
      LAYANG=0.
      NLAY=20
C      print*,'MUSELAYER. NLAYER = ',nlayer
C     looking for keywords in file
2     READ(2,1,END=3)TEXT
1     FORMAT(A)
      CALL UPCASE(TEXT)
      CALL REMSP(TEXT)
      IF(TEXT(1:1).EQ.' ')THEN
        GOTO 3
       ELSE IF(TEXT(1:5).EQ.'NLAYG')THEN
        READ(TEXT(7:),*)NLAYG
       ELSE IF(TEXT(1:7).EQ.'NLAYTOP')THEN
        READ(TEXT(8:),*)NLAYTOP
       ELSE IF(TEXT(1:7).EQ.'NLAYBOT')THEN
        READ(TEXT(8:),*)NLAYBOT
       ELSE IF(TEXT(1:6).EQ.'NCLOUD')THEN
        READ(TEXT(7:),*)NCLOUD
        Hold=0.0
        DO 157 J=1,NCLOUD
         READ(2,*,END=3)A,B,C	! cloud top temp, press(mb), fraction(%)
         B=B/1013.0
	 C=1.0-C*0.01		! Convert cloud fraction (%) to cloud
C				! transmittance

         CALL INTERP_T(T,H,NPRO,H1,A,Hold)
         print*,J,'cloud T,H',A,H1
	 IF(H1.LT.0)THEN
          CALL VERINT(P,H,NPRO,H1,B)
          print*,J,'cloud P,H',B,H1
         ENDIF

         CLBOT(J)=H1-0.1
	 CLTOP(J)=H1

	 NCLAY(J)=1
	 CDENS(J,1)=-LOG(C)	! Convert cloud transmittance to
C				! optical depth
157     CONTINUE
       ELSE IF(TEXT(1:6).EQ.'LAYANG')THEN
        READ(TEXT(7:),*)LAYANG
       ELSE IF(TEXT(1:5).EQ.'LAYHT')THEN
        READ(TEXT(6:),*)LAYHT
       ELSE IF(TEXT(1:6).EQ.'LAYINT')THEN
        READ(TEXT(7:),*)LAYINT
       ELSE IF(TEXT(1:11).EQ.'GROUND TEMP')THEN
        READ(TEXT(12:),*)TGR
        PRINT*,'Ground Temperature = ',TGR
       ELSE
        CALL WTEXT('invalid muselayer keyword')
        PRINT*,TEXT
        STOP
        END IF
      GOTO 2
3     CONTINUE
C      WRITE(*,42)NLAY,LAYHT,LAYANG
42    FORMAT(' computing',I4,' layers from height',F7.1,
     1' at',F7.3,' degrees from nadir')
C      WRITE(*,43)LAYTYP,LAYINT
43    FORMAT(' layer type =',I3,'   integration over layer type =',I3)
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

C
C     computing the bases of each layer
      CALL VERINT(H,P,NPRO,PBOT,LAYHT)
C      print*,PBOT
      SMAX=SQRT((RADIUS+H(NPRO))**2-(SIN2A*(Z0)**2))
     1-COSA*(Z0)
C      print*,SMAX
C     Assume LAYTYP=2. splitting by equal height
      DO I=1,MAXCON
       DO J=1,MAXLAY
        CONT(I,J)=0
       END DO
      END DO
      I=1
      BASEH(I)=LAYHT-1.0		! Ground layer
      DO 104 J=1,NLAYBOT
       I=I+1
       BASEH(I)=LAYHT+FLOAT(J-1)*(CLBOT(1)-LAYHT)/FLOAT(NLAYBOT)
104   CONTINUE
      DO 1005 K = 1,NCLOUD-1
        DO 655 J=1,NCLAY(K)
          I=I+1
          BASEH(I)=CLBOT(K) + FLOAT(J-1)*(CLTOP(K)-CLBOT(K))/NCLAY(K)
          CONT(K,I)=CDENS(K,J)
655      CONTINUE
         DO 666 J=1,NLAYG
          I=I+1
          BASEH(I)=CLTOP(K) + FLOAT(J-1)*(CLBOT(K+1)-CLTOP(K))/NLAYG
666      CONTINUE
1005    CONTINUE
        DO 667 J=1,NCLAY(NCLOUD)
          I=I+1
          BASEH(I)=CLBOT(K) + FLOAT(J-1)*(CLTOP(K)-CLBOT(K))/NCLAY(K)
          CONT(NCLOUD,I)=CDENS(NCLOUD,J)
667     CONTINUE
        DO 668 J=1,NLAYTOP
         I=I+1
         BASEH(I)=CLTOP(NCLOUD) + 
     &    FLOAT(J-1)*(H(NPRO)-CLTOP(NCLOUD))/NLAYTOP
668     CONTINUE
      NLAY=I
      DO 200 I=1,NLAY-1
        DELH(I)=BASEH(I+1)-BASEH(I)
200   CONTINUE
      DELH(NLAY)=H(NPRO)-BASEH(NLAY)

C
C     computing details of each layer
      DO 110 L=1,NLAY
      I=L+NLAYER
C       print*,L,I
C     find the pressure and temperature associated with the base of this layer
      CALL VERINT(H,P,NPRO,BASEP(I),BASEH(I))
      CALL VERINT(H,T,NPRO,BASET(I),BASEH(I))
      BASEP(1)=BASEP(2)			! Keep 'ground layer' simple
      BASET(1)=BASET(2)
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
      NCONT=NCLOUD+1

C     now find the layer details depending upon type specified LAYINT
      IF(LAYINT.EQ.0)THEN
C       just using the values at the centre of the layer
        S=0.5*(S0+S1)
        HEIGHT=SQRT((S+(Z0)*COSA)**2+SIN2A*(Z0)**2)
     1  -RADIUS
        CALL VERINT(H,P,NPRO,PRESS(I),HEIGHT)
        CALL VERINT(H,T,NPRO,TEMP(I),HEIGHT)

        IF(INFATM.AND.L.EQ.NLAY)THEN
C         fudge so that top layer represents whole of rest of  
          CALL WTEXT('INFATM not implemented yet')
          STOP
         END IF
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

C       calculating the number of molecules per km per cm2
        DUDS=MODBOLTZ*PNOW/TNOW
        TOTAM(I)=TOTAM(I)+DUDS*W(K)
        TEMP(I)=TEMP(I)+TNOW*DUDS*W(K)
        PRESS(I)=PRESS(I)+PNOW*DUDS*W(K)
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
C        DO 145 J=1,NCONT
C         CONT(J,I)=CONT(J,I)*DELS/3
C145     CONTINUE
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
C     Set up a pseudo-ground later of near-infinite optical depth.
      TEMP(1)=TGR
      CONT(NCONT,1)=1e20      



      FSTLAY=NLAYER+1
      NLAYER=NLAYER+NLAY
      LSTLAY=NLAYER
      LAYERS=.TRUE.
C
C     Ensure bottom dummy layer has positive gas pressures!
      DO 131 J=1,NVMR
       F = PP(2,J)/PRESS(2)
       PP(1,J)=F*PRESS(1)
       AMOUNT(1,J)=F*TOTAM(1)
131   CONTINUE

      RETURN
      END



      SUBROUTINE INTERP_T(T,H,NPRO,H1,A,Hold)
      INTEGER NPRO
      REAL T(NPRO),H(NPRO),H1,A,Hold

      H1=-1.0
      DO I=1,NPRO-1
       IF((T(I).LE.A.AND.T(I+1).GT.A).OR.
     &			(T(I).GT.A.AND.T(I+1).LE.A))THEN
         F=(A-T(I))/(T(I+1)-T(I))
         H1 = (1-F)*H(I)+F*H(I+1)
         IF(H1.GT.Hold)GOTO 20
       ENDIF
      ENDDO
20    Hold=H1

      RETURN

      END
