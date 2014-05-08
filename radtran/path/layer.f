      SUBROUTINE LAYER(TEXT)
C     $Id: layer.f,v 1.6 2011-06-17 14:46:16 irwin Exp $
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
      INCLUDE '../includes/planrad.f'
C     note that laycom uses parameters defined in pathcom
C--------------------------------------------------------------
C     miscellaneous variables used in code
      INTEGER NINT,ICHECK
      DOUBLE PRECISION SUM
      PARAMETER (NINT=101)
C     integrations are performed using simpson's rule, W holds the weights
      REAL W(NINT)
      REAL DTR,TMP(500)
      PARAMETER (DTR=3.1415926/180.)
C     DTR is conversion factor for degrees to radians
      REAL PBOT,LNPBOT,PNOW,LNPNOW,TNOW,S,S0,S1,SMAX,SIN2A,COSA,Z0,DUDS
      REAL DTDS,DUSCON,CONOP,D(MAXPRO),DNOW
C     PBOT is the pressure at the bottom of the lowest layer
C     LNPBOT is ln(PBOT)
C     PNOW is the current working pressure
C     LNPNOW = ln(PNOW)
C     TNOW = the temperature at current height
C     S is the distance along the atmospheric path
C     S0 and S1 are the values of S at the bottom and top of the layer
C     SMAX is the distance S along the path from the bottom to the top of the 
C          model atmosphere
C     SIN2A = square of the sine of the angle from the nadir
C     COSA = the cosine
C     Z0 = the distance of the start of the path from the centre of the planet
C     DUDS = DU/DS = the number of molecules per cm2 per km along the path
C     DTDS = DT/DS = the dust optical depth per km
      INTEGER I,J,K,k1,L,LAYINT
      REAL HEIGHT,VMR1,X,DELS,stmp,XICNOW(MAXCON),XIFC(MAXCON,MAXPAT)
      REAL XMOLWT,CALCMOLWT,XVMR(MAXGAS)
C--------------------------------------------------------------
C
C     setting defaults for the layer parameters defined in laycom.f


      LAYTYP=1
      LAYINT=1
      LAYHT=H(1)
      LAYANG=0.
      NLAY=20
C     If radius is being retrieved then we need to use this radius when
C     calculating the layers
      if(jradf.gt.0)radius=radius2
C     looking for keywords in file
2     READ(2,1,END=3)TEXT
1     FORMAT(A)
      CALL UPCASE(TEXT)
      CALL REMSP(TEXT)
      IF(TEXT(1:1).EQ.' ')THEN
        GOTO 3
       ELSE IF(TEXT(1:4).EQ.'NLAY')THEN
        READ(TEXT(5:),*)NLAY
       ELSE IF(TEXT(1:6).EQ.'LAYANG')THEN
        READ(TEXT(7:),*)LAYANG
       ELSE IF(TEXT(1:5).EQ.'LAYHT')THEN
        READ(TEXT(6:),*)LAYHT
       ELSE IF(TEXT(1:6).EQ.'LAYINT')THEN
        READ(TEXT(7:),*)LAYINT
       ELSE IF(TEXT(1:6).EQ.'LAYTYP')THEN
        READ(TEXT(7:),*)LAYTYP
       ELSE
        CALL WTEXT('invalid layer keyword')
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

      IF(NLAY.GT.MAXLAY)THEN
       PRINT*,'Error in layer.f NLAY > MAXLAY',NLAY,MAXLAY
       STOP
      ENDIF
      IF(NLAY.GT.NPRO)THEN
       PRINT*,'WARNING: layer.f: NLAY > NPRO',NLAY,NPRO
       PRINT*,'Code may possibly be unstable.'
      ENDIF
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

c      print*,'layer radius = ',radius
C
C     computing the bases of each layer
      CALL VERINT(H,P,NPRO,PBOT,LAYHT)
C      CALL CUBINT(H,P,NPRO,LAYHT,PBOT)
      LNPBOT=ALOG(PBOT)
      SMAX=SQRT((RADIUS+H(NPRO))**2-(SIN2A*(Z0)**2))
     1-COSA*(Z0)
      IF(LAYTYP.EQ.0)THEN
C       splitting by equal pressure
        DO 101 I=1,NLAY
        PNOW=PBOT+FLOAT(I-1)*(P(NPRO)-PBOT)/FLOAT(NLAY)
        CALL VERINT(P,H,NPRO,BASEH(I),PNOW)
C        CALL CUBINT(P,H,NPRO,PNOW,BASEH(I))
101     CONTINUE
       ELSE IF(LAYTYP.EQ.1)THEN
C       splitting by equal log pressure
        DO 103 I=1,NLAY
        LNPNOW=ALOG(P(NPRO))
        LNPNOW=LNPBOT+FLOAT(I-1)*(LNPNOW-LNPBOT)/FLOAT(NLAY)
        PNOW=EXP(LNPNOW)
        CALL VERINT(P,H,NPRO,BASEH(I),PNOW)
C        CALL CUBINT(P,H,NPRO,PNOW,BASEH(I))
103     CONTINUE
       ELSE IF(LAYTYP.EQ.2)THEN
C       splitting by equal height
        print*,'LAYER',npro,layht,h(npro),nlay
        DO 104 I=1,NLAY
        BASEH(I)=LAYHT+FLOAT(I-1)*(H(NPRO)-LAYHT)/FLOAT(NLAY)
        print*,i,baseh(i)
104     CONTINUE
       ELSE IF(LAYTYP.EQ.3)THEN
C       splitting by equal distance
        DO 105 I=1,NLAY
        S=FLOAT(I-1)*SMAX/FLOAT(NLAY)
        BASEH(I)=SQRT(SIN2A*(Z0)**2+(S+COSA*(Z0))**2)
     1  -RADIUS
105     CONTINUE
       ELSE
        WRITE(*,102)
102     FORMAT(' no code for tangent layer type')
        STOP
      END IF
      DO 200 I=1,NLAY-1
      DELH(I)=BASEH(I+1)-BASEH(I)
200   CONTINUE
      DELH(NLAY)=H(NPRO)-BASEH(NLAY)
C
C     computing details of each layer
      CONOP=0.
      DO 110 L=1,NLAY
      I=L+NLAYER
C     find the pressure and temperature associated with the base of this layer
      CALL VERINT(H,P,NPRO,BASEP(I),BASEH(I))
      CALL VERINT(H,T,NPRO,BASET(I),BASEH(I))
C     and the distance from the begining of the path and the ratio of the
C     distance of the path through to the vertical extent
      stmp = (RADIUS+BASEH(I))**2-SIN2A*(Z0)**2
      if (stmp.lt.0) stmp=0
      S0=SQRT(stmp) -(Z0)*COSA
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
      NCONT=NDUST
      FLAGH2P = JFP
      FLAGC = JFC
      DO 122 J=1,NCONT
       CONT(J,I)=0.0
       XIFC(J,I)=0.0
122   CONTINUE
      HFP(I)=0.0
      HFC(I)=0.0
  
C     now find the layer details depending upon type specified LAYINT
      IF(LAYINT.EQ.0)THEN
C       just using the values at the centre of the layer
        S=0.5*(S0+S1)
        HEIGHT=SQRT((S+(Z0)*COSA)**2+SIN2A*(Z0)**2)
     1  -RADIUS
        CALL VERINT(H,P,NPRO,PRESS(I),HEIGHT)
        CALL VERINT(H,T,NPRO,TEMP(I),HEIGHT)
C        CALL CUBINT(H,P,NPRO,HEIGHT,PRESS(I))
C        CALL CUBINT(H,T,NPRO,HEIGHT,TEMP(I))
        DO 883 J=1,NDUST
         DO M=1,NPRO
          D(M)=DUST(J,M)
         END DO
         CALL VERINT(H,D,NPRO,CONT(J,I),HEIGHT)
883     CONTINUE
        CALL VERINT(H,FPH2I,NPRO,HFP(I),HEIGHT)		! para-H2
        CALL VERINT(H,FCLOUDI,NPRO,HFC(I),HEIGHT)	! fcloud

        DO J=1,NDUST
         DO M=1,NPRO
          TMP(M)=1.0*ICLOUDI(J,M)
         ENDDO
         CALL VERINT(H,TMP,NPRO,XOUT,HEIGHT)
         IFC(J,I)=INT(XOUT+0.5)
C         print*,'direct',J,I,IFC(J,I)
        ENDDO

        DUDS=MODBOLTZ*PRESS(I)/TEMP(I)
        TOTAM(I)=DUDS*(S1-S0)

        IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
        ELSE
            DO J=1,NVMR
C             print*,'J VMR =',J,VMR(1,J)
             XVMR(J)=VMR(I,J)
            ENDDO
            XMOLWT=CALCMOLWT(NVMR,XVMR,ID,ISO)
        ENDIF

        DO 127 J=1,NCONT
         CONT(J,I)=CONT(J,I)*TOTAM(I)*XMOLWT/AVOGAD
127     CONTINUE
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
        CALL VERINT(H,FPH2I,NPRO,FPNOW,HEIGHT)	! para-H2
        CALL VERINT(H,FCLOUDI,NPRO,FCNOW,HEIGHT) ! fcloud

        DO J=1,NDUST
         DO M=1,NPRO
          TMP(M)=1.0*ICLOUDI(J,M)
         ENDDO
         CALL VERINT(H,TMP,NPRO,XICNOW(J),HEIGHT)
        ENDDO
 
C       calculating the number of molecules per km per cm2
        DUDS=MODBOLTZ*PNOW/TNOW
        TOTAM(I)=TOTAM(I)+DUDS*W(K)
        TEMP(I)=TEMP(I)+TNOW*DUDS*W(K)
C        if (i.eq.1.and.k.eq.1) then
C           write (*, *) 's = ', s, '  s0 = ', s0
C           write (*, *) 'dels = ', dels
C           write (*, *) 'cosa = ', cosa, '  sin2a = ', sin2a
C           write (*, *) 'z0 = ', z0
C        end if
        PRESS(I)=PRESS(I)+PNOW*DUDS*W(K)
        HFP(I)=HFP(I)+FPNOW*DUDS*W(K)
        HFC(I)=HFC(I)+FCNOW*DUDS*W(K)

C        print*,'amform = ',amform
        IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
        ELSE
            DO J=1,NVMR
C             print*,'J VMR =',J,VMR(1,J)
             XVMR(J)=VMR(I,J)
C             print*,J,XVMR(J),ID(J),ISO(J)
            ENDDO
            XMOLWT=CALCMOLWT(NVMR,XVMR,ID,ISO)
        ENDIF

        DO 124 J=1,NCONT
         XIFC(J,I)=XIFC(J,I)+XICNOW(J)*DUDS*W(K)
         DO M=1,NPRO
          D(M)=DUST(J,M)
         END DO
         CALL VERINT(H,D,NPRO,DNOW,HEIGHT)
         CONT(J,I)=CONT(J,I)+DNOW*(XMOLWT/AVOGAD)*DUDS*W(K)
124     CONTINUE
        DO 123 J=1,NVMR
         CALL VERINT(H,VMR(1,J),NPRO,VMR1,HEIGHT)
         AMOUNT(I,J)=AMOUNT(I,J)+VMR1*DUDS*W(K)
         PP(I,J)=PP(I,J)+VMR1*PNOW*DUDS*W(K)
123     CONTINUE
120     CONTINUE

        TEMP(I)=TEMP(I)/TOTAM(I)
        PRESS(I)=PRESS(I)/TOTAM(I)
        HFP(I)=HFP(I)/TOTAM(I)
        HFC(I)=HFC(I)/TOTAM(I)
        DO 173 J=1,NVMR
        AMOUNT(I,J)=AMOUNT(I,J)*DELS/3.
        PP(I,J)=PP(I,J)/TOTAM(I)
173     CONTINUE
        DO 145 J=1,NCONT
         IFC(J,I)=INT(XIFC(J,I)/TOTAM(I)+0.5)
         CONT(J,I)=CONT(J,I)*DELS/3
C         print*,'CG, J,I,',J,I,IFC(J,I)
145     CONTINUE
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

      DO 139 J=1,NCONT
       CONT(J,I)=CONT(J,I)/LAYSF(I)
139   CONTINUE

C     Now check that layer properties are realistic
      ICHECK=0
      IF(PRESS(I).LT.0.0)ICHECK=1
      IF(TEMP(I).LT.0.0)ICHECK=1
      IF(TOTAM(I).LT.0.0)ICHECK=1
      IF(HFC(I).LT.0.0)ICHECK=1
      IF(TOTAM(I).LT.0.0)ICHECK=1
      DO J=1,NVMR
       IF(AMOUNT(I,J).LT.0.0)ICHECK=1
      ENDDO
      DO J=1,NVMR
       IF(AMOUNT(I,J).LT.0.0)ICHECK=1
       IF(PP(I,J).LT.0.0)ICHECK=1
      ENDDO
      DO J=1,NCONT
       IF(CONT(J,I).LT.0.0)ICHECK=1
      ENDDO

      IF(ICHECK.EQ.1)THEN
       PRINT*,'************ WARNING ******************'
       PRINT*,'LAYER.F: LAYER PROPERTIES GONE NEGATIVE'
       PRINT*,'LAYER = ',I
       PRINT*,'BASEH(I),PRESS(I),TEMP(I),TOTAM(I),HFC(I)',
     1   BASEH(I),PRESS(I),TEMP(I),TOTAM(I),HFC(I)
       PRINT*,'AMOUNTS = ',(AMOUNT(I,J),J=1,NVMR) 
       PRINT*,'Partial Pressures = ',(PP(I,J),J=1,NVMR) 
       PRINT*,'Dust opacities = ',(CONT(J,I),J=1,NCONT) 
       PRINT*,' '
       PRINT*,'You should check to see if you have set'
       PRINT*,'the lowest layer height (LAYHT) in either'
       PRINT*,'the .set or .pat file to be less than the'
       PRINT*,'lowest height listed in the .prf file'
       PRINT*,' '
       PRINT*,'Aborting...'
       PRINT*,'***************************************'
       STOP
      ENDIF
110   CONTINUE
      
C      Dust (in CONT(J,I) has now been multiplied by the density (in g cm-3)
C      and multiplied by the path length (km). CONT(J,I) has then been scaled
C      by the scaling factor to give the equivalent abundance in the vertical
C      path through the layer. Hence values are number of dust particles per 
C      cm2 in vertical path for each layer.
C
C      Similarly TOTAM(I) is the number of molecules per cm2 for each layer
C      in the vertical path. 

      PRINT*,'Number of dust types = ',NCONT
      IF(NCONT.GT.0)THEN
C       PRINT*,' '
C       PRINT*,'Number of dust particles / cm2 for each layer : '
C       PRINT*,'(Vertical path)'
       DO J=1,NCONT
        PRINT*,'Dust type = ',J
        SUM=0.
        DO L=1,NLAY
         I=L+NLAYER
         SUM=SUM+CONT(J,I)
        END DO
        PRINT*,'Total number/cm2 for vertical path = ',SUM
       END DO
      END IF

      FSTLAY=NLAYER+1
      NLAYER=NLAYER+NLAY
      LSTLAY=NLAYER
      LAYERS=.TRUE.
C
      RETURN
      END
