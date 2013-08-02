      SUBROUTINE LAYERG(TEXT)
C     $Id: layerg.f,v 1.6 2011-06-17 14:46:16 irwin Exp $
C***********************************************************************
C_TITL:	LAYERG: calculates atmospheric layers for path.f
C
C_DESC:	Calculates atmospheric layers for path.f
C
C_ARGS:	Input variable:
C	text	CHARACTER*(*)	Text string containing keywords needed for
C				this subroutine. CHARACTER*(*) declares an
C				incoming character variable whose length
C				is unknown.
C
C_FILE:	No files openned.
C
C_CALL:	remsp	Removes leading spaces from a character string.  
C	upcase	Capitolises text string.
C	verint	Performs vertical interpolation and integration of
C		atmospheric profiles.
C	verintg	Performs vertical interpolation and integration of
C		atmospheric profiles.
C
C_HIST:	26feb93	SBC	Original version.
C***********************************************************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/pathcom.f'
C ../includes/pathcom.f holds the variables used by the software when
C calculating atmospheric paths (e.g. NPATH, IMOD and LINKEY).
      INCLUDE '../includes/laycom.f'
C ../includes/laycom.f holds variables used only by the path software
C parameters are passed between routines mostly using common blocks
C because of the extensive use of large arrays. NOTE: laycom uses
C parameters defined in pathcom.
      INCLUDE '../includes/laygrad.f'
      INCLUDE '../includes/constdef.f'
C ../includes/laygrad.f holds the variables for use in gradient
C calculation.
      INCLUDE '../includes/planrad.f'

      INTEGER I,J,K,L,JJ,M,ICHECK
      INTEGER NINT,LAYINT
      PARAMETER (NINT=101)
      REAL SUM,F,F1

      REAL DTR,W(NINT),XOUT
      PARAMETER (DTR=3.1415926/180.)
C DTR: Conversion factor for degrees to radians.
C W: The weights used in the simpson's rule integrations.

      REAL PBOT,LNPBOT,PNOW,LNPNOW,TNOW,FPNOW,FCNOW
C PBOT: Pressure at the bottom of the lowest layer.
C LNPBOT: =ln(PBOT).
C PNOW: Pressure at current height.
C LNPNOW: =ln(PNOW).
C TNOW: Temperature at current height.

      REAL D(MAXPRO),DNOW,XICNOW(MAXCON),XIFC(MAXCON,MAXPAT)
      REAL TMP(500)

      REAL S,S0,S1,SMAX,SIN2A,COSA,Z0,DUDS,conttest
C S: Distance along the atmospheric path.
C S0: S at the bottom of the layer.
C S1: S at the top of the layer.
C SMAX: The distance S along the path to the top of the model atmosphere.
C SIN2A: Square of the sine of the angle from the nadir.
C COSA: Cosine of the angle from the nadir.
C Z0: Distance of the start of the path from the centre of the planet.
C DUDS: Number of molecules per cm2 per km along the path

      REAL HEIGHT,VMR1,X,DELS,stmp,XMOLWT,CALCMOLWT,XVMR(MAXGAS)

      CHARACTER*(*) TEXT

C********************************* CODE ********************************

C Setting defaults for the layer parameters defined in laycom.f
      LAYTYP = 1
      LAYINT = 1
      LAYHT = H(1)
      LAYANG = 0.0
      NLAY = 20
      
C     If radius is being retrieved then we need to use this radius when
C     calculating the layers 
      if(jradf.gt.0)radius=radius2

      print*,'radius layer=',radius

C Looking for keywords in file
2     READ(2,1,END=3)TEXT
1     FORMAT(A)
      CALL REMSP(TEXT)
      CALL UPCASE(TEXT)
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
        WRITE(*,*)' LAYERG.f :: Error: invalid layer keyword.'
        WRITE(*,*)' Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)' Invalid keyword: ',text
        STOP
      ENDIF
      GOTO 2
3     CONTINUE

      IF(NLAY.GT.MAXLAY)THEN
        WRITE(*,*)' LAYERG.f :: Error: NLAY > MAXLAY'
        WRITE(*,*)' Stopping program.'
        WRITE(*,*)' '
        WRITE(*,*)' NLAY, MAXLAY: ',NLAY,MAXLAY
        STOP
      ENDIF

      IF(NLAY.GT.NPRO)THEN
       PRINT*,'WARNING: layerg.f: NLAY > NPRO',NLAY,NPRO 
       PRINT*,'Code may possibly be unstable.'
      ENDIF

      IF(LAYHT.LT.H(1))THEN
        PRINT*,'**** WARNING: LAYERG ******'
        PRINT*,'LAYHT < H(1)',LAYHT,H(1)
        PRINT*,'***************************'
      ENDIF
      
C      print*,'nlay = ',nlay

C Simpsons rule weights
      DO 130 I=1,NINT
        W(I) = 2.0
        IF(I.EQ.2*(I/2))W(I) = 4.0
130   CONTINUE
      W(1) = 1.0
      W(NINT) = 1.0

      SIN2A = SIN(DTR*LAYANG)**2
      COSA = COS(DTR*LAYANG)
      Z0 = RADIUS + LAYHT
C            print*,'layer radius = ',radius

C Computing the bases of each layer
      CALL VERINT(H,P,NPRO,PBOT,LAYHT)
      LNPBOT = ALOG(PBOT)
      SMAX = SQRT((RADIUS + H(NPRO))**2 - (SIN2A*(Z0)**2)) - COSA*(Z0)

      IF(LAYTYP.EQ.0)THEN
C Splitting by equal pressure
        DO 101 I=1,NLAY
          PNOW = PBOT + FLOAT(I-1)*(P(NPRO)-PBOT)/FLOAT(NLAY)
          CALL VERINT(P,H,NPRO,BASEH(I),PNOW)
101     CONTINUE

      ELSE IF(LAYTYP.EQ.1)THEN
C Splitting by equal log pressure
        DO 103 I=1,NLAY
          LNPNOW = ALOG(P(NPRO))
          LNPNOW = LNPBOT + FLOAT(I-1)*(LNPNOW - LNPBOT)/FLOAT(NLAY)
          PNOW = EXP(LNPNOW)
          CALL VERINT(P,H,NPRO,BASEH(I),PNOW)
103     CONTINUE

      ELSE IF(LAYTYP.EQ.2)THEN
C Splitting by equal height
        print*,npro,h(npro),layht,nlay
        DO 104 I=1,NLAY
          BASEH(I) = LAYHT + FLOAT(I-1)*(H(NPRO) - LAYHT)/FLOAT(NLAY)
          print*,i,baseh(i)
104     CONTINUE

      ELSE IF(LAYTYP.EQ.3)THEN
C Splitting by equal distance
        DO 105 I=1,NLAY
          S = FLOAT(I-1)*SMAX/FLOAT(NLAY)
          BASEH(I) = SQRT(SIN2A*(Z0)**2 + (S + COSA*(Z0))**2) -RADIUS
105     CONTINUE
      ELSE

        WRITE(*,*)' LAYERG.f :: Error: no code for tangent layer type.'
        WRITE(*,*)' Stopping program.'
        WRITE(*,*)' '
        STOP
      ENDIF

      DO 200 I=1,NLAY-1
        DELH(I) = BASEH(I+1) - BASEH(I)
200   CONTINUE
      DELH(NLAY) = H(NPRO) - BASEH(NLAY)

C-----------------------------------------------------------------------
C
C	Computing details of each layer
C
C-----------------------------------------------------------------------
      DO 110 L=1,NLAY
        I = L + NLAYER
C Find the pressure and temperature associated with the base of this layer
        CALL VERINT(H,P,NPRO,BASEP(I),BASEH(I))
        CALL VERINT(H,T,NPRO,BASET(I),BASEH(I))
C and the distance from the begining of the path and the ratio of the
C distance of the path through to the vertical extent
        STMP = (RADIUS + BASEH(I))**2 - SIN2A*(Z0)**2
        IF(STMP.LT.0)STMP = 0
        S0 = SQRT(STMP) - (Z0)*COSA
        IF(L.LT.NLAY)THEN
          S1 = SQRT((RADIUS+BASEH(I+1))**2 - SIN2A*(Z0)**2) - (Z0)*COSA
          LAYSF(I) = (S1 - S0)/(BASEH(I+1) - BASEH(I))
        ELSE
          S1 = SMAX
          LAYSF(I) = (S1 - S0)/(H(NPRO) - BASEH(I))
        ENDIF

C Initialise the variables describing the layer
        PRESS(I) = 0.0
        TEMP(I) = 0.0
        TOTAM(I) = 0.0
        HFP(I) = 0.0
        HFC(I) = 0.0
        DO JJ=1,NPRO
          DTE(I,JJ) = 0.0
          DFP(I,JJ) = 0.0
          DFC(I,JJ) = 0.0
        ENDDO
        DO JJ=1,NPRO
          DAM(I,JJ) = 0.0
        ENDDO
        DO 121 J=1,NVMR
          AMOUNT(I,J) = 0.0
          PP(I,J) = 0.0
121     CONTINUE
        NCONT = NDUST
        FLAGH2P = JFP
        DO 122 J=1,NCONT
          CONT(J,I) = 0.0
          XIFC(J,I) = 0.0
122     CONTINUE
        DO JJ=1,NPRO
          DCO(I,JJ) = 0.0
        ENDDO

C Now find the layer details depending upon type specified LAYINT

        IF(LAYINT.EQ.0)THEN
C Just using the values at the centre of the layer
          S = 0.5*(S0 + S1)
          HEIGHT = SQRT((S + (Z0)*COSA)**2 + SIN2A*(Z0)**2) - RADIUS
          CALL VERINT(H,P,NPRO,PRESS(I),HEIGHT)
          CALL VERINTG(H,T,NPRO,TEMP(I),HEIGHT,JJ,F)
          DTE(I,JJ) = DTE(I,JJ) + (1 - F)
          DTE(I,JJ+1) = DTE(I,JJ+1) + F
          DO 883 J=1,NDUST
            DO M=1,NPRO
              D(M) = DUST(J,M)
            ENDDO
            CALL VERINTG(H,D,NPRO,CONT(J,I),HEIGHT,JJ,F)
883       CONTINUE
          DCO(I,JJ) = DCO(I,JJ) + (1.0 - F)
          DCO(I,JJ+1) = DCO(I,JJ+1) + F
          CALL VERINTG(H,FPH2I,NPRO,HFP(I),HEIGHT,JJ,F)	! para-H2
          CALL VERINTG(H,FCLOUDI,NPRO,HFC(I),HEIGHT,JJ,F1) !fcloud

          DO J=1,NDUST
           DO M=1,NPRO
            TMP(M)=1.0*ICLOUDI(J,M)
           ENDDO
           CALL VERINT(H,TMP,NPRO,XOUT,HEIGHT)
           IFC(J,I)=INT(XOUT+0.5)
          ENDDO

          DFP(I,JJ) = DFP(I,JJ) + (1.0 - F)
          DFP(I,JJ+1) = DFP(I,JJ+1) + F
          DFC(I,JJ) = DFC(I,JJ) + (1.0 - F1)
          DFC(I,JJ+1) = DFC(I,JJ+1) + F1

C         Not doing gradients for partial cloud identifiers yet
C         Pat Irwin   2/4/07

          DUDS = MODBOLTZ*PRESS(I)/TEMP(I)
          TOTAM(I) = DUDS*(S1 - S0)

          IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
          ELSE
            DO J=1,NVMR
C             print*,'J, I, VMR =',J,I,VMR(I,J)
             XVMR(J)=VMR(I,J)
            ENDDO
C            print*,'NVMR,XVMR,ID,ISO',NVMR,XVMR,ID,ISO
            XMOLWT=CALCMOLWT(NVMR,XVMR,ID,ISO)
          ENDIF
C          print*,'AMFORM,XMOLWT',AMFORM,XMOLWT

          DO 127 J=1,NCONT
            CONT(J,I) = CONT(J,I)*TOTAM(I)*XMOLWT/AVOGAD
127       CONTINUE
	  DO JJ=1,NPRO
            DCO(I,JJ) = DCO(I,JJ)*TOTAM(I)*XMOLWT/AVOGAD
          ENDDO
          DO 128 J=1,NVMR
            CALL VERINTG(H,VMR(1,J),NPRO,VMR1,HEIGHT,JJ,F)
            AMOUNT(I,J) = AMOUNT(I,J) + VMR1*TOTAM(I)
            PP(I,J) = PP(I,J) + VMR1*PRESS(I)
128       CONTINUE
          DAM(I,JJ) = DAM(I,JJ) + (1.0 - F)*TOTAM(I)
 	  DAM(I,JJ+1) = DAM(I,JJ+1) + F*TOTAM(I)

        ELSE IF(LAYINT.EQ.1)THEN
C Computing a Curtis-Godson equivalent path for a gas with constant mixing
C ratio
          DELS = (S1-S0)/FLOAT(NINT-1)
          DO 120 K=1,NINT
            S = S0 + FLOAT(K-1)*DELS
            HEIGHT = SQRT((S + (Z0)*COSA)**2 + SIN2A*(Z0)**2) - RADIUS

            CALL VERINT(H,P,NPRO,PNOW,HEIGHT)
            CALL VERINT(H,T,NPRO,TNOW,HEIGHT)
            CALL VERINTG(H,FPH2I,NPRO,FPNOW,HEIGHT,JJ,F)	! para-H2
            CALL VERINTG(H,FCLOUDI,NPRO,FCNOW,HEIGHT,JJ,F1)	! fcloud

            DO J=1,NDUST
             DO M=1,NPRO
              TMP(M)=1.0*ICLOUDI(J,M)
             ENDDO
             CALL VERINT(H,TMP,NPRO,XICNOW(J),HEIGHT)
            ENDDO 

C Calculating the number of molecules per km per cm2
            DUDS = MODBOLTZ*PNOW/TNOW
            TOTAM(I) = TOTAM(I) + DUDS*W(K)
            TEMP(I) = TEMP(I) + TNOW*DUDS*W(K)
            DTE(I,JJ) = DTE(I,JJ) + (1 - F)*DUDS*W(K)
            DTE(I,JJ+1) = DTE(I,JJ+1) + F*DUDS*W(K)

            PRESS(I) = PRESS(I) + PNOW*DUDS*W(K)
            HFP(I) = HFP(I) + FPNOW*DUDS*W(K)
            HFC(I) = HFC(I) + FCNOW*DUDS*W(K)

            DFP(I,JJ) = DFP(I,JJ) + (1 - F)*DUDS*W(K)
            DFP(I,JJ+1) = DFP(I,JJ+1) + F*DUDS*W(K)
            DFC(I,JJ) = DFC(I,JJ) + (1 - F1)*DUDS*W(K)
            DFC(I,JJ+1) = DFC(I,JJ+1) + F1*DUDS*W(K)


            IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
            ELSE
              DO J=1,NVMR
C                print*,'J, I, VMR =',J,I,VMR(I,J)
               XVMR(J)=VMR(I,J)
C               print*,'NVMR,XVMR,ID,ISO',NVMR,XVMR,ID,ISO
              ENDDO
              XMOLWT=CALCMOLWT(NVMR,XVMR,ID,ISO)
            ENDIF
C            print*,'AMFORM,XMOLWT',AMFORM,XMOLWT

            DO 124 J=1,NCONT
              XIFC(J,I)=XIFC(J,I)+XICNOW(J)*DUDS*W(K)
              DO M=1,NPRO
                D(M) = DUST(J,M)
              ENDDO
              CALL VERINT(H,D,NPRO,DNOW,HEIGHT)
            conttest=CONT(1,I)+conttest
              CONT(J,I) = CONT(J,I) + DNOW*DUDS*W(K)*XMOLWT/AVOGAD
124         CONTINUE
            DCO(I,JJ) = DCO(I,JJ) + (1 - F)*DUDS*W(K)*XMOLWT/AVOGAD
            DCO(I,JJ+1) = DCO(I,JJ+1) + F*DUDS*W(K)*XMOLWT/AVOGAD
            DO 123 J=1,NVMR
              CALL VERINT(H,VMR(1,J),NPRO,VMR1,HEIGHT)
              AMOUNT(I,J) = AMOUNT(I,J) + VMR1*DUDS*W(K)
              PP(I,J) = PP(I,J) + VMR1*PNOW*DUDS*W(K)
123         CONTINUE
            DAM(I,JJ) = DAM(I,JJ) + (1 - F)*DUDS*W(K)
            DAM(I,JJ+1) = DAM(I,JJ+1) + F*DUDS*W(K)
120       CONTINUE
C			print*,'MOLWT,AVOGAD=',XMOLWT,AVOGAD
          TEMP(I) = TEMP(I)/TOTAM(I)
          PRESS(I) = PRESS(I)/TOTAM(I)
          HFP(I) = HFP(I)/TOTAM(I)
          HFC(I) = HFP(I)/TOTAM(I)
          DO JJ=1,NPRO
            DTE(I,JJ) = DTE(I,JJ)/TOTAM(I)
            DFP(I,JJ) = DFP(I,JJ)/TOTAM(I)
            DFC(I,JJ) = DFC(I,JJ)/TOTAM(I)
          ENDDO
          DO 173 J=1,NVMR
            AMOUNT(I,J) = AMOUNT(I,J)*DELS/3.0
            PP(I,J) = PP(I,J)/TOTAM(I)
173       CONTINUE
C      print*,'CONTTEST=',CONTTEST
	  DO JJ=1,NPRO
            DAM(I,JJ) = DAM(I,JJ)*DELS/3.0
          ENDDO
          TOTAM(I) = TOTAM(I)*DELS/3.0
          DO 145 J=1,NCONT
            IFC(J,I)=INT(XIFC(J,I)/TOTAM(I)+0.5)
            CONT(J,I) = CONT(J,I)*DELS/3.0
C            print*,'CONT(J,I)2=',CONT(J,I)
145       CONTINUE
          DO JJ=1,NPRO
            DCO(I,JJ) = DCO(I,JJ)*DELS/3.0
          ENDDO

        ELSE
          WRITE(*,*)' LAYERG.f :: Error: unidentified layer-integration'
          WRITE(*,*)' type. Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' layint: ',layint
          STOP

        ENDIF

C-----------------------------------------------------------------------
C
C	Now scaling all the amounts back to a vertical path.
C
C-----------------------------------------------------------------------
        TOTAM(I) = TOTAM(I)/LAYSF(I)
        DO 129 J=1,NVMR
          AMOUNT(I,J) = AMOUNT(I,J)/LAYSF(I)
129     CONTINUE
        DO JJ=1,NPRO
          DAM(I,JJ) = DAM(I,JJ)/LAYSF(I)
        ENDDO

        DO 139 J=1,NCONT
          CONT(J,I) = CONT(J,I)/LAYSF(I)
C         print*,'CONT(J,I)3=',CONT(J,I)
139     CONTINUE
        DO JJ=1,NPRO
          DCO(I,JJ) = DCO(I,JJ)/LAYSF(I)
        ENDDO

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
      
C-----------------------------------------------------------------------
C
C	Dust (in CONT(J,I) has now been multiplied by the density [g 
C	cm^-3] and multiplied by the path length (km). CONT(J,I) has then
C	been scaled by the scaling factor to give the equivalent abundance
C	in the vertical path through the layer. Hence values are number of
C	dust particles per cm2 in vertical path for each layer.
C
C	Similarly TOTAM(I) is the number of molecules per cm2 for each
C	layer in the vertical path.
C
C-----------------------------------------------------------------------
      WRITE(*,*)' LAYERG.f :: Number of dust types: ',NCONT
      IF(NCONT.GT.0)THEN
        DO J=1,NCONT
          WRITE(*,*)' LAYERG.f :: Dust type = ',J
          SUM = 0.0
          DO L=1,NLAY
            I = L + NLAYER
            SUM = SUM + CONT(J,I)
          ENDDO
          WRITE(*,*)' LAYERG.f :: Total number/cm2 for path= ',SUM
        ENDDO
      ENDIF

      FSTLAY = NLAYER + 1
      NLAYER = NLAYER + NLAY
      LSTLAY = NLAYER
      LAYERS = .TRUE.

C      OPEN(47,FILE='layerg.dat',STATUS='unknown')
C      WRITE(47,*)NLAY,NPRO
C      DO I=1,NLAY
C       WRITE(47,*)(DTE(I,JJ),JJ=1,NPRO)
C      ENDDO
C      DO I=1,NLAY
C       WRITE(47,*)(DFP(I,JJ),JJ=1,NPRO)
C      ENDDO
C      DO I=1,NLAY
C       WRITE(47,*)(DAM(I,JJ),JJ=1,NPRO)
C      ENDDO
C      DO I=1,NLAY
C       WRITE(47,*)(DCO(I,JJ),JJ=1,NPRO)
C      ENDDO
C
C      CLOSE(47)

      RETURN

      END
