      SUBROUTINE ACKERMANMARLEYX1(IPLANET,LATITUDE,AMFORM,NPROIN,NVMR,
     1 IDGAS,ISOGAS,PIN,TIN,HIN,VMRIN,XMOLIN,NCONT,CONTIN,FLUX,IMODEL,
     2 FRAIN,JVMR,XDEEP,DENSCOND,RADCOND,MWCOND,X1,X2,QCOUT)
C     **************************************************************   
C     Subroutine to condense clouds as per the formalism of 
C     Ackerman and Marley (2001). 

C     Input variables
C	IPLANET	INTEGER		Planet number (needed to compute gravity)
C	LATITUDE  REAL		Latitude (needed to compute gravity)
C	AMFORM	INTEGER		Profile type
C	NPROIN	INTEGER		Number of vertical levels
C	NVMR	INTEGER		Number of gas vmr profiles
C	IDGAS(NVMR) INTEGER	Gas IDs
C	ISOGAS(NVMR) INTEGER	Gas ISOs
C	PIN(MAXPRO)	REAL	Pressure(atm)
C	TIN(MAXPRO)	REAL	Temperaures(K)
C	HIN(MAXPRO)	REAL	Heights (km)
C	VMRIN(MAXPRO,MAXGAS) REAL	Vmrs
C	XMOLIN(MAXPRO)	REAL	Molecular weight profile
C	NCONT	INTEGER		Number of cloud types
C	CONTIN(MAXCON,MAXPRO) REAL	Specific densities
C       FLUX	REAL		Convective heat flux (can set to
C				 STEF_BOLTZ*T_eff**4). Assume units of W m-2
C       IMODEL  INTEGER		Model. 0 = Lewis, 1=Ackerman
C       FRAIN	REAL		f_rain parameter of Ackerman model
C       JVMR	INTEGER		Gas IVMR to consider applying this model to
C	XDEEP	REAL		Deep vmr of gas IVMR
C       DENSCOND REAL		Density of condensate material (g/cm3)
C       RADCOND	REAL		Radius of condensate (microns)
C       MWCOND	REAL		Molecular weight (g) of condensed phase
C
C     Output variables
C	X1(MAXPRO)	REAL	Output vmr profile
C	X2(MAXPRO)	REAL	Output cloud density profile
C	QCOUT(MAXPRO)	REAL	Condensed mole fraction profile
C
C     	Pat Irwin	4/1/16	Original
C	Pat Irwin	15/1/20	Revised to apply to one gas only.	
C			
C     **************************************************************   
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE '../radtran/includes/constdef.f'

      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS)
      REAL CONT(MAXCON,MAXPRO),MOLWT,LATITUDE,R,DTDZ(MAXPRO)
      REAL QCONT(MAXCON,MAXPRO),FLUX,CP,RADIUS,DALR(MAXPRO)
      REAL CALCMOLWT,XVMR(MAXGAS),GASDATA(20,5),DELH,XDEEP
      REAL XMOLWT,CPSPEC,EDDY(MAXPRO),PIN(MAXPRO),TIN(MAXPRO)
      REAL HIN(MAXPRO),VMRIN(MAXPRO,MAXGAS),CONTIN(MAXCON,MAXPRO)
      REAL XMOLIN(MAXPRO),QCOUT(MAXPRO)
      REAL A,B,C,D,SVP,QS,DPEXP,SCALE(MAXPRO),G,X,L,LH
      REAL QV(MAXPRO),QT(MAXPRO),QC(MAXPRO),RHO(MAXPRO),XLAMBDA,KBOLTZ
      REAL WS(MAXPRO),FRAIN,DELQT,QT1,KMIN,X1(MAXPRO),X2(MAXPRO)
      REAL DENSCOND,RADCOND,MWCOND,N(MAXPRO),V,MP,M1,NM,NC
      REAL NAVAGADRO,NP,XMOL(MAXPRO),X1A(MAXPRO),X2A(MAXPRO)
      REAL CLOUD(MAXPRO),XOD(MAXPRO),XCOL1(MAXPRO),XCOL2(MAXPRO)
      REAL AA,BB,AREA,VOLUME
      INTEGER IDIST
      INTEGER K,AMFORM,NPRO,NVMR,I,IERR,NGAS,IPLANET,JVMR
      INTEGER IDGAS(MAXGAS),ISOGAS(MAXGAS),ISCALE(MAXGAS),J
      INTEGER JGAS,IVMR,NCONT,ICONDENSE,IMODEL,NPROIN
      CHARACTER*100 IPFILE,BUFFER,OPFILE
      CHARACTER*100 AEFILE,QCFILE,ANAME
      CHARACTER*8 PNAME
      PARAMETER (XLAMBDA = 0.1,KMIN=1e5,KBOLTZ=1.38064852E-23)
      PARAMETER (NAVAGADRO=6.022E23)
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

C      print*,IPLANET,LATITUDE,AMFORM,NPROIN,NVMR
C      print*,(IDGAS(J),J=1,NVMR)
C      print*,(ISOGAS(J),J=1,NVMR)
C      DO I=1,NPROIN
C        print*,PIN(I),TIN(I),HIN(I),XMOLIN(I),(VMRIN(I,J),J=1,NVMR)
C      ENDDO

C     RGAS read in from constdef.f, so set R accordingly and in correct
C     units of J mol-1 K-1
      R=RGAS*0.001

C     Set molar heat capacity at constant pressure (ideal gas) for a polyatomic
C     gas with translational and rotational degrees of freedom activated.
C     CP is J K-1 mol-1
      CP = 4*R



C     Calculate density of atmosphere (g/cm3) and DALR (K/km)


C     First interpolate profile on to a grid with step of ~1 km

      DELH = HIN(NPROIN)-HIN(1)
      NPRO = 1+INT(DELH/1.0)

      IF(NPRO.GT.MAXPRO)NPRO=MAXPRO

      DELH=(HIN(NPROIN)-HIN(1))/FLOAT(NPRO-1)


      if(idiag.gt.0)print*,'ack',NPROIN,HIN(1),HIN(NPROIN),DELH
      if(idiag.gt.0)print*,'ack NPRO = ',NPRO

      DO 101 J=1,NPRO
       H(J)=HIN(1)+(J-1)*DELH
       CALL VERINT(HIN,PIN,NPROIN,P(J),H(J))
       CALL VERINT(HIN,TIN,NPROIN,T(J),H(J))
       CALL VERINT(HIN,XMOLIN,NPROIN,XMOL(J),H(J))
C       print*,J,P(J),T(J),XMOL(J)
       DO I=1,NVMR
        DO K=1,NPROIN
         X1(K)=VMRIN(K,I)
        ENDDO
        CALL VERINT(HIN,X1,NPROIN,VMR(J,I),H(J))      
       ENDDO
       DO I=1,NCONT
        DO K=1,NPROIN
         X1(K)=CONTIN(I,K)
        ENDDO
        CALL VERINT(HIN,X1,NPROIN,CONT(I,J),H(J))      
       ENDDO

101   CONTINUE
 

      DO 201 J=1,NPRO
        XMOLWT=XMOL(J)
C       Compute density in g/cm3
        RHO(J) = P(J)*0.1013*XMOLWT/(R*T(J))
C       Compute number density in molecule/cm3
        N(J)=P(J)*0.1013/(KBOLTZ*T(J))
        CALL NEWGRAV(IPLANET,LATITUDE,H(J),RADIUS,G,PNAME) 
C       CPSPEC is J K-1 Kg-1
        CPSPEC = CP/(XMOLWT*1e-3)
C       DALR is in units of K/km
        DALR(J)=1000.*G/CPSPEC
C       Compute scale height (km)
        SCALE(J)=R*T(J)/(XMOLWT*G)
201   CONTINUE

C     Calculate DT_DZ
      CALL FIRSTDERIV(T,-H,NPRO,DTDZ)

      if(idiag.gt.0)then
       open(12,file='test.prf',status='unknown')
       write(12,*)NPRO
       write(12,*)'H(km) ,P(atm),T(K),SCALE(km),RHO(g/cm3),L(km),
     &EDDY(cm2/s),WS(cm/s),WS*FRAIN*DELH/EDDY'
      endif

      DO 202 J=1,NPRO
C      Calculate eddy diffusion coefficient
C      Compute x = R*F/mu*rho*Cp. Since Cp is ideal this simplifies to
C      F/(rho*4). Factor of 1e3 below converts g/cm3 to kg/m3
       X = FLUX/(RHO(J)*1e3*CP/R)
       LH = MAX(XLAMBDA,DTDZ(J)/DALR(J))
C      L is mixing length in km
       L=LH*SCALE(J)
C      EDDY(J) is eddy diffusion coefficient in cm2 s-1
       EDDY(J)=1e4*(1e3*SCALE(J)/3.0)*(LH**1.33333)*(X**0.333333)
       EDDY(J)=MAX(EDDY(J),KMIN)
C      WS is convective velocity scale in cm s-1
       WS(J) = EDDY(J)/(L*1e5)
C       DELH is height of layer in cm
       IF(J.LT.NPRO)THEN 
         DELH=1e5*(H(J+1)-H(J))
       ENDIF
C       WS(J) = 1000.
C       EDDY(J)=2e8
        if(idiag.gt.0)then
         write(12,*)H(J),P(J),T(J),SCALE(J),RHO(J),L,EDDY(J),WS(J),
     1    WS(J)*FRAIN*DELH/EDDY(J)
        endif
202   CONTINUE

      if(idiag.gt.0)close(12)

1     FORMAT(A)
      ANAME='SVP.dat'
      CALL DATARCHIVE(ANAME)
      OPEN(13,FILE=ANAME,STATUS='OLD')
C     First skip header
57    READ(13,1)BUFFER
      IF(BUFFER(1:1).EQ.'#')GOTO 57
      READ(BUFFER,*)NGAS
      DO 21 I=1,NGAS
       READ(13,*)(GASDATA(I,J),J=1,5)
21    CONTINUE
      CLOSE(13)


      if(idiag.gt.0)then
       print*,'JVMR, IDGAS(JVMR), ISOGAS(JVMR) = ',JVMR, IDGAS(JVMR),
     1  ISOGAS(JVMR)
      endif

C      CALL PROMPT('Enter IDIST,AA,BB : ')
C      READ*,IDIST,AA,BB

C      print*,IDIST,AA,BB
C      CALL SIZEDIST(IDIST,AA,BB,AREA,VOLUME)
C      print*,area,volume

C     Volume of condensed particles (cm3)
      V = 1.3333*PI*(RADCOND*1e-4)**3
C      V = VOLUME

      if(idiag.gt.0)print*,'R (micron),D (g/cm3) = ',RADCOND,DENSCOND
      if(idiag.gt.0)print*,'V (cm3) = ',V
      if(idiag.gt.0)print*,'old V (cm3) = ',1.3333*PI*(RADCOND*1e-4)**3

C     Mass of condensed particle (g) (DENSCOND is g/cm3)
      MP = DENSCOND*V
      if(idiag.gt.0)print*,'MP (g) = ',MP

C     Mass of one molecule (g)
      M1 = MWCOND/NAVAGADRO
      if(idiag.gt.0)then
       print*,'M1(g),mwcond (g),NAVAGADRO',M1,mwcond,NAVAGADRO
      endif

C     Number of molecules per condensed particle
      NM = MP/M1
      if(idiag.gt.0)print*,'Number of molecules/particle = ',NM


      DO 99 JGAS=1,NGAS
        IF(GASDATA(JGAS,1).EQ.IDGAS(JVMR))THEN
C        Target constituent will condense. Need to modify mole fraction and cloud
         A = GASDATA(JGAS,2)
         B = GASDATA(JGAS,3)
         C = GASDATA(JGAS,4)
         D = GASDATA(JGAS,5)
 
       
         CLOUD(1)=0.
     
         QV(1)=XDEEP
         QC(1)=0.
         QT(1)=QV(1)

         X1A(1)=QV(1)

         DO 97 J=2,NPRO
C         Calculate saturated vapour mole fraction
          SVP=DPEXP(A+B/T(J)+C*T(J)+D*T(J)*T(J))
          QS = SVP/P(J)

C         Calculate DELH in cm
          DELH = 1e5*(H(J)-H(J-1))

          IF(IMODEL.EQ.0)THEN

           QV(J)=MIN(QV(J-1),QS)
           QC(J)=MAX(0.0,QV(J-1)-QS)
           QT(J)=QV(J)+QC(J)

          ELSE

           DELQT = -DELH*FRAIN*WS(J-1)*QC(J-1)/EDDY(J-1)
           QT(J)=QT(J-1)+DELQT
           QC(J)=MAX(0.0,QT(J)-QS)
           QV(J)=MIN(QT(J),QS)
           QV(J)=MIN(QV(J),QV(J-1))

          ENDIF

C         X1A(J) is mole fraction of gas remaining
          X1A(J)=QV(J)

C         So QV is mole fraction of gas remaining and QC is equivalent mole 
C         fraction of condensed gas, QT is equivalent mole fraction of gas and condensate remaining . 
C         To convert QC into particles per gram we need to do the following:

C         NC is number of condensed molecules per cm3, since N(J) is number density of air (molecule/cm3)
          NC = QC(J)*N(J)

C         NP is number of condensed particles per cm3 since NM is number of molecules per particle
          NP = NC/NM

C         CLOUD is number of condensed particles/cm2 in layer
          CLOUD(J)=NP*DELH
  
C         X2A(J) is number of condensed particles/gram of atmosphere          
          X2A(J)=NP/RHO(J)

97       CONTINUE

        ENDIF

99    CONTINUE

      XOD(NPRO)=0.
      XCOL1(NPRO)=0.
      XCOL2(NPRO)=0.
      if(idiag.gt.0)then
       print*,'I, P(ATM), part./cm2, integ/cm2, cloud g/cm2, vap g/cm2,
     & deepvap g/cm2'
      endif   
      DO I=NPRO-1,1,-1
       DELH=1E5*(H(I+1)-H(I))
       XOD(I) = XOD(I+1)+CLOUD(I)
       XCOL1(I) = XCOL1(I+1)+RHO(I)*QV(I)*DELH		
       XCOL2(I) = XCOL2(I+1)+RHO(I)*QV(1)*DELH		
       if(idiag.gt.0)print*,I,P(I),CLOUD(I),XOD(I),XOD(I)*MP*1000.0,
     &  XCOL1(I),XCOL2(I)
      ENDDO

C     Interpolate answers back on to original grid
      DO 401 I=1,NPROIN
       CALL VERINT(H,X1A,NPRO,X1(I),HIN(I))
       CALL VERINT(H,X2A,NPRO,X2(I),HIN(I))
       CALL VERINT(H,QC,NPRO,QCOUT(I),HIN(I))
401   CONTINUE


      RETURN

      END


      SUBROUTINE FIRSTDERIV(Y,X,NX,DYDX)      
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      REAL X(MAXPRO),Y(MAXPRO),DYDX(MAXPRO)
      INTEGER I,NX
      DO 10 I=1,NX
       IF(I.EQ.1) THEN
        DYDX(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
       ELSE
        IF(I.EQ.NX)THEN
         DYDX(I)=(Y(I-1)-Y(I))/(X(I-1)-X(I))
        ELSE
         DYDX(I)=(Y(I+1)-Y(I-1))/(X(I+1)-X(I-1))
        ENDIF
       ENDIF
10    CONTINUE

      RETURN

      END


      SUBROUTINE SIZEDIST(IDIST,AA,BB,AREA,VOLUME)
C     ************************************************************
C     Computes size distribution and returns mean area and volume of
C     particles
C
C     Input variables:
C	IDIST	INTEGER	Distribution: 0=Gamma, 1=log-normal
C	AA	REAL	Mean radius (micron)
C	BB	REAL 	Variance
C
C     Returned variables:
C	AREA	REAL	Mean area per particle (m2)
C	VOLUME  REAL	Mean volume per particle (m3)
C     ************************************************************
      IMPLICIT NONE
      INTEGER IDIST,NINT,I
      REAL AA,BB,AREA,VOLUME,XNORM,PI
      PARAMETER (NINT=1000,PI=3.1415927)
      DOUBLE PRECISION R(NINT),R1,R2,DR,RX,N(NINT),NMAX,WT(NINT),F

C     Find necessary integration limits
      R1 = DBLE(AA/2.0)
      R2 = DBLE(AA*2.0)

103   NMAX=0.0
      DO I=1,NINT
       RX=R1+(R2-R1)*DBLE(I-1)/DBLE(NINT-1)
       R(I)=RX
       IF(IDIST.EQ.0)THEN
        N(I)=(RX**((1.0-3*BB)/BB))* EXP(-RX/(AA*BB))
       ELSE
        N(I)=(1.0/RX)*EXP(-(LOG(RX)-LOG(AA))**2/(2*BB**2))
       ENDIF
       IF(N(I).GT.NMAX)NMAX=N(I)
      ENDDO

      IF(N(1).GT.0.01*NMAX)THEN
       R1=R1/3.0
       GOTO 103
      ENDIF
      IF(N(NINT).GT.0.001*NMAX)THEN
       R2=R2*2.0
       GOTO 103
      ENDIF

      DO I=1,NINT
       WT(I)=1.0
      ENDDO
      WT(1)=0.5
      WT(NINT)=0.5
      DR=(R2-R1)/FLOAT(NINT-1)

      XNORM=0.0
      AREA=0.0
      VOLUME=0.0

      DO I=1,NINT
       F=WT(I)*N(I)*DR
       XNORM=XNORM+SNGL(F)
       AREA=AREA+SNGL(PI*(R(I)**2)*F)
       VOLUME=VOLUME+SNGL((4.0*PI/3.0)*(R(I)**3)*F)
      ENDDO

      AREA=AREA/XNORM
      VOLUME=VOLUME/XNORM
   
      AREA=AREA*1E-12
      VOLUME=VOLUME*1E-18

      RETURN

      END


