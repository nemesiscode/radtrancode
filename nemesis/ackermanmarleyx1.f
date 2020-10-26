      SUBROUTINE ACKERMANMARLEYX1(IPLANET,LATITUDE,AMFORM,NPROIN,NVMR,
     1 IDGAS,ISOGAS,PIN,TIN,HIN,VMRIN,XMOLIN,NCONT,CONTIN,FLUX,IMODEL,
     2 FRAIN,JVMR,DENSCOND,RADCOND,MWCOND,X1,X2)
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
C       DENSCOND REAL		Density of condensate material (g/cm3)
C       RADCOND	REAL		Radius of condensate (microns)
C       MWCOND	REAL		Molecular weight (g) of condensed phase
C
C     Output variables
C	X1(MAXPRO)	REAL	Output vmr profile
C	X2(MAXPRO)	REAL	Output cloud density profile
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
      REAL CALCMOLWT,XVMR(MAXGAS),GASDATA(20,5),DELH
      REAL XMOLWT,CPSPEC,EDDY(MAXPRO),PIN(MAXPRO),TIN(MAXPRO)
      REAL HIN(MAXPRO),VMRIN(MAXPRO,MAXGAS),CONTIN(MAXCON,MAXPRO)
      REAL XMOLIN(MAXPRO)
      REAL A,B,C,D,SVP,QS,DPEXP,SCALE(MAXPRO),G,X,L,LH
      REAL QV(MAXPRO),QT(MAXPRO),QC(MAXPRO),RHO(MAXPRO),XLAMBDA,KBOLTZ
      REAL WS(MAXPRO),FRAIN,DELQT,QT1,KMIN,X1(MAXPRO),X2(MAXPRO)
      REAL DENSCOND,RADCOND,MWCOND,N(MAXPRO),V,MP,M1,NM,NC
      REAL NAVAGADRO,NP,XMOL(MAXPRO),X1A(MAXPRO),X2A(MAXPRO)
      INTEGER K,AMFORM,NPRO,NVMR,I,IERR,NGAS,IPLANET,JVMR
      INTEGER IDGAS(MAXGAS),ISOGAS(MAXGAS),ISCALE(MAXGAS),J
      INTEGER JGAS,IVMR,NCONT,ICONDENSE,IMODEL,NPROIN
      CHARACTER*100 IPFILE,BUFFER,OPFILE
      CHARACTER*100 AEFILE,QCFILE,ANAME
      CHARACTER*8 PNAME
      PARAMETER (XLAMBDA = 0.1,KMIN=1e5,KBOLTZ=1.38064852E-23)
      PARAMETER (NAVAGADRO=6.022e23)
C     RGAS read in from constdef.f, so set R accordingly and in correct
C     units of J mol-1 K-1
      R=RGAS*0.001

C     Set molar heat capacity at constant pressure (ideal gas) for a polyatomic
C     gas with translational and rotational degrees of freedom activated.
C     CP is J K-1 mol-1
      CP = 4*R

C      print*,'Test ',1.0*0.1013*28.0/(R*273.0)
C     Calculate density of atmosphere (g/cm3) and DALR (K/km)


C     First interpolate profile on to a grid with step of ~1 km

      DELH = HIN(NPROIN)-HIN(1)
      NPRO = 1+INT(DELH/1.0)

      IF(NPRO.GT.MAXPRO)NPRO=MAXPRO

      DELH=(HIN(NPROIN)-HIN(1))/FLOAT(NPRO-1)

      print*,NPROIN,HIN(1),HIN(NPROIN),DELH
      print*,'NPRO = ',NPRO

      DO 101 J=1,NPRO
       H(J)=HIN(1)+(J-1)*DELH
       print*,J,H(J)
       CALL VERINT(HIN,PIN,NPROIN,P(J),H(J))
       CALL VERINT(HIN,TIN,NPROIN,T(J),H(J))
       CALL VERINT(HIN,XMOLIN,NPROIN,XMOL(J),H(J))
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
       print*,J,H(J)

101   CONTINUE
 
      DO 102 J=1,NPRO
       print*,H(J),P(J),T(J),(VMR(J,K),K=1,NVMR)
102   CONTINUE
      DO 103 J=1,NPRO
       print*,H(J),P(J),(CONT(K,J),K=1,NCONT)
103   CONTINUE

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
C        print*,J,P(J),T(J),RHO(J),SCALE(J),XMOLWT
201   CONTINUE

C     Calculate DT_DZ
      CALL FIRSTDERIV(T,-H,NPRO,DTDZ)

      open(12,file='test.prf',status='unknown')
      write(12,*)NPRO
      write(12,*)'H(km) ,P(atm),T(K),SCALE(km),RHO(g/cm3),L(km),
     &EDDY(cm2/s),WS(cm/s)'

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
C       WS(J) = 1000.
C       EDDY(J)=2e8
C       print*,J,H(J),SCALE(J),RHO(J),L,EDDY(J),LH,X,WS(J)
        write(12,*)H(J),P(J),T(J),SCALE(J),RHO(J),L,EDDY(J),WS(J)
202   CONTINUE

      close(12)

1     FORMAT(A)
      ANAME='SVP.dat'
      CALL DATARCHIVE(ANAME)
      OPEN(13,FILE=ANAME,STATUS='OLD')
C      WRITE(*,*)' '
C      WRITE(*,*)'ackermanmarleyx1: reading saturated-vapour-pressure'
C      WRITE(*,*)'  data from ',ANAME
C     First skip header
57    READ(13,1)BUFFER
      IF(BUFFER(1:1).EQ.'#')GOTO 57
      READ(BUFFER,*)NGAS
      DO 21 I=1,NGAS
       READ(13,*)(GASDATA(I,J),J=1,5)
21    CONTINUE
      CLOSE(13)


C      PRINT*,'JVMR, IDGAS(JVMR), ISOGAS(JVMR) = ',JVMR, IDGAS(JVMR),
C     1  ISOGAS(JVMR)


C     Volume of condensed particles (m3)
      V = 1.3333*PI*(RADCOND*1e-6)**3
C      print*,'R,D = ',RADCOND,DENSCOND
C      print*,'V = ',V
C     Mass of condensed particle (kg)
      MP = DENSCOND*V
C      print*,'MP = ',MP
C     Mass of one molecule
      M1 = 1E-3*MWCOND/NAVAGADRO
C     Number of molecules per condensed particle
C      print*,'M1,mwcond,NAVAGADRO',M1,mwcond,NAVAGADRO
      NM = MP/M1

      DO 99 JGAS=1,NGAS
        IF(GASDATA(JGAS,1).EQ.IDGAS(JVMR))THEN
C        Constituent may condense. May need to modify mole fraction and cloud
         A = GASDATA(JGAS,2)
         B = GASDATA(JGAS,3)
         C = GASDATA(JGAS,4)
         D = GASDATA(JGAS,5)
 
       
     
         QV(1)=VMR(1,JVMR)
         QC(1)=0.
         QT(1)=QV(1)

         X1A(1)=QV(1)

C         IF(IMODEL.EQ.0)THEN
C           PRINT*,'J, QT(J), QV(J), QC(J)'
C         ELSE
C           PRINT*,'J, QT(J), QV(J), QC(J), DELQT'
C         ENDIF

         DO 97 J=2,NPRO
C         Calculate saturated vapour mole fraction
          SVP=DPEXP(A+B/T(J)+C*T(J)+D*T(J)*T(J))
          QS = SVP/P(J)

          IF(IMODEL.EQ.0)THEN

           QV(J)=MIN(QV(J-1),QS)
           QC(J)=MAX(0.0,QV(J-1)-QS)
           QT(J)=QV(J)+QC(J)

C           print*,J,QT(J),QV(J),QC(J)

          ELSE

C          Calculate DELH in cm
           DELH = 1e5*(H(J)-H(J-1))
           DELQT = -DELH*FRAIN*WS(J-1)*QC(J-1)/EDDY(J-1)
           QT(J)=QT(J-1)+DELQT
           QC(J)=MAX(0.0,QT(J)-QS)
           QV(J)=MIN(QT(J),QS)
           QV(J)=MIN(QV(J),QV(J-1))

C           print*,J,QT(J),QV(J),QC(J),DELQT

          ENDIF

          X1A(J)=QV(J)
C         So QV is mole fraction of gas remaining and QC is equivalent mole 
C         fraction of condensed gas, QT is equivalent mole fraction of gas and condensate remaining . 
C         To convert QC into particles per gram we need to do the following:

C         NC is number of condensed molecules per cm3, since N(J) is number density of air (molecule/cm3)
          NC = QC(J)*N(J)
C          print*,J,QC(J),N(J),NC
C         NP is number of condensed particles per cm3 since NM is number of molecules per particle
          NP = NC/NM
C          print*,NM,NP
C         X2A(J) is number of condensed particles/gram of atmosphere
          
          X2A(J)=NP/RHO(J)
C          print*,J,NP,RHO(J),X2A(J)
97       CONTINUE

        ENDIF

99    CONTINUE

C     Interpolate answers back on to original grid
      DO 401 I=1,NPROIN
       CALL VERINT(H,X1A,NPRO,X1(I),HIN(I))
       CALL VERINT(H,X2A,NPRO,X2(I),HIN(I))
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
