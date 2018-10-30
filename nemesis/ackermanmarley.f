      SUBROUTINE ACKERMANMARLEY(IPFILE,OPFILE,AEFILE,QCFILE,FLUX,
     1 IMODEL,FRAIN)
C     **************************************************************   
C     Subroutine to condense clouds as per the formalism of 
C     Ackerman and Marley (2001). 

C     Input variables
C    	IPFILE	CHARACTER*100	Name of input .prf file
C    	OPFILE	CHARACTER*100	Name of output .prf file
C    	AEFILE	CHARACTER*100	Name of output condensate profile
C				(particles/gram)
C    	QCFILE	CHARACTER*100	Name of output condensate profile
C				(particles/cm3)
C       FLUX	REAL		Convective heat flux (can set to
C				 STEF_BOLTZ*T_eff**4). Assume units of W m-2
C       IMODEL  INTEGER		Model. 0 = Lewis, 1=Ackerman
C       FRAIN	REAL		f_rain parameter of Ackerman model
C       
C     Pat Irwin	4/1/16
C
C     **************************************************************   
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE '../radtran/includes/constdef.f'

      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS)
      REAL CONT(MAXCON,MAXPRO),MOLWT,LATITUDE,R,DTDZ(MAXPRO)
      REAL QCONT(MAXCON,MAXPRO),FLUX,CP,RADIUS,DALR(MAXPRO)
      REAL CALCMOLWT,XVMR(MAXGAS),GASDATA(20,5),DELH
      REAL XMOLWT,XXMOLWT(MAXPRO),CPSPEC,EDDY(MAXPRO)
      REAL A,B,C,D,SVP,QS,DPEXP,SCALE(MAXPRO),G,X,L,LH
      REAL QV(MAXPRO),QT(MAXPRO),QC(MAXPRO),RHO(MAXPRO),XLAMBDA
      REAL WS(MAXPRO),FRAIN,DELQT,QT1,KMIN
      INTEGER K,AMFORM,IPLANET,NPRO,NVMR,I,IERR,NGAS
      INTEGER IDGAS(MAXGAS),ISOGAS(MAXGAS),ISCALE(MAXGAS),J
      INTEGER JGAS,IVMR,NCONT,ICONDENSE,IMODEL
      CHARACTER*100 IPFILE,BUFFER,OPFILE
      CHARACTER*100 AEFILE,QCFILE,ANAME
      CHARACTER*8 PNAME
      PARAMETER (XLAMBDA = 0.1,KMIN=1e5)
      CALL FILE(IPFILE,IPFILE,'prf')
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
      MOLWT=-1.
C     First skip header
54    READ(1,1)BUFFER
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)AMFORM
1     FORMAT(A)

      IF(AMFORM.EQ.0)THEN
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
      ELSE
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR
      ENDIF
      print*,IPLANET,LATITUDE,NPRO,NVMR,MOLWT

      IF(NPRO.GT.MAXPRO)THEN
          PRINT*,'Error in subprofretg. NPRO>MAXPRO ',NPRO,MAXPRO
          STOP
      ENDIF

      DO 20 I=1,NVMR
       READ(1,*)IDGAS(I),ISOGAS(I)
       ISCALE(I)=1
20    CONTINUE

C     Skip header
      READ(1,*)
      DO 30 I=1,NPRO
        READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,NVMR)
30    CONTINUE

      CLOSE(UNIT=1)

C     Make sure that vmrs add up to 1 if AMFORM=1
      IF(AMFORM.EQ.1)THEN
        CALL ADJUSTVMR(NPRO,NVMR,VMR,ISCALE,IERR)

        DO 301 I=1,NPRO
         DO K=1,NVMR
          XVMR(K)=VMR(I,K)
         ENDDO
         XXMOLWT(I)=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
301     CONTINUE
C     Calculate molecular weight but not add vmrs to 1 if AMFORM=2
      ELSEIF(AMFORM.EQ.2)THEN
        DO I=1,NPRO
         DO K=1,NVMR
          XVMR(K)=VMR(I,K)
         ENDDO
         XXMOLWT(I)=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
        ENDDO
      ENDIF

C     Make sure profile is hydrostatically balanced and compute
C     atmospheric scale height.
      CALL XHYDROSTATH(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1  IDGAS,ISOGAS,H,P,T,VMR,SCALE)

C     RGAS read in from constdef.f, so set R accordingly and in correct
C     units of J mol-1 K-1
      R=RGAS*0.001

C     Set molar heat capacity at constant pressure (ideal gas) for a polyatomic
C     gas with translational and rotational degrees of freedom activated.
C     CP is J K-1 mol-1
      CP = 4*R

C     Calculate density of atmosphere (g/cm3) and DALR (K/km)
      DO 201 J=1,NPRO
        IF(AMFORM.EQ.0)THEN
          XMOLWT=MOLWT
        ELSE
          XMOLWT=XXMOLWT(J)
        ENDIF
        RHO(J) = P(J)*0.1013*XMOLWT/(R*T(J))
C        print*,P(J),XMOLWT,R,T(J),RHO(J)
        CALL NEWGRAV(IPLANET,LATITUDE,H(J),RADIUS,G,PNAME) 
C       CPSPEC is J K-1 Kg-1
        CPSPEC = CP/(XMOLWT*1e-3)
C       DALR is in units of K/km
        DALR(J)=1000.*G/CPSPEC
201   CONTINUE

C     Calculate DT_DZ
      CALL FIRSTDERIV(T,-H,NPRO,DTDZ)

C     Calculate eddy diffusion coefficient

C     Compute R*F/mu*rho*Cp. Since Cp is ideal this simplifies to
C     F/(rho*4) )

      open(12,file='test.prf',status='unknown')
      write(12,*)NPRO
      write(12,*)'H(J),P(J),T(J),SCALE(J),RHO(J),L,EDDY(J),WS(J)'
      DO 202 J=1,NPRO
C      Factor of 1e3 below converts g/cm3 to kg/m3
       X = FLUX/(RHO(J)*1e3*CP/R)
       LH = MAX(XLAMBDA,DTDZ(J)/DALR(J))
C      L is mixing length in km
       L=LH*SCALE(J)
C      EDDY(J) is eddy diffusion coefficeint in cm2 s-1
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

      ANAME='SVP.dat'
      CALL DATARCHIVE(ANAME)
      OPEN(13,FILE=ANAME,STATUS='OLD')
      WRITE(*,*)' '
      WRITE(*,*)'SUBPROFRETG: reading saturated-vapour-pressure'
      WRITE(*,*)'  data from ',ANAME
C     First skip header
57    READ(13,1)BUFFER
      IF(BUFFER(1:1).EQ.'#')GOTO 57
      READ(BUFFER,*)NGAS
      DO 21 I=1,NGAS
       READ(13,*)(GASDATA(I,J),J=1,5)
21    CONTINUE
      CLOSE(13)


      NCONT=0

      DO 101 IVMR=1,NVMR
       DO 99 JGAS=1,NGAS
        IF(GASDATA(JGAS,1).EQ.IDGAS(IVMR))THEN
C        Constituent may condense. May need to modify mole fraction and cloud
         A = GASDATA(JGAS,2)
         B = GASDATA(JGAS,3)
         C = GASDATA(JGAS,4)
         D = GASDATA(JGAS,5)
 
     
         QV(1)=VMR(1,IVMR)
         QC(1)=0.
         QT(1)=QV(1)
         ICONDENSE=0

         DO 97 J=2,NPRO
C         Calculate saturated vapour mole fraction
          SVP=DPEXP(A+B/T(J)+C*T(J)+D*T(J)*T(J))
          QS = SVP/P(J)

          IF(IMODEL.EQ.0)THEN

           QV(J)=MIN(QV(J-1),QS)
           QC(J)=MAX(0.0,QV(J-1)-QS)
           QT(J)=QV(J)+QC(J)

           print*,J,QT(J),QV(J),QC(J)

          ELSE

           DELH = 1e3*(H(J)-H(J-1))
           DELQT = -DELH*FRAIN*WS(J-1)*QC(J-1)/EDDY(J-1)
           QT(J)=QT(J-1)+DELQT
           QC(J)=MAX(0.0,QT(J)-QS)
           QV(J)=MIN(QT(J),QS)

           print*,J,QT(J),QV(J),QC(J),DELQT

          ENDIF


          IF(QC(J).GT.0.0)THEN
            ICONDENSE=1
          ENDIF

97       CONTINUE

         IF(ICONDENSE.GT.0)THEN
           print*,'Gas Condenses : ',IVMR
           print*,'ID,ISO : ',IDGAS(IVMR),ISOGAS(IVMR)
C          Adjust vmr of gas in question
           print*,IVMR,NPRO
           DO 95 J=1,NPRO
            VMR(J,IVMR)=QV(J)
C            print*,J,QV(J)
95         CONTINUE

C          Check to see if vmrs should add up to 1.0 If so scale the other
C          gases
           print*,'A',AMFORM
           IF(AMFORM.EQ.1)THEN
            DO I=1,NVMR
             ISCALE(I)=1
            ENDDO
            ISCALE(IVMR)=0
            CALL ADJUSTVMR(NPRO,NVMR,VMR,ISCALE,IERR)

            DO 302 I=1,NPRO
             DO K=1,NVMR
              XVMR(K)=VMR(I,K)
             ENDDO
             XXMOLWT(I)=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
302         CONTINUE
C          Calculate molecular weight but not add vmrs to 1 if AMFORM=2
           ELSEIF(AMFORM.EQ.2)THEN
            DO I=1,NPRO
             DO K=1,NVMR
              XVMR(K)=VMR(I,K)
             ENDDO
             XXMOLWT(I)=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
            ENDDO
           ENDIF

           NCONT=NCONT+1
           DO 96 J=1,NPRO
            CONT(NCONT,J)=QC(J)/RHO(J)
            QCONT(NCONT,J)=QC(J)
96         CONTINUE

         ENDIF

        ENDIF

99     CONTINUE
101   CONTINUE


C     ************* Write out modified profiles *********

      CALL FILE(OPFILE,OPFILE,'prf')
      OPEN(UNIT=2,FILE=OPFILE,STATUS='UNKNOWN',ERR=52)
      WRITE(2,*)AMFORM
      IF(AMFORM.EQ.0)THEN
       WRITE(2,501)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
      ELSE
       WRITE(2,501)IPLANET,LATITUDE,NPRO,NVMR
      ENDIF
501   FORMAT(1X,I3,F7.2,1X,I3,I3,F8.3)
      DO 503 I=1,NVMR
        WRITE(2,502)IDGAS(I),ISOGAS(I)
502     FORMAT(1X,I3,I5)
503   CONTINUE
      WRITE(2,504)(I,I=1,NVMR)
504   FORMAT(1X,' height (km) ',' press (atm) ','  temp (K)   ',
     1  40(' VMR gas',I3,2X))
      DO 505 I=1,NPRO
          WRITE(2,506)H(I),P(I),T(I),(VMR(I,J),J=1,NVMR)
506       FORMAT(1X,F13.3,E13.5,F13.4,40(E13.5))
505   CONTINUE

52    CONTINUE
      CLOSE(2)

      CALL FILE(AEFILE,AEFILE,'prf')
      OPEN(UNIT=2,FILE=AEFILE,STATUS='UNKNOWN')
      BUFFER='# '//AEFILE
      WRITE(2,10)BUFFER
      WRITE(2,*)NPRO, NCONT
      DO 41 I=1,NPRO
        WRITE(2,*) H(I),(CONT(J,I),J=1,NCONT)
41    CONTINUE
      CLOSE(2)


      CALL FILE(QCFILE,QCFILE,'prf')
      OPEN(UNIT=2,FILE=QCFILE,STATUS='UNKNOWN')
      BUFFER='# '//QCFILE
      WRITE(2,10)BUFFER
      WRITE(2,*)NPRO, NCONT
      DO 42 I=1,NPRO
        WRITE(2,*) H(I),(QCONT(J,I),J=1,NCONT)
42    CONTINUE
      CLOSE(2)

10    FORMAT(A)

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
