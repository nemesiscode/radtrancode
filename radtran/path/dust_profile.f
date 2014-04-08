      PROGRAM DUST_PROFILE
C     $Id: dust_profile.f,v 1.13 2010-02-11 15:43:05 fletcher Exp $
C     *********************************************************************
C     Generates dust specific concentration profile file for use by GENBAND
C     where specific concentration is the number of dust particles per gram 
C     of atmosphere.
C
C     Pat Irwin    13/12/93
C
C     *********************************************************************
      IMPLICIT NONE
      INCLUDE '../includes/constdef.f'
      INCLUDE '../includes/dbcom.f'
C      CHARACTER*100 IPFILE
      CHARACTER*200 BUFFER
      CHARACTER*8 PNAME
      INTEGER I,J,NPRO,NCONT,NVMR,ID(15),ISO(15),N,AMFORM,IPLANET
      REAL CONT(100,1000),VMR(1000,15),CP0,CXS,CSUM(1000)
      REAL H(1000),P(1000),T(1000),MOLWT,HC,RADIUS,H1,P1(1000)
      INTEGER IDUST,K,M,M1,M2,IP
      REAL DQNU,DELH,R1,R2,CALCMOLWT,XMOLWT
      REAL XVMR(20)
      REAL SCALEH(1000),G
      REAL CLOUDFSH,CLOUDP,FACT,Q0,SUM,S(1000),CLOUDP1,CLOUDP2
      REAL D0,XM1,D,RHO(1000),LATITUDE
C     *******************************************

      CALL RESERVEGAS

      CALL PROMPT('Enter name of temp/press profile file?')
      READ(*,10)IPFILE
10    FORMAT(A)
      CALL FILE(IPFILE,IPFILE,'prf')
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
54    READ(1,10)BUFFER
      IF(BUFFER(1:1).EQ.'#') GOTO 54
      READ(BUFFER,*)AMFORM
      IF(AMFORM.EQ.1)THEN
       READ(1,*)IPLANET,LATITUDE,NPRO,NVMR
      ELSE
       READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
      ENDIF
      DO 20 I=1,NVMR
      READ(1,*)ID(I),ISO(I)
20    CONTINUE
C     Skip Header
      READ(1,*)
      DO 30 I=1,NPRO
        READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,NVMR)
30    CONTINUE
      CLOSE(UNIT=1)
C     ******************************************

      CP0 = P(1)

      CALL PROMPT('Enter name of dust file : ')
      READ(*,10)IPFILE

      CALL FILE(IPFILE,IPFILE,'prf')


      CALL PROMPT('Enter number of particle types : ')
      READ(5,*) NCONT

      DO 1000 J=1,NCONT
    
      WRITE(*,*)' For dust/aerosol particle type ',J,', choose model:'
      WRITE(*,*)'   (1) Homogeneous (constant concentration with H)'
      WRITE(*,*)'   (2) Steady state (developed for Mars: constant'
      WRITE(*,*)'       at low altitudes, thereafter falls off at a' 
      WRITE(*,*)'       defined rate)'
      WRITE(*,*)'   (3) Cloud layers (generate simple cloud model via'
      WRITE(*,*)'       the definition of base pressure, fractional' 
      WRITE(*,*)'       scale height, cloud density and cloud-particle'
      WRITE(*,*)'       mass)'
      WRITE(*,*)'   (4) New cloud layers (generate simple cloud model'
      WRITE(*,*)'       via the definition of cloudbase pressure,'
      WRITE(*,*)'       cloudtop pressure, fractional scale height,'
      WRITE(*,*)'       required OD, and reference pressure and cloud'
      WRITE(*,*)'       cross-section. This model *is* different to'
      WRITE(*,*)'       that defined in option #3)'
      WRITE(*,*)'   (5) As model 3, but using actual cloud scale height'
      WRITE(*,*)'       (in km)'
      WRITE(*,*)' Model selection: '
      READ*,IDUST

      IF(IDUST.LT.3)THEN
       PRINT*,'Enter q0, the specific limiting aerosol concentration '
       CALL PROMPT('of this distribution : ')
       READ*,Q0
      END IF

      IF(IDUST.EQ.2)THEN
       CALL PROMPT('Enter inflection coefficient (0.007 Mars) : ')
       READ*,DQNU
      END IF

      IF(IDUST.EQ.3)THEN
       PRINT*,'Enter base pressure(atm) and fractional scale height'
       CALL PROMPT('of cloud : ')
       READ*,CLOUDP,CLOUDFSH
       PRINT*,'Enter Cloud density at bottom (g/cm3)'
       READ*,D0
       PRINT*,'Enter mass of single cloud particle'
       READ*,XM1
      END IF

      IF(IDUST.EQ.4)THEN
       PRINT*,'Enter base P(atm) top P and frac. scale height'
       CALL PROMPT('of cloud : ')
       READ*,CLOUDP1,CLOUDP2,CLOUDFSH
       PRINT*,'Enter required OD and base pressure above which you'
       CALL PROMPT('want the cloud density to be integrated : ')
       READ*,D0,CP0
       PRINT*,'Enter the cloud x-section in the .xsc file at the '
       PRINT*,'wavelength where you want the cloud to have the' 
       CALL PROMPT('required optical depth : ')
       READ*,CXS
      END IF

      IF(IDUST.EQ.5)THEN
       PRINT*,'Enter base pressure(atm) and scale height (km)'
       CALL PROMPT('of cloud : ')
       READ*,CLOUDP,CLOUDFSH
       PRINT*,'Enter Cloud density at bottom (g/cm3)'
       READ*,D0
       PRINT*,'Enter mass of single cloud particle'
       READ*,XM1
      END IF

C     Compute variation of scale heights and densities in atmosphere
      DO 50 I=1,NPRO

          CALL NEWGRAV(IPLANET,LATITUDE,H(I),RADIUS,G,PNAME)

          IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
          ELSE
	    GASFIL='gasinfo04.dat'
	    CALL DATARCHIVE(GASFIL)
	    DBLUN=12
	    CALL RDGAS
            DO K=1,NVMR
             XVMR(K)=VMR(I,K)
            ENDDO
            XMOLWT=CALCMOLWT(NVMR,XVMR,ID,ISO)
          ENDIF

          SCALEH(I)=RGAS*T(I)/(XMOLWT*G*1000.0)
          RHO(I)=P(I)*1.013E5*1E-6*XMOLWT/(8.31*T(I))
          CONT(J,I)=0.
 	  S(I)=0.0
50    CONTINUE

C     *****************************************************
C     Note: In subsequent code the following units apply:
C      CONT() dust particles/gram
C      S()    dust particles/cm3
C     *****************************************************
      IF(IDUST.GE.0.AND.IDUST.LE.2)THEN
c                                 ***   Column tau above tan height
         DO 55 I=1,NPRO
         
          IF(IDUST.EQ.1)THEN
           CONT(J,I)=Q0
          ELSE
           CONT(J,I)=Q0*EXP(DQNU*(1 - EXP(H(I)/SCALEH(I))))
          END IF
 	  S(I)=CONT(J,I)*RHO(I)		! particles/cm3

55       CONTINUE

       ELSE IF(IDUST.EQ.3)THEN
        DO I=1,NPRO
         CONT(J,I)=0.
         IF(P(I).GE.CLOUDP)THEN
          M=I
         ENDIF
        ENDDO
        PRINT*,'Setting cloud base to pressure (atm) ',P(M)
        HC = CLOUDFSH*SCALEH(I)
        PRINT*,'Scale height (km) = ',HC
        PRINT*,'I,S(I),CONT(J,I)'
        DO I=M,NPRO
         FACT=(H(I)-H(M))/HC
         D=D0*EXP(-FACT)
         S(I) = D/XM1
         CONT(J,I)=S(I)/RHO(I)
        END DO

       ELSE IF(IDUST.EQ.4)THEN
        M1=-1
        M2=-1
	print*,'CLOUDP1,CLOUDP2,CLOUDFSH',CLOUDP1,CLOUDP2,CLOUDFSH
        DO I=1,NPRO
         IF(P(I).GE.CLOUDP1)THEN
          M1=I
         ENDIF
	 IF(P(I).GE.CLOUDP2)THEN
	  M2=I+1
         ENDIF
        ENDDO
        IF(M1.LT.0)M1=1
        IF(M2.LT.0.OR.M2.GT.NPRO)M2=NPRO
        print*,'M1,M2 = ',M1,M2
        PRINT*,'Setting cloud base to pressure (atm) ',P(M1)
        PRINT*,'Setting cloud top to pressure (atm) ',P(M2)

C       Calculate air density in g/cm3
        S(M1) = 1.0		 ! number of particles/cm3
        CONT(J,M1)=S(M1)/RHO(M1) ! number of particles/gram of atmosphere

        print*,'Number of particles/cm3 at cloud base = ',S(M1)
        print*,'Density at cloud base (g/cm3) = ',RHO(M1)
	print*,'Pressure, temperature = ',P(M1),T(M1)
        print*,'S(M1),CONT(J,M1) = ',S(M1),CONT(J,M1)

        DO I=M1+1,M2
         H1 = 0.5*(SCALEH(I-1)+SCALEH(I))
         S(I) = S(I-1)*EXP( (H(I-1)-H(I)) / (CLOUDFSH*H1))
         CONT(J,I)=S(I)/RHO(I)
        END DO

       ELSE IF(IDUST.EQ.5)THEN
        DO I=1,NPRO
         CONT(J,I)=0.
         IF(P(I).GE.CLOUDP)THEN
          M=I
         ENDIF
        ENDDO
        PRINT*,'Setting cloud base to pressure (atm) ',P(M)
        HC = CLOUDFSH
        PRINT*,'Scale height (km) = ',HC
        PRINT*,'I,S(I),CONT(J,I)'
        DO I=M,NPRO
         FACT=(H(I)-H(M))/HC
         D=D0*EXP(-FACT)
         S(I) = D/XM1
         CONT(J,I)=S(I)/RHO(I)
        END DO
       ELSE
          WRITE(6,*) 'Dust model specified as something other than 1-5'
          WRITE(6,*) 'Stopped in Dust_profile'
          STOP
       END IF

C      Now integrate the specific dust concentrations to estimate the total
C      optical depth. This section used to use a nifty integration scheme which
C      was mathematically more accurate than the straight trapezium rule 
C      integration. Unfortunately, it was not very compatible with the subsequent
C      layer property calculations which just read the aerosol.prf file and assumed
C      a linear interpolation between layers. Hence, the trapezium rule calculation
C      is actually more appropriate here! 
C      Since we know that the the number density of particles drops as
C      S = S0exp(-z/H), then integrating from z1 to z2 we find that 
C      the integral is H(S(z1)-S(z2)) where H is the scale height. 

       SUM=0
       DO 59 I=NPRO,1,-1

 	 IF(I.EQ.NPRO)THEN
C          Find optical depth above top layer, using the scale height
           SUM = SUM + S(I)*SCALEH(I)*1E5

         ELSE

           DELH = H(I+1)-H(I)
           SUM = SUM + DELH*(S(I)+S(I+1))*0.5*1E5

         ENDIF

         CSUM(NPRO-I+1)=SUM
         P1(NPRO-I+1)=P(I)

59     CONTINUE

       PRINT*,'Integration from space to bottom of atmosphere '
       PRINT*,'Number of particles / cm2 in nadir for type : ',J 
       PRINT*,'( = Optical Depth / X-sectional area) = ',SUM

       CALl VERINT(P1,CSUM,NPRO,SUM,CP0)

       PRINT*,'Integration to pressure level : ',CP0,' is ',SUM

       IF(IDUST.EQ.4)THEN
        PRINT*,'Modifying CONT: to get right OD'
        DO I =1,NPRO
         CONT(J,I)=CONT(J,I)*D0/(SUM*CXS)
        ENDDO
       ENDIF

1000  CONTINUE

      OPEN(UNIT=1,FILE=IPFILE,STATUS='UNKNOWN')
      BUFFER='# '
      DO I=1,LEN(IPFILE)
       BUFFER(I+2:I+2)=IPFILE(I:I)
      END DO
      WRITE(1,10)BUFFER
      WRITE(1,*)NPRO, NCONT
      DO 31 I=1,NPRO
        WRITE(1,*) H(I),(CONT(J,I),J=1,NCONT)
31    CONTINUE
      CLOSE(UNIT=1)

      PRINT*,'*************** PLEASE NOTE **********************'
      PRINT*,'Please note that the normalisation scheme used to' 
      PRINT*,'calculate integrated opacities in this program are'
      PRINT*,'not entirely consistent with the the final mean'
      PRINT*,'layer properties calculated in the radiative transfer'
      PRINT*,'programs - Radtrans and Nemesis.'
      PRINT*,'   You should check that the aerosol column amounts'
      PRINT*,'reported by these programs (or Path) match your'
      PRINT*,'expectations! If not, then you will need to scale the'
      PRINT*,'aerosol profile accordingly'
      PRINT*,'**************************************************'

      END
