      SUBROUTINE SUBPROFRETGDISC(IPFILE,ISPACE,ISCAT,GASGIANT,XLAT,NVAR,
     1 VARIDENT,VARPARAM,NX,XN,JPRE,NCONT,FLAGH2P,XMAP)
C     $Id:
C     ***********************************************************************
C     Subroutine to modify an existing ipfile.ref T/P/vmr profile and 
C     aerosol.ref cloud density file according to the contents of the state
C     vector XN. New profiles are written to ipfile.prf and aerosol.prf 
C     files respectively.
C
C     Code also warns if vapour pressures exceeding SVP but does not
C     reset vmrs since the temperatures could also be in error.
C 
C     Code also shifts .ref file to the required latitude and recalulates
C     the heights from the hydrostatic equation.
C
C     Routine also returns the XMAP matrix which relates the functional
C     derivatives calculated by CIRSADG with the elements of the state
C     vector  
C     
C     Input variables
C	IPFILE	CHARACTER*100	Root run name
C	ISPACE	INTEGER		0=cm-1, 1=microns
C	ISCAT	INTEGER		Scattering indicator
C	GASGIANT LOGICAL	Indicates if planet is a Gas Giant
C	XLAT	REAL		Latitude of spectrum to be simulated
C	NVAR	INTEGER		Number of variable profiles
C	VARIDENT(MVAR,3) INTEGER Identity of profiles and parameterisation
C					scheme
C	VARPARAM(MVAR,MPARAM) REAL	Additional parameterisation
C	NX	INTEGER		Number if elements in state vector
C	XN(MX)	REAL		State vector
C     Output variables
C	NCONT	INTEGER		Number of cloud particle types
C	FLAGH2P INTEGER		Set to 1 if para-H2 profile is variable
C	XMAP(MAXV,MAXGAS+2+MAXCON,MAXPRO) REAL Matrix relating functional
C				derivatives calculated by CIRSRADG to the 
C				elements of the state vector
C
C     Pat Irwin	29/7/96		Original
C     Pat Irwin 17/10/03	Revised for Nemesis
C
C     ***********************************************************************
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

      REAL XN(MX),DPEXP,DELH,XFAC,DXFAC,XTMP,SUM,PMIN
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS)
      REAL CONT(MAXGAS,MAXPRO),XLAT,X,XREF(MAXPRO),X1(MAXPRO)
      REAL PKNEE,HKNEE,XDEEP,XFSH,PARAH2(MAXPRO),XH,XKEEP,PSHLD,HSHLD
      REAL OD(MAXPRO),ND(MAXPRO),Q(MAXPRO),RHO,F,XOD,DQDX(MAXPRO)
      REAL DNDH(MAXPRO),DQDH(MAXPRO),FCLOUD(MAXPRO)
      INTEGER ICLOUD(MAXGAS,MAXPRO),NCONT1
      INTEGER NPRO,NPRO1,NVMR,JZERO,IV,IP,IVAR,JCONT,JVMR
      INTEGER IDGAS(MAXGAS),ISOGAS(MAXGAS),IPAR,IVMR,NP
      INTEGER IDAT,NCONT,FLAGH2P,JFSH,JPRE,JHYDRO,ISCAT
      REAL HTAN,PTAN
      REAL G,R
      PARAMETER (R=8.31)
      CHARACTER*8 PNAME   
      INTEGER JSWITCH
      REAL GASDATA(20,5)
      CHARACTER*100 IPFILE,BUFFER
      CHARACTER*100 ANAME
      INTEGER I,J,K,N,IFORM,IPLANET,NGAS,IGAS,NX,NXTEMP,IX,ISPACE
      REAL TEMP
      REAL A,B,C,D,SVP,PP,LATITUDE
      REAL MOLWT,SCALE(MAXPRO)
      REAL XMAP(MAXV,MAXGAS+2+MAXCON,MAXPRO)

      INTEGER NVAR,VARIDENT(MVAR,3)
      REAL VARPARAM(MVAR,MPARAM)
      LOGICAL GASGIANT


C----------------------------------------------------------------------------
C
C     First read in reference ATMOSPHERIC profile

10    FORMAT(A)

      CALL FILE(IPFILE,IPFILE,'ref')
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
C     First skip header
54     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 54
       READ(BUFFER,*)IFORM
1      FORMAT(A)
       READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
       IF(NPRO.GT.MAXPRO)THEN
        PRINT*,'Error in subprofretg. NPRO>MAXPRO ',NPRO,MAXPRO
        STOP
       ENDIF

C      reset latitude to required input value. Will need when recomputing
C      scale heights
       LATITUDE = XLAT

       DO 20 I=1,NVMR
       READ(1,*)IDGAS(I),ISOGAS(I)
20     CONTINUE

C      reading the first block of profiles
       READ(1,*)
       N=MIN(NVMR,3)
C      N is the maximum VMR which can be read in from the next block
       DO 30 I=1,NPRO
       IF(IFORM.EQ.0)THEN
         READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,N)
        ELSE IF(IFORM.EQ.1)THEN
         READ(1,*)H(I),T(I),(VMR(I,J),J=1,N)
        ELSE
         CALL WTEXT('invalid format')
         STOP
         END IF
30     CONTINUE
C      reading in additional blocks if any
C      N VMR profiles have been read in so far
33     IF(NVMR.GT.N)THEN 
        READ(1,*)
C       profiles up to VMR(?,K) to be read from this block
        K=MIN(NVMR,(N+6))
        DO 32 I=1,NPRO
        READ(1,*)(VMR(I,J),J=N+1,K)
32      CONTINUE
        N=K
        GOTO 33
       END IF
      CLOSE(UNIT=1)

C     all processing below assumes that heights are in ascending order
C     so sorting just in case
      DO 12 I=1,NPRO-1
        IF(ABS(H(I)-H(I+1)).LT.0.01)THEN
 	 WRITE(*,14)
14       FORMAT(' identical height values found')
C	 do j=1,npro
C          print*,j,h(j),p(j),t(j)
C         enddo
	 STOP
        END IF
       IF(H(I).GT.H(I+1))THEN
	TEMP=H(I+1)
	H(I+1)=H(I)
	H(I)=TEMP
	TEMP=P(I+1)
	P(I+1)=P(I)
	P(I)=TEMP
	TEMP=T(I+1)
	T(I+1)=T(I)
	T(I)=TEMP
	DO 15 K=1,NVMR
 	 TEMP=VMR(I+1,K)
	 VMR(I+1,K)=VMR(I,K)
	 VMR(I,K)=TEMP
15      CONTINUE
       END IF
12    CONTINUE


C **************** Modify profile via hydrostatic equation ********
      JHYDRO=0
      DO I=1,NVAR
C       print*,I,VARIDENT(I,1)
       IF(VARIDENT(I,1).EQ.666)THEN
        JHYDRO=1
        HTAN = VARPARAM(I,1)
C        print*,'jpre,xn(jpre)',jpre,xn(jpre)
        PTAN = EXP(XN(JPRE))
       ENDIF
      ENDDO
C      print*,'JHYDRO',JHYDRO
      IF(JHYDRO.EQ.0)THEN
       CALL HYDROSTATH(NPRO,H,P,T,MOLWT,IPLANET,LATITUDE,SCALE)
      ELSE
c       CALL HYDROSTATP(NPRO,H,P,T,MOLWT,IPLANET,LATITUDE,
c     1  HTAN,PTAN,SCALE)
      ENDIF
	
C     Read in reference AEROSOL profile
      OPEN(UNIT=1,FILE='aerosol.ref',STATUS='OLD')
C     First skip header
55     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 55
       READ(BUFFER,*)NPRO1,NCONT
       IF (NPRO1.NE.NPRO)THEN
        print*,'subprofretg: NPRO <> NPRO1'
        stop
       ENDIF
       DO 31,I=1,NPRO
        READ(1,*)X,(CONT(J,I),J=1,NCONT)
31     CONTINUE
      CLOSE(1)

C     See if planet is a Giant Planet and we're working in wavenumbers. 
C        If so then then read in reference para-H2 fraction file
      FLAGH2P = 0
      IF(GASGIANT.AND.ISPACE.EQ.0.AND.
     1 (IPLANET.LE.8.OR.IPLANET.GE.12))THEN
       FLAGH2P=1
       OPEN(1,FILE='parah2.ref',STATUS='OLD')
C      First skip header
56     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 56
       READ(BUFFER,*)NPRO1
       IF(NPRO1.NE.NPRO)THEN
          PRINT*,'Error in subprofretg'
          PRINT*,'Para-H2 profile has wrong number of levels'
          STOP
       ENDIF
       DO I=1,NPRO
        READ(1,*)XH,PARAH2(I)
       ENDDO
       CLOSE(1)
      ENDIF


C     See if this is a scattering calculation. If so, read in the 
C     fractional cloud cover file.
      IF(ISCAT.EQ.1)THEN
       OPEN(1,FILE='fcloud.ref',STATUS='OLD')
C      First skip header
58     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 58
       READ(BUFFER,*)NPRO1,NCONT1
       IF(NPRO1.NE.NPRO)THEN
          PRINT*,'Error in subprofretg'
          PRINT*,'fcloud profile has wrong number of levels'
          STOP
       ENDIF
       DO I=1,NPRO
        READ(1,*)XH,FCLOUD(I),(ICLOUD(J,I),J=1,NCONT)
       ENDDO
       CLOSE(1)
      ENDIF

      DO IV=1,MAXV
       DO IP=1,NVMR+3+NCONT
        DO I=1,NPRO
         XMAP(IV,IP,I)=0.0
        ENDDO
       ENDDO
      ENDDO

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


      NXTEMP=0
      DO 1000 IVAR = 1,NVAR
C       PRINT*,'SUBPROFRETG: IVAR = ',IVAR
       JCONT=-1
       JVMR=-1
       IPAR=-1
C       print*,(varident(ivar,j),j=1,3)      
       IF(VARIDENT(IVAR,1).LE.100)THEN
        IF(VARIDENT(IVAR,1).EQ.0)THEN
C        variable is Temperature
C         PRINT*,'Temperature'
         DO I=1,NPRO
          XREF(I)=T(I)
         ENDDO
         IPAR = NVMR+1

        ELSEIF(VARIDENT(IVAR,1).LT.0)THEN
C        variable is aerosol amount
         JCONT = -VARIDENT(IVAR,1)
         IF(JCONT.GT.NCONT+2)THEN
          PRINT*,'Error in subprofretg, JCONT > NCONT+2',JCONT,NCONT+2
          STOP
         ENDIF
C        Note if JCONT = NCONT+1 then profile contains para-H2 fraction
C        Note if JCONT = NCONT+2 then profile contains fraction cloud  cover
         IF(JCONT.EQ.NCONT+1)THEN
C          print*,'para-H2'
CCCCCCCCCCCC*****************CCCCCCCCCCCCCCC********C*C*C*C*C*C*C*C*C*C*C*C*
c          IF(FLAGH2P.EQ.1)THEN
c           DO I=1,NPRO
c            XREF(I)=PARAH2(I)
c           ENDDO
c          ELSE
c           PRINT*,'Error in subprofretg, para-H2 fraction declared as'
c           PRINT*,'variable but atmosphere is not Giant Planet.'
c           STOP
c          ENDIF
         ELSEIF(JCONT.EQ.NCONT+2)THEN
C          print*,'fractional cloud cover'
          DO I=1,NPRO
            XREF(I)=FCLOUD(I)
          ENDDO
         ELSE
C          print*,'Aerosol : ',JCONT
          DO I=1,NPRO
           XREF(I)=CONT(JCONT,I)
          ENDDO
         ENDIF
         IPAR = NVMR+1+JCONT
        ELSE
C        Must be gas v.m.r., find which one
         DO IVMR=1,NVMR
          IF(VARIDENT(IVAR,1).EQ.IDGAS(IVMR).AND.
     1     VARIDENT(IVAR,2).EQ.ISOGAS(IVMR))JVMR = IVMR
         ENDDO
         IF(JVMR.LT.1)THEN
          PRINT*,'Subprofretg: Gas could not be found'
          STOP
         ENDIF
C         print*,'Gas : ',IDGAS(JVMR),ISOGAS(JVMR)
         DO I=1,NPRO
          XREF(I)=VMR(I,JVMR)
         ENDDO
         IPAR = JVMR

        ENDIF

C        print*,'VARIDENT : ',VARIDENT(IVAR,1),VARIDENT(IVAR,2),
C     1    VARIDENT(IVAR,3)

        IF(VARIDENT(IVAR,3).EQ.0)THEN
         NP = NPRO
         DO I=1,NPRO
           IF(VARIDENT(IVAR,1).EQ.0)THEN
            X1(I) = XN(NXTEMP+I)
            XMAP(NXTEMP+I,IPAR,I)=1.0
           ELSE
             IF(XN(NXTEMP+I).GT.-82.8931)THEN
               X1(I) = EXP(XN(NXTEMP+I))
             ELSE
               X1(I) = 1.0E-36
             ENDIF
             IF(XN(NXTEMP+I).LT.80.0)THEN
               X1(I) = EXP(XN(NXTEMP+I))
             ELSE
               X1(I) = EXP(80.0)
             ENDIF
             XMAP(NXTEMP+I,IPAR,I)=X1(I)
           ENDIF
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.1)THEN
         NP = 2
         PKNEE = VARPARAM(IVAR,1)
         PSHLD = VARPARAM(IVAR,2)
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
         ENDIF
         XFSH  = EXP(XN(NXTEMP+2))
         XFAC = (1.0-XFSH)/XFSH
C         DXFAC = -1.0/(XFSH*XFSH)
C        New gradient correction if fsh is held as logs
         DXFAC = -1.0/XFSH

         CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
         CALL VERINT(P,H,NPRO,HSHLD,PSHLD)
         JFSH = 0

         DO J=1,NPRO
          X1(J)=XDEEP
          IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,J)=1.0
          ELSE
            XMAP(NXTEMP+1,IPAR,J)=X1(J)
          ENDIF
          IF(P(J).LT.PKNEE.AND.P(J).GT.PSHLD)THEN  

             IF(JFSH.EQ.0)THEN
               DELH=H(J)-HKNEE
             ELSE
               DELH=H(J)-H(J-1)
             ENDIF
             
             X1(J)=X1(J-1)*EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1          X1(J-1)*EXP(-DELH*XFAC/SCALE(J)) +
     2          XMAP(NXTEMP+2,IPAR,J-1)*EXP(-DELH*XFAC/SCALE(J))

             JFSH = 1

             IF(X1(J).LT.1e-36)X1(J)=1e-36

          ELSEIF(P(J).LT.PSHLD)THEN

             X1(J)=X1(J-1)
             XMAP(NXTEMP+1,IPAR,J)=X1(J)
         
          ENDIF

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.2)THEN
         NP = 1
         DO J = 1,NPRO
          X1(J) = XREF(J)*XN(NXTEMP+1)
          XMAP(NXTEMP+1,IPAR,J)=XREF(J)
         ENDDO
	
        ELSEIF(VARIDENT(IVAR,3).EQ.3)THEN
         NP = 1
         DO J = 1,NPRO
          X1(J) = XREF(J)*EXP(XN(NXTEMP+1))
          XMAP(NXTEMP+1,IPAR,J)=X1(J)
         ENDDO
	
        ELSEIF(VARIDENT(IVAR,3).EQ.4)THEN
         NP = 3
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
         ENDIF

         XFSH  = EXP(XN(NXTEMP+2))
         XFAC = (1.0-XFSH)/XFSH
C         DXFAC = -1.0/(XFSH*XFSH)
C        New gradient correction if fsh is held as logs
         DXFAC = -1.0/XFSH
         PKNEE = EXP(XN(NXTEMP+3))

         CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
         JFSH = 0

         DO J=1,NPRO
          X1(J)=XDEEP
          IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,J)=1.0
          ELSE
            XMAP(NXTEMP+1,IPAR,J)=X1(J)
          ENDIF
          IF(P(J).LT.PKNEE)THEN  
             IF(JFSH.EQ.0)THEN
               DELH=H(J)-HKNEE
             ELSE
               DELH=H(J)-H(J-1)
             ENDIF
             X1(J)=X1(J-1)*EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1          X1(J-1)*EXP(-DELH*XFAC/SCALE(J)) +
     2          XMAP(NXTEMP+2,IPAR,J-1)*EXP(-DELH*XFAC/SCALE(J))

             IF(JFSH.EQ.0)THEN
               XMAP(NXTEMP+3,IPAR,J)=PKNEE*(XFAC/P(J))*
     1             X1(J-1)*EXP(-DELH*XFAC/SCALE(J))
             ENDIF
             XMAP(NXTEMP+3,IPAR,J)=XMAP(NXTEMP+3,IPAR,J)+
     1          XMAP(NXTEMP+2,IPAR,J-1)*EXP(-DELH*XFAC/SCALE(J))

             JFSH = 1

             IF(X1(J).LT.1e-36)X1(J)=1e-36

           END IF
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.6)THEN

         NP = 2
         HKNEE = VARPARAM(IVAR,1)
         IF(VARIDENT(IVAR,1).GE.0)THEN
          PRINT*,'Warning from SUBPROFRETG. You are using a Venusian'
          PRINT*,'cloud profile parameterisation for a non-cloud'
          PRINT*,'variable'         
          STOP 
         ENDIF

         XDEEP = EXP(XN(NXTEMP+1))
         XFAC  = EXP(XN(NXTEMP+2))

C         print*,HKNEE,XDEEP,XFAC

C        Calculate density of atmosphere (g/cm3)
         RHO = P(NPRO)*0.1013*MOLWT/(R*T(NPRO))
C         print*,P(NPRO),T(NPRO),MOLWT,R,RHO

C        Start ND(NPRO) at a random value. Will be rescaled anyway
         ND(NPRO)=1e-5
C        OD is in units of particles/cm2 = particles/cm3 x length(cm)
         OD(NPRO)=ND(NPRO)*SCALE(NPRO)*1E5
C        Q is specific density = particles/gram = particles/cm3 x g/cm3
         Q(NPRO)=ND(NPRO)/RHO         
         DNDH(NPRO)=0.0
         DQDH(NPRO)=0.0


         JFSH=-1
         DO J=NPRO-1,1,-1
          DELH = H(J+1)-H(J)
C         Calculate density of atmosphere (g/cm3)
C          print*,MOLWT,R,P(J),T(J),H(J)
          RHO = (0.1013*MOLWT/R)*(P(J)/T(J))
          ND(J)=ND(J+1)*EXP(DELH/XFAC)
          DNDH(J)=-ND(J)*DELH/(XFAC**2)+EXP(DELH/XFAC)*DNDH(J+1)
          OD(J)=OD(J+1)+(ND(J) - ND(J+1))*XFAC*1E5
          Q(J)=ND(J)/RHO
          DQDH(J) = DNDH(J)/RHO

C          print*,J,DELH,RHO,ND(J),DNDH(J),OD(J),Q(J),DQDH(J)

          IF(H(J).LE.HKNEE.AND.JFSH.LT.0)THEN
           F = (HKNEE-H(J))/DELH
           XOD = (1.-F)*OD(J) + F*OD(J+1)
C           PRINT*,'J,F,XOD = ',J,F,XOD 
           JFSH=1
          ENDIF
         ENDDO

         DO J=1,NPRO
          OD(J)=XDEEP*OD(J)/XOD
          ND(J)=XDEEP*ND(J)/XOD
          Q(J)=XDEEP*Q(J)/XOD

          IF(Q(J).GT.1e36)Q(J)=1e36
          IF(Q(J).LT.1e-36)Q(J)=1e-36

          X1(J)=Q(J)

          DNDH(J)=DNDH(J)*XDEEP/XOD
          DQDH(J)=DQDH(J)*XDEEP/XOD
C          DQDX(J)=Q(J)/XOD
          DQDX(J)=Q(J)

          XMAP(NXTEMP+1,IPAR,J)=DQDX(J)*XDEEP
          XMAP(NXTEMP+2,IPAR,J)=DQDH(J)*XFAC

         ENDDO

C         print*,NPRO,H(NPRO),ND(NPRO),OD(NPRO),Q(NPRO)
C         print*,XMAP(NXTEMP+1,IPAR,NPRO),XMAP(NXTEMP+2,IPAR,NPRO)
         DO J=NPRO-1,1,-1
          SUM=SUM+OD(J)
C          print*,J,H(J),ND(J),OD(J),Q(J)
C          print*,XMAP(NXTEMP+1,IPAR,J),XMAP(NXTEMP+2,IPAR,J)
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.7)THEN
         NP = 2
         PKNEE = VARPARAM(IVAR,1)
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
         ENDIF
         XFSH  = EXP(XN(NXTEMP+2))
         XFAC = (1.0-XFSH)/XFSH
C         DXFAC = -1.0/(XFSH*XFSH)
C        New gradient correction if fsh is held as logs
         DXFAC = -1.0/XFSH

         CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
         JFSH = 0       

         DELH=1000.0
         DO J=1,NPRO
          PMIN = ABS(P(J)-PKNEE)
          IF(PMIN.LT.DELH)THEN
           JFSH = J
           DELH=PMIN
          ENDIF
         ENDDO

         IF(JFSH.LT.2.OR.JFSH.GT.NPRO-1)THEN
          PRINT*,'SUBPROFRETG. Must choose pressure level'
          PRINT*,'within range of profile'
          STOP
         ENDIF

         print*,'Requested knee pressure : ',PKNEE
         print*,'Snapping to : ',P(JFSH)

         X1(JFSH)=XDEEP

         IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,JFSH)=1.0
         ELSE
            XMAP(NXTEMP+1,IPAR,JFSH)=X1(JFSH)
         ENDIF


         DO J=JFSH+1,NPRO 
             DELH = H(J)-H(J-1)         
             X1(J)=X1(J-1)*EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1          X1(J-1)*EXP(-DELH*XFAC/SCALE(J)) +
     2          XMAP(NXTEMP+2,IPAR,J-1)*EXP(-DELH*XFAC/SCALE(J))

             IF(X1(J).LT.1e-36)X1(J)=1e-36

         ENDDO

         DO J=JFSH-1,1,-1             
             DELH = H(J)-H(J+1)
             X1(J)=X1(J+1)*EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J+1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1          X1(J+1)*EXP(-DELH*XFAC/SCALE(J)) +
     2          XMAP(NXTEMP+2,IPAR,J+1)*EXP(-DELH*XFAC/SCALE(J))

             IF(X1(J).GT.1e36)X1(J)=1e36

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.8)THEN
         NP = 3
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
         ENDIF
         XFSH  = EXP(XN(NXTEMP+2))
         XFAC = (1.0-XFSH)/XFSH
C         DXFAC = -1.0/(XFSH*XFSH)
C        New gradient correction if fsh is held as logs
         DXFAC = -1.0/XFSH
         PKNEE = EXP(XN(NXTEMP+3))

         CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
         JFSH = 0
      
         DO J=1,NPRO
          X1(J)=1.0e-36
          IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,J)=1.0
          ELSE
            XMAP(NXTEMP+1,IPAR,J)=X1(J)
          ENDIF

          IF(P(J).LT.PKNEE)THEN  
             IF(JFSH.EQ.0)THEN
               DELH=H(J)-HKNEE
               XKEEP = XDEEP
             ELSE
               DELH=H(J)-H(J-1)
               XKEEP = X1(J-1)
             ENDIF

             X1(J)=XKEEP*EXP(-DELH*XFAC/SCALE(J))

             XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1          XKEEP*EXP(-DELH*XFAC/SCALE(J)) +
     2          XMAP(NXTEMP+2,IPAR,J-1)*EXP(-DELH*XFAC/SCALE(J))

             IF(JFSH.EQ.0)THEN
               XMAP(NXTEMP+3,IPAR,J)=PKNEE*(XFAC/P(J))*
     1             XKEEP*EXP(-DELH*XFAC/SCALE(J))
             ENDIF
             XMAP(NXTEMP+3,IPAR,J)=XMAP(NXTEMP+3,IPAR,J)+
     1          XMAP(NXTEMP+2,IPAR,J-1)*EXP(-DELH*XFAC/SCALE(J))

             JFSH = 1

             IF(X1(J).LT.1e-36)X1(J)=1e-36

           END IF

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.9)THEN

         NP = 3
         IF(VARIDENT(IVAR,1).GE.0)THEN
          PRINT*,'Warning from SUBPROFRETG. You are using a'
          PRINT*,'cloud profile parameterisation for a non-cloud'
          PRINT*,'variable'         
          STOP 
         ENDIF

         XDEEP = EXP(XN(NXTEMP+1))
         XFSH  = EXP(XN(NXTEMP+2))
         HKNEE = XN(NXTEMP+3)

C         print*,HKNEE,XDEEP,XFSH

C        Calculate density of atmosphere (g/cm3)
         RHO = P(NPRO)*0.1013*MOLWT/(R*T(NPRO))
C         print*,P(NPRO),T(NPRO),MOLWT,R,RHO

C        Start ND(NPRO) at a random value. Will be rescaled anyway
         ND(NPRO)=1e-35
C        OD is in units of particles/cm2 = particles/cm3 x length(cm)
         OD(NPRO)=ND(NPRO)*SCALE(NPRO)*1E5
C        Q is specific density = particles/gram = particles/cm3 x g/cm3
         Q(NPRO)=ND(NPRO)/RHO         
         DNDH(NPRO)=0.0
         DQDH(NPRO)=0.0

         
         JFSH=-1
         DO J=NPRO-1,1,-1
          DELH = H(J+1)-H(J)
          XFAC = SCALE(J)*XFSH
C         Calculate density of atmosphere (g/cm3)
C          print*,MOLWT,R,P(J),T(J),H(J)
          RHO = (0.1013*MOLWT/R)*(P(J)/T(J))
          ND(J)=ND(J+1)*EXP(DELH/XFAC)
          DNDH(J)=-ND(J)*DELH/(XFAC**2)+EXP(DELH/XFAC)*DNDH(J+1)
          OD(J)=OD(J+1)+(ND(J) - ND(J+1))*XFAC*1E5
          Q(J)=ND(J)/RHO
          DQDH(J) = DNDH(J)/RHO

C          print*,J,DELH,RHO,ND(J),DNDH(J),OD(J),Q(J),DQDH(J)

          IF(H(J).LE.HKNEE.AND.JFSH.LT.0)THEN
           F = (HKNEE-H(J))/DELH
           XOD = (1.-F)*OD(J) + F*OD(J+1)
C           PRINT*,'J,F,XOD = ',J,F,XOD 
           JFSH=1
          ENDIF
         ENDDO

         DO J=1,NPRO
          OD(J)=XDEEP*OD(J)/XOD
          ND(J)=XDEEP*ND(J)/XOD
          Q(J)=XDEEP*Q(J)/XOD
          IF(H(J).LE.HKNEE)Q(J) = 0.0

          IF(Q(J).GT.1e36)Q(J)=1e36
          IF(Q(J).LT.1e-36)Q(J)=1e-36

          X1(J)=Q(J)
C          print*,J,x1(J)
          DNDH(J)=DNDH(J)*XDEEP/XOD
          DQDH(J)=DQDH(J)*XDEEP/XOD
C          DQDX(J)=Q(J)/XOD
          DQDX(J)=Q(J)

          XMAP(NXTEMP+1,IPAR,J)=DQDX(J)*XDEEP
          XMAP(NXTEMP+2,IPAR,J)=DQDH(J)*XFAC

         ENDDO

C         print*,NPRO,H(NPRO),ND(NPRO),OD(NPRO),Q(NPRO)
C         print*,XMAP(NXTEMP+1,IPAR,NPRO),XMAP(NXTEMP+2,IPAR,NPRO)
         DO J=NPRO-1,1,-1
          SUM=SUM+OD(J)
C          print*,J,H(J),ND(J),OD(J),Q(J)
C          print*,XMAP(NXTEMP+1,IPAR,J),XMAP(NXTEMP+2,IPAR,J)
         ENDDO
        
        ENDIF

       ELSE

C       Must hold non-atmospheric parameter - find which.
        IF(VARIDENT(IVAR,1).EQ.999)THEN
C         print*,'Surface temperature'
C        no atmospheric mapping
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.888)THEN
C         print*,'Surface albedo spectrum'
         IPAR = -1
         NP = INT(VARPARAM(IVAR,1))
        ELSEIF(VARIDENT(IVAR,1).EQ.777)THEN
C         print*,'Tangent height correction'
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.666)THEN
C         print*,'Tangent pressure'
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.555)THEN
C         print*,'Radius of Planet'
         IPAR = -1
         NP = 1
c         DO J=1,NPRO
c          XMAP(NXTEMP+1,1,J)=EXP(XN(NXTEMP+1))
c	 ENDDO
c	print*,nxtemp+1,xn(nxtemp+1),'subprofretg'
        ELSE
         PRINT*,'SUBPROFRETG: VARTYPE NOT RECOGNISED'
         STOP
        ENDIF

       ENDIF
       
       IF(IPAR.GT.0)THEN
        IF(IPAR.LE.NVMR)THEN
         DO I=1,NPRO
          VMR(I,IPAR)=X1(I)
         ENDDO
        ELSEIF(IPAR.EQ.NVMR+1)THEN
         DO I=1,NPRO
          T(I)=X1(I)
         ENDDO
        ELSE
         IF(JCONT.LE.NCONT)THEN
          DO I=1,NPRO
           CONT(JCONT,I)=X1(I)
          ENDDO
         ELSEIF(JCONT.EQ.NCONT+2)THEN
          DO I=1,NPRO
           FCLOUD(I)=X1(I)
          ENDDO
         ELSE
          DO I=1,NPRO
           PARAH2(I)=X1(I)
          ENDDO
         ENDIF 
        ENDIF
       ENDIF
 
       NXTEMP = NXTEMP+NP

1000  CONTINUE

C     ********  Modify profile with hydrostatic equation ********
      IF(JHYDRO.EQ.0)THEN
       CALL HYDROSTATH(NPRO,H,P,T,MOLWT,IPLANET,LATITUDE,SCALE)
      ELSE
c       CALL HYDROSTATP(NPRO,H,P,T,MOLWT,IPLANET,LATITUDE,
c     1  HTAN,PTAN,SCALE)
      ENDIF


C     ********* Make sure nothing saturates *************
      DO 233 IGAS=1,NVMR

         IDAT=0
         DO I=1,NGAS
          IF(GASDATA(I,1).EQ.IDGAS(IGAS))THEN
            A = GASDATA(I,2)
            B = GASDATA(I,3)
            C = GASDATA(I,4)
            D = GASDATA(I,5)
            IDAT=1
          ENDIF
	 ENDDO

         IF(IDAT.EQ.1)THEN
 	   JSWITCH=0
           DO J=1,NPRO
            SVP=DPEXP(A+B/T(J)+C*T(J)+D*T(J)*T(J))
	    PP=VMR(J,IGAS)*P(J)

            IF(PP.GT.SVP)THEN
             IF(JSWITCH.EQ.0)THEN
              PRINT*,'Subprofretg: following gas predicted to condense'
              PRINT*,IDGAS(IGAS),ISOGAS(IGAS)
              PRINT*,'However, vmr curve left unchanged'
             ENDIF 
C             VMR(J,IGAS)=SVP/P(J)
             JSWITCH=1
            ENDIF

C            IF(JSWITCH.EQ.1.AND.VMR(J,IGAS).GT.VMR(J-1,IGAS))THEN
C             VMR(J,IGAS)=VMR(J-1,IGAS)
C            ENDIF
C            IF(JSWITCH.EQ.1)THEN
C             DO IX=1,NX
C              XMAP(IX,IGAS,J)=0.0
C             ENDDO
C            ENDIF

           ENDDO

         ENDIF

233   ENDDO
      
C     ************* Write out modified profiles *********

      IFORM=0
      CALL FILE(IPFILE,IPFILE,'prf')
      OPEN(UNIT=2,FILE=IPFILE,STATUS='UNKNOWN',ERR=52)
      WRITE(2,*)IFORM
      WRITE(2,501)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
501   FORMAT(1X,I3,F7.2,1X,I3,I3,F8.3)
      DO 503 I=1,NVMR
	WRITE(2,502)IDGAS(I),ISOGAS(I)
502     FORMAT(1X,I3,I5)
503   CONTINUE
C     writing the first block of profiles
      N=MIN(NVMR,3) 
      WRITE(2,504)(I,I=1,N)
504   FORMAT(1X,' height (km) ',' press (atm) ','  temp (K)   ',
     1  3(' VMR gas',I3,2X))
      DO 505 I=1,NPRO
        IF(IFORM.EQ.0)THEN
	  WRITE(2,506)H(I),P(I),T(I),(VMR(I,J),J=1,N)
506       FORMAT(1X,F13.3,E13.5,F13.4,3(E13.5))
         ELSE IF(IFORM.EQ.1)THEN
	  WRITE(2,511)H(I),T(I),(VMR(I,J),J=1,N)
511       FORMAT(1X,F13.3,13X,F13.4,3(E13.5))
          END IF
505   CONTINUE
C     writing additional blocks if any
C     N VMR profiles have been written so far
507   IF(NVMR.GT.N)THEN
C         profiles up to VMR(?,K) to be written to this block
	  K=MIN(NVMR,(N+6))
	  WRITE(2,508)(I,I=N+1,K)
508       FORMAT(1X,6(' VMR gas',I3,2X))
	  DO 509 I=1,NPRO
	  WRITE(2,510)(VMR(I,J),J=N+1,K)
510       FORMAT(1X,6E13.5)
509       CONTINUE
	  N=K
	  GOTO 507
	  END IF
C
52    CONTINUE
      CLOSE(2)


      OPEN(UNIT=2,FILE='aerosol.prf',STATUS='UNKNOWN')
C      print*,'subprofretg. Writing aerosol.prf'
      BUFFER='# aerosol.prf'   
      WRITE(2,10)BUFFER
      WRITE(2,*)NPRO, NCONT
      DO 41 I=1,NPRO
        WRITE(2,*) H(I),(CONT(J,I),J=1,NCONT)
C        print*,H(I),(CONT(J,I),J=1,NCONT)
41    CONTINUE
      CLOSE(2)

CCCCCCCCCCC********************CCCCCCCCCCCCC******CC***CC*CCCC*C
c      IF(FLAGH2P.EQ.1)THEN
c       OPEN(UNIT=2,FILE='parah2.prf',STATUS='UNKNOWN')
c        BUFFER = '# parah2.prf'
c        WRITE(2,*)NPRO
c        DO I=1,NPRO
c         WRITE(2,*)H(I),PARAH2(I)
c        ENDDO
c       CLOSE(2)
c      ENDIF

      IF(ISCAT.EQ.1)THEN
       OPEN(UNIT=2,FILE='fcloud.prf',STATUS='UNKNOWN')
        BUFFER = '# fcloud.prf'
        WRITE(2,*)NPRO,NCONT
        DO I=1,NPRO
         WRITE(2,*)H(I),FCLOUD(I),(ICLOUD(J,I),J=1,NCONT)
        ENDDO
       CLOSE(2)
      ENDIF

      RETURN

      END

