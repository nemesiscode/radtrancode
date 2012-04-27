      SUBROUTINE SUBPROFRETGPTX(IPFILE,ISPACE,ISCAT,GASGIANT,NVAR,
     1 VARIDENT,VARPARAM,NX,XN,JPRE,FLAGH2P)
C     $Id:
C     ***********************************************************************
C     Subroutine to modify an existing ipfile.prf T/P/vmr profile and 
C     aerosol.prf cloud density file according to the contents of the state
C     vector XN retrieved in a PREVIOUS run.  New profiles are written to
C     ipfile.prf and aerosol.prf files respectively.
C
C     Code also warns if vapour pressures exceeding SVP but does not
C     reset vmrs since the temperatures could also be in error.
C 
C     Code recalulates the heights from the hydrostatic equation.
C
C     Input variables
C	IPFILE	CHARACTER*100	Root run name
C	ISPACE	INTEGER		Wavenumber/Wavelength space
C	ISCAT	INTEGER		Scattering flag
C	GASGIANT LOGICAL	Indicates if planet is a Gas Giant
C	NVAR	INTEGER		Number of variable profiles
C	VARIDENT(MVAR,3) INTEGER Identity of profiles and parameterisation
C					scheme
C	VARPARAM(MVAR,MPARAM) REAL	Additional parameterisation
C	NX	INTEGER		Number if elements in state vector
C	XN(MX)	REAL		State vector
C	JPRE	INTEGER		Location of tangent pressure (if defined)
C     Output variables
C	None
C     Pat Irwin	29/7/96		Original
C     Pat Irwin 17/10/03	Revised for Nemesis
C
C     ***********************************************************************
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'

      REAL XN(MX),DPEXP,DELH,XFAC,DXFAC,XTMP,PMIN
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS)
      REAL CONT(MAXGAS,MAXPRO),X,XREF(MAXPRO),X1(MAXPRO)
      REAL PKNEE,HKNEE,XDEEP,XFSH,PARAH2(MAXPRO),XH,XKEEP
      REAL FCLOUD(MAXPRO)
      INTEGER ICLOUD(MAXGAS,MAXPRO),ISCAT,NCONT1
      INTEGER NPRO,NPRO1,NVMR,JZERO,IV,IP,IVAR,JCONT,JVMR
      INTEGER IDGAS(MAXGAS),ISOGAS(MAXGAS),IPAR,IVMR,NP
      INTEGER IDAT,FLAGH2P,NCONT,JFSH,JPRE,JHYDRO
      REAL HTAN,PTAN
      CHARACTER*8 PNAME   
      INTEGER JSWITCH,ISPACE
      REAL GASDATA(20,5)
      CHARACTER*100 IPFILE,BUFFER
      CHARACTER*100 ANAME
      INTEGER I,J,K,N,IFORM,IPLANET,NGAS,IGAS,NX,NXTEMP,IX
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

      CALL FILE(IPFILE,IPFILE,'prf')
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
C     First skip header
54     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 54
       READ(BUFFER,*)IFORM
1      FORMAT(A)
       READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
       IF(NPRO.GT.MAXPRO)THEN
        PRINT*,'Error in subprofretgx. NPRO>MAXPRO ',NPRO,MAXPRO
        STOP
       ENDIF

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
	 do j=1,npro
          print*,j,h(j),p(j),t(j)
         enddo
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
       IF(VARIDENT(I,1).EQ.666)THEN
        JHYDRO=1
        HTAN = VARPARAM(I,1)
        PTAN = EXP(XN(JPRE))
       ENDIF
      ENDDO
      IF(JHYDRO.EQ.0)THEN
       CALL HYDROSTATH(NPRO,H,P,T,MOLWT,IPLANET,LATITUDE,SCALE)
      ELSE
c       CALL HYDROSTATP(NPRO,H,P,T,MOLWT,IPLANET,LATITUDE,
c    1   HTAN,PTAN,SCALE)
      ENDIF

C     Read in reference AEROSOL profile
      OPEN(UNIT=1,FILE='aerosol.prf',STATUS='OLD')
C     First skip header
55     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 55
       READ(BUFFER,*)NPRO1,NCONT
       IF (NPRO1.NE.NPRO)THEN
        print*,'subprofretgx: NPRO <> NPRO1'
        stop
       ENDIF
       DO 31,I=1,NPRO
        READ(1,*)X,(CONT(J,I),J=1,NCONT)
31     CONTINUE
      CLOSE(1)

C     See if we need to read in a para-H2 fraction file
      IF(FLAGH2P.EQ.1)THEN
       OPEN(1,FILE='parah2.prf',STATUS='OLD')
C      First skip header
56     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 56
       READ(BUFFER,*)NPRO1
       IF(NPRO1.NE.NPRO)THEN
          PRINT*,'Error in subprofretgx'
          PRINT*,'Para-H2 profile has wrong number of levels',NPRO1,NPRO
          STOP
       ENDIF
       DO I=1,NPRO
        READ(1,*)XH,PARAH2(I)
       ENDDO
       CLOSE(1)
      ENDIF


      IF(ISCAT.EQ.1)THEN
       OPEN(1,FILE='fcloud.prf',STATUS='OLD')
C      First skip header
58     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 58
       READ(BUFFER,*)NPRO1,NCONT1
       IF(NPRO1.NE.NPRO)THEN
          PRINT*,'Error in subprofretgx'
          PRINT*,'fcloud profile has wrong number of levels'
          STOP
       ENDIF
       DO I=1,NPRO
        READ(1,*)XH,FCLOUD(I),(ICLOUD(J,I),J=1,NCONT)
       ENDDO
       CLOSE(1)
      ENDIF

      DO IV=1,MAXV
       DO IP=1,NVMR+2+NCONT
        DO I=1,NPRO
         XMAP(IV,IP,I)=0.0
        ENDDO
       ENDDO
      ENDDO

      ANAME='SVP.dat'
      CALL DATARCHIVE(ANAME)
      OPEN(13,FILE=ANAME,STATUS='OLD')
      WRITE(*,*)' '
      WRITE(*,*)'SUBPROFRETGX: reading saturated-vapour-pressure'
      WRITE(*,*)'  data from ',aname
C     First skip header
57    READ(13,1)BUFFER
      IF(BUFFER(1:1).EQ.'#')GOTO 57
      READ(BUFFER,*)NGAS
      DO 21 I=1,NGAS
       READ(13,*)(GASDATA(I,J),J=1,5)
21    CONTINUE
      CLOSE(13)

c      OPEN(13,FILE='/home/jupiter/plan2/oxpln5/radtran/source/'//
c     1 'release2/data/SVP.DAT',
c     1 STATUS='OLD')
c      READ(13,*)NGAS
c      DO 21 I=1,NGAS
c       READ(13,*)(GASDATA(I,J),J=1,5)
c21    CONTINUE
c      CLOSE(13)

      NXTEMP=0
      DO 1000 IVAR = 1,NVAR

C       PRINT*,'SUBPROFRETGX: IVAR = ',IVAR
       JCONT=-1
       JVMR=-1
       IPAR=-1
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
C         print*,'JCONT = ',JCONT
         IF(JCONT.GT.NCONT+2)THEN
          PRINT*,'Error in subprofretgx, JCONT > NCONT+2',JCONT,NCONT+2
          STOP
         ENDIF
C        Note if JCONT = NCONT+1 then profile contains para-H2 fraction
C        Note if JCONT = NCONT+2 then profile contains fractional cloud cover
         IF(JCONT.EQ.NCONT+1)THEN
C          print*,'para-H2'
          IF(FLAGH2P.EQ.1)THEN
           DO I=1,NPRO
            XREF(I)=PARAH2(I)
           ENDDO
          ELSE
           PRINT*,'Error in subprofretgx'
           PRINT*,'Para-H2 fraction declared as variable but'
           PRINT*,'Atmosphere is not Giant Planet'
           STOP
          ENDIF
         ELSEIF(JCONT.EQ.NCONT+2)THEN
C          print*,'fcloud'
           DO I=1,NPRO
            XREF(I)=FCLOUD(I)
           ENDDO
         ELSE
C          print*,'Aerosol : ',JCONT
          DO I=1,NPRO
           XREF(I)=CONT(JCONT,I)
C           print*,i,xref(i)
          ENDDO
         ENDIF
         IPAR = NVMR+1+JCONT
        ELSE
C        Must hold gas v.m.r. - find which one.
         DO IVMR=1,NVMR
          IF(VARIDENT(IVAR,1).EQ.IDGAS(IVMR).AND.
     1     VARIDENT(IVAR,2).EQ.ISOGAS(IVMR))JVMR = IVMR
         ENDDO
         IF(JVMR.LT.1)THEN
          PRINT*,'Subprofretgx: Gas could not be found'
          STOP
         ENDIF
C         print*,'Gas : ',IDGAS(JVMR),ISOGAS(JVMR)
         DO I=1,NPRO
          XREF(I)=VMR(I,JVMR)
         ENDDO
         IPAR = JVMR

        ENDIF

C        print*,'Subprofretgx: NVMR,NCONT,IPAR = ',NVMR,NCONT,IPAR

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
             IF(XN(NXTEMP+I).LT.20.0)THEN
               X1(I) = EXP(XN(NXTEMP+I))
             ELSE
               X1(I) = EXP(20.0)
             ENDIF
             XMAP(NXTEMP+I,IPAR,I)=X1(I)
           ENDIF
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.1)THEN
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
         DXFAC = -1.0/XFSH


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

             JFSH = 1

          END IF
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.2)THEN
         NP = 1
         DO J=1,NPRO
          X1(J) = XREF(J)*XN(NXTEMP+1)
          XMAP(NXTEMP+1,IPAR,J)=XREF(J)
         ENDDO
	
        ELSEIF(VARIDENT(IVAR,3).EQ.3)THEN
         NP = 1
         DO J=1,NPRO
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

          END IF
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.6)THEN

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

         DO J=1,NPRO
          IF(P(J).LT.PKNEE)THEN
           JFSH=J
           GOTO 201
          ENDIF
         ENDDO
201      CONTINUE

         IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,JFSH)=1.0
         ELSE
            XMAP(NXTEMP+1,IPAR,JFSH)=XDEEP
         ENDIF
         X1(JFSH)=XDEEP

         DO J=JFSH+1,NPRO
          DELH = H(J)-H(J-1)
          XTMP = X1(J-1)
          X1(J)=XTMP*EXP(-DELH*XFAC/SCALE(J))
          XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
          XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1          X1(J-1)*EXP(-DELH*XFAC/SCALE(J)) +
     2          XMAP(NXTEMP+2,IPAR,J-1)*EXP(-DELH*XFAC/SCALE(J))

          IF(X1(J).LT.1e-36)X1(J)=1e-36

         ENDDO

         DO J=JFSH-1,1,-1
          DELH = H(J)-H(J+1)
          XTMP = X1(J+1)
          X1(J)=XTMP*EXP(-DELH*XFAC/SCALE(J))
          XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J+1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
           XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1          X1(J+1)*EXP(-DELH*XFAC/SCALE(J)) +
     2          XMAP(NXTEMP+2,IPAR,J+1)*EXP(-DELH*XFAC/SCALE(J))

          IF(X1(J).LT.1e-36)X1(J)=1e-36

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
          PRINT*,'SUBPROFRETGX. Must choose pressure level'
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
         DXFAC = -1.0/XFSH
         PKNEE = EXP(XN(NXTEMP+3))

         CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
         JFSH = 0

         DO J=1,NPRO
          X1(J)=1e-36
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

          END IF
         ENDDO

        ELSE
         print*,'vartype profile parametrisation not recognised'
         stop
        ENDIF

       ELSE

         IF(VARIDENT(IVAR,1).EQ.999)THEN
C         Surface temperature element. Hence no atmospheric mapping
C          print*,'Surface temperature'
          NP = 1
          IPAR = -1
         ELSEIF(VARIDENT(IVAR,1).EQ.888)THEN
C         Surface albedoelement. Hence no atmospheric mapping
C          print*,'Surface albedo'
          NP = INT(VARPARAM(IVAR,1))
          IPAR = -1
         ELSEIF(VARIDENT(IVAR,1).EQ.777)THEN
C         Tangent height correc. Hence no atmospheric mapping
C          print*,'Tangent height correction'
          NP = 1
          IPAR = -1
         ELSEIF(VARIDENT(IVAR,1).EQ.666)THEN
C         Tangent pressure correc. Hence no atmospheric mapping
C          print*,'Tangent pressure correction'
          NP = 1
          IPAR = -1
        ELSEIF(VARIDENT(IVAR,1).EQ.555)THEN
C         Radius of planet correc. Hence no atmospheric mapping
          print*,'Radius of Planet'
          NP = 1
          IPAR = -1
         ELSE
          PRINT*,'SUBPROFRETGX: VARTYPE NOT RECOGNISED : ',
     1      VARIDENT(IVAR,1)
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
C           print*,jcont,I,X1(I)
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
              PRINT*,'Subprofretgx, following gas predicted to condense'
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

      CALL FILE(IPFILE,IPFILE,'prf')
      OPEN(UNIT=2,FILE=IPFILE,STATUS='UNKNOWN',ERR=52)
      IFORM=0
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
C      print*,'Subprofretgx. Writing aerosol.prf'
      BUFFER='# aerosol.prf'   
      WRITE(2,10)BUFFER
      WRITE(2,*)NPRO, NCONT
      DO 41 I=1,NPRO
        WRITE(2,*) H(I),(CONT(J,I),J=1,NCONT)
C        print*,H(I),(CONT(J,I),J=1,NCONT)
41    CONTINUE
      CLOSE(2)

      IF(FLAGH2P.EQ.1)THEN
       OPEN(UNIT=2,FILE='parah2.prf',STATUS='UNKNOWN')
        BUFFER = '# parah2.prf'
        WRITE(2,*)NPRO
        DO I=1,NPRO
         WRITE(2,*)H(I),PARAH2(I)
        ENDDO
       CLOSE(2)
      ENDIF

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
