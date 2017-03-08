      SUBROUTINE SUBPROFRETG(XFLAG,IPFILE,ISPACE,ISCAT,GASGIANT,XLAT,
     1 NVAR,VARIDENT,VARPARAM,NX,XN,JPRE,NCONT,FLAGH2P,XMAP,IERR)
C     $Id:
C     ***********************************************************************
C     Subroutine to modify an existing ipfile.ref T/P/vmr profile and 
C     aerosol.ref cloud density file according to the contents of the state
C     vector XN. New profiles are written to ipfile.prf and aerosol.prf 
C     files respectively.
C
C     If XFLAG=1, then it is assumed that the .ref files have already been
C     updated to .prf files for the variables currently being retrieved,
C     but need to be further updated for variables retrieved in a previous
C     Nemesis run, i.e. LIN=1 or 3
C
C     Code looks to see if an ipfile.vpf file is present and if so reads in
C     a list of which gas abundances are to be limited to a specified fraction
C     of the SVP, together with whether this limit applies at all altiitudes
C     or only in the troposphere (internal source and tropopause cold trap) or 
C     only in the stratosphere (external source and tropopause cold trap).
C 
C     Code also shifts .ref file to the required latitude and recalulates
C     the heights from the hydrostatic equation.
C
C     Routine also returns the XMAP matrix which relates the functional
C     derivatives calculated by CIRSRADG with the elements of the state
C     vector  
C     
C     Input variables
C	XFLAG	INTEGER		Updates .ref file if XFLAG=0, 
C				updates .prf file if XFLAG=1.
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
C	JPRE	INTEGER		Level of tanget pressure for limb pressure 
C				retrieval.
C	FLAGH2P INTEGER		Set to 1 if para-H2 profile is variable
C     Output variables
C	NCONT	INTEGER		Number of cloud particle types
C	XMAP(MAXV,MAXGAS+2+MAXCON,MAXPRO) REAL Matrix relating functional
C				derivatives calculated by CIRSRADG to the 
C				elements of the state vector. 
C				Elements of XMAP are the rate of change of
C				the profile vectors (i.e. temperature, vmr prf
C				files) with respect to the change in the state
C				vector elements. So if X1(J) is the modified 
C				temperature,vmr,clouds at level J to be 
C				written out to runname.prf or aerosol.prf then
C				XMAP(K,L,J) is d(X1(J))/d(XN(K)) and where
C				 L is the identifier (1 to NGAS+1+2*NCONT)
C	IERR	INTEGER		Error reporting flag set to 1 if a vmr has
C				gone negative
C
C     Pat Irwin	29/7/96		Original
C     Pat Irwin 17/10/03	Revised for Nemesis
C     Pat Irwin 9/5/12		Updated for new Nemesis
C
C     ***********************************************************************
      IMPLICIT NONE
      INCLUDE '../radtran/includes/arrdef.f'
      INCLUDE '../radtran/includes/constdef.f'
      INCLUDE '../radtran/includes/emcee.f'
      INCLUDE 'arraylen.f'

      REAL XN(MX),DPEXP,DELH,XFAC,DXFAC,XTMP,SUM,PMIN
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS)
      REAL CONT(MAXCON,MAXPRO),XLAT,X,XREF(MAXPRO),X1(MAXPRO)
      REAL PKNEE,HKNEE,XDEEP,XFSH,PARAH2(MAXPRO),XH,XKEEP,X2(MAXPRO)
      REAL RHO,F,DQDX(MAXPRO),DX,PLIM,XFACP,CWID,PTROP, NEWF
      REAL DNDH(MAXPRO),DQDH(MAXPRO),FCLOUD(MAXPRO)
      DOUBLE PRECISION Q(MAXPRO),OD(MAXPRO),ND(MAXPRO),XOD
      INTEGER ISCALE(MAXGAS),XFLAG,NPVAR,MAXLAT,IERR
      INTEGER NLATREF,ILATREF,JLAT,KLAT,ICUT,JX
      PARAMETER (MAXLAT=20)
      REAL HREF(MAXLAT,MAXPRO),TREF(MAXLAT,MAXPRO),FLAT
      REAL PREF(MAXLAT,MAXPRO),VMRREF(MAXLAT,MAXPRO,MAXGAS)
      REAL LATREF(MAXLAT),MOLWTREF(MAXLAT),XFSHREF
      REAL XRH,XCDEEP,P1,PS,PS1,PH,Y1,Y2,YY1,YY2
      REAL XCH4,PCH4,PCUT,GETRADIUS,RPARTICLE
      INTEGER ICLOUD(MAXCON,MAXPRO),NCONT1,JSPEC,IFLA,I1
      INTEGER NPRO,NPRO1,NVMR,JZERO,IV,IP,IVAR,JCONT,JVMR
      INTEGER IDGAS(MAXGAS),ISOGAS(MAXGAS),IPAR,JPAR,IVMR,NP
      INTEGER IDAT,NCONT,FLAGH2P,JFSH,JPRE,JHYDRO,ISCAT,ICOND
      REAL HTAN,PTAN,R,REFRADIUS
      REAL GASDATA(20,5),XVMR(MAXGAS),XMOLWT,XXMOLWT(MAXPRO)
      REAL CALCMOLWT
      INTEGER JSWITCH,ITEST,ICL
      CHARACTER*100 IPFILE,BUFFER
      CHARACTER*100 ANAME
      INTEGER I,J,K,N,AMFORM,IPLANET,NGAS,IGAS,NX,NXTEMP,IX,ISPACE
      REAL TEMP
      REAL A,B,C,D,SVP,PP,LATITUDE
      REAL MOLWT,SCALE(MAXPRO),XMAP1,SCALEH
      REAL XMAP(MAXV,MAXGAS+2+MAXCON,MAXPRO)

      INTEGER NVAR,VARIDENT(MVAR,3)
      REAL VARPARAM(MVAR,MPARAM)
      REAL XWID,Y,Y0,XX
      LOGICAL GASGIANT,FEXIST,VPEXIST,NTEST,ISNAN
      REAL VP(MAXGAS),VP1
      INTEGER SVPFLAG(MAXGAS),SVPFLAG1
      INTEGER NVP,ISWITCH(MAXGAS),IP1,IP2,JKNEE
      REAL XLDEEP,XLHIGH
      COMMON /SROM223/PCUT

C----------------------------------------------------------------------------
C
C     First zero-fill HREF, PREF, TREF and VMRREF arrays
      do i=1,MAXLAT
       do j=1,MAXPRO
        HREF(i,j)=0.
        PREF(i,j)=0.
        TREF(i,j)=0.
        do k=1,MAXGAS
         VMRREF(i,j,k)=0.
        enddo
       enddo
      enddo

C     First read in reference ATMOSPHERIC profile

10    FORMAT(A)

C     RGAS now read in from constdef.f, so set R accordingly and in correct
C     units
      R=RGAS*0.001

C     If XFLAG = 0, then we're generating a .prf file from the .ref file
C     If XFLAG = 1, then we further updating the .prf file with previously
C                   retrieved parameters (LIN=1 or 3)
      IF(XFLAG.EQ.0)THEN
       CALL FILE(IPFILE,IPFILE,'ref')
      ELSE
       CALL FILE(IPFILE,IPFILE,'prf')
      ENDIF
      OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
      MOLWT=-1.
C     First skip header
54     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 54
       READ(BUFFER,*)AMFORM
1      FORMAT(A)
       IF(XFLAG.EQ.0)THEN
         READ(1,*)NLATREF
         print*,'NLATREF = ',NLATREF
         IF(NLATREF.GT.MAXLAT)THEN
          PRINT*,'MAXLAT in subprofretg is too small. Increase'
          PRINT*,'and recompile'
          PRINT*,'MAXLAT, NLATREF = ',MAXLAT,NLATREF
          STOP
         ENDIF

         DO 601 ILATREF=1,NLATREF
          IF(AMFORM.EQ.1)THEN
           READ(1,*)IPLANET,LATREF(ILATREF),NPRO,NVMR
          ELSE
           READ(1,*)IPLANET,LATREF(ILATREF),NPRO,NVMR,
     1      MOLWTREF(ILATREF)
          ENDIF

          IF(NPRO.GT.MAXPRO)THEN
           PRINT*,'Error in subprofretg. NPRO>MAXPRO ',NPRO,MAXPRO
           STOP
          ENDIF

          DO 23 I=1,NVMR
           READ(1,*)IDGAS(I),ISOGAS(I)
           ISCALE(I)=1
23        CONTINUE

C         Skip header
          READ(1,*)
          DO 33 I=1,NPRO
            READ(1,*)HREF(ILATREF,I),PREF(ILATREF,I),
     & TREF(ILATREF,I),(VMRREF(ILATREF,I,J),J=1,NVMR)
            IF(VMRflag.eq.1)THEN
             HREF(ILATREF,I)=MCMCheight(I)
             TREF(ILATREF,I)=MCMCtemp(I)
             DO J=1, NVMR
              VMRREF(ILATREF,I,J)=MCMCvmr(J)
             ENDDO
            ENDIF
33        CONTINUE
601      CONTINUE

         

C        Now interpolate to correct latitude
         IF(NLATREF.EQ.1)THEN
          JLAT=1
          FLAT=0.
          PRINT*,'Snapping to first latitude'
         ELSE
          KLAT=-1
          DO ILATREF=1,NLATREF
           IF(XLAT.GE.LATREF(ILATREF))KLAT=ILATREF
          ENDDO
          IF(KLAT.LT.1)THEN
           PRINT*,'Requested latitude is less than range given'
           PRINT*,'Using lowest latitude available'
           PRINT*,'Requested : ',XLAT
           PRINT*,'Lowest available : ',LATREF(1)
           JLAT=1
           FLAT=0.
          ELSEIF(KLAT.EQ.NLATREF)THEN
           PRINT*,'Requested latitude is greater than range given'
           PRINT*,'Using highest latitude available'
           PRINT*,'Requested : ',XLAT
           PRINT*,'Highest available : ',LATREF(NLATREF)
           JLAT=NLATREF-1
           FLAT=1.0
          ELSE
           JLAT=KLAT
           FLAT=(XLAT-LATREF(JLAT))/
     &			(LATREF(JLAT+1)-LATREF(JLAT))
           PRINT*,'JLAT,FLAT',JLAT,FLAT
           PRINT*,'LATREF(JLAT),LATREF(JLAT+1)',
     &		LATREF(JLAT),LATREF(JLAT+1)
          ENDIF
         ENDIF

         MOLWT=(1.0-FLAT)*MOLWTREF(JLAT)+FLAT*MOLWTREF(JLAT+1)
         DO I=1,NPRO
          H(I)=(1.0-FLAT)*HREF(JLAT,I)+FLAT*HREF(JLAT+1,I)
          P(I)=(1.0-FLAT)*PREF(JLAT,I)+FLAT*PREF(JLAT+1,I)
          T(I)=(1.0-FLAT)*TREF(JLAT,I)+FLAT*TREF(JLAT+1,I)
          DO J=1,NVMR
           VMR(I,J)=(1.0-FLAT)*VMRREF(JLAT,I,J)+
     &		FLAT*VMRREF(JLAT+1,I,J)
          ENDDO
         ENDDO
 
         LATITUDE=XLAT

       ELSE

         IF(AMFORM.EQ.1)THEN
          READ(1,*)IPLANET,LATITUDE,NPRO,NVMR
         ELSE
          READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
         ENDIF
         print*,IPLANET,LATITUDE,NPRO,NVMR,MOLWT
         IF(NPRO.GT.MAXPRO)THEN
          PRINT*,'Error in subprofretg. NPRO>MAXPRO ',NPRO,MAXPRO
          STOP
         ENDIF

C        reset latitude to required input value. Will need when recomputing
C        scale heights
         LATITUDE = XLAT

         DO 20 I=1,NVMR
          READ(1,*)IDGAS(I),ISOGAS(I)
          ISCALE(I)=1
20       CONTINUE

C        Skip header
         READ(1,*)
         DO 30 I=1,NPRO
           READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,NVMR)
            IF(VMRflag.eq.1)THEN
             H(I)=MCMCheight(I)
             T(I)=MCMCtemp(I)
             DO J=1, NVMR
              VMR(I,J)=MCMCvmr(J)
             ENDDO
            ENDIF
30       CONTINUE

        CLOSE(UNIT=1)

       ENDIF


C      Make sure that vmrs add up to 1 if AMFORM=1
       IF(AMFORM.EQ.1)THEN
        CALL ADJUSTVMR(NPRO,NVMR,VMR,ISCALE,IERR)
        
        DO 301 I=1,NPRO
         DO K=1,NVMR
          XVMR(K)=VMR(I,K)
c          print*,I,K,XVMR(K)
         ENDDO
         XXMOLWT(I)=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
301     CONTINUE
       ENDIF


C     all processing below assumes that heights are in ascending order
C     so sorting just in case
      DO 12 I=1,NPRO-1
        IF(ABS(H(I)-H(I+1)).LT.0.01)THEN
         WRITE(*,14)
14       FORMAT(' identical height values found')
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
       CALL XHYDROSTATH(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1  IDGAS,ISOGAS,H,P,T,VMR,SCALE)
      ELSE
       CALL XHYDROSTATP(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1  IDGAS,ISOGAS,H,P,T,VMR,HTAN,PTAN,SCALE)
      ENDIF

C     Read in reference AEROSOL profile
      IF(XFLAG.EQ.0)THEN
       OPEN(UNIT=1,FILE='aerosol.ref',STATUS='OLD')
      ELSE
       OPEN(UNIT=1,FILE='aerosol.prf',STATUS='OLD')
      ENDIF
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

C     See if we need to deal with a variable para-H2 fraction
      IF(FLAGH2P.EQ.1)THEN

       inquire(file='parah2.ref',exist=fexist)
       if(.not.fexist) then
        print*,'Error in subprofretg.f - you have a variable para-H2'      
        print*,'para-H2 CIA file, but no parah2.ref file'
        stop
       endif

       IF(XFLAG.EQ.0)THEN
        OPEN(1,FILE='parah2.ref',STATUS='OLD')
       ELSE
        OPEN(1,FILE='parah2.prf',STATUS='OLD')
       ENDIF
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
      IF(ISCAT.GT.0)THEN
       IF(XFLAG.EQ.0)THEN
        OPEN(1,FILE='fcloud.ref',STATUS='OLD')
       ELSE
        OPEN(1,FILE='fcloud.prf',STATUS='OLD')
       ENDIF
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
       JSPEC=-1
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
          IF(FLAGH2P.EQ.1)THEN
           DO I=1,NPRO
            XREF(I)=PARAH2(I)
           ENDDO
          ELSE
           PRINT*,'Error in subprofretg, para-H2 fraction declared as'
           PRINT*,'variable but atmosphere is not Giant Planet.'
           STOP
          ENDIF
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
          PRINT*,'Subprofretg: Gas could not be found',
     1     VARIDENT(IVAR,1),VARIDENT(IVAR,2)
          STOP
         ENDIF
C         print*,'Gas : ',IDGAS(JVMR),ISOGAS(JVMR)
C        Set ISCALE=0 for this gas to prevent vmr being scaled to give a 
C        total sum or vmrs=1 for AMFORM=1 format profile
         ISCALE(JVMR)=0
         DO I=1,NPRO
          XREF(I)=VMR(I,JVMR)
         ENDDO
         IPAR = JVMR

        ENDIF

C        print*,'VARIDENT : ',VARIDENT(IVAR,1),VARIDENT(IVAR,2),
C     1    VARIDENT(IVAR,3)

C       Look up number of parameters needed to define this type of profile
        NP = NPVAR(VARIDENT(IVAR,3),NPRO)
C        print*,'IVAR,VARIDENT,NP',IVAR,(VARIDENT(IVAR,I),I=1,3),NP
        IF(VARIDENT(IVAR,3).EQ.0)THEN
C        Model 0. Continuous profile
C        ***************************************************************
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
C        Model 1. Variable deep abundance, fixed knee pressure and variable
C        fractional scale height.
C        ***************************************************************
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

             IF(X1(J).LT.1e-36)X1(J)=1e-36

          END IF
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.2)THEN
C        Model 2. Profile is scaled fraction of reference profile
C        ***************************************************************
         DO J = 1,NPRO
          X1(J) = XREF(J)*XN(NXTEMP+1)
          XMAP(NXTEMP+1,IPAR,J)=XREF(J)
         ENDDO
	
        ELSEIF(VARIDENT(IVAR,3).EQ.3)THEN
C        Model 3. Profile is scaled fraction of reference profile, but 
C         code uses log scaling, which is more robust
C        ***************************************************************
         DO J = 1,NPRO
          IF(XREF(J).GT.(1.0E37/EXP(XN(NXTEMP+1))).AND.XREF(J).GT.0.0)
     1   	   XN(NXTEMP+1)=LOG(1.0E37/XREF(J))
          X1(J) = XREF(J)*EXP(XN(NXTEMP+1))
          XMAP(NXTEMP+1,IPAR,J)=X1(J)
         ENDDO
	
        ELSEIF(VARIDENT(IVAR,3).EQ.4)THEN
C        Model 4. Variable deep abundance, variable knee pressure and variable
C        fractional scale height.
C        ***************************************************************
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

        ELSEIF(VARIDENT(IVAR,3).EQ.5)THEN
C        Model 5. Continuous profile, but variable with latitude.
C        Option is defunct.
C        ***************************************************************

         PRINT*,'VARIDENT(IVAR,3).EQ.5 - Model no longer supported'
         STOP

        ELSEIF(VARIDENT(IVAR,3).EQ.6)THEN
C        Model 6. Venus-type cloud profile represented by a fixed base height, 
C        variable integrated optical depth above that height and variable cloud
C        scale height (km)
C        ***************************************************************
         HKNEE = VARPARAM(IVAR,1)
         IF(VARIDENT(IVAR,1).GE.0)THEN
          PRINT*,'Warning from SUBPROFRETG. You are using a Venusian'
          PRINT*,'cloud profile parameterisation for a non-cloud'
          PRINT*,'variable'         
          STOP 
         ENDIF

         XDEEP = EXP(XN(NXTEMP+1))
         XFAC  = EXP(XN(NXTEMP+2))

C        Calculate density of atmosphere (g/cm3)
         IF(AMFORM.EQ.0)THEN
          XMOLWT=MOLWT
         ELSE
          XMOLWT=XXMOLWT(NPRO)
         ENDIF
         RHO = P(NPRO)*0.1013*XMOLWT/(R*T(NPRO))

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

          IF(AMFORM.EQ.0)THEN
           XMOLWT=MOLWT
          ELSE
           XMOLWT=XXMOLWT(J)
          ENDIF

          RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))
          ND(J)=ND(J+1)*EXP(DELH/XFAC)
          DNDH(J)=-ND(J)*DELH/(XFAC**2)+EXP(DELH/XFAC)*DNDH(J+1)
          OD(J)=OD(J+1)+(ND(J) - ND(J+1))*XFAC*1E5
          Q(J)=ND(J)/RHO
          DQDH(J) = DNDH(J)/RHO

          IF(H(J).LE.HKNEE.AND.JFSH.LT.0)THEN
           F = (HKNEE-H(J))/DELH
           XOD = (1.-F)*OD(J) + F*OD(J+1)
           JFSH=1
          ENDIF
         ENDDO

         DO J=1,NPRO
          OD(J)=XDEEP*OD(J)/XOD
          ND(J)=XDEEP*ND(J)/XOD
          Q(J)=XDEEP*Q(J)/XOD

          IF(Q(J).GT.1e10)Q(J)=1e10
          IF(Q(J).LT.1e-36)Q(J)=1e-36

          DNDH(J)=DNDH(J)*XDEEP/XOD
          DQDH(J)=DQDH(J)*XDEEP/XOD
          DQDX(J)=Q(J)

          XMAP(NXTEMP+1,IPAR,J)=DQDX(J)*XDEEP
          XMAP(NXTEMP+2,IPAR,J)=DQDH(J)*XFAC

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.7)THEN
C        Model 7. Profile modelled with a variable abundance at a fixed 
C        pressure level and variable fractional scale height above and below
C        that level.
C        ***************************************************************
         PKNEE = VARPARAM(IVAR,1)
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
         ENDIF
         XFSH  = EXP(XN(NXTEMP+2))
         XFAC = (1.0-XFSH)/XFSH
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

C         print*,'Requested knee pressure : ',PKNEE
C         print*,'Snapping to : ',P(JFSH)

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

             IF(X1(J).GT.1e10)X1(J)=1e10

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.8)THEN
C        Model 8. Profile is represented by a value at a variable pressure level
C        plus a fractional scale height. Below the knee pressure the profile is 
C        set to zero - a simple cloud in other words!
C        ***************************************************************
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
         ENDIF
         XFSH  = EXP(XN(NXTEMP+2))
         XFAC = (1.0-XFSH)/XFSH
         DXFAC = -1.0/XFSH
         PKNEE = EXP(XN(NXTEMP+3))

         CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
         JFSH = 0
      
         DO J=1,NPRO

          X1(J)=1e-36

          IF(P(J).LT.PKNEE)THEN  
             IF(JFSH.EQ.0)THEN
               DELH = H(J)-HKNEE
               XKEEP = XDEEP

               IF(VARIDENT(IVAR,1).EQ.0)THEN
                XMAP1=1.
               ELSE
                XMAP1=XDEEP
               ENDIF

               XMAP(NXTEMP+1,IPAR,J)=XMAP1*
     1                  EXP(-DELH*XFAC/SCALE(J))
             ELSE

               DELH=H(J)-H(J-1)
               XKEEP = X1(J-1)
               XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
             ENDIF

             X1(J)=XKEEP*EXP(-DELH*XFAC/SCALE(J))

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
C        Model 9. Similar to model 8. Profile is represented by a value 
C        at a variable HEIGHT, rather than PRESSURE plus a fractional scale height. 
C        Below the reference height the profile is set to zero. In addition,
C        this model scales the profile to give the requested integrated cloud 
C        optical depth.
C        ***************************************************************
         IF(VARIDENT(IVAR,1).GE.0)THEN
          PRINT*,'Warning from SUBPROFRETG. You are using a'
          PRINT*,'cloud profile parameterisation for a non-cloud'
          PRINT*,'variable'         
          STOP 
         ENDIF


C        Calculate gradient numerically as it's just too hard otherwise
         DO 22 ITEST=1,4

          XDEEP = EXP(XN(NXTEMP+1))
          XFSH  = EXP(XN(NXTEMP+2))
          HKNEE = XN(NXTEMP+3)

          DX=0.05*XN(NXTEMP+ITEST-1)
          IF(DX.EQ.0.)DX=0.1

          IF(ITEST.GT.1)THEN
C            print*,'ITEST,IPAR,XN,DX = ',ITEST,IPAR,
C     1		XN(NXTEMP+ITEST-1),DX
C            print*,'XDEEP,XFSH,HKNEE',XDEEP,XFSH,HKNEE
          ENDIF
          IF(ITEST.EQ.2)THEN
            XDEEP=EXP(XN(NXTEMP+1)+DX)
          ENDIF
          IF(ITEST.EQ.3)THEN
            XFSH  = EXP(XN(NXTEMP+2)+DX)
          ENDIF
          IF(ITEST.EQ.4)THEN
            HKNEE = XN(NXTEMP+3)+DX
          ENDIF
C          print*,'Mod: ITEST,XDEEP,XFSH,HKNEE',ITEST,XDEEP,
C     &		XFSH,HKNEE




C         Start ND,Q,OD at zero
C         OD is in units of particles/cm2 = particles/cm3 x length(cm)
C         Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
          DO J=1,NPRO
           ND(J)=0.
           OD(J)=0
           Q(J)=0.
           DNDH(J)=0.
           DQDH(J)=0.
          ENDDO

          JFSH=-1

          if(VMRflag.eq.1)then
           if(HKNEE.ge.H(NPRO))THEN
            HKNEE = H(NPRO)
           endif
          endif

          IF(H(1).GE.HKNEE)THEN
            IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
            ELSE
             XMOLWT=XXMOLWT(1)
            ENDIF
C           Calculate density of atmosphere  (g/cm3)
            ND(1)=1.0
            RHO = P(1)*0.1013*XMOLWT/(R*T(1))
            Q(1)=ND(1)/RHO
            JFSH=1
          ENDIF

          DO J=2,NPRO
           DELH = H(J)-H(J-1)
           XFAC = SCALE(J)*XFSH


           IF(H(J).GE.HKNEE)THEN
            IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
            ELSE
             XMOLWT=XXMOLWT(J)
            ENDIF

C           Calculate density of atmosphere  (g/cm3)
            RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))

            IF(JFSH.LT.0)THEN
             ND(J)=1.0
             JFSH=1
            ELSE
             ND(J)=ND(J-1)*DPEXP(-DELH/XFAC)
             DNDH(J)=ND(J)*DELH/(XFAC**2)+DPEXP(-DELH/XFAC)*DNDH(J-1)
            ENDIF
            Q(J)=ND(J)/RHO
            DQDH(J) = DNDH(J)/RHO
           ENDIF

c           print*,'BLAH:',J,H(J), HKNEE, ND(J),JFSH

          ENDDO


C         Integrate optical thickness
          OD(NPRO)=ND(NPRO)*SCALE(NPRO)*XFSH*1E5
          JFSH=-1
          DO J=NPRO-1,1,-1
           DELH = H(J+1) - H(J)
           XFAC = SCALE(J)*XFSH         
           OD(J)=OD(J+1)+0.5*(ND(J) + ND(J+1))*XFAC*1E5
           IF(H(J).LE.HKNEE.AND.JFSH.LT.0)THEN
            F = (HKNEE-H(J))/DELH
            XOD = (1.-F)*OD(J) + F*OD(J+1)
            JFSH=1
           ENDIF
c           print*, 'OD F XOD JFSH',OD(J),F,XOD,JFSH,ND(J),ND(J+1)
          ENDDO

C         The following section was found not to be as accurate as
C         desired due to misalignments at boundaries and so needs some 
C         post-processing in gsetrad.f
          DO J=1,NPRO
           OD(J)=XDEEP*OD(J)/XOD
           ND(J)=XDEEP*ND(J)/XOD
           Q(J)=XDEEP*Q(J)/XOD
           IF(H(J).LT.HKNEE)THEN
            IF(H(J+1).GE.HKNEE)THEN
             Q(J)=Q(J)*(1.0 - (HKNEE-H(J))/(H(J+1)-H(J)))
            ELSE
             Q(J) = 0.0
            ENDIF
           ENDIF
c           print*, 'Q: ', Q(J), XDEEP, XOD
           IF(Q(J).GT.1e10)Q(J)=1e10
           IF(Q(J).LT.1e-36)Q(J)=1e-36
           NTEST=ISNAN(Q(J))
           IF(NTEST)THEN
            print*,'Error in subprofretg.f, cloud density is NAN'
	    STOP
           ENDIF

           IF(ITEST.EQ.1)THEN
            X1(J)=Q(J)
C            print*,'J,X1(J)',J,X1(J)
           ELSE
            XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
C            PRINT*,'ITEST,J,XMAP',ITEST,J,Q(J),X1(J),(Q(J)-X1(J))/DX
           ENDIF

           DNDH(J)=DNDH(J)*XDEEP/XOD
           DQDH(J)=DQDH(J)*XDEEP/XOD
           DQDX(J)=Q(J)/XOD

          ENDDO

22       CONTINUE
        
        ELSEIF(VARIDENT(IVAR,3).EQ.10)THEN
C        Model 10. Profile is a condensing cloud, parameterisation contains
C        the deep gas vmr, the required relative humidity above the condensation
C        level, the required optical depth of the condensed cloud and the 
C        fractional scale height of the condensed cloud. Cloud will condense
C        in cloud profile VARPARAM(IVAR,1).
C        ***************************************************************

         XDEEP = EXP(XN(NXTEMP+1))
         XRH  = EXP(XN(NXTEMP+2))
         XCDEEP = EXP(XN(NXTEMP+3))
         XFSH  = EXP(XN(NXTEMP+4))

         IDAT=0
         DO I=1,NGAS
          IF(GASDATA(I,1).EQ.IDGAS(IPAR))THEN
            A = GASDATA(I,2)
            B = GASDATA(I,3)
            C = GASDATA(I,4)
            D = GASDATA(I,5)
            IDAT=1
          ENDIF
	 ENDDO    

         IF(IDAT.EQ.0)THEN
          print*,'Subprofretg: Gas SVP data cannot be found'
          print*,IPAR,IDGAS(IPAR)
         ENDIF

C        Find where the gas will condense.
         JSPEC=INT(VARPARAM(IVAR,1))

         JPAR = NVMR+1+JSPEC

         IFLA=0
         HKNEE=0.
         DO I=1,NPRO
          P1=P(I)*XDEEP
          PS=DPEXP(A+B/T(I)+C*T(I)+D*T(I)*T(I))
          PH = PS*XRH
          IF(P1.LT.PS)THEN
           X1(I)=XDEEP
           XMAP(NXTEMP+1,IPAR,I)=X1(I)
          ELSE

           IF(IFLA.EQ.0)THEN
            Y1=ALOG(P(I-1)*XDEEP)
            Y2=ALOG(P(I)*XDEEP)
            I1=I-1
            PS1=DPEXP(A+B/T(I1)+C*T(I1)+D*T(I1)*T(I1))

            YY1=ALOG(PS1)
            YY2=ALOG(PS)

            F = (YY1-Y1)/((Y2-Y1)-(YY2-YY1))

            HKNEE = H(I-1)+F*(H(I)-H(I-1))
            IFLA=1
           ENDIF

           X1(I)=PH/P(I)
           XMAP(NXTEMP+2,IPAR,I)=PS/P(I)
          ENDIF

C         This section limits the vmr to relative RH at
C         all altitudes, not just above the condensation level. HKNEE has 
C         already been set

          IF(P1.GT.PH)THEN
            X1(I)=PH/P(I)
            XMAP(NXTEMP+2,IPAR,I)=PS/P(I)
          ELSE
            X1(I)=XDEEP
            XMAP(NXTEMP+1,IPAR,I)=X1(I)
          ENDIF

C         Now make sure that vmr does not rise again once condensation has
C         begun. i.e. freeze vmr at the cold trap.

          IF(IFLA.EQ.1.AND.X1(I).GT.X1(I-1))THEN
           X1(I)=X1(I-1)
           XMAP(NXTEMP+2,IPAR,I)=XMAP(NXTEMP+2,IPAR,I-1)
          ENDIF

         ENDDO

C        Now put a cloud at the condensation level

         IF(AMFORM.EQ.0)THEN
          XMOLWT=MOLWT
         ELSE
          XMOLWT=XXMOLWT(NPRO)
         ENDIF

C        Calculate density of atmosphere (g/cm3)
         RHO = P(NPRO)*0.1013*XMOLWT/(R*T(NPRO))

C        Start ND(NPRO) at a random value. Will be rescaled anyway
         ND(NPRO)=1e-35
C        OD is in units of particles/cm2 = particles/cm3 x length(cm)
C        In this case this is the scale height at the top level.
         OD(NPRO)=ND(NPRO)*SCALE(NPRO)*1E5
C        Q is specific density = particles/gram = particles/cm3 x g/cm3
         Q(NPRO)=ND(NPRO)/RHO         
         DNDH(NPRO)=0.0
         DQDH(NPRO)=0.0
         
         JFSH=-1
         DO J=NPRO-1,1,-1
          DELH = H(J+1)-H(J)
          XFAC = SCALE(J)*XFSH

          IF(AMFORM.EQ.0)THEN
           XMOLWT=MOLWT
          ELSE
           XMOLWT=XXMOLWT(J)
          ENDIF

C         Calculate density of atmosphere (g/cm3)
          RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))
          ND(J)=ND(J+1)*EXP(DELH/XFAC)
          DNDH(J)=-ND(J)*DELH/(XFAC**2)+EXP(DELH/XFAC)*DNDH(J+1)

          OD(J)=OD(J+1)+(ND(J) - ND(J+1))*XFAC*1E5
          Q(J)=ND(J)/RHO
          DQDH(J) = DNDH(J)/RHO

          IF(H(J).LE.HKNEE.AND.JFSH.LT.0)THEN
           F = (HKNEE-H(J))/DELH
           XOD = (1.-F)*OD(J) + F*OD(J+1)
           JFSH=1
          ENDIF
         ENDDO

C        The following section was found not to be as accurate as desired
C        due to misalignments at boundaries and so needs some post-processing in 
C        gsetrad.f
         DO J=1,NPRO
          OD(J)=XCDEEP*OD(J)/XOD
          ND(J)=XCDEEP*ND(J)/XOD
          Q(J)=XCDEEP*Q(J)/XOD
          IF(H(J).LT.HKNEE)THEN
           IF(H(J+1).GE.HKNEE)THEN
            Q(J)=Q(J)*(1.0 - (HKNEE-H(J))/(H(J+1)-H(J)))
           ELSE
            Q(J) = 0.0
           ENDIF
          ENDIF
          IF(Q(J).GT.1e10)Q(J)=1e10
          IF(Q(J).LT.1e-36)Q(J)=1e-36

          X2(J)=Q(J)
          DNDH(J)=DNDH(J)*XCDEEP/XOD
          DQDH(J)=DQDH(J)*XCDEEP/XOD
          DQDX(J)=Q(J)/XOD

          IF(H(J).LT.HKNEE)THEN
           XMAP(NXTEMP+3,JPAR,J)=DQDX(J)*XCDEEP
           XMAP(NXTEMP+4,JPAR,J)=DQDH(J)*XFAC
          ENDIF
         ENDDO


        ELSEIF(VARIDENT(IVAR,3).EQ.11)THEN
C        Model 11. Condensing gas, but no associated cloud. Model requires
C        the deep gas abundance and the desired relative humidity above the 
C        condensation level only (VARPARAM(IVAR,1))=0) or at all levels
C        (VARPARAM(IVAR,1))=1)
C        ***************************************************************

         XDEEP = EXP(XN(NXTEMP+1))
         XRH  = EXP(XN(NXTEMP+2))
         ICOND = NINT(VARPARAM(IVAR,1))

         IDAT=0
         DO I=1,NGAS
          IF(GASDATA(I,1).EQ.IDGAS(IPAR))THEN
            A = GASDATA(I,2)
            B = GASDATA(I,3)
            C = GASDATA(I,4)
            D = GASDATA(I,5)
            IDAT=1
          ENDIF
	 ENDDO
    

         IF(IDAT.EQ.0)THEN
          print*,'Subprofretg: Gas SVP data cannot be found'
          print*,IPAR,IDGAS(IPAR)
         ENDIF

         IFLA=0
         HKNEE=0.
         DO I=1,NPRO
          P1=P(I)*XDEEP
          PS=DPEXP(A+B/T(I)+C*T(I)+D*T(I)*T(I))
          PH = PS*XRH
          IF(P1.LT.PS)THEN
           X1(I)=XDEEP
           XMAP(NXTEMP+1,IPAR,I)=X1(I)
          ELSE

           IF(IFLA.EQ.0)THEN
            Y1=ALOG(P(I-1)*XDEEP)
            Y2=ALOG(P(I)*XDEEP)
            I1=I-1
            PS1=DPEXP(A+B/T(I1)+C*T(I1)+D*T(I1)*T(I1))

            YY1=ALOG(PS1)
            YY2=ALOG(PS)

            F = (YY1-Y1)/((Y2-Y1)-(YY2-YY1))

            HKNEE = H(I-1)+F*(H(I)-H(I-1))
            IFLA=1
           ENDIF

           X1(I)=PH/P(I)
           XMAP(NXTEMP+2,IPAR,I)=PS/P(I)

          ENDIF

C         Section to determine if RH is to apply at all levels (ICOND = 0)
C         or only above the condensation level (ICOND = 1)
C         ---------------------------------------------------------------
          IF(ICOND.EQ.0)THEN
           IF(P1.GT.PH)THEN
             X1(I)=PH/P(I)
             XMAP(NXTEMP+2,IPAR,I)=PS/P(I)
           ELSE
             X1(I)=XDEEP
             XMAP(NXTEMP+1,IPAR,I)=X1(I)
           ENDIF
          ENDIF
C         ---------------------------------------------------------------


C         Now make sure that vmr does not rise again once condensation has
C         begun. i.e. freeze vmr at the cold trap.

          IF(IFLA.EQ.1.AND.X1(I).GT.X1(I-1))THEN
           X1(I)=X1(I-1)
           XMAP(NXTEMP+2,IPAR,I)=XMAP(NXTEMP+2,IPAR,I-1)
          ENDIF

         ENDDO


        ELSEIF(VARIDENT(IVAR,3).EQ.12)THEN
C        Model 12. Profile is represented a Gaussian with a peak value
C        at a variable pressure level plus a variable FWHM in log pressure.
C        ***************************************************************
         XDEEP = EXP(XN(NXTEMP+1))
         PKNEE = EXP(XN(NXTEMP+2))
         XWID  = EXP(XN(NXTEMP+3))


         Y0=ALOG(PKNEE)

         DO J=1,NPRO

          Y=ALOG(P(J))          
          
          X1(J)=XDEEP*EXP(-((Y-Y0)/XWID)**2)
          IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,J)=X1(J)/XDEEP
          ELSE
            XMAP(NXTEMP+1,IPAR,J)=X1(J)
          ENDIF

          XMAP(NXTEMP+2,IPAR,J)=Y0*2.*(Y-Y0)*X1(J)/XWID**2
          XMAP(NXTEMP+3,IPAR,J)=-2.0*((Y-Y0)**2)*X1(J)/XWID**2

         ENDDO


        ELSEIF(VARIDENT(IVAR,3).EQ.13)THEN
C        Model 13. Profile is represented a Lorentzian with a peak value
C        at a variable pressure level plus a variable FWHM in log pressure.
C        ***************************************************************
         XDEEP = EXP(XN(NXTEMP+1))
         PKNEE = EXP(XN(NXTEMP+2))
         XWID  = EXP(XN(NXTEMP+3))

         Y0=ALOG(PKNEE)

         DO J=1,NPRO

          Y=ALOG(P(J))
         
          XX = (Y-Y0)**2 + XWID**2 
          
          X1(J)=XDEEP*(XWID**2)/XX

          IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,J)=X1(J)/XDEEP
          ELSE
            XMAP(NXTEMP+1,IPAR,J)=X1(J)
          ENDIF
         
          XMAP(NXTEMP+2,IPAR,J)=Y0*XDEEP*2.*(Y-Y0)*(XWID**2)/XX**2
          XMAP(NXTEMP+3,IPAR,J)=XWID*2*XDEEP*XWID*(XX - XWID**2)/XX**2

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.14)THEN
C        Model 14. Profile is represented a Gaussian with a specified optical
C        thickness centred at a variable altitude level plus a variable FWHM in 
C        height.
C        FHWM is also folded into total opacity, but this gets renormalised
C        by gsetrad.f anyway.
C        ***************************************************************

         XDEEP = EXP(XN(NXTEMP+1))
         HKNEE = XN(NXTEMP+2)
         XWID  = EXP(XN(NXTEMP+3))

         Y0=HKNEE

C         print*,'Xdeep, hknee, xwid',XDEEP,HKNEE,XWID

C        **** Want to normalise to get optical depth right. ***
C        ND is the particles per cm3 (but will be rescaled)
C        OD is in units of particles/cm2 = particles/cm3 x length(cm)
C        OD(J)=ND(J)*SCALE(J)*1E5
C        Q is specific density = particles/gram = particles/cm3 x g/cm3
C        Q(J)=ND(J)/RHO         
         
         XOD=0.
         DO J=1,NPRO

          IF(AMFORM.EQ.0)THEN
           XMOLWT=MOLWT
          ELSE
           XMOLWT=XXMOLWT(J)
          ENDIF

          RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))

          Y=H(J)
                
          Q(J) = 1./(XWID*SQRT(PI))*EXP(-((Y-Y0)/XWID)**2)
          ND(J) = Q(J)*RHO 
          OD(J) = ND(J)*SCALE(J)*1e5

          XOD=XOD+OD(J)
   
          X1(J)=Q(J)

         ENDDO

C        Empirical correction to XOD
         XOD = XOD*0.25

C         print*,'XOD = ',XOD


         DO J=1,NPRO

          X1(J)=Q(J)*XDEEP/XOD   

          IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,J)=X1(J)/XDEEP
          ELSE
            XMAP(NXTEMP+1,IPAR,J)=X1(J)
          ENDIF

          XMAP(NXTEMP+2,IPAR,J)=2.*(Y-Y0)*X1(J)/XWID**2
          XMAP(NXTEMP+3,IPAR,J)=-2.0*((Y-Y0)**2)*X1(J)/XWID**3
     &             -  X1(J)/XWID

         ENDDO

C        *** This renormalisation is pretty accurate, but not quite accurate
C        *** enough and so it gets updated in gsetrad.f

        ELSEIF(VARIDENT(IVAR,3).EQ.15)THEN
C        Model 15. Profile is represented a Lorentzian with a specified optical
C        thickness centred at a variable altitude level plus a variable FWHM in 
C        height.
C        FHWM is also folded into total opacity, but this gets renormalised
C        by gsetrad.f anyway.
C        ***************************************************************
         XDEEP = EXP(XN(NXTEMP+1))
         HKNEE = XN(NXTEMP+2)
         XWID  = EXP(XN(NXTEMP+3))

         Y0=HKNEE


C        ND is the particles per cm3 (but will be rescaled)
C        OD is in units of particles/cm2 = particles/cm3 x length(cm)
C        OD(J)=ND(J)*SCALE(J)*1E5
C        Q is specific density = particles/gram = particles/cm3 x g/cm3
C        Q(J)=ND(J)/RHO         
         
         XOD=0.

         DO J=1,NPRO

          IF(AMFORM.EQ.0)THEN
           XMOLWT=MOLWT
          ELSE
           XMOLWT=XXMOLWT(J)
          ENDIF

          RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))

          Y=H(J)

          XX = (Y-Y0)**2 + XWID**2 
          
          Q(J) = (1./PI)*XWID/XX
          ND(J) = Q(J)*RHO 
          OD(J) = ND(J)*SCALE(J)*1e5

          XOD=XOD+OD(J)
   
          X1(J)=Q(J)

         ENDDO

C        Empirical correction to XOD
         XOD = XOD*0.25

C         print*,'XOD = ',XOD

         DO J=1,NPRO
         
          X1(J)=Q(J)*XDEEP/XOD   

          IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,J)=X1(J)/XDEEP
          ELSE
            XMAP(NXTEMP+1,IPAR,J)=X1(J)
          ENDIF
         
          XMAP(NXTEMP+2,IPAR,J)=(XDEEP/(XOD*PI))*2.*
     &                                  (Y-Y0)*XWID/XX**2
          XMAP(NXTEMP+3,IPAR,J)=XWID*(XDEEP/(XOD*PI))*
     &					(XX - 2*XWID**2)/XX**2

         ENDDO

C        *** This renormalisation is pretty accurate, but not quite accurate
C        *** enough and so it gets updated in gsetrad.f

        ELSEIF(VARIDENT(IVAR,3).EQ.16)THEN
C        Model 16. Profile modelled with a variable abundance at a fixed 
C        pressure level and a lapse rate above and lapse rate below
C        that level.
C        ***************************************************************
         PKNEE = VARPARAM(IVAR,1)
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
         ENDIF
         PKNEE  = EXP(XN(NXTEMP+2))
         XLDEEP  = EXP(XN(NXTEMP+3))
         XLHIGH  = EXP(XN(NXTEMP+4))

         CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
         JFSH = 0       

         JKNEE=-1
         DO J=1,NPRO
          IF(H(J).GT.HKNEE)THEN
           JKNEE=J-1
           GOTO 211
          ENDIF
         ENDDO
211      CONTINUE

         IF(JKNEE.LT.1.OR.JKNEE.GT.NPRO-1)THEN
          PRINT*,'SUBPROFRETG. Must choose pressure level'
          PRINT*,'within range of profile'
          PRINT*,'Model = 16'
          STOP
         ENDIF

        
         DELH = HKNEE-H(JKNEE)
         X1(JKNEE)=XDEEP + DELH*XLDEEP
         IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,JKNEE)=1.0
         ELSE
            XMAP(NXTEMP+1,IPAR,JKNEE)=XDEEP
         ENDIF
         XMAP(NXTEMP+2,IPAR,JKNEE)=-XLDEEP*SCALE(JKNEE)
         XMAP(NXTEMP+3,IPAR,JKNEE)=DELH*XLDEEP
         XMAP(NXTEMP+4,IPAR,JKNEE)=0.


         DELH = H(JKNEE+1)-HKNEE
         X1(JKNEE+1)=XDEEP + DELH*XLHIGH
         IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,JKNEE+1)=1.0
         ELSE
            XMAP(NXTEMP+1,IPAR,JKNEE+1)=XDEEP
         ENDIF
         XMAP(NXTEMP+2,IPAR,JKNEE+1)=XLHIGH*SCALE(JKNEE)
         XMAP(NXTEMP+3,IPAR,JKNEE+1)=0.
         XMAP(NXTEMP+4,IPAR,JKNEE+1)=DELH*XLHIGH


         DO J=JKNEE+2,NPRO 
             
             DELH = H(J)-H(J-1)         
             X1(J)=X1(J-1)+DELH*XLHIGH

             XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)
C             XMAP(NXTEMP+2,IPAR,J)=XMAP(NXTEMP+2,IPAR,J-1)+
C     &					XLHIGH*SCALE(J)
             XMAP(NXTEMP+2,IPAR,J)=XLHIGH*SCALE(J)

             XMAP(NXTEMP+3,IPAR,J)=0.
             XMAP(NXTEMP+4,IPAR,J)=XMAP(NXTEMP+4,IPAR,J-1)+
     &					DELH*XLHIGH

             IF(X1(J).LT.1e-36)X1(J)=1e-36
             IF(X1(J).GT.1e10)X1(J)=1e10

         ENDDO

         DO J=JKNEE-1,1,-1             
             DELH = H(J+1)-H(J)
             X1(J)=X1(J+1)+DELH*XLDEEP

             XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J+1)
C             XMAP(NXTEMP+2,IPAR,J)=XMAP(NXTEMP+2,IPAR,J+1)-
C     &				       XLDEEP*SCALE(J)
             XMAP(NXTEMP+2,IPAR,J)=-XLDEEP*SCALE(J)
             XMAP(NXTEMP+3,IPAR,J)=XMAP(NXTEMP+3,IPAR,J+1)+
     &				       DELH*XLDEEP
             XMAP(NXTEMP+4,IPAR,J)=0.

             IF(X1(J).GT.1e10)X1(J)=1e10
             IF(X1(J).LT.1e-36)X1(J)=1e-36

         ENDDO

C         open(12,file='test.dump',status='unknown')
C         write(12,*)npro
C         do j=1,npro
C          write(12,*)x1(j),(xmap(k,ipar,j),k=nxtemp+1,nxtemp+4)
C         enddo
C         close(12)


        ELSEIF(VARIDENT(IVAR,3).EQ.17)THEN
C        Model 17. Profile modelled with a variable abundance at a fixed 
C        pressure level and variable fractional scale height above and below
C        that level. Profile is limited at pressures less than a specified 
C        level
C        ***************************************************************
         PKNEE = VARPARAM(IVAR,1)
         PLIM = VARPARAM(IVAR,2)
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
         ENDIF
         XFSH  = EXP(XN(NXTEMP+2))
         XFAC = (1.0-XFSH)/XFSH
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

C         print*,'Requested knee pressure : ',PKNEE
C         print*,'Snapping to : ',P(JFSH)

         X1(JFSH)=XDEEP

         IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,JFSH)=1.0
         ELSE
            XMAP(NXTEMP+1,IPAR,JFSH)=X1(JFSH)
         ENDIF


         DO J=JFSH+1,NPRO 
             DELH = H(J)-H(J-1)
             IF(P(J).GT.PLIM)THEN          
              X1(J)=X1(J-1)*EXP(-DELH*XFAC/SCALE(J))
              XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
              XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1          X1(J-1)*EXP(-DELH*XFAC/SCALE(J)) +
     2          XMAP(NXTEMP+2,IPAR,J-1)*EXP(-DELH*XFAC/SCALE(J))
             ELSE
              X1(J)=X1(J-1)
              XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)
              XMAP(NXTEMP+2,IPAR,J)=XMAP(NXTEMP+2,IPAR,J-1)
             ENDIF

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

             IF(X1(J).GT.1e10)X1(J)=1e10

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.18)THEN
C        Model 7. Profile modelled with a variable abundance at a fixed 
C        pressure level and variable fractional scale height below
C        that level. Abundance above is fixed
C        ***************************************************************
         PKNEE = VARPARAM(IVAR,1)
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
         ENDIF
         XFSH  = EXP(XN(NXTEMP+2))
         XFAC = (1.0-XFSH)/XFSH
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

C         print*,'Requested knee pressure : ',PKNEE
C         print*,'Snapping to : ',P(JFSH)

         X1(JFSH)=XDEEP

         IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,JFSH)=1.0
         ELSE
            XMAP(NXTEMP+1,IPAR,JFSH)=X1(JFSH)
         ENDIF


         DO J=JFSH+1,NPRO 
             X1(J)=X1(J-1)
             XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)
             XMAP(NXTEMP+2,IPAR,J)=0.

         ENDDO

         DO J=JFSH-1,1,-1             
             DELH = H(J)-H(J+1)
             X1(J)=X1(J+1)*EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J+1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
             XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1          X1(J+1)*EXP(-DELH*XFAC/SCALE(J)) +
     2          XMAP(NXTEMP+2,IPAR,J+1)*EXP(-DELH*XFAC/SCALE(J))

             IF(X1(J).GT.1e10)X1(J)=1e10
             IF(X1(J).LT.1e-36)X1(J)=1e-36

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.19)THEN
C        Model 19. Identical to model 9, but pinches off abudance near tropopause.
C        ***************************************************************
         IF(VARIDENT(IVAR,1).GE.0)THEN
          PRINT*,'Warning from SUBPROFRETG. You are using a'
          PRINT*,'cloud profile parameterisation for a non-cloud'
          PRINT*,'variable'         
          STOP 
         ENDIF

         IF(AMFORM.EQ.0)THEN
          XMOLWT=MOLWT
         ELSE
          XMOLWT=XXMOLWT(NPRO)
         ENDIF

C        Calculate density of atmosphere (g/cm3)
         RHO = P(NPRO)*0.1013*XMOLWT/(R*T(NPRO))

C        Start ND(NPRO) at a random value. Will be rescaled anyway
         ND(NPRO)=1e-35
C        OD is in units of particles/cm2 = particles/cm3 x length(cm)
C        In this case this is the scale height at the top level.
         OD(NPRO)=ND(NPRO)*SCALE(NPRO)*1E5
C        Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
         Q(NPRO)=ND(NPRO)/RHO         
         DNDH(NPRO)=0.0
         DQDH(NPRO)=0.0


C        Calculate gradient numerically as it's just too hard otherwise
         DO 32 ITEST=1,5

          XDEEP = EXP(XN(NXTEMP+1))
          XFSH  = EXP(XN(NXTEMP+2))
          HKNEE = XN(NXTEMP+3)
          CWID = XN(NXTEMP+4)

          DX=0.05*XN(NXTEMP+ITEST-1)
          IF(DX.EQ.0.)DX=0.1

          IF(ITEST.GT.1)THEN
C            print*,'ITEST,IPAR,XN,DX = ',ITEST,IPAR,
C     1		XN(NXTEMP+ITEST-1),DX
C            print*,'XDEEP,XFSH,HKNEE',XDEEP,XFSH,HKNEE
          ENDIF
          IF(ITEST.EQ.2)THEN
            XDEEP=EXP(XN(NXTEMP+1)+DX)
          ENDIF
          IF(ITEST.EQ.3)THEN
            XFSH  = EXP(XN(NXTEMP+2)+DX)
          ENDIF
          IF(ITEST.EQ.4)THEN
            HKNEE = XN(NXTEMP+3)+DX
          ENDIF
          IF(ITEST.EQ.4)THEN
            HKNEE = XN(NXTEMP+4)+DX
          ENDIF
C          print*,'Mod: XDEEP,XFSH,HKNEE',XDEEP,XFSH,HKNEE

          JFSH=-1
          DO J=NPRO-1,1,-1
           DELH = H(J+1)-H(J)
           XFAC = SCALE(J)*XFSH
C          Calculate density of atmosphere (g/cm3)
         
           IF(AMFORM.EQ.0)THEN
            XMOLWT=MOLWT
           ELSE
            XMOLWT=XXMOLWT(J)
           ENDIF

           RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))

           ND(J)=ND(J+1)*EXP(DELH/XFAC)

           IF(P(J).GT.0.1)THEN
             XFACP=1.0-EXP(-((LOG(P(J))-LOG(0.1))/CWID)**2)
             ND(J)=ND(J)*XFACP
             IF(ND(J).LT.1E-35)ND(J)=1e-25
           ELSE
             ND(J)=1e-35
           ENDIF 

           DNDH(J)=-ND(J)*DELH/(XFAC**2)+EXP(DELH/XFAC)*DNDH(J+1)

           OD(J)=OD(J+1)+(ND(J) - ND(J+1))*XFAC*1E5
           Q(J)=ND(J)/RHO
           DQDH(J) = DNDH(J)/RHO
           IF(H(J).LE.HKNEE.AND.JFSH.LT.0)THEN
            F = (HKNEE-H(J))/DELH
            XOD = (1.-F)*OD(J) + F*OD(J+1)
            JFSH=1
           ENDIF
          ENDDO

C         The following section was found not to be as accurate as
C         desired due to misalignments at boundaries and so needs some 
C         post-processing in gsetrad.f
          DO J=1,NPRO
           OD(J)=XDEEP*OD(J)/XOD
           ND(J)=XDEEP*ND(J)/XOD
           Q(J)=XDEEP*Q(J)/XOD
           IF(H(J).LT.HKNEE)THEN
            IF(H(J+1).GE.HKNEE)THEN
             Q(J)=Q(J)*(1.0 - (HKNEE-H(J))/(H(J+1)-H(J)))
            ELSE
             Q(J) = 0.0
            ENDIF
           ENDIF
           IF(Q(J).GT.1e10)Q(J)=1e10
           IF(Q(J).LT.1e-36)Q(J)=1e-36
           NTEST=ISNAN(Q(J))
           IF(NTEST)THEN
            print*,'Error in subprofretg.f, cloud density is NAN'
            print*,'XDEEP,XOD,HKNEE = ',XDEEP,XOD,HKNEE
            DO JX =1,NPRO
             OD(JX)=XDEEP*OD(JX)/XOD
             ND(JX)=XDEEP*ND(JX)/XOD
             Q(JX)=XDEEP*Q(JX)/XOD
             print*,'JX,OD,ND,Q = ',JX,OD(JX),ND(JX),Q(JX)
             IF(H(JX).LT.HKNEE)THEN
              IF(H(JX+1).GE.HKNEE)THEN
               Q(JX)=Q(JX)*(1.0 - (HKNEE-H(JX))/(H(JX+1)-H(JX)))
              ELSE
               Q(JX) = 0.0
              ENDIF
              print*,'H(JX),H(JX+1),HKNEE',H(JX),H(JX+1),HKNEE
              print*,HKNEE-H(JX),H(JX+1)-H(JX)
              print*,(1.0 - (HKNEE-H(JX))/(H(JX+1)-H(JX)))
             ENDIF
             IF(Q(JX).GT.1e10)Q(JX)=1e10
             IF(Q(JX).LT.1e-36)Q(JX)=1e-36
            ENDDO
	    STOP
           ENDIF

           IF(ITEST.EQ.1)THEN
            X1(J)=Q(J)
C            print*,'J,X1(J)',J,X1(J)
           ELSE
            XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
C            PRINT*,'ITEST,J,XMAP',ITEST,J,Q(J),X1(J),(Q(J)-X1(J))/DX
           ENDIF

           DNDH(J)=DNDH(J)*XDEEP/XOD
           DQDH(J)=DQDH(J)*XDEEP/XOD
           DQDX(J)=Q(J)/XOD

          ENDDO

32       CONTINUE
        

        ELSEIF(VARIDENT(IVAR,3).EQ.20)THEN
C        Model 1. Variable deep abundance, fixed knee pressure, fixed
C        tropopause temperature and variable fractional scale height.
C        ***************************************************************
         PKNEE = VARPARAM(IVAR,1)
         PTROP = VARPARAM(IVAR,2)
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

             IF(X1(J).LT.1e-36)X1(J)=1e-36

          END IF
          IF(P(J).LT.PTROP)THEN
           X1(J)=1e-36
          ENDIF
         ENDDO


        ELSEIF(VARIDENT(IVAR,3).EQ.21)THEN
C        Model 21. Similar to model 9. Profile is represented by a value 
C        at a variable HEIGHT, rather than PRESSURE plus a fractional scale height. 
C        Below the reference height the profile is set to zero. In addition,
C        this model scales the profile to give the requested integrated cloud 
C        optical depth.
C        ***************************************************************
         IF(VARIDENT(IVAR,1).GE.0)THEN
          PRINT*,'Warning from SUBPROFRETG. You are using a'
          PRINT*,'cloud profile parameterisation for a non-cloud'
          PRINT*,'variable'         
          STOP 
         ENDIF


C        Calculate gradient numerically as it's just too hard otherwise
         DO 24 ITEST=1,3

          XDEEP = EXP(XN(NXTEMP+1))
          HKNEE = XN(NXTEMP+2)
          XFSH  = VARPARAM(IVAR,1)

C         Need radius from associated 444 particle parameterisation
          ICL=-VARIDENT(IVAR,1)
C          print*,'ICLOUD = ',ICL
          RPARTICLE = GETRADIUS(ICL,NVAR,VARIDENT,VARPARAM,XN,NPRO)
          IF(RPARTICLE.LT.0)THEN
           PRINT*,'Particle size cannot be found'
           STOP
          ENDIF
          REFRADIUS=VARPARAM(IVAR,2)
          XFSH=XFSH*REFRADIUS/RPARTICLE

C          print*,'Test rparticle,refradius,xfsh,icl = ',
C     1     RPARTICLE,REFRADIUS,XFSH,ICL
          DX=0.05*XN(NXTEMP+ITEST-1)
          IF(DX.EQ.0.)DX=0.1

          IF(ITEST.GT.1)THEN
C            print*,'ITEST,IPAR,XN,DX = ',ITEST,IPAR,
C     1		XN(NXTEMP+ITEST-1),DX
C            print*,'XDEEP,XFSH,HKNEE',XDEEP,XFSH,HKNEE
          ENDIF
          IF(ITEST.EQ.2)THEN
            XDEEP=EXP(XN(NXTEMP+1)+DX)
          ENDIF
          IF(ITEST.EQ.3)THEN
            HKNEE = XN(NXTEMP+3)+DX
          ENDIF
C          print*,'Mod: ITEST,XDEEP,XFSH,HKNEE',ITEST,XDEEP,
C     &		XFSH,HKNEE




C         Start ND,Q,OD at zero
C         OD is in units of particles/cm2 = particles/cm3 x length(cm)
C         Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
          DO J=1,NPRO
           ND(J)=0.
           OD(J)=0
           Q(J)=0.
           DNDH(J)=0.
           DQDH(J)=0.
          ENDDO

          JFSH=-1
          IF(H(1).GE.HKNEE)THEN
            IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
            ELSE
             XMOLWT=XXMOLWT(1)
            ENDIF
C           Calculate density of atmosphere  (g/cm3)
            ND(1)=1.0
            RHO = P(1)*0.1013*XMOLWT/(R*T(1))
            Q(1)=ND(1)/RHO
            JFSH=1
          ENDIF

          DO J=2,NPRO
           DELH = H(J)-H(J-1)
           XFAC = SCALE(J)*XFSH


           IF(H(J).GE.HKNEE)THEN
            IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
            ELSE
             XMOLWT=XXMOLWT(J)
            ENDIF

C           Calculate density of atmosphere  (g/cm3)
            RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))

            IF(JFSH.LT.0)THEN
             ND(J)=1.0
             JFSH=1
            ELSE
             ND(J)=ND(J-1)*DPEXP(-DELH/XFAC)
             DNDH(J)=ND(J)*DELH/(XFAC**2)+DPEXP(-DELH/XFAC)*DNDH(J-1)
            ENDIF
            Q(J)=ND(J)/RHO
            DQDH(J) = DNDH(J)/RHO
           ENDIF

C           print*,J,H(J),ND(J),Q(J)

          ENDDO


C         Integrate optical thickness
          OD(NPRO)=ND(NPRO)*SCALE(NPRO)*XFSH*1E5
          JFSH=-1
          DO J=NPRO-1,1,-1
           XFAC = SCALE(J)*XFSH         
           OD(J)=OD(J+1)+0.5*(ND(J) + ND(J+1))*XFAC*1E5
           IF(H(J).LE.HKNEE.AND.JFSH.LT.0)THEN
            F = (HKNEE-H(J))/DELH
            XOD = (1.-F)*OD(J) + F*OD(J+1)
            JFSH=1
           ENDIF
C           print*,'h',J,OD(J)
          ENDDO

C         The following section was found not to be as accurate as
C         desired due to misalignments at boundaries and so needs some 
C         post-processing in gsetrad.f
          DO J=1,NPRO
           OD(J)=XDEEP*OD(J)/XOD
           ND(J)=XDEEP*ND(J)/XOD
           Q(J)=XDEEP*Q(J)/XOD
           IF(H(J).LT.HKNEE)THEN
            IF(H(J+1).GE.HKNEE)THEN
             Q(J)=Q(J)*(1.0 - (HKNEE-H(J))/(H(J+1)-H(J)))
            ELSE
             Q(J) = 0.0
            ENDIF
           ENDIF
           IF(Q(J).GT.1e10)Q(J)=1e10
           IF(Q(J).LT.1e-36)Q(J)=1e-36
           NTEST=ISNAN(Q(J))
           IF(NTEST)THEN
            print*,'Error in subprofretg.f, cloud density is NAN'
	    STOP
           ENDIF

           IF(ITEST.EQ.1)THEN
            X1(J)=Q(J)
C            print*,'J,X1(J)',J,X1(J)
           ELSE
            XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
C            PRINT*,'ITEST,J,XMAP',ITEST,J,Q(J),X1(J),(Q(J)-X1(J))/DX
           ENDIF

           DNDH(J)=DNDH(J)*XDEEP/XOD
           DQDH(J)=DQDH(J)*XDEEP/XOD
           DQDX(J)=Q(J)/XOD

          ENDDO

24       CONTINUE


        ELSE

         PRINT*,'Subprofretg: Model parametrisation code is not defined'
         PRINT*,(VARIDENT(IVAR,J),J=1,3)
         STOP

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
        ELSEIF(VARIDENT(IVAR,1).EQ.887)THEN
C         print*,'Cloud x-section spectrum'
         IPAR = -1
         NP = INT(VARPARAM(IVAR,1))
        ELSEIF(VARIDENT(IVAR,1).EQ.889)THEN
C         print*,'Surface albedo spectrum multiplier'
         IPAR = -1
         NP = 1
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
        ELSEIF(VARIDENT(IVAR,1).EQ.444)THEN
C         print*,'Variable size and RI'
C        See if there is an associated IMOD=21 cloud. In which case
C        modifying the radius will affect the vertical cloud distribution.
         IPAR = -1
         ICL=VARIDENT(IVAR,2)
         RPARTICLE=EXP(XN(NXTEMP+1))
         CALL GETCLOUD21(ICLOUD,NVAR,VARIDENT,VARPARAM,XN,NPRO,
     1 XDEEP,HKNEE,XFSHREF,REFRADIUS)

         IF(XDEEP.GT.0)THEN
          IPAR=NVMR+1+ICL
          XFSH=XFSHREF*REFRADIUS/RPARTICLE

C         Calculate gradient numerically as it's just too hard otherwise
          DO 25 ITEST=1,2

           DX=0.05*XN(NXTEMP+ITEST-1)
           IF(DX.EQ.0.)DX=0.1

           IF(ITEST.EQ.2)THEN
            XFSH=XFSHREF*REFRADIUS/EXP(XN(NXTEMP+1)+DX)
           ENDIF


C          Start ND,Q,OD at zero
C          OD is in units of particles/cm2 = particles/cm3 x length(cm)
C          Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
           DO J=1,NPRO
            ND(J)=0.
            OD(J)=0
            Q(J)=0.
           ENDDO

           JFSH=-1
           IF(H(1).GE.HKNEE)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(1)
             ENDIF
C            Calculate density of atmosphere  (g/cm3)
             ND(1)=1.0
             RHO = P(1)*0.1013*XMOLWT/(R*T(1))
             Q(1)=ND(1)/RHO
             JFSH=1
           ENDIF

           DO J=2,NPRO
            DELH = H(J)-H(J-1)
            XFAC = SCALE(J)*XFSH


            IF(H(J).GE.HKNEE)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(J)
             ENDIF

C            Calculate density of atmosphere  (g/cm3)
             RHO = (0.1013*XMOLWT/R)*(P(J)/T(J)) 

             IF(JFSH.LT.0)THEN
              ND(J)=1.0
              JFSH=1
             ELSE
              ND(J)=ND(J-1)*DPEXP(-DELH/XFAC)
             ENDIF
             Q(J)=ND(J)/RHO
            ENDIF


           ENDDO


C          Integrate optical thickness
           OD(NPRO)=ND(NPRO)*SCALE(NPRO)*XFSH*1E5
           JFSH=-1
           DO J=NPRO-1,1,-1
            XFAC = SCALE(J)*XFSH         
            OD(J)=OD(J+1)+0.5*(ND(J) + ND(J+1))*XFAC*1E5
            IF(H(J).LE.HKNEE.AND.JFSH.LT.0)THEN
             F = (HKNEE-H(J))/DELH
             XOD = (1.-F)*OD(J) + F*OD(J+1)
             JFSH=1
            ENDIF
C            print*,'h',J,OD(J)
           ENDDO

C          The following section was found not to be as accurate as
C          desired due to misalignments at boundaries and so needs some 
C          post-processing in gsetrad.f
           DO J=1,NPRO
            OD(J)=XDEEP*OD(J)/XOD
            Q(J)=XDEEP*Q(J)/XOD
            IF(H(J).LT.HKNEE)THEN
             IF(H(J+1).GE.HKNEE)THEN
              Q(J)=Q(J)*(1.0 - (HKNEE-H(J))/(H(J+1)-H(J)))
             ELSE
              Q(J) = 0.0
             ENDIF
            ENDIF
            IF(Q(J).GT.1e10)Q(J)=1e10
            IF(Q(J).LT.1e-36)Q(J)=1e-36
            NTEST=ISNAN(Q(J))
            IF(NTEST)THEN
             print*,'Error in subprofretg.f, cloud density is NAN'
	     STOP
            ENDIF

            IF(ITEST.EQ.1)THEN
             X2(J)=Q(J)
            ELSE
             XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X2(J))/DX
            ENDIF

           ENDDO

25       CONTINUE


         ENDIF
         NP = 2+INT(VARPARAM(IVAR,1))
        ELSEIF(VARIDENT(IVAR,1).EQ.333)THEN
C         print*,'Surface gravity (log10(g))'
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.222)THEN
C         print*,'Sromovsky cloud layering'
         IPAR = -1
         NP = 8
        ELSEIF(VARIDENT(IVAR,1).EQ.223)THEN
C         print*,'Sromovsky cloud layering with methane'
         IPAR = -1
         DO I=1,NVMR
           IF(IDGAS(I).EQ.6)IPAR=I
         ENDDO
         IF(IPAR.LT.0)THEN
           PRINT*,'Error in subprofretg. Model 223 defined, but'
           PRINT*,'no CH4 in .ref file'
           STOP
         ENDIF
         DO I=1,NPRO
           X1(I)=VMR(I,IPAR)
         ENDDO
         PCH4 = EXP(XN(NXTEMP+5))/1.013
         XFAC = EXP(XN(NXTEMP+6))
         IF(XFAC.GT.1.0)THEN
          PRINT*,'Error in subprofretg, model 223. XFAC > 1.'
          PRINT*,'Limiting to 1.0'
          XFAC=1.
         ENDIF
         XCH4 = X1(1)*XFAC
C        Set PCUT to where CH4 starts reducing in cases where XFAC=1
         ICUT=0
         DO I=1,NPRO
          IF(X1(I).LT.X1(1).AND.ICUT.EQ.0)THEN
            PCUT=P(I)
            ICUT=1
          ENDIF
         ENDDO
         DO I=1,NPRO
          IF(P(I).LT.PCH4)THEN
           IF(XCH4.LT.X1(I))THEN
            X1(I)=XCH4
            PCUT=P(I)
           ENDIF
          ENDIF
         ENDDO
         NP = 9
        ELSEIF(VARIDENT(IVAR,1).EQ.224)THEN
C         print*,'Sromovsky cloud layering with extended UTC'
         IPAR = -1
         NP = 9
        ELSEIF(VARIDENT(IVAR,1).EQ.225)THEN
C         print*,'Sromovsky cloud layering with extended UTC, cut-off'
C         print*,'and variable methane'
         IPAR = -1
         DO I=1,NVMR
           IF(IDGAS(I).EQ.6)IPAR=I
         ENDDO
         IF(IPAR.LT.0)THEN
           PRINT*,'Error in subprofretg. Model 225 defined, but'
           PRINT*,'no CH4 in .ref file'
           STOP
         ENDIF
         DO I=1,NPRO
           X1(I)=VMR(I,IPAR)
         ENDDO
         PCH4 = EXP(XN(NXTEMP+5))/1.013
         XFAC = EXP(XN(NXTEMP+11))
C         IF(XFAC.GT.1.0)THEN
C          PRINT*,'Error in subprofretg, model 223. XFAC > 1.'
C          PRINT*,'Limiting to 1.0'
C          XFAC=1.
C         ENDIF
C         XCH4 = X1(1)*XFAC
C        Set PCUT to where CH4 starts reducing in cases where XFAC=1
C         ICUT=0
C         DO I=1,NPRO
C          IF(X1(I).LT.X1(1).AND.ICUT.EQ.0)THEN
C            PCUT=P(I)
C            ICUT=1
C          ENDIF
C         ENDDO
C         DO I=1,NPRO
C          IF(P(I).LT.PCH4)THEN
C           IF(XCH4.LT.X1(I))THEN
C            X1(I)=XCH4
C            PCUT=P(I)
C           ENDIF
C          ENDIF
C         ENDDO

         NP = 11
        ELSEIF(VARIDENT(IVAR,1).EQ.226)THEN
C         print*,'Two cloud layering'
         IPAR = -1
         NP = 8

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

C        Extra section for combined cloud/gas profile - Model 10.
C        **********************************************************
         IF(JSPEC.GT.0.AND.JSPEC.LE.NCONT)THEN
          DO I=1,NPRO
           CONT(JSPEC,I)=X2(I)
          ENDDO
         ENDIF
C        **********************************************************

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


C     Now make sure the resulting VMRs add up to 1.0 for an
C     AMFORM=1 profile
      IF(AMFORM.EQ.1)THEN
        CALL ADJUSTVMR(NPRO,NVMR,VMR,ISCALE,IERR)
      ENDIF


C     ********  Modify profile with hydrostatic equation ********
      IF(JHYDRO.EQ.0)THEN
       CALL XHYDROSTATH(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1  IDGAS,ISOGAS,H,P,T,VMR,SCALE)
      ELSE
       CALL XHYDROSTATP(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1  IDGAS,ISOGAS,H,P,T,VMR,HTAN,PTAN,SCALE)
      ENDIF


C     ********* Check to see if anything should saturate *************

C     Default is that nothing condenses, even if it should. Set flags
C     accordingly
      DO J=1,NGAS
       ISWITCH(J) = 0
       VP(J)      = 1.0
       SVPFLAG(J) = 0
      ENDDO

C     First see if a list of gases to be forced to condense exists
      CALL FILE(IPFILE,IPFILE,'vpf')
      INQUIRE(FILE=IPFILE,EXIST=VPEXIST)

      IF(VPEXIST)THEN

       PRINT*,'Reading in svp flags from : ',IPFILE
       OPEN(12,FILE=IPFILE,STATUS='OLD')
        READ(12,*)NVP
        DO I=1,NVP
         READ(12,*)IP1,IP2,VP1,SVPFLAG1
         DO J=1,NVMR
           IF(IDGAS(J).EQ.IP1.AND.ISOGAS(J).EQ.IP2) THEN
             IF (SVPFLAG1.GT.0) THEN
               ISWITCH(J) = 1
             ENDIF
             VP(J)      = VP1
             SVPFLAG(J) = SVPFLAG1
           ENDIF
         ENDDO

        ENDDO
       CLOSE(12)
       
      ENDIF

c  ** loop around gases **
      DO 233 IGAS=1,NVMR

C      Can we condense this gas if we want to?
       IF(ISWITCH(IGAS).EQ.1)THEN

c        look for SVP data on this gas
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

c        if gas SVP data exists then go ahead
         IF(IDAT.EQ.1)THEN
 	     JSWITCH=0
c          loop over all pressure levels
           DO J=1,NPRO
C            Calculate SVP multiplied by required relative humidity VP
C		 Also calculate partial pressure of gas in question PP
             SVP=VP(IGAS)*DPEXP(A+B/T(J)+C*T(J)+D*T(J)*T(J))
	       PP=VMR(J,IGAS)*P(J)
	      
             IF(PP.GT.SVP)THEN
               IF(JSWITCH.EQ.0)THEN
                 PRINT*,'Subprofretg: gas predicted to condense'
                 PRINT*,'setting to SVP x VP'
                 PRINT*,IDGAS(IGAS),ISOGAS(IGAS)
               ENDIF 
               VMR(J,IGAS)=SVP/P(J)
               JSWITCH=1
             ENDIF

c ** NB technically xmap should be set=0.0 if condensation occurs
c however, this hard limit leads to undesirable retrieval behaviour such as:
c 1: if gas drops just below condensation on one iteration it can never return
c 2: leads to a sharp edge and erratic retrieval behaviour
c therefore, to  solve this XMAP is scaled by svp/pp. this also has a steep dropoff
c but gives a more gentle response and more desirable retrieval behaviour.
	        IF(JSWITCH.EQ.1)THEN
	          DO IX=1,NX
cc Neither of these work very well. best to leave unaltered (NAT)
cc		      XMAP(IX,IGAS,J)=0.0
cc		      XMAP(IX,IGAS,J)=XMAP(IX,IGAS,J)*SVP/PP
	          ENDDO
	        ENDIF

           ENDDO
           
c ** inspect the resulting profile for this gas **
c ** allow for presence of a cold trap if condensation occurs anywhere on the profile **
	     IF (JSWITCH.EQ.1) THEN
	       IF (SVPFLAG(IGAS).EQ.1) THEN
c		 * don't bother applying a cold trap, use simple condensation only *
	       ELSE IF (SVPFLAG(IGAS).EQ.2) THEN
c		 * assume gas is sourced from interior and has no local minima 
c		   ie the gas can only decrease with increasing altitude *
		   DO J=2,NPRO
		     IF (VMR(J,IGAS).GT.VMR(J-1,IGAS))THEN
		       VMR(J,IGAS)=VMR(J-1,IGAS)
		     ENDIF
		   ENDDO
	       ELSE IF (SVPFLAG(IGAS).EQ.3) THEN
c		 * assume gas is sourced from upper atmosphere 
c		   and has no local minima for levels deeper than the 0.05atm pressure level *
c		   ie gas can only decrease with decreasing altitude for pressures higher than 0.05atm *
		   DO J=NPRO-1,1,-1
		     IF (VMR(J,IGAS).GT.VMR(J+1,IGAS)) THEN
c                ** NB don't apply above far above tropopause as some photochemical profiles
c			  have local minima and maxima at high altitude**		     
		       IF ( P(J) .GE. 0.05 ) THEN
		         VMR(J,IGAS)=VMR(J+1,IGAS)
		       ENDIF
		     ENDIF
		   ENDDO
		 ENDIF 
	     ENDIF
	     
	     
         ENDIF

       ENDIF

233   ENDDO
c  ** end of loop around gases **
      
C     ************* Write out modified profiles *********

      CALL FILE(IPFILE,IPFILE,'prf')
      OPEN(UNIT=2,FILE=IPFILE,STATUS='UNKNOWN',ERR=52)
      WRITE(2,*)AMFORM
      IF(AMFORM.EQ.1)THEN
       WRITE(2,501)IPLANET,LATITUDE,NPRO,NVMR
      ELSE
       WRITE(2,501)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
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
c        print*, 'AERO:', H(I), CONT(1,I)
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

      IF(ISCAT.GT.0)THEN
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

