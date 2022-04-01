      SUBROUTINE SUBPROFRETG(XFLAG,IPFILE,ISPACE,ISCAT,GASGIANT,XLAT,
     1 XLON,NVAR,VARIDENT,VARPARAM,NX,XN,JPRE,NCONT,FLAGH2P,XMAP,IERR)
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
C	XLON	REAL		Longitude of spectrum to be simulated
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
C				the profile vectors (i.e., temperature, vmr prf
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

      integer, intent(in) :: xflag,ispace,iscat,nvar,nx,jpre
      integer, intent(in) :: varident(mvar,3),flagh2p
      real, intent(in) ::varparam(mvar,mparam),xlat,xlon
      character (len=100) :: ipfile,aname,buffer
      logical, intent(in) :: gasgiant

      real, intent(out) :: xmap(maxv,maxgas+2+maxcon,maxpro)
      integer, intent(out) :: ncont,ierr

      REAL XN(MX),DPEXP,DELH,XFAC,DXFAC,XTMP,SUM,PMIN,XSTEP
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS)
      REAL CONT(MAXCON,MAXPRO),X,XREF(MAXPRO),X1(MAXPRO)
      REAL PKNEE,HKNEE,XDEEP,XFSH,PARAH2(MAXPRO),XH,XKEEP,X2(MAXPRO)
      REAL HKNEE1,XDEEP1,XWID1,XT(MAXPRO),XTC(MAXPRO)
      REAL RHO,F,DQDX(MAXPRO),DX,PLIM,XFACP,CWID,PTROP, NEWF
      REAL DNDH(MAXPRO),DQDH(MAXPRO),FCLOUD(MAXPRO),HTOP,PTOP
      REAL dtempdx(MAXPRO,5),T0,Teff,alpha,ntemp,tau0,QC(MAXPRO)
      REAL XP1(MAXPRO),LP1(MAXPRO),XP2(MAXPRO),GRAD1
      REAL LPMIN,LPMAX,DLP,XPS(MAXPRO),XP2S(MAXPRO)
      REAL CTAU(MAXPRO),SETXLAPSE,XLAPSE,U,XF,TB,XCORR
      DOUBLE PRECISION Q(MAXPRO),OD(MAXPRO),ND(MAXPRO),XOD
      INTEGER ISCALE(MAXGAS),NPVAR,MAXLAT,ILAPSE
      INTEGER NLATREF,ILATREF,JLAT,KLAT,ICUT,JX,JLEV,ILEV
      PARAMETER (MAXLAT=20)
      REAL HREF(MAXLAT,MAXPRO),TREF(MAXLAT,MAXPRO),FLAT
      REAL PREF(MAXLAT,MAXPRO),VMRREF(MAXLAT,MAXPRO,MAXGAS)
      REAL LATREF(MAXLAT),MOLWTREF(MAXLAT),XFSHREF
      REAL XRH,XCDEEP,P1,PS,PS1,PH,Y1,Y2,YY1,YY2,ODX
      real plog,p1log,p2log,p2,v1log,v2log,grad,p3
      REAL XCH4,PCH4,PCUT,GETRADIUS,RPARTICLE,SHAPE
      INTEGER ICLOUD(MAXCON,MAXPRO),NCONT1,JSPEC,IFLA,I1,I2
      INTEGER NPRO,NPRO1,NVMR,JZERO,IV,IP,IVAR,JCONT,JVMR
      INTEGER IDGAS(MAXGAS),ISOGAS(MAXGAS),IPAR,JPAR,IVMR,NP
      INTEGER IDAT,JFSH,JHYDRO,ICOND,JFSH1
      REAL HTAN,PTAN,R,REFRADIUS,XCOL,ETOP
      REAL GASDATA(20,5),XVMR(MAXGAS),XMOLWT,XXMOLWT(MAXPRO)
      REAL CALCMOLWT,XRHO(MAXPRO)
      REAL RADCOND,DENSCOND,MWCOND
      INTEGER JSWITCH,ITEST,ICL,KJ
      INTEGER I,J,K,N,AMFORM,IPLANET,NGAS,IGAS,NXTEMP,IX
      REAL TEMP,RADIUS,G
      CHARACTER*8 PNAME
      REAL A,B,C,D,SVP,PP,LATITUDE,LONGITUDE,LAT1,LON1
      REAL V1(3),V0(3),XP,DLONG,FLONG,LONGITUDE1
      REAL DLAT,GRADTOUT(MAXPRO,5),TINT
      INTEGER ILONG,JLONG,ILAT
      REAL MOLWT,SCALE(MAXPRO),XMAP1,SCALEH,DTR
      PARAMETER (DTR=PI/180.)

      INTEGER NLOCATE,J1,MLON,MTHET,MLAT,NLAT
      PARAMETER(MLON=20,MTHET=3*MLON+1,MLAT=20)
      REAL YTH,YYTH(MTHET),YYTH2(MTHET),FTH,FI,FJ
      REAL CPHI1,CPHI2
      REAL YLONG(MLON),YLAT(MLAT)
      REAL XWID,Y,Y0,XX,L1,THETA(MTHET),GRADTMP(MAXPRO,MX)
      LOGICAL FEXIST,VPEXIST,NTEST,ISNAN,FVIVIEN
      REAL VP(MAXGAS),VP1,XS,GRADL(MAXPRO,MX)
      REAL GRADLINE(MAXPRO,MAXPRO)
      INTEGER SVPFLAG(MAXGAS),SVPFLAG1,NLONG,NTHETA,N1,N2
      INTEGER NVP,ISWITCH(MAXGAS),IP1,IP2,JKNEE,NLEVEL,JTOP
      REAL XLDEEP,XLHIGH,HVS,dlogp,XPC
      COMMON /SROM223/PCUT
      REAL XC1,XC2,RH,SLOPE,xch4new(maxpro),xch4newgrad(maxpro)
      REAL PD,RHC,RHM,VX,PT,FLUX,FRAIN
      INTEGER IMODEL
      REAL ALPHAP,BETAP,KIRP,GAMMAV1,GAMMAV2,TSTAR,RSTAR,SDIST
      REAL XLINE(4,MAXPRO),GLINE(4,MAXPRO,5),DXD1,DXD2,DXD3,DXD4
      REAL LONZ(4),TOUT(MAXPRO)
      INTEGER IX1,IX2,IX3,IX4,ILON
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      CALL RESERVEGAS

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
         if(idiag.gt.0)print*,'NLATREF = ',NLATREF
         IF(NLATREF.GT.MAXLAT)THEN
           print*,'MAXLAT in subprofretg is too small. Increase'
           print*,'and recompile'
           print*,'MAXLAT, NLATREF = ',MAXLAT,NLATREF
          STOP
         ENDIF

         DO 601 ILATREF=1,NLATREF
          IF(AMFORM.EQ.0)THEN
           READ(1,*)IPLANET,LATREF(ILATREF),NPRO,NVMR,MOLWTREF(ILATREF)
          ELSE
           READ(1,*)IPLANET,LATREF(ILATREF),NPRO,NVMR
          ENDIF

          IF(NPRO.GT.MAXPRO)THEN
            print*,'Error in subprofretg. NPRO>MAXPRO ',NPRO,MAXPRO
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
            IF(MCMCflag.eq.1)THEN
             HREF(ILATREF,I)=MCMCheight(I)
             TREF(ILATREF,I)=MCMCtemp(I)
             PREF(ILATREF,I)=MCMCpres(I)
             DO J=1, NVMR
              VMRREF(ILATREF,I,J)=MCMCvmr(I,J)
             ENDDO
            ENDIF
33        CONTINUE
601      CONTINUE

         

C        Now interpolate to correct latitude
         IF(NLATREF.EQ.1)THEN
          JLAT=1
          FLAT=0.
          if(idiag.gt.0)print*,'Snapping to first latitude'
         ELSE
          KLAT=-1
          DO ILATREF=1,NLATREF
           IF(XLAT.GE.LATREF(ILATREF))KLAT=ILATREF
          ENDDO
          IF(KLAT.LT.1)THEN
           if(idiag.gt.0)then
            print*,'Requested latitude is less than range given'
            print*,'Using lowest latitude available'
            print*,'Requested : ',XLAT
            print*,'Lowest available : ',LATREF(1)
           endif
           JLAT=1
           FLAT=0.
          ELSEIF(KLAT.EQ.NLATREF)THEN
           if(idiag.gt.0)then
            print*,'Requested latitude is greater than range given'
            print*,'Using highest latitude available'
            print*,'Requested : ',XLAT
            print*,'Highest available : ',LATREF(NLATREF)
           endif
           JLAT=NLATREF-1
           FLAT=1.0
          ELSE
           JLAT=KLAT
           FLAT=(XLAT-LATREF(JLAT))/
     &			(LATREF(JLAT+1)-LATREF(JLAT))
           if(idiag.gt.0)then
            print*,'JLAT,FLAT',JLAT,FLAT
            print*,'LATREF(JLAT),LATREF(JLAT+1)',
     &		LATREF(JLAT),LATREF(JLAT+1)
           endif
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
         LONGITUDE=XLON

         FVIVIEN=.FALSE.
C         FVIVIEN=.TRUE.

         IF(FVIVIEN)THEN

          CALL INTERPVIVIEN(XLAT,XLON,NPRO,NVMR,P,H,T,VMR)

          CALL XHYDROSTATH(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,
     1  MOLWT,IDGAS,ISOGAS,H,P,T,VMR,SCALE)

C          OPEN(12,FILE='TEST.OUT',STATUS='UNKNOWN')
C          WRITE(12,*)NPRO,NVMR
C          DO I=1,NPRO
C           WRITE(12,*)P(I),H(I),T(I),(VMR(I,J),J=1,NVMR)
C          ENDDO
C          CLOSE(12)

C          STOP


         ENDIF

       ELSE

         IF(AMFORM.EQ.0)THEN
          READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
         ELSE
          READ(1,*)IPLANET,LATITUDE,NPRO,NVMR
         ENDIF
         if(idiag.gt.0)print*,IPLANET,LATITUDE,NPRO,NVMR,MOLWT
         IF(NPRO.GT.MAXPRO)THEN
          if(idiag.gt.0)then
           print*,'Error in subprofretg. NPRO>MAXPRO ',NPRO,MAXPRO
          endif
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
            IF(MCMCflag.eq.1)THEN
             H(I)=MCMCheight(I)
             T(I)=MCMCtemp(I)
             P(I)=MCMCpres(I)
             DO J=1, NVMR
              VMR(I,J)=MCMCvmr(I,J)
             ENDDO
            ENDIF
30       CONTINUE

        CLOSE(UNIT=1)

       ENDIF


C      Make sure that vmrs add up to 1 if AMFORM=1
       IF(AMFORM.EQ.1)THEN
        if(idiag.gt.0)print*,'XX. ISCALE = ',(ISCALE(J),J=1,NVMR) 
        CALL ADJUSTVMR(NPRO,NVMR,VMR,ISCALE,IERR)
        
        IF(IERR.EQ.1)THEN
         if(idiag.gt.0)then
          print*,'XX. Warning from Adjustvmr: IERR = ',IERR
          print*,'Warning from subprofretg. VMRS do not add to 1'
          print*,'Resetting to reference'      
         endif 
         DO I=1,NPRO
          DO J=1,NVMR
           VMR(I,J)=(1.0-FLAT)*VMRREF(JLAT,I,J)+
     &		FLAT*VMRREF(JLAT+1,I,J)
          ENDDO
         ENDDO
         CALL ADJUSTVMR(NPRO,NVMR,VMR,ISCALE,IERR)
         if(idiag.gt.0)print*,'XX. IERRX = ',IERR
        ENDIF

        DO 301 I=1,NPRO
         DO K=1,NVMR
          XVMR(K)=VMR(I,K)
c          if(idiag.gt.0)print*,I,K,XVMR(K)
         ENDDO
         XXMOLWT(I)=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
301     CONTINUE
c      Calculate MOLWT but dont add VMRs to 1 if AMFORM=2

       ELSEIF(AMFORM.EQ.2)THEN

        DO I=1,NPRO
         DO K=1,NVMR
          XVMR(K)=VMR(I,K)
c          if(idiag.gt.0)print*,I,K,XVMR(K)
         ENDDO
         XXMOLWT(I)=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
        ENDDO
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
       if(idiag.gt.0)print*,'Calling xhydrostath'
       if(idiag.gt.0)then
        do i=1,npro
         print*,'i, P(i),H(i) = ',i,P(i),H(i)
        enddo
       endif
C       print*,AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT
       CALL XHYDROSTATH(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1  IDGAS,ISOGAS,H,P,T,VMR,SCALE)
       if(idiag.gt.0)then
        do i=1,npro
         print*,'Mod: i, P(i),H(i) = ',i,P(i),H(i)
        enddo
       endif
      ELSE
       if(idiag.gt.0)print*,'Calling xhydrostatp'
       if(idiag.gt.0)print*,'P(1),H(1) = ',P(1),H(1)
       CALL XHYDROSTATP(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1  IDGAS,ISOGAS,H,P,T,VMR,HTAN,PTAN,SCALE)
       if(idiag.gt.0)print*,'Mod: P(1),H(1) = ',P(1),H(1)
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

C     Calculate atmospheric density
      DO 307 I = 1,NPRO
C      Calculate density of atmosphere (g/cm3)
       IF(AMFORM.EQ.0)THEN
          XMOLWT=MOLWT
       ELSE
          XMOLWT=XXMOLWT(I)
       ENDIF
       XRHO(I) = P(I)*0.1013*XMOLWT/(R*T(I))
307   CONTINUE


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
          print*,'Error in subprofretg'
          print*,'Para-H2 profile has wrong number of levels'
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
          print*,'Error in subprofretg'
          print*,'fcloud profile has wrong number of levels'
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
C      WRITE(*,*)'SUBPROFRETG: reading saturated-vapour-pressure'
C      WRITE(*,*)'  data from ',ANAME
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
       JCONT=-1
       JSPEC=-1
       JVMR=-1
       IPAR=-1
       IF(VARIDENT(IVAR,1).LE.100)THEN
        IF(VARIDENT(IVAR,1).EQ.0)THEN
C        variable is Temperature
         DO I=1,NPRO
          XREF(I)=T(I)
         ENDDO
         IPAR = NVMR+1

        ELSEIF(VARIDENT(IVAR,1).LT.0)THEN
C        variable is aerosol amount
         JCONT = -VARIDENT(IVAR,1)
         IF(JCONT.GT.NCONT+2)THEN
          print*,'Error in subprofretg, JCONT > NCONT+2',JCONT,NCONT+2
          STOP
         ENDIF
C        Note if JCONT = NCONT+1 then profile contains para-H2 fraction
C        Note if JCONT = NCONT+2 then profile contains fraction cloud  cover
         IF(JCONT.EQ.NCONT+1)THEN
          IF(FLAGH2P.EQ.1)THEN
           DO I=1,NPRO
            XREF(I)=PARAH2(I)
           ENDDO
          ELSE
           print*,'Error in subprofretg, para-H2 fraction declared as'
           print*,'variable but atmosphere is not Giant Planet.'
           STOP
          ENDIF
         ELSEIF(JCONT.EQ.NCONT+2)THEN
          DO I=1,NPRO
            XREF(I)=FCLOUD(I)
          ENDDO
         ELSE
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
          print*,'Subprofretg: Gas could not be found',
     1     VARIDENT(IVAR,1),VARIDENT(IVAR,2)
          STOP
         ENDIF
C        Set ISCALE=0 for this gas to prevent vmr being scaled to give a 
C        total sum or vmrs=1 for AMFORM=1 format profile
         ISCALE(JVMR)=0
         DO I=1,NPRO
          XREF(I)=VMR(I,JVMR)
         ENDDO
         IPAR = JVMR

        ENDIF


C       Look up number of parameters needed to define this type of profile
        NP = NPVAR(VARIDENT(IVAR,3),NPRO,VARPARAM(IVAR,1))
        if(idiag.gt.0)print*,'IPAR = ',IPAR
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

        ELSEIF(VARIDENT(IVAR,3).EQ.-1)THEN
C        Model -1. Continuous particles/cm3 profile
C        ***************************************************************
         DO I=1,NPRO
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
             XMAP(NXTEMP+I,IPAR,I)=X1(I)/XRHO(I)
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
         if(idiag.gt.0)print*,NXTEMP+1,XN(NXTEMP+1)
         DO J = 1,NPRO
          X1(J) = XREF(J)*XN(NXTEMP+1)
          if(idiag.gt.0)print*,x1(j)
          XMAP(NXTEMP+1,IPAR,J)=XREF(J)
         ENDDO
	
        ELSEIF(VARIDENT(IVAR,3).EQ.3)THEN
C        Model 3. Profile is scaled fraction of reference profile, but 
C         code uses log scaling, which is more robust
C        ***************************************************************


         DO J = 1,NPRO
C          IF(XREF(J).GT.(1.0E37/EXP(XN(NXTEMP+1))).AND.XREF(J).GT.0.0)
C     1   	   XN(NXTEMP+1)=LOG(1.0E37/XREF(J))
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

         print*,'VARIDENT(IVAR,3).EQ.5 - Model no longer supported'
         STOP

        ELSEIF(VARIDENT(IVAR,3).EQ.6)THEN
C        Model 6. Venus-type cloud profile represented by a fixed base height, 
C        variable integrated optical depth above that height and variable cloud
C        scale height (km)
C        ***************************************************************
         HKNEE = VARPARAM(IVAR,1)
         IF(VARIDENT(IVAR,1).GE.0)THEN
          print*,'Warning from SUBPROFRETG. You are using a Venusian'
          print*,'cloud profile parameterisation for a non-cloud'
          print*,'variable'         
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

C         print*,'xdeep,xfsh = ',xdeep,xfsh

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

C         print*,pknee,jfsh,p(jfsh)

         IF(JFSH.LT.2.OR.JFSH.GT.NPRO-1)THEN
          print*,'SUBPROFRETG. Must choose pressure level'
          print*,'within range of profile'
          STOP
         ENDIF


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


         IF(VARIDENT(IVAR,1).GE.0)THEN
          print*,'Warning from SUBPROFRETG. You are using a'
          print*,'cloud profile parameterisation for a non-cloud'
          print*,'variable'         
          STOP 
         ENDIF

C        Calculate gradient numerically as it's just too hard otherwise
         DO 27 ITEST=1,4


          XDEEP = EXP(XN(NXTEMP+1))
          XFSH  = EXP(XN(NXTEMP+2))
          PKNEE = EXP(XN(NXTEMP+3))

          DX=0.05*XN(NXTEMP+ITEST-1)
          IF(DX.EQ.0.)DX=0.1

          IF(ITEST.EQ.2)THEN
            XDEEP=EXP(XN(NXTEMP+1)+DX)
          ENDIF
          IF(ITEST.EQ.3)THEN
            XFSH  = EXP(XN(NXTEMP+2)+DX)
          ENDIF
          IF(ITEST.EQ.4)THEN
            PKNEE = EXP(XN(NXTEMP+3)+DX)
          ENDIF


C         Start ND,Q,OD at zero
C         N is in units of particles/cm3
C         OD is in units of particles/cm2 = particles/cm3 x length(cm)
C         Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
          DO J=1,NPRO
           ND(J)=0.
           OD(J)=0
           Q(J)=0.
          ENDDO

          JFSH=-1

          CALL VERINT(P,H,NPRO,HKNEE,PKNEE)

          if(MCMCflag.eq.1)then
           if(HKNEE.ge.H(NPRO))THEN
            PKNEE = P(NPRO)
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
           ENDIF

          ENDDO

C         Integrate optical thickness
          OD(NPRO)=ND(NPRO)*SCALE(NPRO)*XFSH*1E5
          JFSH=-1
          DO J=NPRO-1,1,-1
           DELH = H(J+1) - H(J)
           XFAC = SCALE(J)*XFSH         
           OD(J)=OD(J+1)+(ND(J) - ND(J+1))*XFAC*1E5
           IF(H(J).LT.HKNEE.AND.JFSH.LT.0)THEN
              OD(J)=OD(J+1)+XFAC*1E5
              JFSH=1.
           ENDIF
          ENDDO
          XOD=OD(1)

C         The following section was found not to be as accurate as
C         desired due to misalignments at boundaries and so needs some 
C         post-processing in gsetrad.f
          DO J=1,NPRO
           OD(J)=XDEEP*OD(J)/XOD
           ND(J)=XDEEP*ND(J)/XOD
           Q(J)=XDEEP*Q(J)/XOD
           IF(Q(J).GT.1e10)Q(J)=1e10
           IF(Q(J).LT.1e-36)Q(J)=1e-36
           NTEST=ISNAN(Q(J))
           IF(NTEST)THEN
            if(idiag.gt.0)then
             print*,'Error in subprofretg.f, cloud density is NAN'
             print*,'Setting to 1e-36'
            endif
	    Q(J)=1e-36
           ENDIF

           IF(ITEST.EQ.1)THEN
            X1(J)=Q(J)
           ELSE
            XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
           ENDIF

          ENDDO

27       CONTINUE

        ELSEIF(VARIDENT(IVAR,3).EQ.9)THEN
C        Model 9. Similar to model 8. Profile is represented by a value 
C        at a variable HEIGHT, rather than PRESSURE plus a fractional scale height. 
C        Below the reference height the profile is set to zero. In addition,
C        this model scales the profile to give the requested integrated cloud 
C        optical depth.
C        ***************************************************************
         IF(VARIDENT(IVAR,1).GE.0)THEN
          print*,'Warning from SUBPROFRETG. You are using a'
          print*,'cloud profile parameterisation for a non-cloud'
          print*,'variable'         
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
C            if(idiag.gt.0)print*,'ITEST,IPAR,XN,DX = ',ITEST,IPAR,
C     1		XN(NXTEMP+ITEST-1),DX
C            if(idiag.gt.0)print*,'XDEEP,XFSH,HKNEE',XDEEP,XFSH,HKNEE
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

          if(MCMCflag.eq.1)then
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
            if(idiag.gt.0)then
             print*,'Error in subprofretg.f, cloud density is NAN'
             print*,'Setting to 1e-36'
            endif
            Q(J)=1e-36
           ENDIF

           IF(ITEST.EQ.1)THEN
            X1(J)=Q(J)
           ELSE
            XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
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

         if(idiag.gt.0) print*,'Mod10, IPAR, IDGAS(IPAR) = ',
     & IPAR, IDGAS(IPAR)
         IF(IDAT.EQ.0)THEN
          if(idiag.gt.0)then
           print*,'Subprofretg: Gas SVP data cannot be found'
           print*,IPAR,IDGAS(IPAR)
          endif
         ENDIF


C        Find where the gas will condense and if so put cloud there with required FSH.
C        Set JSPEC to required cloud ID
         JSPEC=ABS(INT(VARPARAM(IVAR,1)))

         JPAR = NVMR+1+JSPEC
         if(idiag.gt.0) print*,'Mod10, JPAR = ', JPAR

         DO I=1,NPRO
C         Preset profiles and gradients to prevent numerical instability later
          X1(I)=1e-36
          XMAP(NXTEMP+1,IPAR,I)=X1(I)
          XMAP(NXTEMP+2,IPAR,I)=X1(I)
          XMAP(NXTEMP+3,JPAR,I)=X1(I)
          XMAP(NXTEMP+4,JPAR,I)=X1(I)
         ENDDO 

C        Calculate gradient numerically as it's just too hard otherwise
         DO 302 ITEST=1,5

           XDEEP = EXP(XN(NXTEMP+1))
           XRH  = EXP(XN(NXTEMP+2))
           XCDEEP = EXP(XN(NXTEMP+3))
           XFSH  = EXP(XN(NXTEMP+4))

C           DX=0.05*XN(NXTEMP+ITEST-1)
C           IF(DX.EQ.0.)DX=0.1
           DX=0.1

C           IF(ITEST.GT.1)THEN
C            if(idiag.gt.0)print*,'ITEST,IPAR,XN,DX = ',ITEST,IPAR,
C     1         XN(NXTEMP+ITEST-1),DX
C            if(idiag.gt.0)print*,'XDEEP,XFSH,HKNEE',XDEEP,XFSH,HKNEE
C           ENDIF

           IF(ITEST.EQ.2)THEN
             XDEEP=EXP(XN(NXTEMP+1)+DX)
           ENDIF
           IF(ITEST.EQ.3)THEN
             XRH  = EXP(XN(NXTEMP+2)+DX)
           ENDIF
           IF(ITEST.EQ.4)THEN
             XCDEEP = EXP(XN(NXTEMP+3)+DX)
           ENDIF
           IF(ITEST.EQ.5)THEN
             XFSH = EXP(XN(NXTEMP+4)+DX)
           ENDIF

           IFLA=0
           HKNEE=0.
           DO I=1,NPRO
            XT(I)=1e-36

            P1=P(I)*XDEEP
            PS=DPEXP(A+B/T(I)+C*T(I)+D*T(I)*T(I))
            PH = PS*XRH
            IF(P1.LT.PS)THEN
             XT(I)=XDEEP
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

             XT(I)=PH/P(I)
             IF(XT(I).LT.1E-36)XT(I)=1E-36
            ENDIF

C           This section limits the vmr to relative RH at
C           all altitudes, not just above the condensation level. HKNEE has 
C           already been set

            IF(P1.GT.PH)THEN
             XT(I)=PH/P(I)
             IF(XT(I).LT.1E-36)XT(I)=1E-36
            ELSE
             XT(I)=XDEEP
            ENDIF

C           Now make sure that vmr does not rise again once condensation has
C           begun. i.e. freeze vmr at the cold trap.

            IF(IFLA.EQ.1.AND.XT(I).GT.XT(I-1))THEN
             XT(I)=XT(I-1)
            ENDIF

            IF(ITEST.EQ.1)THEN
             X1(I)=XT(I)
            ELSE
             IF(ITEST.LT.4)XMAP(NXTEMP+ITEST-1,IPAR,I)=(XT(I)-X1(I))/DX
            ENDIF

           ENDDO

C          Now put a cloud at the condensation level

           IF(AMFORM.EQ.0)THEN
            XMOLWT=MOLWT
           ELSE
            XMOLWT=XXMOLWT(NPRO)
           ENDIF

C          Calculate density of atmosphere (g/cm3)
           RHO = P(NPRO)*0.1013*XMOLWT/(R*T(NPRO))

C          Start ND(NPRO) at a random value. Will be rescaled anyway
           ND(NPRO)=1e-35
C          OD is in units of particles/cm2 = particles/cm3 x length(cm)
C          In this case this is the scale height at the top level.
           OD(NPRO)=ND(NPRO)*SCALE(NPRO)*1E5
C          Q is specific density = particles/gram = particles/cm3 x g/cm3
           XTC(NPRO)=ND(NPRO)/RHO         
         
C           print*,'HKNEE = ',HKNEE
C           print*,H(NPRO),ND(NPRO),RHO,OD(NPRO),XTC(NPRO)

           DO J=NPRO-1,1,-1
            DELH = H(J+1)-H(J)
            XFAC = SCALE(J)*XFSH
             IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
            ELSE
             XMOLWT=XXMOLWT(J)
            ENDIF

C           Calculate density of atmosphere (g/cm3)
            RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))
            ND(J)=ND(J+1)*EXP(DELH/XFAC)

            IF(H(J).GE.HKNEE)THEN
             OD(J)=OD(J+1)+(ND(J) - ND(J+1))*XFAC*1E5
            ELSE
             ND(J)=0.0
             OD(J)=OD(J+1)
            ENDIF
            XTC(J)=ND(J)/RHO
 
C            print*,H(J),ND(J),RHO,OD(J),XTC(J)

           ENDDO
           XOD=OD(1)

C           print*,'XCDEEP,XOD,Ratio',XCDEEP,XOD,XCDEEP/XOD

C          The following section was found not to be as accurate as desired
C          due to misalignments at boundaries and so needs some post-processing in 
C          gsetrad.f
           DO J=1,NPRO
            OD(J)=XCDEEP*OD(J)/XOD
            ND(J)=XCDEEP*ND(J)/XOD
            XTC(J)=XCDEEP*XTC(J)/XOD
            IF(H(J).LT.HKNEE)THEN
             IF(H(J+1).GE.HKNEE)THEN
              XTC(J)=XTC(J)*(1.0 - (HKNEE-H(J))/(H(J+1)-H(J)))
             ELSE
              XTC(J) = 0.0
             ENDIF
            ENDIF
            IF(XTC(J).GT.1e10)XT(J)=1e10
            IF(XTC(J).LT.1e-36)XT(J)=1e-36

            IF(ITEST.EQ.1)THEN
             X2(J)=XTC(J)
C             print*,J,H(J),XTC(J)
            ELSE
             XMAP(NXTEMP+ITEST-1,JPAR,J)=(XTC(J)-X2(J))/DX
            ENDIF

           ENDDO

302      CONTINUE

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
          if(idiag.gt.0)then
           print*,'Subprofretg: Gas SVP data cannot be found'
           print*,IPAR,IDGAS(IPAR)
          endif
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
          print*,'SUBPROFRETG. Must choose pressure level'
          print*,'within range of profile'
          print*,'Model = 16'
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
          print*,'SUBPROFRETG. Must choose pressure level'
          print*,'within range of profile'
          STOP
         ENDIF

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
          print*,'SUBPROFRETG. Must choose pressure level'
          print*,'within range of profile'
          STOP
         ENDIF

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
          print*,'Warning from SUBPROFRETG. You are using a'
          print*,'cloud profile parameterisation for a non-cloud'
          print*,'variable'         
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
C            if(idiag.gt.0)print*,'ITEST,IPAR,XN,DX = ',ITEST,IPAR,
C     1		XN(NXTEMP+ITEST-1),DX
C            if(idiag.gt.0)print*,'XDEEP,XFSH,HKNEE',XDEEP,XFSH,HKNEE
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
            if(idiag.gt.0)then
             print*,'Error in subprofretg.f, cloud density is NAN'
             print*,'XDEEP,XOD,HKNEE = ',XDEEP,XOD,HKNEE
            endif
            DO JX =1,NPRO
             OD(JX)=XDEEP*OD(JX)/XOD
             ND(JX)=XDEEP*ND(JX)/XOD
             Q(JX)=XDEEP*Q(JX)/XOD
             if(idiag.gt.0)then
              print*,'JX,OD,ND,Q = ',JX,OD(JX),ND(JX),Q(JX)
             endif
             IF(H(JX).LT.HKNEE)THEN
              IF(H(JX+1).GE.HKNEE)THEN
               Q(JX)=Q(JX)*(1.0 - (HKNEE-H(JX))/(H(JX+1)-H(JX)))
              ELSE
               Q(JX) = 0.0
              ENDIF
              if(idiag.gt.0)then
               print*,'H(JX),H(JX+1),HKNEE',H(JX),H(JX+1),HKNEE
               print*,HKNEE-H(JX),H(JX+1)-H(JX)
               print*,(1.0 - (HKNEE-H(JX))/(H(JX+1)-H(JX)))
              endif
             ENDIF
             IF(Q(JX).GT.1e10)Q(JX)=1e10
             IF(Q(JX).LT.1e-36)Q(JX)=1e-36
            ENDDO
	    STOP
           ENDIF

           IF(ITEST.EQ.1)THEN
            X1(J)=Q(J)
           ELSE
            XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
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
          print*,'Warning from SUBPROFRETG. You are using a'
          print*,'cloud profile parameterisation for a non-cloud'
          print*,'variable'         
          STOP 
         ENDIF


C        Calculate gradient numerically as it's just too hard otherwise
         DO 24 ITEST=1,3

          XDEEP = EXP(XN(NXTEMP+1))
          HKNEE = XN(NXTEMP+2)
          XFSH  = VARPARAM(IVAR,1)

C         Need radius from associated 444/445 particle parameterisation
          ICL=-VARIDENT(IVAR,1)
          RPARTICLE = GETRADIUS(ICL,NVAR,VARIDENT,VARPARAM,XN,NPRO)
          IF(RPARTICLE.LT.0)THEN
           print*,'Particle size cannot be found'
           STOP
          ENDIF
          REFRADIUS=VARPARAM(IVAR,2)
          XFSH=XFSH*REFRADIUS/RPARTICLE

          DX=0.05*XN(NXTEMP+ITEST-1)
          IF(DX.EQ.0.)DX=0.1

          IF(ITEST.GT.1)THEN
C            if(idiag.gt.0)print*,'ITEST,IPAR,XN,DX = ',ITEST,IPAR,
C     1		XN(NXTEMP+ITEST-1),DX
C            if(idiag.gt.0)print*,'XDEEP,XFSH,HKNEE',XDEEP,XFSH,HKNEE
          ENDIF
          IF(ITEST.EQ.2)THEN
            XDEEP=EXP(XN(NXTEMP+1)+DX)
          ENDIF
          IF(ITEST.EQ.3)THEN
            HKNEE = XN(NXTEMP+3)+DX
          ENDIF




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
            if(idiag.gt.0)then
             print*,'Error in subprofretg.f, cloud density is NAN'
             print*,'Setting to 1e-36'
            endif
            Q(J)=1e-36
           ENDIF

           IF(ITEST.EQ.1)THEN
            X1(J)=Q(J)
           ELSE
            XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
           ENDIF

           DNDH(J)=DNDH(J)*XDEEP/XOD
           DQDH(J)=DQDH(J)*XDEEP/XOD
           DQDX(J)=Q(J)/XOD

          ENDDO

24       CONTINUE

        ELSEIF(VARIDENT(IVAR,3).EQ.22)THEN
C        Brown dwarf parameterised temperature profile.

         tau0=EXP(XN(NXTEMP+1))
         ntemp=EXP(XN(NXTEMP+2))
         Teff=EXP(XN(NXTEMP+3))
         alpha=EXP(XN(NXTEMP+4))
         T0=EXP(XN(NXTEMP+5))

         call tbrownrc(npro,p,T0,Teff,tau0,ntemp,alpha,x1,
     &	dtempdx)

        do J=1,npro
         do K=1,np
          XMAP(NXTEMP+K,IPAR,J)=dtempdx(J,K)*EXP(XN(NXTEMP+K))
         enddo
        enddo

        ELSEIF( (VARIDENT(IVAR,3).EQ.23) .or. 
     >          (VARIDENT(IVAR,3).EQ.26) )THEN
C        Model 23/26. 2 point vmr gradient (NAT)
c		Profile is defined by two (p,v) points, with a linear gradient (in log p)
c		in between. The low pressure point is at (p1,v1) and the high pressure point 
c		is at (p2,v2). 
c		23: Profile is constant above/below this gradient region (i.e.
c		p<p1 v=v1 and p>p2 v=v2.)
c		26: Profile is constant above this gradient region and zero below (i.e.
c		p<p1 v=v1 and p>p2 v=0.)
c		All variable are retrieved. 
c		Not yet fully implemented for T	
C        ***************************************************************
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           print*,'ERROR: not coded for temperature yet'
           stop
         ENDIF
         if(xn(nxtemp+2).gt.xn(nxtemp+4)) then
           if(idiag.gt.0)then
            print*,'Warning: p1>p2 non-single valued profile.'
            print*,'       : setting p2=p1'
           endif
           xn(nxtemp+4) = xn(nxtemp+2)
         endif
         v1log = xn(nxtemp+1)
         p1log = xn(nxtemp+2)
         v2log = xn(nxtemp+3)
         p2log = xn(nxtemp+4)
         p1 = exp(p1log)
         p2 = exp(p2log)
         grad = (v2log-v1log)/(p2log-p1log)
         DO J=1,NPRO
          plog = alog(p(j))
          if (p(j).le.p1) then
c         * low pressure, constant continuance of vmr at p1
            x1(j) = exp(v1log)
            XMAP(NXTEMP+1,IPAR,J)=exp(xn(nxtemp+1))
            XMAP(NXTEMP+2,IPAR,J)=0.0
            XMAP(NXTEMP+3,IPAR,J)=0.0
            XMAP(NXTEMP+4,IPAR,J)=0.0
          elseif (p(j).ge.p2) then
            if (VARIDENT(IVAR,3).EQ.23) then
c           * high pressure, constant continuance of vmr at p2
              x1(j) = exp(v2log)
              XMAP(NXTEMP+1,IPAR,J)=0.0
              XMAP(NXTEMP+2,IPAR,J)=0.0
              XMAP(NXTEMP+3,IPAR,J)=exp(xn(nxtemp+3))
              XMAP(NXTEMP+4,IPAR,J)=0.0
            else
c           * high pressure, vmr=0 for p>p2
              x1(j) = 0.0
              XMAP(NXTEMP+1,IPAR,J)=0.0
              XMAP(NXTEMP+2,IPAR,J)=0.0
              XMAP(NXTEMP+3,IPAR,J)=0.0
              XMAP(NXTEMP+4,IPAR,J)=0.0
            endif
          else
c         * linear interpolation in log pressure / log vmr *
            x1(j)=exp( v1log + grad*(plog-p1log) )
c         * d X1 /d log v1 *
            XMAP(NXTEMP+1,IPAR,J)=(1-(plog-p1log)/(p2log-p1log))*x1(j)
c         * d X1 /d log p1 *
            XMAP(NXTEMP+2,IPAR,J)=-grad*x1(j) +
     >        x1(j)*(plog-p1log)*(v2log-v1log)*(p2log-p1log)**(-2)
c         * d X1 /d log v2 *
            XMAP(NXTEMP+3,IPAR,J)=(  (plog-p1log)/(p2log-p1log))*x1(j)
c         * d X1 /d log p2 *
            XMAP(NXTEMP+4,IPAR,J)= -x1(j)*
     >        (plog-p1log)*(v2log-v1log)*(p2log-p1log)**(-2)
          endif
	    if(X1(J).LT.1e-36)X1(J)=1e-36
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.24)THEN
C        Model 24: Variable constant vmr below knee, variable constant vmr above knee, variable knee pressure.
C        ***************************************************************
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
           XSTEP = XN(NXTEMP+2)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
           XSTEP = EXP(XN(NXTEMP+2))
         ENDIF

         PKNEE = EXP(XN(NXTEMP+3))

C        Technically, this profile would be parametrised by a Heaviside step function (X1(P)=(XDEEP-XSTEP)*H(P-PKNEE)+XSTEP).
C        However, this would also have the effect of making XMAP(NXTEMP+3,IPAR,J) = 0 for all values of J,
C        as the derivative of X1(J)/PKNEE is a multiple of a Dirac delta centred on PKNEE.
C        Hence, we use the approximation H(P) = 0.5(1+tanh(HVS*P)) where we choose an arbitrary value of HVS = 20
C        to roughly optimise the Heaviside step function given the orders of magnitude of the other input parameters and the pressure grid.
         CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
         HVS = 20.0
         JFSH = 0
         DO J=1,NPRO
          X1(J)=0.5*(XDEEP-XSTEP)*(1+tanh(HVS*(P(J)-PKNEE)))+XSTEP
          IF(VARIDENT(IVAR,1).EQ.0)THEN
           XMAP(NXTEMP+1,IPAR,J)=0.5*(1+tanh(HVS*(P(J)-PKNEE)))
           XMAP(NXTEMP+2,IPAR,J)=1-0.5*(1+tanh(HVS*(P(J)-PKNEE)))
          ELSE
           XMAP(NXTEMP+1,IPAR,J)=XDEEP*(1+tanh(HVS*(P(J)-PKNEE)))/2
           XMAP(NXTEMP+2,IPAR,J)=XSTEP*(1-tanh(HVS*(P(J)-PKNEE)))/2
          ENDIF
          IF(HVS*(P(J)-PKNEE).GT.40)THEN!get rid of infinity errors (cosh goes to infinity very quickly)
           XMAP(NXTEMP+3,IPAR,J)=0
          ELSE
           XMAP(NXTEMP+3,IPAR,J)=(-HVS*PKNEE*(XDEEP-XSTEP)) / 
     1                            (2*cosh(HVS*(P(J)-PKNEE))**2)
          ENDIF

          IF(X1(J).LT.1e-36)X1(J)=1e-36

         ENDDO


        ELSEIF(VARIDENT(IVAR,3).EQ.25)THEN
C        Model 25: Continuous profile with fewer than NPRO vertical levels
C        ***************************************************************
C         LPMAX=ALOG(P(1))
C         LPMIN=ALOG(P(NPRO))
C         NP = INT(VARPARAM(IVAR,1))
C         DLP = (LPMAX-LPMIN)/FLOAT(NP-1)    
CC        Find tabulated pressures and profile
C         DO J=1,NP
C          LP1(J)=LPMIN + DLP*FLOAT(J-1)
CC         Reverse profile array which by convention goes up in the atmosphere
C          XP1(1+NP-J)=XN(NXTEMP+J)
C         ENDDO

         NP = INT(VARPARAM(IVAR,1))
         DO J=1,NP
C         Reverse profile array which by convention goes up in the atmosphere
          LP1(1+NP-J)=ALOG(VARPARAM(IVAR,J+2))!read in pressure grid
          XP1(1+NP-J)=XN(NXTEMP+J)!read in profile
         ENDDO

C        Fit a cubic spline to the points
C         CALL CSPLINE(LP1,XP1,NP,1e30,1e30,XP2)

C        Fit a cubic spline to the points, dealing with discontinuities where necessary
         k=0
         DO J=2,NP-1
          IF(ABS(xp1(j+1)-xp1(j))/ABS(lp1(j+1)-lp1(j)).GT.10)THEN
           k=j!mark location of discontinuity
          ENDIF
          IF(k.gt.0)EXIT
         ENDDO
         IF(k.gt.0)THEN
          CALL CSPLINE(LP1(1:k),XP1(1:k),k,1e30,1e30,XP2(1:k))
          CALL CSPLINE(LP1(k+1:np),XP1(k+1:np),np-k,1e30,
     1                  1e30,XP2(k+1:np))
C         Use a tanh approximation at the discontinuity to calculate the 2nd derivative (see model 27)
          HVS = 20
          dlogp = LP1(k) - ((LP1(k)+LP1(k+1))/2)
          XP2(k)=HVS*(XP1(k)-XP1(k+1))*tanh(HVS*dlogp)*
     1           (1/cosh(HVS*dlogp)**2)
          dlogp = LP1(k+1) - ((LP1(k)+LP1(k+1))/2)
          XP2(k+1)=HVS*(XP1(k)-XP1(k+1))*tanh(HVS*dlogp)*
     1           (1/cosh(HVS*dlogp)**2)
         ELSE
          CALL CSPLINE(LP1,XP1,NP,1e30,1e30,XP2)
         ENDIF

        
C        If straight line, set second derivatives to 0
         IF(ABS(EXP(XP1(2))-EXP(XP1(1))).LT.1e-36)XP2(1)=0
         IF(ABS(EXP(XP1(NP))-EXP(XP1(NP-1))).LT.1e-36)XP2(NP)=0
         DO J=2,NP-1
          IF(ABS(EXP(XP1(J))-EXP(XP1(J-1))).LT.1e-36)THEN
           IF(ABS(EXP(XP1(J))-EXP(XP1(J+1))).LT.1e-36)XP2(J)=0
          ENDIF
         ENDDO


         DO J=1,NPRO
          L1 = ALOG(P(J))
          CALL CSPLINT(LP1,XP1,XP2,NP,L1,XX)
          IF(VARIDENT(IVAR,1).EQ.0)THEN
           X1(J)=XX
          ELSE
           X1(J)=EXP(XX)
          ENDIF
          IF(X1(J).LT.1e-36)X1(J)=1e-36
         ENDDO

C        Numerical differentiation
         XPS=XP1
         DO I = 1,NP
           IF(VARIDENT(IVAR,1).EQ.0)THEN
            DX=2.0
           ELSE
            DX=0.1
           ENDIF
           XPS(1+NP-I)=XP1(1+NP-I)+DX
C          Fit a cubic spline to the points
           CALL CSPLINE(LP1,XPS,NP,1e30,1e30,XP2S)
           DO J=1,NPRO
            L1 = ALOG(P(J))
            CALL CSPLINT(LP1,XPS,XP2S,NP,L1,XX)
            IF(VARIDENT(IVAR,1).EQ.0)THEN
             X2(J)=XX
            ELSE
             X2(J)=EXP(XX)
            ENDIF
            IF(X2(J).LT.1e-36)X2(J)=1e-36
            XMAP(NXTEMP+I,IPAR,J)=(X2(J)-X1(J))/DX
           ENDDO
           XPS(1+NP-I)=XP1(1+NP-I)
         ENDDO

C        I thought that getting accurate gradients from a cubic spline 
C        interpolated curve was going to be well hard. So I did a linear 
C        interpolation for this part instead. This actually turned out to
C        be too inaccurate, but I have left the code here for reference.

C         DO J=1,NPRO
C          L1 = ALOG(P(J))          
C          I = 1+INT((NP-1)*(L1-LPMIN)/(LPMAX-LPMIN))
C          IF(I.LT.1)I=1
C          IF(I.GE.NP)I=NP-1
C          XX = (L1-LP1(I))/DLP
C          if(idiag.gt.0)print*,'B',J,P(J),L1,I,LP1(I),XX,1.0-XX
C          I1 = NP+1-I
C          I2 = NP+1-(I+1)
C          IF(VARIDENT(IVAR,1).EQ.0)THEN
C           XMAP(NXTEMP+I1,IPAR,J)=XX
C           XMAP(NXTEMP+I2,IPAR,J)=1.-XX
C          ELSE
C           XMAP(NXTEMP+I1,IPAR,J)=X1(J)*XX
C           XMAP(NXTEMP+I2,IPAR,J)=X1(J)*(1.-XX)
C          ENDIF
C         ENDDO
C         stop
C         IF(X1(J).LT.1e-36)X1(J)=1e-36

         IF(VARPARAM(IVAR,3).LT.P(1))THEN
          print*,'Warning subprofretg:'
          print*,'Maximum pressure in interpolated grid lower'
          print*,'than maximum pressure in jupiter.ref.'
          print*,'This could lead to infinity errors'
          stop
         ENDIF
         IF(VARPARAM(IVAR,NP+2).GT.P(NPRO))THEN
          print*,'Warning subprofretg:'
          print*,'Minimum pressure in interpolated grid higher'
          print*,'than minimum pressure in jupiter.ref.'
          print*,'This could lead to infinity errors'
          stop
         ENDIF

        ELSEIF(VARIDENT(IVAR,3).EQ.27)THEN
C        Model 27: Step profile. (Nick Teanby)
C		XN1 = log(xdeep)		[= xdeep if Temperature]
C		XN2 = log(xshallow)	[= xshallow if Temperature]
C		XN3 = log(pknee)
C 
C        similar to model 24 except transition is steeper, more step like, and is not sensitive to pressure level.
C        Implemented by taking difference of log(p) (instead of p in option 24)
C        ***************************************************************
         IF(VARIDENT(IVAR,1).EQ.0)THEN
           XDEEP = XN(NXTEMP+1)
           XSTEP = XN(NXTEMP+2)
         ELSE
           XDEEP = EXP(XN(NXTEMP+1))
           XSTEP = EXP(XN(NXTEMP+2))
         ENDIF
         PKNEE = EXP(XN(NXTEMP+3))

C        Technically, this profile would be parametrised by a Heaviside step function (X1(P)=(XDEEP-XSTEP)*H(P-PKNEE)+XSTEP).
C        However, this would also have the effect of making XMAP(NXTEMP+3,IPAR,J) = 0 for all values of J.
C        Hence, we use the approximation H(P) = 0.5(1+tanh(HVS*delta_LogP)) where we choose an arbitrary value of HVS = 10
C        to roughly approximate the Heaviside step function.
C	   If step is not steep enough then increase HVS

         HVS = 10.0
         DO J=1,NPRO
          dlogp = log(P(J)) - log(PKNEE)
          X1(J)=0.5*(XDEEP-XSTEP)*( 1+tanh(HVS*dlogp) )+XSTEP
          IF(VARIDENT(IVAR,1).EQ.0)THEN
           XMAP(NXTEMP+1,IPAR,J)=  0.5*(1+tanh(HVS*dlogp))
           XMAP(NXTEMP+2,IPAR,J)=1-0.5*(1+tanh(HVS*dlogp))
          ELSE
           XMAP(NXTEMP+1,IPAR,J)=0.5*XDEEP*(1+tanh(HVS*dlogp))
           XMAP(NXTEMP+2,IPAR,J)=XSTEP-0.5*XSTEP*(1+tanh(HVS*dlogp))
          ENDIF
          IF(HVS*dlogp.GT.20)THEN
           XMAP(NXTEMP+3,IPAR,J)= 0.0
          ELSE
           XMAP(NXTEMP+3,IPAR,J)= -0.5*(XDEEP-XSTEP)*HVS*
     >                                 (1.0-tanh(HVS*dlogp)**2)
          ENDIF

          IF(X1(J).LT.1e-36)X1(J)=1e-36

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.28)THEN
C        Model 28: Modify just one element of a continuous profile
C        ***************************************************************

c         Level in profile to be changed
          JLEV = VARPARAM(IVAR,1)

          DO I=1,NPRO

           IF(I.EQ.JLEV)THEN
            IF(VARIDENT(IVAR,1).EQ.0)THEN
             X1(I) = XN(NXTEMP+1)
             XMAP(NXTEMP+1,IPAR,I)=1.0
            ELSE
              IF(XN(NXTEMP+I).GT.-82.8931)THEN
                X1(I) = EXP(XN(NXTEMP+1))
              ELSE
                X1(I) = 1.0E-36
              ENDIF
              IF(XN(NXTEMP+I).LT.80.0)THEN
                X1(I) = EXP(XN(NXTEMP+1))
              ELSE
                X1(I) = EXP(80.0)
              ENDIF
              XMAP(NXTEMP+1,IPAR,I)=X1(I)
            ENDIF
           ELSE
            X1(I) = XREF(I) 
            XMAP(NXTEMP+1,IPAR,I)=0.0 
           ENDIF
          ENDDO
 
        ELSEIF(VARIDENT(IVAR,3).EQ.29)THEN
C        Model 29. Continuous profile at multiple locations
C        ***************************************************************

C        Need to find nearest entry to requested lat/long
         CALL XPROJ(LATITUDE,LONGITUDE,V0)

         NLOCATE=INT(VARPARAM(IVAR,1))
         J=2
         XX = -1000.
         DO I=1,NLOCATE
          LAT1=VARPARAM(IVAR,J)
          LON1=VARPARAM(IVAR,J+1)
          CALL XPROJ(LAT1,LON1,V1)
          XP=0.
          DO K=1,3
           XP=XP+V0(K)*V1(K)
          ENDDO
          IF(XP.GT.XX)THEN
           XX=XP
           I1=I
          ENDIF
          J=J+2 
         ENDDO
         DO I=1, NPRO
           J1 = NXTEMP+(I1-1)*NPRO+I
           IF(VARIDENT(IVAR,1).EQ.0)THEN
            X1(I) = XN(J1)
            XMAP(J1,IPAR,I)=1.0
           ELSE
             IF(XN(J1).GT.-82.8931)THEN
               X1(I) = EXP(XN(J1))
             ELSE
               X1(I) = 1.0E-36
             ENDIF
             IF(XN(J1).LT.80.0)THEN
               X1(I) = EXP(XN(J1))
             ELSE
               X1(I) = EXP(80.0)
             ENDIF
             XMAP(J1,IPAR,I)=X1(I)
           ENDIF
         ENDDO


        ELSEIF(VARIDENT(IVAR,3).EQ.30)THEN
C        Model 30. Inhomogenous disc profile at multiple locations
C        ***************************************************************

C        Need to interpolation in longitude for each vertical level
         NLONG = INT(VARPARAM(IVAR,1)/VARPARAM(IVAR,2)+0.1)
         NLEVEL = INT(VARPARAM(IVAR,2))
         if(idiag.gt.0)then
          print*,'Model 30 - nlong,nlevel,np = ',nlong,nlevel,np
          print*,'Model 30 - latitude,longitude = ',LATITUDE,
     &           LONGITUDE
         endif
         DLONG=360.0/FLOAT(NLONG)
         LONGITUDE1=LONGITUDE
         IF(LONGITUDE.LT.0.0)LONGITUDE1=LONGITUDE+360.
         IF(LONGITUDE.GE.360.0)LONGITUDE1=LONGITUDE-360.
         ILONG=1+INT(LONGITUDE1/DLONG)
         FLONG = (LONGITUDE1 - (ILONG-1)*DLONG)/DLONG
         if(flong.gt.1.0)then 
          print*,'Error: flong > 1.0'
          stop
         endif
         IF(ILONG.LT.NLONG)THEN
          JLONG=ILONG+1
         ELSE
          JLONG=1
         ENDIF

         if(idiag.gt.0)print*,'LONGITUDE1,ILONG,JLONG,FLONG',
     &    LONGITUDE1,ILONG,JLONG,FLONG

         IF(VARIDENT(IVAR,1).EQ.0)THEN
            DX=2.0
         ELSE
            DX=0.1
         ENDIF

C        Read in pressure grid from varparam
         DO J=1,NLEVEL
          LP1(J)=ALOG(VARPARAM(IVAR,J+2))
          DO K=1,NP
           GRADL(J,K)=0.
          ENDDO
         ENDDO
       
C        Set exponent of cos(lat) variation

         XPC = VARPARAM(IVAR,NLEVEL+3) 

         DO J=1,NLEVEL
          SUM=0.
          DO I=1,NLONG
           J1=NXTEMP+(I-1)*NLEVEL+J
           YLONG(I)=XN(J1)
           SUM=SUM+XN(J1)/FLOAT(NLONG)
          ENDDO

C         Alternate - setting pole temperature to average of morning and
C                     afternoon terminators.
          N1=NLONG/4
          N2=N1+2*N1        
          SUM=0.5*(YLONG(N1+1)+YLONG(N2+1))

          YTH = (1.0-FLONG)*YLONG(ILONG)+FLONG*YLONG(JLONG)
          XP1(J) = SUM + (YTH-SUM)*(COS(LATITUDE*DTR))**XPC

C         Still need to incorporate gradient of SUM

          K=(ILONG-1)*NLEVEL+J
          GRADL(J,K)=(1.0-FLONG)*(COS(LATITUDE*DTR))**XPC
          K=(JLONG-1)*NLEVEL+J
          GRADL(J,K)=FLONG*(COS(LATITUDE*DTR))**XPC

         ENDDO

C        Now need to interpolate local NLEVEL profile to NPRO profile
         DO J=1,NPRO

          L1 = ALOG(P(J))

          F=-1.
          DO JLEV=1,NLEVEL-1
           IF(L1.LE.LP1(JLEV).AND.L1.GT.LP1(JLEV+1))THEN
             ILEV=JLEV
             F = (L1-LP1(JLEV))/(LP1(JLEV+1)-LP1(JLEV))
             GOTO 111
            ENDIF
          ENDDO
111       CONTINUE
          IF(F.LT.0)THEN
            IF(L1.GT.LP1(1))THEN
             ILEV=1
             F = (L1-LP1(ILEV))/(LP1(ILEV+1)-LP1(ILEV))
            ELSE
             ILEV=NLEVEL-1
             F = (L1-LP1(ILEV))/(LP1(ILEV+1)-LP1(ILEV))
            ENDIF
            if(idiag.gt.0)then
             print*,'Model 30 warning - pressure out of range'
             print*,'Having to extrapolate'
             print*,L1,LP1(1),LP1(NLEVEL)
             print*,ILEV,LP1(ILEV),LP1(ILEV+1),F
            endif
          ENDIF

          XX = (1.0-F)*XP1(ILEV)+F*XP1(ILEV+1)

          IF(VARIDENT(IVAR,1).EQ.0)THEN
           X1(J)=XX
          ELSE
           X1(J)=EXP(XX)
          ENDIF
          IF(X1(J).LT.1e-36)X1(J)=1e-36
          
          DO K=1,NP
           GRADTMP(J,K)=(1.0-F)*GRADL(ILEV,K)+F*GRADL(ILEV+1,K)         
           IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+K,IPAR,J)=GRADTMP(J,K)
           ELSE
            XMAP(NXTEMP+K,IPAR,J)=X1(J)*GRADTMP(J,K)
           ENDIF
          ENDDO

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.31)THEN
C        Model 31. Inhomogenous disc scaling factor
C        ***************************************************************

C        Need to interpolation in longitude
         NLONG = INT(VARPARAM(IVAR,1))


         if(idiag.gt.0)then
          print*,'Model 31 - nlong,np = ',nlong,np
          print*,'Model 31 - latitude,longitude = ',LATITUDE,
     &           LONGITUDE
         endif
         LONGITUDE1=LONGITUDE
         IF(LONGITUDE.LT.0.0)LONGITUDE1=LONGITUDE+360.
         IF(LONGITUDE.GE.360.0)LONGITUDE1=LONGITUDE-360.
         DLONG=360.0/FLOAT(NLONG)
         ILONG=1+INT(LONGITUDE1/DLONG)
         FLONG = (LONGITUDE1 - (ILONG-1)*DLONG)/DLONG
         if(flong.gt.1.0)then 
          print*,'Error: flong > 1.0'
          stop
         endif
         IF(ILONG.LT.NLONG)THEN
          JLONG=ILONG+1
         ELSE
          JLONG=1
         ENDIF



         if(idiag.gt.0)print*,ILONG,JLONG,FLONG,NP

         IF(VARIDENT(IVAR,1).EQ.0)THEN
            DX=2.0
         ELSE
            DX=0.1
         ENDIF

         DO K=1,NP
           GRADL(1,K)=0.
         ENDDO

         SUM=0.
         DO I=1,NLONG
          J1=NXTEMP+I
          YLONG(I)=XN(J1)
          SUM=SUM+XN(J1)/FLOAT(NLONG)
         ENDDO


C        Set exponent of cos(lat) variation
         XPC=VARPARAM(IVAR,2)
         if(idiag.gt.0)print*,'XPC = ',XPC    


         YTH = (1.0-FLONG)*YLONG(ILONG)+FLONG*YLONG(JLONG)
         XS = SUM + (YTH-SUM)*(COS(LATITUDE*DTR))**XPC

         GRADL(1,ILONG)=(1.0-FLONG)*(COS(LATITUDE*DTR))**XPC
         GRADL(1,JLONG)=FLONG*(COS(LATITUDE*DTR))**XPC

         DO J=1,NPRO
           X1(J) = XREF(J)*EXP(XS)
           XMAP(NXTEMP+ILONG,IPAR,J)=X1(J)*GRADL(1,ILONG)
           XMAP(NXTEMP+JLONG,IPAR,J)=X1(J)*GRADL(1,JLONG)
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.32)THEN
C        ***************************************************************
C        Model 32. Profile is represented by a value at a variable pressure level
C        plus a fractional scale height. Below the knee pressure the profile is 
C        set to drop exponentially. Similar model to model 8.
C        ***************************************************************

         IF(VARIDENT(IVAR,1).GE.0)THEN
          print*,'Warning from SUBPROFRETG. You are using a'
          print*,'cloud profile parameterisation for a non-cloud'
          print*,'variable'         
          STOP 
         ENDIF

C        Calculate gradient numerically as it's just too hard otherwise
         DO 207 ITEST=1,4


          XDEEP = EXP(XN(NXTEMP+1))
          XFSH  = EXP(XN(NXTEMP+2))
          PKNEE = EXP(XN(NXTEMP+3))

          DX=0.05*XN(NXTEMP+ITEST-1)
          IF(DX.EQ.0.)DX=0.1

          IF(ITEST.EQ.2)THEN
            XDEEP=EXP(XN(NXTEMP+1)+DX)
          ENDIF
          IF(ITEST.EQ.3)THEN
            XFSH  = EXP(XN(NXTEMP+2)+DX)
          ENDIF
          IF(ITEST.EQ.4)THEN
            PKNEE = EXP(XN(NXTEMP+3)+DX)
          ENDIF


          CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
          if(idiag.gt.0)print*,pknee,hknee

C         Start ND,Q,OD at zero
C 	  N is in units of particles/cm3
C         OD is in units of particles/cm2 = particles/cm3 x length(cm)
C         Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
          DO J=1,NPRO
           ND(J)=0.
           OD(J)=0
           Q(J)=0.
          ENDDO

          JKNEE=-1
          XF=1.
C         find levels in atmosphere that span pknee
          DO J=1,NPRO-1
           IF(P(J).GE.PKNEE.AND.P(J+1).LT.PKNEE)THEN
            JKNEE=J
           ENDIF
          ENDDO
 
          IF(JKNEE.LT.0)THEN
           print*,'subprofretg: Error in model 32. Stop'
           print*,'IVAR,XDEEP,XFSH,PKNEE'
           print*,IVAR,XDEEP,XFSH,PKNEE
           STOP
          ENDIF          

          DELH=H(JKNEE+1)-HKNEE
          XFAC=0.5*(SCALE(JKNEE)+SCALE(JKNEE+1))*XFSH
          ND(JKNEE+1)=DPEXP(-DELH/XFAC)

          DELH = HKNEE-H(JKNEE)
          XFAC=XF
          ND(JKNEE)=DPEXP(-DELH/XFAC)
         
          DO J=JKNEE+2,NPRO
           DELH = H(J)-H(J-1)
           XFAC = SCALE(J)*XFSH
           ND(J) = ND(J-1)*DPEXP(-DELH/XFAC)
          ENDDO
          
          DO J=1,JKNEE-1
           DELH = H(JKNEE)-H(J)
           XFAC = XF
           ND(J) = DPEXP(-DELH/XFAC)
          ENDDO

          DO J=1,NPRO
            IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
            ELSE
             XMOLWT=XXMOLWT(J)
            ENDIF
C           Calculate density of atmosphere  (g/cm3)
            RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))
C           Calculate initial particles/gram
            Q(J)=ND(J)/RHO
          ENDDO

C         Now integrate optical thickness
          OD(NPRO)=ND(NPRO)*SCALE(NPRO)*XFSH*1E5
          JFSH=-1
          DO J=NPRO-1,1,-1
           IF(J.GT.JKNEE)THEN
             DELH = H(J+1) - H(J)
             XFAC = SCALE(J)*XFSH
             OD(J)=OD(J+1)+(ND(J) - ND(J+1))*XFAC*1E5
           ELSE
             IF(J.EQ.JKNEE)THEN
              DELH = H(J+1)-HKNEE
              XFAC = 0.5*(SCALE(J)+SCALE(J+1))*XFSH         
              OD(J)=OD(J+1)+(1. - ND(J+1))*XFAC*1E5
              DELH = HKNEE-H(J)
              XFAC = XF
              OD(J)=OD(J)+(1. - ND(J))*XFAC*1E5
             ELSE
              DELH = H(J+1)-H(J)
              XFAC = XF
              OD(J)=OD(J+1)+(ND(J+1) - ND(J))*XFAC*1E5
             ENDIF
           ENDIF
          ENDDO

          ODX=OD(1)

C         Now normalise specific density profile.
C         This is also redone in gsetrad.f to make this totally secure.
          DO J=1,NPRO
           OD(J)=OD(J)*XDEEP/ODX
           ND(J)=ND(J)*XDEEP/ODX
           Q(J)=Q(J)*XDEEP/ODX
           IF(Q(J).GT.1e10)Q(J)=1e10
           IF(Q(J).LT.1e-36)Q(J)=1e-36
           NTEST=ISNAN(Q(J))
           IF(NTEST)THEN
            if(idiag.gt.0)then
             print*,'Error in subprofretg.f, cloud density is NAN'
             print*,'Setting to 1e-36'
            endif
	    Q(J)=1e-36
           ENDIF

           IF(ITEST.EQ.1)THEN
            X1(J)=Q(J)
           ELSE
            XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
           ENDIF

          ENDDO

207       CONTINUE

        ELSEIF(VARIDENT(IVAR,3).EQ.33)THEN
C        Model 33. Inhomogenous disc profile at multiple locations
C        ***************************************************************

C        Need to interpolation in longitude for each vertical level
         NLONG = INT(VARPARAM(IVAR,1)/VARPARAM(IVAR,2)+0.1)
         NLEVEL = INT(VARPARAM(IVAR,2))
         if(idiag.gt.0)then
          print*,'Model 33 - nlong,nlevel,np = ',nlong,nlevel,np
          print*,'Model 33 - latitude,longitude = ',LATITUDE,
     &           LONGITUDE
         endif
         DLONG=360.0/FLOAT(NLONG)
         LONGITUDE1=LONGITUDE
         IF(LONGITUDE.LT.0.0)LONGITUDE1=LONGITUDE+360.
         IF(LONGITUDE.GE.360.0)LONGITUDE1=LONGITUDE-360.
         ILONG=1+INT(LONGITUDE1/DLONG)
         FLONG = (LONGITUDE1 - (ILONG-1)*DLONG)/DLONG
         if(flong.gt.1.0)then 
          print*,'Error: flong > 1.0'
          stop
         endif
         IF(ILONG.LT.NLONG)THEN
          JLONG=ILONG+1
         ELSE
          JLONG=1
         ENDIF

         if(idiag.gt.0)print*,'LONGITUDE1,ILONG,JLONG,FLONG',
     &    LONGITUDE1,ILONG,JLONG,FLONG

         IF(VARIDENT(IVAR,1).EQ.0)THEN
            DX=2.0
         ELSE
            DX=0.1
         ENDIF

C        Read in pressure grid from varparam
         DO J=1,NLEVEL
          LP1(J)=ALOG(VARPARAM(IVAR,J+2))
          if(idiag.gt.0)print*,J,VARPARAM(IVAR,J+2),LP1(J),EXP(LP1(J))
          DO K=1,NP
           GRADL(J,K)=0.
          ENDDO
         ENDDO
       
         DO J=1,NLEVEL
          SUM=0.
          DO I=1,NLONG
           J1=NXTEMP+(I-1)*NLEVEL+J
           YLONG(I)=XN(J1)
           SUM=SUM+XN(J1)/FLOAT(NLONG)
          ENDDO

          FI = XN(NXTEMP+NLONG*NLEVEL+ILONG)
          FJ = XN(NXTEMP+NLONG*NLEVEL+JLONG)
        
          YTH = (1.0-FLONG)*YLONG(ILONG)+FLONG*YLONG(JLONG)
          FTH = (1.0-FLONG)*FI+FLONG*FJ

          CPHI1=(COS(LATITUDE*DTR))**0.25
          CPHI2=(COS(LATITUDE*DTR))**2.0

          XP1(J) = SUM + (YTH-SUM)*(FTH*CPHI1+(1.0-FTH)*CPHI2)

          K=(ILONG-1)*NLEVEL+J
          GRADL(J,K)=(1.0-FLONG)*(FTH*CPHI1+(1.0-FTH)*CPHI2)

          K=(JLONG-1)*NLEVEL+J
          GRADL(J,K)=FLONG*(FTH*CPHI1 +(1.0-FTH)*CPHI2)

          GRADL(J,NLONG*NLEVEL+ILONG)=
     &     (YTH-SUM)*(CPHI1 - CPHI2)*(1.0-FLONG)
          GRADL(J,NLONG*NLEVEL+JLONG)=
     &     (YTH-SUM)*(CPHI1 - CPHI2)*FLONG

         ENDDO

C        Now need to interpolate local NLEVEL profile to NPRO profile

         DO J=1,NPRO

          L1 = ALOG(P(J))

          F=-1.
          DO JLEV=1,NLEVEL-1
           IF(L1.LE.LP1(JLEV).AND.L1.GT.LP1(JLEV+1))THEN
             ILEV=JLEV
             F = (L1-LP1(JLEV))/(LP1(JLEV+1)-LP1(JLEV))
             GOTO 151
            ENDIF
          ENDDO
151       CONTINUE
          IF(F.LT.0)THEN
            IF(L1.GT.LP1(1))THEN
             ILEV=1
             F = (L1-LP1(ILEV))/(LP1(ILEV+1)-LP1(ILEV))
            ELSE
             ILEV=NLEVEL-1
             F = (L1-LP1(ILEV))/(LP1(ILEV+1)-LP1(ILEV))
            ENDIF
            if(idiag.gt.0)then
             print*,'Model 33 warning - pressure out of range'
             print*,'Having to extrapolate'
             print*,L1,LP1(1),LP1(NLEVEL)
             print*,ILEV,LP1(ILEV),LP1(ILEV+1),F
            endif
          ENDIF

          XX = (1.0-F)*XP1(ILEV)+F*XP1(ILEV+1)

          IF(VARIDENT(IVAR,1).EQ.0)THEN
           X1(J)=XX
          ELSE
           X1(J)=EXP(XX)
          ENDIF
          IF(X1(J).LT.1e-36)X1(J)=1e-36
          
          DO K=1,NP
           GRADTMP(J,K)=(1.0-F)*GRADL(ILEV,K)+F*GRADL(ILEV+1,K)         
           IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+K,IPAR,J)=GRADTMP(J,K)
           ELSE
            XMAP(NXTEMP+K,IPAR,J)=X1(J)*GRADTMP(J,K)
           ENDIF
          ENDDO

         ENDDO


        ELSEIF(VARIDENT(IVAR,3).EQ.34)THEN
C        Model 34. Milne-Eddington Temperature Profile
C        ***************************************************************
         SETXLAPSE = VARPARAM(IVAR,1)

         DO J=1,NPRO
          CTAU(I)=0.
         ENDDO

         TB=XN(NXTEMP+1)
         XF=EXP(XN(NXTEMP+2))

         CTAU(NPRO)=P(NPRO)
         X2(NPRO)=TB*((2.0+XF*CTAU(NPRO))/4.0)**0.25 

         DO J=NPRO-1,1,-1
          CTAU(J)=CTAU(J+1)+P(J)
          X2(J)=TB*((2.0+XF*CTAU(J))/4.0)**0.25
         ENDDO

         U = (2.0+XF*CTAU(NPRO))/4.0
         X1(NPRO)=TB*U**0.25
         XMAP(NXTEMP+1,IPAR,NPRO)=X1(NPRO)/TB
         XMAP(NXTEMP+2,IPAR,NPRO)=0.0625*TB*(U**(-0.75))*
     &    CTAU(NPRO)*XF

         ILAPSE=1
         DO J=NPRO-1,1,-1
          DELH = H(J+1)-H(J)
          XLAPSE = (X2(J)-X2(J+1))/DELH
          IF(XLAPSE.GT.SETXLAPSE.OR.ILAPSE.EQ.0)THEN
           X1(J)=X1(J+1)+SETXLAPSE*DELH
           ILAPSE=1
           XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J+1)
           XMAP(NXTEMP+2,IPAR,J)=XMAP(NXTEMP+2,IPAR,J+1)            
          ELSE
           X1(J)=X2(J)
           U = (2.0+XF*CTAU(J))/4.0
           XMAP(NXTEMP+1,IPAR,J)=X1(J)/TB
           XMAP(NXTEMP+2,IPAR,J)=0.0625*TB*(U**(-0.75))*
     &    CTAU(J)*XF
          ENDIF
         ENDDO


        ELSEIF(VARIDENT(IVAR,3).EQ.35)THEN
C        Model 35. Axisymmetric model continuous
C        ***************************************************************

C        Need to interpolation in longitude for each vertical level
         NLAT = INT(VARPARAM(IVAR,1)/VARPARAM(IVAR,2)+0.1)
         NLEVEL = INT(VARPARAM(IVAR,2))
         if(idiag.gt.0)then
          print*,'Model 35 - nlat,nlevel,np = ',nlat,nlevel,np
          print*,'Model 35 - latitude,longitude = ',LATITUDE,
     &           LONGITUDE
         endif
         DLAT=180.0/FLOAT(NLAT-1)
         ILAT=1+INT((90.0+LATITUDE)/DLAT)
         FLAT = (90+LATITUDE - (ILAT-1)*DLAT)/DLAT
         if(flat.gt.1.0)then 
          if(idiag.gt.0)print*,'Error: flat > 1.0'
          stop
         endif
         IF(ILAT.LT.NLAT)THEN
          JLAT=ILAT+1
         ELSE
          stop
         ENDIF

         IF(VARIDENT(IVAR,1).EQ.0)THEN
            DX=2.0
         ELSE
            DX=0.1
         ENDIF

C        Read in pressure grid from varparam
         DO J=1,NLEVEL
          LP1(J)=ALOG(VARPARAM(IVAR,J+2))
          if(idiag.gt.0)print*,J,VARPARAM(IVAR,J+2),LP1(J),EXP(LP1(J))
          DO K=1,NP
           GRADL(J,K)=0.
          ENDDO
         ENDDO
       
         DO J=1,NLEVEL
          DO I=1,NLAT
           J1=NXTEMP+(I-1)*NLEVEL+J
           YLAT(I)=XN(J1)
          ENDDO

          XP1(J) = (1.0-FLAT)*YLAT(ILAT)+FLAT*YLAT(JLAT)

          K=(ILAT-1)*NLEVEL+J
          GRADL(J,K)=(1.0-FLAT)
          K=(JLAT-1)*NLEVEL+J
          GRADL(J,K)=FLAT

         ENDDO

         DO J=1,NPRO

          L1 = ALOG(P(J))

          F=-1.
          DO JLEV=1,NLEVEL-1
           IF(L1.LE.LP1(JLEV).AND.L1.GT.LP1(JLEV+1))THEN
             ILEV=JLEV
             F = (L1-LP1(JLEV))/(LP1(JLEV+1)-LP1(JLEV))
             GOTO 1181
            ENDIF
          ENDDO
1181      CONTINUE
          IF(F.LT.0)THEN
            IF(L1.GT.LP1(1))THEN
             ILEV=1
             F = (L1-LP1(ILEV))/(LP1(ILEV+1)-LP1(ILEV))
            ELSE
             ILEV=NLEVEL-1
             F = (L1-LP1(ILEV))/(LP1(ILEV+1)-LP1(ILEV))
            ENDIF
            if(idiag.gt.0)then
             print*,'Model 35 warning - pressure out of range'
             print*,'Having to extrapolate'
             print*,L1,LP1(1),LP1(NLEVEL)
             print*,ILEV,LP1(ILEV),LP1(ILEV+1),F
            endif
          ENDIF

          XX = (1.0-F)*XP1(ILEV)+F*XP1(ILEV+1)

          IF(VARIDENT(IVAR,1).EQ.0)THEN
           X1(J)=XX
          ELSE
           X1(J)=EXP(XX)
          ENDIF
          IF(X1(J).LT.1e-36)X1(J)=1e-36
          
          DO K=1,NP
           GRADTMP(J,K)=(1.0-F)*GRADL(ILEV,K)+F*GRADL(ILEV+1,K)         
           IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+K,IPAR,J)=GRADTMP(J,K)
           ELSE
            XMAP(NXTEMP+K,IPAR,J)=X1(J)*GRADTMP(J,K)
           ENDIF
          ENDDO


         ENDDO


        ELSEIF(VARIDENT(IVAR,3).EQ.36)THEN
C        Model 36. Axisymmetric scaling
C        ***************************************************************

C        Need to interpolation in longitude
         NLAT = INT(VARPARAM(IVAR,1))


         if(idiag.gt.0)then
          print*,'Model 36 - nlat,np = ',nlat,np
          print*,'Model 36 - latitude,longitude = ',LATITUDE,
     &           LONGITUDE
         endif

         DLAT=180.0/FLOAT(NLAT-1)
         ILAT=1+INT((90+LATITUDE)/DLAT)
         FLAT = (90+LATITUDE - (ILAT-1)*DLAT)/DLAT
         if(flat.gt.1.0)then 
          print*,'Error: flong > 1.0'
          stop
         endif
         IF(ILAT.LT.NLAT)THEN
          JLAT=ILAT+1
         ELSE
          stop
         ENDIF

         IF(VARIDENT(IVAR,1).EQ.0)THEN
            DX=2.0
         ELSE
            DX=0.1
         ENDIF

         DO K=1,NP
           GRADL(1,K)=0.
         ENDDO

         SUM=0.
         DO I=1,NLAT
          J1=NXTEMP+I
          YLAT(I)=XN(J1)
         ENDDO

         XS = (1.0-FLAT)*YLONG(ILAT)+FLAT*YLAT(JLAT)

         GRADL(1,ILAT)=(1.0-FLAT)
         GRADL(1,JLAT)=FLAT

         DO J=1,NPRO
           X1(J) = XREF(J)*EXP(XS)
           XMAP(NXTEMP+ILONG,IPAR,J)=X1(J)*GRADL(1,ILONG)
           XMAP(NXTEMP+JLONG,IPAR,J)=X1(J)*GRADL(1,JLONG)
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.37)THEN
C        Model 37. constant opacity/bar model
C        ***************************************************************
         P1=VARPARAM(IVAR,1)
         P2=VARPARAM(IVAR,2)
         XDEEP = EXP(XN(NXTEMP+1))

         DO J=1,NPRO
          X1(J)=0.0
C         Need to convert opacity/bar to particles/gram

          IF(P(J)*1.013.LE.P1.AND.P(J)*1.013.GE.P2)THEN
           CALL NEWGRAV(IPLANET,XLAT,H(J),RADIUS,G,PNAME)
           X1(J)=10*XDEEP*G/1e5
           XMAP(NXTEMP+1,IPAR,J)=X1(J)
          ENDIF
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.38)THEN
C        Model 38. Karkoschka CH4 model
C        ***************************************************************
         XC1=VARPARAM(IVAR,1)
         XC2=VARPARAM(IVAR,2)
         RH=VARPARAM(IVAR,3)

         SLOPE = EXP(XN(NXTEMP+1))
         
         CALL modifych4kark(npro,P,T,xc1,xc2,RH,slope,
     1    xch4new,xch4newgrad)
        
         DO J=1,NPRO
          X1(J)=xch4new(J)
          XMAP(NXTEMP+1,IPAR,J)=X1(J)*xch4newgrad(J)
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.39)THEN
C        Model 39. Irwin CH4 model
C        ***************************************************************
         XC1 = EXP(XN(NXTEMP+1))
         XC2 = VARPARAM(IVAR,1)
         RH = VARPARAM(IVAR,2)

         
         CALL modifych4irwin(npro,P,T,xc1,xc2,RH,xch4new,xch4newgrad)
        
         DO J=1,NPRO
          X1(J)=xch4new(J)
          XMAP(NXTEMP+1,IPAR,J)=X1(J)*xch4newgrad(J)
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.40)THEN
C        Model 40. toledo haze model
C        ***************************************************************
         XCOL = EXP(XN(NXTEMP+1))
         XFSH = EXP(XN(NXTEMP+2))
         XFAC = (1.0-XFSH)/XFSH
C        New gradient correction if fsh is held as logs
         DXFAC = -1.0/XFSH

         P1 = VARPARAM(IVAR,1)
         P2 = VARPARAM(IVAR,2)
         P3 = VARPARAM(IVAR,3)


         DO J=1,NPRO
          X1(J)=1E-36
          X2(J)=1E-36
         ENDDO


         JFSH=0
         JFSH1=0     

         DO J=1,NPRO
          IF(P(J).GT.P1)THEN
            X1(J) = 1e-36
          ELSE
            IF(J.GT.1)THEN 
              DELH=H(J)-H(J-1)
            ELSE
              DELH=0
            ENDIF
            IF(JFSH.EQ.0)THEN
	      X1(J)=XCOL
              JFSH=1
            ELSE
              IF(P(J).GT.P2)THEN
                X1(J)=X1(J-1)*EXP(-DELH*XFAC/SCALE(J))
                XMAP(NXTEMP+1,IPAR,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
                XMAP(NXTEMP+2,IPAR,J)=(-DELH/SCALE(J))*DXFAC*
     1            X1(J-1)*EXP(-DELH*XFAC/SCALE(J)) +
     2            XMAP(NXTEMP+2,IPAR,J-1)*EXP(-DELH*XFAC/SCALE(J))
              ELSE
               IF(JFSH1.EQ.0)THEN
                 X2(J)=X1(J-1)*EXP(-DELH*XFAC/SCALE(J))
                 JFSH1=1
                 XMAP(NXTEMP+1,IPAR+1,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
                 XMAP(NXTEMP+2,IPAR+1,J)=(-DELH/SCALE(J))*DXFAC*
     1            X1(J-1)*EXP(-DELH*XFAC/SCALE(J)) +
     2            XMAP(NXTEMP+2,IPAR+1,J-1)*EXP(-DELH*XFAC/SCALE(J))
               ELSE
                 IF(P(J).GE.P3)THEN
                  X2(J)=X2(J-1)*EXP(-DELH*XFAC/SCALE(J))
                  XMAP(NXTEMP+1,IPAR+1,J)=XMAP(NXTEMP+1,IPAR,J-1)*
     1                  EXP(-DELH*XFAC/SCALE(J))
                  XMAP(NXTEMP+2,IPAR+1,J)=(-DELH/SCALE(J))*DXFAC*
     1             X1(J-1)*EXP(-DELH*XFAC/SCALE(J)) +
     2             XMAP(NXTEMP+2,IPAR+1,J-1)*EXP(-DELH*XFAC/SCALE(J))
                 ENDIF
               ENDIF
              ENDIF
            ENDIF
          ENDIF

         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.41)THEN
C        Model 41. Sromovsky descended CH4 model
C        ***************************************************************
         XC1 = EXP(XN(NXTEMP+1))
         PD = EXP(XN(NXTEMP+2))
         RHC = EXP(XN(NXTEMP+3))
         RHM = EXP(XN(NXTEMP+4))
         VX = EXP(XN(NXTEMP+5))
      
         XC2 = VARPARAM(IVAR,1)
         PT = VARPARAM(IVAR,2)
         
         CALL modifych4sromovsky(npro,P,T,xc1,xc2,PD,PT,RHC,RHM,VX,
     1    xch4new)
        
         DO J=1,NPRO
          X1(J)=xch4new(J)
C         Calculating gradients from the Sromovsky model is too hard so set to zero
          XMAP(NXTEMP+1,IPAR,J)=0.0
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.42)THEN
C        Model 42. Ackerman and Marley model
C        ***************************************************************
         XDEEP = EXP(XN(NXTEMP+1))
         FLUX = EXP(XN(NXTEMP+2))
         FRAIN = EXP(XN(NXTEMP+3))
      
         IMODEL = INT(VARPARAM(IVAR,1))
         JSPEC = INT(VARPARAM(IVAR,2))
         DENSCOND = VARPARAM(IVAR,3)
         RADCOND = VARPARAM(IVAR,4)
         MWCOND = VARPARAM(IVAR,5)
         XCORR = VARPARAM(IVAR,6)

C         print*,'VARPARAM',(VARPARAM(IVAR,I),I=1,6)
C         print*,IMODEL,JSPEC,DENSCOND,RADCOND,MWCOND,XCORR

         JVMR=-1
         DO I=1,NVMR
          IF((IDGAS(I).EQ.VARIDENT(IVAR,1)).AND.
     1     (ISOGAS(I).EQ.VARIDENT(IVAR,2)))THEN
           JVMR=I
          ENDIF
         ENDDO

         IF(JVMR.LT.0)THEN
          print*,'Error in subprofretg.f for Model=42'
          print*,'Condensing gas not present'
          stop
         ENDIF

C         print*,IPLANET,LATITUDE,AMFORM,NPRO,NVMR
C         print*,(IDGAS(J),J=1,NVMR)
C         print*,(ISOGAS(J),J=1,NVMR)
C         DO I=1,NPRO
C          print*,P(I),T(I),H(I),XXMOLWT(I),(VMR(I,J),J=1,NVMR)
C         ENDDO

         CALL ACKERMANMARLEYX1(IPLANET,LATITUDE,AMFORM,NPRO,NVMR,
     1    IDGAS,ISOGAS,P,T,H,VMR,XXMOLWT,NCONT,CONT,FLUX,IMODEL,
     2    FRAIN,JVMR,XDEEP,DENSCOND,RADCOND,MWCOND,X1,X2,QC)

C        X1 is set directly by subroutine
C        X2 is associated specific density profile for condensate
C         type JSPEC. This is assigned to the aerosol.prf file later
C         We also need to correct this for the extinction x-section at the 
C         a reverence wavelength as the x-section for other particle types
C         will normally have been normalised.
        
         DO J=1,NPRO
C         Calculating gradients from the Ackerman and Marley  model is too 
C         hard so set to zero
          XMAP(NXTEMP+1,IPAR,J)=0.0
          X2(J)=X2(J)*XCORR
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.43)THEN
C        Model 43. Parmentier and Guillot (2014) and Line et al. (2013)
C        double grey analytic TP profile
C        ***************************************************************
         ALPHAP = EXP(XN(NXTEMP+1))
         BETAP = EXP(XN(NXTEMP+2))
         KIRP = EXP(XN(NXTEMP+3))
         GAMMAV1 = EXP(XN(NXTEMP+4))
         GAMMAV2 = EXP(XN(NXTEMP+5))
         TSTAR = VARPARAM(IVAR,1)
         RSTAR = VARPARAM(IVAR,2)
         SDIST = VARPARAM(IVAR,3)
         TINT = VARPARAM(IVAR,4)

         CALL PARMENTIERGUILLOT1(IPLANET,LATITUDE,NPRO,
     1    P,H,ALPHAP,BETAP,KIRP,GAMMAV1,GAMMAV2,TSTAR,RSTAR,
     2    SDIST,TINT,X1,GRADTOUT)
        
         DO J=1,NPRO
          XMAP(NXTEMP+1,IPAR,J)=ALPHAP*GRADTOUT(J,1)
          XMAP(NXTEMP+2,IPAR,J)=BETAP*GRADTOUT(J,2)
          XMAP(NXTEMP+3,IPAR,J)=KIRP*GRADTOUT(J,3)
          XMAP(NXTEMP+4,IPAR,J)=GAMMAV1*GRADTOUT(J,4)
          XMAP(NXTEMP+5,IPAR,J)=GAMMAV2*GRADTOUT(J,5)
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.44)THEN
C        Model 44. Equivalent to Model 30 (i.e., Inhomogenous disc 
C        profile at multiple locations, but using Parmentier and Guillot 
C        (2014) and Line et al. (2013) double grey analytic TP profile
C        ***************************************************************

C        Need to interpolation in longitude for each parameter
         NLONG = INT(VARPARAM(IVAR,1))
         if(idiag.gt.0)then
          print*,'Model 44 - nlong = ',nlong
          print*,'Model 44 - latitude,longitude = ',LATITUDE,
     &           LONGITUDE
         endif
         DLONG=360.0/FLOAT(NLONG)
         LONGITUDE1=LONGITUDE
         IF(LONGITUDE.LT.0.0)LONGITUDE1=LONGITUDE+360.
         IF(LONGITUDE.GE.360.0)LONGITUDE1=LONGITUDE-360.
         ILONG=1+INT(LONGITUDE1/DLONG)
         FLONG = (LONGITUDE1 - (ILONG-1)*DLONG)/DLONG
         if(flong.gt.1.0)then 
          print*,'Error: flong > 1.0'
          stop
         endif
         IF(ILONG.LT.NLONG)THEN
          JLONG=ILONG+1
         ELSE
          JLONG=1
         ENDIF

C        Set exponent of cos(lat) variation
         XPC = VARPARAM(IVAR,2) 
C        Set other TP parameters
         TSTAR = VARPARAM(IVAR,3)
         RSTAR = VARPARAM(IVAR,4)
         SDIST = VARPARAM(IVAR,5)
         TINT = VARPARAM(IVAR,6)

C        Setting pole temperature to average of morning and
C                     afternoon terminators.
         N1=NLONG/4
         N2=N1+2*N1        

         N1=N1+1
         N2=N2+1
  
         LONZ(1)=ILONG
         LONZ(2)=JLONG
         LONZ(3)=N1
         LONZ(4)=N2

         DO J=1,4
          I=1
          ILON=LONZ(J)
          IX=NXTEMP + (ILON-1)*5+1
          ALPHAP=EXP(XN(IX))
          IX=NXTEMP + (ILON-1)*5+2
          BETAP=EXP(XN(IX))
          IX=NXTEMP + (ILON-1)*5+3
          KIRP=EXP(XN(IX))
          IX=NXTEMP + (ILON-1)*5+4
          GAMMAV1=EXP(XN(IX))
          IX=NXTEMP + (ILON-1)*5+5
          GAMMAV2=EXP(XN(IX))

          CALL PARMENTIERGUILLOT1(IPLANET,LATITUDE,NPRO,
     1    P,H,ALPHAP,BETAP,KIRP,GAMMAV1,GAMMAV2,TSTAR,RSTAR,
     2    SDIST,TINT,TOUT,GRADTOUT)

          DO I=1,NPRO
           XLINE(J,I)=TOUT(I)
           GLINE(J,I,1)=ALPHAP*GRADTOUT(I,1)
           GLINE(J,I,2)=BETAP*GRADTOUT(I,2)
           GLINE(J,I,3)=KIRP*GRADTOUT(I,3)
           GLINE(J,I,4)=GAMMAV1*GRADTOUT(I,4)
           GLINE(J,I,5)=GAMMAV2*GRADTOUT(I,5)
          ENDDO
         ENDDO
   
         DO I=1,NPRO

C        Setting pole temperature to average of morning and
C                     afternoon terminators.
          SUM=0.5*(XLINE(3,I)+XLINE(4,I))

          YTH = (1.0-FLONG)*XLINE(1,I)+FLONG*XLINE(2,I)
          X1(J) = SUM + (YTH-SUM)*(COS(LATITUDE*DTR))**XPC
          
          DO K=1,5
           DXD1=(1.0-FLONG)*GLINE(1,I,K)*(COS(LATITUDE*DTR))**XPC
           DXD2=FLONG*GLINE(2,I,K)*(COS(LATITUDE*DTR))**XPC
           DXD3=(1-(COS(LATITUDE*DTR))**XPC)*0.5*GLINE(3,I,K)
           DXD4=(1-(COS(LATITUDE*DTR))**XPC)*0.5*GLINE(4,I,K)
           IX1 = NXTEMP+LONZ(1)*5+K
           IX2 = NXTEMP+LONZ(2)*5+K
           IX3 = NXTEMP+LONZ(3)*5+K
           IX4 = NXTEMP+LONZ(4)*5+K
           XMAP(IX1,IPAR,I)=DXD1
           XMAP(IX2,IPAR,I)=DXD2
           XMAP(IX3,IPAR,I)=DXD3
           XMAP(IX4,IPAR,I)=DXD4

          ENDDO
    
         ENDDO

        
         DO J=1,NPRO
          XMAP(NXTEMP+1,IPAR,J)=GRADTOUT(J,1)
          XMAP(NXTEMP+2,IPAR,J)=GRADTOUT(J,2)
          XMAP(NXTEMP+3,IPAR,J)=GRADTOUT(J,3)
          XMAP(NXTEMP+4,IPAR,J)=GRADTOUT(J,4)
          XMAP(NXTEMP+5,IPAR,J)=GRADTOUT(J,5)
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.45)THEN
C        Model 39. Irwin CH4 model generalised to have all elements variable
C        ***************************************************************
         XC1 = EXP(XN(NXTEMP+1))
         RH = EXP(XN(NXTEMP+2))
         XC2 = EXP(XN(NXTEMP+3))

         
         CALL modifych4irwin(npro,P,T,xc1,xc2,RH,xch4new,xch4newgrad)

C        Need to update to return gradients of other variables at some point.

         DO J=1,NPRO
          X1(J)=xch4new(J)
          XMAP(NXTEMP+1,IPAR,J)=X1(J)*xch4newgrad(J)
         ENDDO

        ELSEIF(VARIDENT(IVAR,3).EQ.46)THEN
C        Model 46. Profile is represented a double Gaussian with specified 
C        optical thicknesses centred at two variable altitude levels plus 
C        variable FWHM in height. Based on model 14
C        ***************************************************************

         XDEEP = EXP(XN(NXTEMP+1))
         HKNEE = XN(NXTEMP+2)
         XWID  = EXP(XN(NXTEMP+3))
         XDEEP1 = EXP(XN(NXTEMP+4))
         HKNEE1 = XN(NXTEMP+5)
         XWID1  = EXP(XN(NXTEMP+6))


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
                
          Q(J) = XDEEP/(XWID*SQRT(PI))*EXP(-((Y-HKNEE)/XWID)**2)
          Q(J) = Q(J)+
     &      XDEEP1/(XWID1*SQRT(PI))*EXP(-((Y-HKNEE1)/XWID1)**2)
          Q(J)=Q(J)/(XDEEP+XDEEP1)
          ND(J) = Q(J)*RHO 
          OD(J) = ND(J)*SCALE(J)*1e5

          XOD=XOD+OD(J)
   
         ENDDO

C        Empirical correction to XOD
         XOD = XOD*0.25

         DO J=1,NPRO

          X1(J)=Q(J)*(XDEEP+XDEEP1)/XOD   

C         These gradients are not quite correct.
          IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,J)=X1(J)/XDEEP
            XMAP(NXTEMP+4,IPAR,J)=X1(J)/XDEEP1
          ELSE
            XMAP(NXTEMP+1,IPAR,J)=X1(J)
            XMAP(NXTEMP+4,IPAR,J)=X1(J)
          ENDIF

          XMAP(NXTEMP+2,IPAR,J)=2.*(Y-HKNEE)*X1(J)/XWID**2
          XMAP(NXTEMP+3,IPAR,J)=-2.0*((Y-HKNEE)**2)*X1(J)/XWID**3
     &             -  X1(J)/XWID
          XMAP(NXTEMP+5,IPAR,J)=2.*(Y-HKNEE1)*X1(J)/XWID1**2
          XMAP(NXTEMP+6,IPAR,J)=-2.0*((Y-HKNEE1)**2)*X1(J)/XWID1**3
     &             -  X1(J)/XWID1

         ENDDO

C        *** This renormalisation is pretty accurate, but not quite accurate
C        *** enough and so it gets updated in gsetrad.f


        ELSEIF(VARIDENT(IVAR,3).EQ.47)THEN
C        Model 14. Profile is represented a Gaussian with a specified optical
C        thickness centred at a variable pressure level plus a variable FWHM (log press) in 
C        height.
C        FHWM is also folded into total opacity, but this gets renormalised
C        by gsetrad.f anyway.
C        ***************************************************************

         XDEEP = EXP(XN(NXTEMP+1))
         PKNEE = EXP(XN(NXTEMP+2))
         XWID  = EXP(XN(NXTEMP+3))


         Y0=ALOG(PKNEE)


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

          Y=ALOG(P(J))          

          Q(J) = 1./(XWID*SQRT(PI))*EXP(-((Y-Y0)/XWID)**2)
          ND(J) = Q(J)*RHO 
          OD(J) = ND(J)*SCALE(J)*1e5

          XOD=XOD+OD(J)
   
          X1(J)=Q(J)

         ENDDO

C        Empirical correction to XOD
         XOD = XOD*0.25

         DO J=1,NPRO

          X1(J)=Q(J)*XDEEP/XOD

          Y=ALOG(P(J))          
          
          IF(VARIDENT(IVAR,1).EQ.0)THEN
            XMAP(NXTEMP+1,IPAR,J)=X1(J)/XDEEP
          ELSE
            XMAP(NXTEMP+1,IPAR,J)=X1(J)
          ENDIF

          XMAP(NXTEMP+2,IPAR,J)=2.*(Y-Y0)*X1(J)/XWID**2
          XMAP(NXTEMP+3,IPAR,J)=-2.0*((Y-Y0)**2)*X1(J)/XWID**3

          XMAP(NXTEMP+2,IPAR,J)=Y0*2.*(Y-Y0)*X1(J)/XWID**2
          XMAP(NXTEMP+3,IPAR,J)=-2.0*((Y-Y0)**2)*X1(J)/XWID**2
     &             -  X1(J)/XWID

         ENDDO

C        *** This renormalisation is pretty accurate, but not quite accurate
C        *** enough and so it gets updated in gsetrad.f


        ELSEIF(VARIDENT(IVAR,3).EQ.48)THEN
C        ***************************************************************
C        Model 32. Profile is represented by a value at a variable pressure level
C        plus a fractional scale height, plus top pressure. Below the knee pressure 
C        and above ptop the profile is 
C        set to drop exponentially. Similar model to model 32.
C        ***************************************************************

         IF(VARIDENT(IVAR,1).GE.0)THEN
          print*,'Warning from SUBPROFRETG. You are using a'
          print*,'cloud profile parameterisation for a non-cloud'
          print*,'variable'         
          STOP 
         ENDIF

C        Calculate gradient numerically as it's just too hard otherwise
         DO 208 ITEST=1,5


          XDEEP = EXP(XN(NXTEMP+1))
          XFSH  = EXP(XN(NXTEMP+2))
          PKNEE = EXP(XN(NXTEMP+3))
          PTOP = EXP(XN(NXTEMP+4))

          DX=0.05*XN(NXTEMP+ITEST-1)
          IF(DX.EQ.0.)DX=0.1

          IF(ITEST.EQ.2)THEN
            XDEEP=EXP(XN(NXTEMP+1)+DX)
          ENDIF
          IF(ITEST.EQ.3)THEN
            XFSH  = EXP(XN(NXTEMP+2)+DX)
          ENDIF
          IF(ITEST.EQ.4)THEN
            PKNEE = EXP(XN(NXTEMP+3)+DX)
          ENDIF
          IF(ITEST.EQ.5)THEN
            PTOP = EXP(XN(NXTEMP+4)+DX)
          ENDIF


          CALL VERINT(P,H,NPRO,HKNEE,PKNEE)
          CALL VERINT(P,H,NPRO,HTOP,PTOP)
          if(idiag.gt.0)print*,pknee,hknee


C          print*,'ptop,htop = ',ptop,htop
C          print*,'pknee,hknee = ',pknee,hknee

C         Start ND,Q,OD at zero
C 	  N is in units of particles/cm3
C         OD is in units of particles/cm2 = particles/cm3 x length(cm)
C         Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
          DO J=1,NPRO
           ND(J)=0.
           OD(J)=0
           Q(J)=0.
          ENDDO

          JKNEE=-1
C         find levels in atmosphere that span pknee
          DO J=1,NPRO-1
           IF(P(J).GE.PKNEE.AND.P(J+1).LT.PKNEE)THEN
            JKNEE=J
           ENDIF
          ENDDO
 
          IF(JKNEE.LT.0)THEN
           print*,'subprofretg: Error in model 48. Stop'
           print*,'IVAR,XDEEP,XFSH,PKNEE'
           print*,IVAR,XDEEP,XFSH,PKNEE
           STOP
          ENDIF          

          JTOP=-1
C         find levels in atmosphere that span ptop
          DO J=1,NPRO-1
           IF(P(J).GE.PTOP.AND.P(J+1).LT.PTOP)THEN
            JTOP=J
           ENDIF
          ENDDO
 
          IF(JTOP.LT.0)THEN
           print*,'subprofretg: Error in model 48. Stop'
           print*,'IVAR,XDEEP,XFSH,PKNEE,PTOP'
           print*,IVAR,XDEEP,XFSH,PKNEE,PTOP
           STOP
          ENDIF          


C          print*,JKNEE,JTOP

          DELH=H(JKNEE+1)-HKNEE
          XFAC=SCALE(JKNEE)*XFSH
          ND(JKNEE+1)=DPEXP(-DELH/XFAC)

          DELH = HKNEE-H(JKNEE)
          XFAC=0.05*SCALE(JKNEE)
          ND(JKNEE)=DPEXP(-DELH/XFAC)
           
          DO J=1,JKNEE-1
           DELH = H(JKNEE)-H(J)
           ND(J) = DPEXP(-DELH/XFAC)
          ENDDO
         
          DO J=JKNEE+2,JTOP
           DELH = H(J)-H(J-1)
           XFAC = SCALE(J)*XFSH
           ND(J) = ND(J-1)*DPEXP(-DELH/XFAC)
          ENDDO

          XFAC=0.05*SCALE(JTOP)

          DO J=JTOP+1,NPRO
           DELH = H(J)-HTOP
           ND(J) = ND(JTOP)*DPEXP(-DELH/XFAC)
          ENDDO


C          DO J=1,NPRO
C           print*,J,H(J),P(I),ND(J)
C          ENDDO          

          DO J=1,NPRO
            IF(AMFORM.EQ.0)THEN
             XMOLWT=MOLWT
            ELSE
             XMOLWT=XXMOLWT(J)
            ENDIF
C           Calculate density of atmosphere  (g/cm3)
            RHO = (0.1013*XMOLWT/R)*(P(J)/T(J))
C           Calculate initial particles/gram
            Q(J)=ND(J)/RHO
          ENDDO

C         Now integrate optical thickness
          OD(NPRO)=ND(NPRO)*SCALE(NPRO)*XFSH*1E5
          JFSH=-1
          DO J=NPRO-1,1,-1
           IF(J.GT.JKNEE)THEN
             DELH = H(J+1) - H(J)
             XFAC = SCALE(J)*XFSH
             OD(J)=OD(J+1)+(ND(J) - ND(J+1))*XFAC*1E5
           ELSE
             IF(J.EQ.JKNEE)THEN
              DELH = H(J+1)-HKNEE
              XFAC = 0.5*(SCALE(J)+SCALE(J+1))*XFSH         
              OD(J)=OD(J+1)+(1. - ND(J+1))*XFAC*1E5
              DELH = HKNEE-H(J)
              XFAC = XF
              OD(J)=OD(J)+(1. - ND(J))*XFAC*1E5
             ELSE
              DELH = H(J+1)-H(J)
              XFAC = XF
              OD(J)=OD(J+1)+(ND(J+1) - ND(J))*XFAC*1E5
             ENDIF
           ENDIF
          ENDDO

          ODX=OD(1)

C         Now normalise specific density profile.
C         This is also redone in gsetrad.f to make this totally secure.
          DO J=1,NPRO
           OD(J)=OD(J)*XDEEP/ODX
           ND(J)=ND(J)*XDEEP/ODX
           Q(J)=Q(J)*XDEEP/ODX
           IF(Q(J).GT.1e10)Q(J)=1e10
           IF(Q(J).LT.1e-36)Q(J)=1e-36
           NTEST=ISNAN(Q(J))
           IF(NTEST)THEN
            if(idiag.gt.0)then
             print*,'Error in subprofretg.f, cloud density is NAN'
             print*,'Setting to 1e-36'
            endif
	    Q(J)=1e-36
           ENDIF

           IF(ITEST.EQ.1)THEN
            X1(J)=Q(J)
           ELSE
            XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
           ENDIF

          ENDDO

208       CONTINUE


        ELSE

         print*,'Subprofretg: Model parametrisation code is not defined'
         print*,(VARIDENT(IVAR,J),J=1,3)
         STOP

        ENDIF

       ELSE

C       Must hold non-atmospheric parameter - find which.
        IF(VARIDENT(IVAR,1).EQ.999)THEN
C         if(idiag.gt.0)print*,'Surface temperature'
C        no atmospheric mapping
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.888)THEN
C         if(idiag.gt.0)print*,'Surface albedo spectrum'
         IPAR = -1
         NP = INT(VARPARAM(IVAR,1))
        ELSEIF(VARIDENT(IVAR,1).EQ.887)THEN
C         if(idiag.gt.0)print*,'Cloud x-section spectrum'
         IPAR = -1
         NP = INT(VARPARAM(IVAR,1))
        ELSEIF(VARIDENT(IVAR,1).EQ.889)THEN
C         if(idiag.gt.0)print*,'Surface albedo spectrum multiplier'
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.777)THEN
C         if(idiag.gt.0)print*,'Tangent height correction'
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.666)THEN
C         if(idiag.gt.0)print*,'Tangent pressure'
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.555)THEN
C         if(idiag.gt.0)print*,'Radius of Planet'
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.102)THEN
C         if(idiag.gt.0)print*,'Variable profile fraction'
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.444.OR.VARIDENT(IVAR,1).EQ.445.
     &OR.VARIDENT(IVAR,1).EQ.446)THEN
C         if(idiag.gt.0)print*,'Variable size and RI'
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
C            if(idiag.gt.0)print*,'h',J,OD(J)
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
             if(idiag.gt.0)then
              print*,'Error in subprofretg.f, cloud density is NAN'
              print*,'Setting to 1e-36'
             endif
             Q(J)=1e-36
            ENDIF

            IF(ITEST.EQ.1)THEN
             X2(J)=Q(J)
            ELSE
             XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X2(J))/DX
            ENDIF

           ENDDO

25       CONTINUE


         ENDIF
         IF(VARIDENT(IVAR,1).EQ.445)THEN
          NP = 3+(2*INT(VARPARAM(IVAR,1)))
         ELSE
          IF(VARPARAM(IVAR,2).GT.0.0)THEN
           IF(VARIDENT(IVAR,1).EQ.444)THEN
            NP = 2+INT(VARPARAM(IVAR,1))
           ELSE
            NP = 3+2*INT(VARPARAM(IVAR,1))
           ENDIF
          ELSE
           IF(VARIDENT(IVAR,1).EQ.444)THEN
            NP = 3
           ELSE
            NP=5
           ENDIF
          ENDIF
         ENDIF



      ELSEIF(VARIDENT(IVAR,1).EQ.443)THEN
C            ** Cloud with variable top pressure, deep value and power 
C		law cross section with variable index. Currently only
C		works for MultiNest
         IPAR=NVMR+2
         JCONT=1
		 if(MCMCflag.eq.1)then
           hknee = MCMCtop
           if(idiag.gt.0)print*,'Hknee'
           if(idiag.gt.0)print*,hknee
           xdeep = MCMCdeep
           xfsh = MCMCscat
         endif
         IF(XDEEP.GT.0)THEN
C         Calculate gradient numerically as it's just too hard otherwise
C	       DO 26 ITEST=1,2

           DX=0.05*XN(NXTEMP+ITEST-1)
           IF(DX.EQ.0.)DX=0.1

C          Start ND,Q,OD at zero
C          OD is in units of particles/cm2 = particles/cm3 x length(cm)
C          Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
           DO J=1,NPRO
            ND(J)=0.
            OD(J)=0
            Q(J)=0.
           ENDDO

           JFSH=-1
           IF(H(NPRO).LE.HKNEE)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(1)
             ENDIF
C            Calculate density of atmosphere  (g/cm3)
             RHO = P(NPRO)*0.1013*XMOLWT/(R*T(NPRO))
             Q(NPRO)=1.0
             ND(NPRO)=Q(NPRO)*RHO
             JFSH=1
           ENDIF

           DO J=NPRO-1,1,-1
            DELH = H(J+1)-H(J)
            XFAC = SCALE(J)


            IF(H(J).LE.HKNEE)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(J)
             ENDIF

C            Calculate density of atmosphere  (g/cm3)
             RHO = (0.1013*XMOLWT/R)*(P(J)/T(J)) 
             Q(J)=1.0
             ND(J)=Q(J)*RHO
            ENDIF


           ENDDO

           if(idiag.gt.0)print*,'Q1'
           if(idiag.gt.0)print*,Q(1)
C          Integrate optical thickness
           OD(1)=ND(1)*SCALE(1)*1E5
           JFSH=-1
           DO J=2,NPRO
            XFAC = SCALE(J)       
            OD(J)=OD(J-1)+0.5*(ND(J) + ND(J-1))*XFAC*1E5
            IF(H(J).GE.HKNEE.AND.JFSH.LT.0)THEN
             F = (H(J)-HKNEE)/DELH
             XOD = (1.-F)*OD(J) + F*OD(J-1)
             JFSH=1
            ENDIF
C            if(idiag.gt.0)print*,'h',J,OD(J)
           ENDDO

C          The following section was found not to be as accurate as
C          desired due to misalignments at boundaries and so needs some 
C          post-processing in gsetrad.f
           DO J=1,NPRO
            OD(J)=XDEEP*OD(J)/XOD
            Q(J)=XDEEP*Q(J)/XOD
            IF(H(J).GT.HKNEE)THEN
             IF(H(J-1).LE.HKNEE)THEN
              Q(J)=Q(J)*(1.0 - (H(J)-HKNEE))/(H(J)-H(J-1))
             ELSE
              Q(J) = 1.0e-36
             ENDIF
            ENDIF
            IF(Q(J).GT.1e10)Q(J)=1e10
            IF(Q(J).LT.1e-36)Q(J)=1e-36
            NTEST=ISNAN(Q(J))
            IF(NTEST)THEN
             if(idiag.gt.0)then
              print*,'Error in subprofretg.f, cloud density is NAN'
              print*,'Setting Q(J) to 1e-36'
             endif
             Q(J)=1e-36
            ENDIF

C            IF(ITEST.EQ.1)THEN
             X1(J)=Q(J)
C            ELSE
C             XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
C            ENDIF

           ENDDO
           if(idiag.gt.0)print*,'Q1 2'
           if(idiag.gt.0)print*,Q(1)

C26       CONTINUE


         ENDIF
         NP = 3

      ELSEIF(VARIDENT(IVAR,1).EQ.442)THEN
C            ** Cloud with variable top and base pressures, deep value and power 
C    law cross section with variable index. Currently only
C	works for MultiNest***         
         IPAR=NVMR+3
         JCONT=1
		 if(MCMCflag.eq.1)then
           hknee = MCMChknee
           if(idiag.gt.0)print*,'Hknee'
           if(idiag.gt.0)print*,hknee
           htop = MCMCtop
           xdeep = MCMCdeep
           xfsh = MCMCscat
         endif
         IF(XDEEP.GT.0)THEN
C         Calculate gradient numerically as it's just too hard otherwise
C	       DO 26 ITEST=1,2

           DX=0.05*XN(NXTEMP+ITEST-1)
           IF(DX.EQ.0.)DX=0.1

C          Start ND,Q,OD at zero
C          OD is in units of particles/cm2 = particles/cm3 x length(cm)
C          Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
           DO J=1,NPRO
            ND(J)=0.
            OD(J)=0
            Q(J)=0.
         ENDDO

         IF(HTOP.LE.HKNEE)THEN
            DO J=1,NPRO
               Q(J)=1.e-36
               X1(J)=Q(J)
            ENDDO
         ELSE   
           JFSH=-1
           IF(H(NPRO).LE.HTOP.AND.H(NPRO).GT.HKNEE)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(1)
             ENDIF
C            Calculate density of atmosphere  (g/cm3)
             RHO = P(NPRO)*0.1013*XMOLWT/(R*T(NPRO))
             Q(NPRO)=1.0
             ND(NPRO)=Q(NPRO)*RHO
             JFSH=1
          ENDIF


           DO J=NPRO-1,1,-1
            DELH = H(J+1)-H(J)
            XFAC = SCALE(J)

            IF(H(J).GT.HKNEE)THEN
            IF(H(J).LE.HTOP)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(J)
             ENDIF

C            Calculate density of atmosphere  (g/cm3)
             RHO = (0.1013*XMOLWT/R)*(P(J)/T(J)) 
             Q(J)=1.0
             ND(J)=Q(J)*RHO
          ENDIF
          ENDIF


           ENDDO

           if(idiag.gt.0)print*,'Q1'
           if(idiag.gt.0)print*,Q(1)
C          Integrate optical thickness
           OD(1)=ND(1)*SCALE(1)*1E5
           JFSH=-1
           DO J=2,NPRO
            XFAC = SCALE(J)       
            OD(J)=OD(J-1)+0.5*(ND(J) + ND(J-1))*XFAC*1E5
C            IF(H(J).GT.HTOP.AND.JFSH.LT.0)THEN
C             F = (H(J)-HTOP)/DELH
C             XOD = (1.-F)*OD(J) + F*OD(J-1)
C             JFSH=1
C            ENDIF
           ENDDO
           XOD=OD(NPRO)
           if(idiag.gt.0)print*, 'OD(NPRO) ='
           if(idiag.gt.0)print*, OD(NPRO)
           if(idiag.gt.0)print*, 'HKNEE, HTOP'
           if(idiag.gt.0)print*, HKNEE, HTOP
           
C          The following section was found not to be as accurate as
C          desired due to misalignments at boundaries and so needs some 
C          post-processing in gsetrad.f
           DO J=1,NPRO
            OD(J)=XDEEP*OD(J)/XOD
            Q(J)=XDEEP*Q(J)/XOD
C            IF(H(J).GT.HKNEE)THEN
C             IF(H(J-1).LE.HKNEE)THEN
C              Q(J)=Q(J)*(1.0 - (HKNEE-H(J-1)))/(H(J)-H(J-1))
C             ELSE
C              Q(J) = 1.0e-36
C             ENDIF
C          ENDIF
C           IF(H(J).GT.HTOP)THEN
C             IF(H(J-1).LE.HTOP)THEN
C              Q(J)=Q(J)*(1.0 - (H(J)-HTOP))/(H(J)-H(J-1))
C             ELSE
C              Q(J) = 1.0e-36
C             ENDIF
C            ENDIF
            IF(Q(J).GT.1e10)Q(J)=1e10
            IF(Q(J).LT.1e-36)Q(J)=1e-36
            NTEST=ISNAN(Q(J))
            IF(NTEST)THEN
             if(idiag.gt.0)then
              print*,'Error in subprofretg.f, cloud density is NAN'
	      print*, 'Setting Q(J) to 1e-36'
             endif
             Q(J)=1.0e-36
            ENDIF
C            IF(ITEST.EQ.1)THEN
             X1(J)=Q(J)
C            ELSE
C             XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
C            ENDIF

          ENDDO
          ENDIF
           if(idiag.gt.0)print*,'Q1 2'
           if(idiag.gt.0)print*,Q(1)

C26       CONTINUE


         ENDIF
         NP = 4


      ELSEIF(VARIDENT(IVAR,1).EQ.441)THEN
C            ** Haze with variable base pressure & opacity,  and power 
C     law cross section with variable index, with opaque grey cloud beneath.
C     After MacDonald & Madhusudhan 2017
                              
C		 Currently only works for MultiNest
         
         IPAR=NVMR+2
         JCONT=1
		 if(MCMCflag.eq.1)then
           hknee = MCMChknee
           if(idiag.gt.0)print*,'Hknee'
           if(idiag.gt.0)print*,hknee
           xdeep = MCMCdeep
           xfsh = MCMCscat
         endif
         IF(XDEEP.GT.0)THEN
C         Calculate gradient numerically as it's just too hard otherwise
C	       DO 26 ITEST=1,2

           DX=0.05*XN(NXTEMP+ITEST-1)
           IF(DX.EQ.0.)DX=0.1

C          Start ND,Q,OD at zero
C          OD is in units of particles/cm2 = particles/cm3 x length(cm)
C          Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
           DO J=1,NPRO
            ND(J)=0.
            OD(J)=0
            Q(J)=0.
         ENDDO
  
           JFSH=-1
           IF(H(NPRO).GT.HKNEE)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(1)
             ENDIF
C            Calculate density of atmosphere  (g/cm3)
             RHO = P(NPRO)*0.1013*XMOLWT/(R*T(NPRO))
             Q(NPRO)=1.0
             ND(NPRO)=Q(NPRO)*RHO
             JFSH=1
          ENDIF


           DO J=NPRO-1,1,-1
            DELH = H(J+1)-H(J)
            XFAC = SCALE(J)

            IF(H(J).GT.HKNEE)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(J)
             ENDIF

C            Calculate density of atmosphere  (g/cm3)
             RHO = (0.1013*XMOLWT/R)*(P(J)/T(J)) 
             Q(J)=1.0
             ND(J)=Q(J)*RHO
          ELSE
             RHO = (0.1013*XMOLWT/R)*(P(J)/T(J)) 
             Q(J)=1.0e30
             ND(J)=Q(J)*RHO
          ENDIF


           ENDDO

           if(idiag.gt.0)print*,'Q1'
           if(idiag.gt.0)print*,Q(1)
C     Integrate optical thickness above knee pressure
           JFSH=-1
           OD(1)=ND(1)*SCALE(1)*1E5
           JFSH=-1
           DO J=1,NPRO
              IF(H(J).GT.HKNEE)THEN
                 IF(JFSH.EQ.-1) THEN
                 OD(J)=ND(J)*SCALE(J)*1E5
                 JFSH=1
              ELSE
                 XFAC = SCALE(J)       
                 OD(J)=OD(J-1)+0.5*(ND(J) + ND(J-1))*XFAC*1E5
              ENDIF
              ENDIF
C            IF(H(J).GT.HTOP.AND.JFSH.LT.0)THEN
C             F = (H(J)-HTOP)/DELH
C             XOD = (1.-F)*OD(J) + F*OD(J-1)
C             JFSH=1
C            ENDIF
C     if(idiag.gt.0)print*,'h',J,OD(J)
           ENDDO
           XOD=OD(NPRO)
           if(idiag.gt.0)print*, 'OD(NPRO) ='
           if(idiag.gt.0)print*, OD(NPRO)
           if(idiag.gt.0)print*, 'HKNEE'
           if(idiag.gt.0)print*, HKNEE
           
C          The following section was found not to be as accurate as
C          desired due to misalignments at boundaries and so needs some 
C     post-processing in gsetrad.f
  
           DO J=1,NPRO
            IF(H(J).GT.HKNEE)THEN
               OD(J)=XDEEP*OD(J)/XOD
               Q(J)=XDEEP*Q(J)/XOD
            ENDIF
C            IF(H(J).GT.HKNEE)THEN
C             IF(H(J-1).LE.HKNEE)THEN
C              Q(J)=Q(J)*(1.0 - (HKNEE-H(J-1)))/(H(J)-H(J-1))
C             ELSE
C              Q(J) = 1.0e-36
C             ENDIF
C          ENDIF
C           IF(H(J).GT.HTOP)THEN
C             IF(H(J-1).LE.HTOP)THEN
C              Q(J)=Q(J)*(1.0 - (H(J)-HTOP))/(H(J)-H(J-1))
C             ELSE
C              Q(J) = 1.0e-36
C             ENDIF
C            ENDIF
            IF(Q(J).GT.1.0e30)Q(J)=1.0e30
            IF(Q(J).LT.1.0e-36)Q(J)=1.0e-36
            NTEST=ISNAN(Q(J))
            IF(NTEST)THEN
             if(idiag.gt.0)then
              print*,'Error in subprofretg.f, cloud density is NAN'
 	      print*,'Setting Q(J) to 1e-36'
             endif
	     Q(J)=1.0e-36
            ENDIF
C            IF(ITEST.EQ.1)THEN
             X1(J)=Q(J)
C            ELSE
C             XMAP(NXTEMP+ITEST-1,IPAR,J)=(Q(J)-X1(J))/DX
C            ENDIF

          ENDDO

           if(idiag.gt.0)print*,'Q1 2'
           if(idiag.gt.0)print*,Q(1)

C26       CONTINUE


         ENDIF
         NP = 4

        ELSEIF(VARIDENT(IVAR,1).EQ.440)THEN
C     Benneke et al. 2015 cloud. Variable particle size,
C     fixed refractive index Benneke vertical profile (3 params).
C     Still a work in progress - very slow!
           IF(MCMCflag.eq.1)then
              RPARTICLE=MCMCpr
              PKNEE = 10**(MCMChknee)
              XDEEP = MCMCdeep
              SHAPE = MCMCshape

              if(idiag.gt.0)print*,'PKNEE, SHAPE, XDEEP'
              if(idiag.gt.0)print*,PKNEE,SHAPE,XDEEP
           ELSE
              print*,'Error: Using Multinest parameterisation'
              print*, 'with normal NEMESIS'  
              STOP
           ENDIF 

         IF(XDEEP.GT.0)THEN
          IPAR=NVMR+2
C          XFSH=XFSHREF*REFRADIUS/RPARTICLE

C         Calculate gradient numerically as it's just too hard otherwise
C          DO 25 ITEST=1,2


           DX=0.05*XN(NXTEMP+ITEST-1)
           IF(DX.EQ.0.)DX=0.1
C
C           IF(ITEST.EQ.2)THEN
            XFSH=XFSHREF*REFRADIUS/EXP(XN(NXTEMP+1)+DX)
C           ENDIF


C          Start ND,Q,OD at zero
C          OD is in units of particles/cm2 = particles/cm3 x length(cm)
C          Q is specific density = particles/gram = (particles/cm3) / (g/cm3)
           DO J=1,NPRO
            ND(J)=0.
            OD(J)=0
            Q(J)=0.
           ENDDO

           JFSH=-1
           IF(P(1).LT.1.0.AND.P(1).GE.PKNEE)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(1)
           ENDIF
C            Calculate density of atmosphere  (g/cm3)
             ND(1)=XDEEP*(LOG(P(1))-LOG(PKNEE))**SHAPE
             RHO = P(1)*0.1013*XMOLWT/(R*T(1))
             Q(1)=ND(1)/RHO
             JFSH=1
           ENDIF

           DO J=2,NPRO
            DELH = H(J)-H(J-1)
            XFAC = SCALE(J)*XFSH


            IF(P(J).LT.1.0.AND.P(J).GE.PKNEE)THEN
             IF(AMFORM.EQ.0)THEN
              XMOLWT=MOLWT
             ELSE
              XMOLWT=XXMOLWT(J)
             ENDIF
              if(idiag.gt.0)print*,'PJ, XMOLWT'
              if(idiag.gt.0)print*,P(J),XMOLWT
C            Calculate density of atmosphere  (g/cm3)
             ND(J)=XDEEP*(LOG(P(J))-LOG(PKNEE))**SHAPE
             RHO = P(J)*0.1013*XMOLWT/(R*T(J))
             if(idiag.gt.0)print*,'RHO,ND(J)'
             if(idiag.gt.0)print*,RHO,ND(J)
             Q(J)=ND(J)/RHO
           ENDIF


           ENDDO


C          Integrate optical thickness
C           OD(NPRO)=ND(NPRO)*SCALE(NPRO)*XFSH*1E5
C           JFSH=-1
C           DO J=NPRO-1,1,-1
C            XFAC = SCALE(J)*XFSH         
C            OD(J)=OD(J+1)+0.5*(ND(J) + ND(J+1))*XFAC*1E5
C            IF(H(J).LE.HKNEE.AND.JFSH.LT.0)THEN
C             F = (HKNEE-H(J))/DELH
C             XOD = (1.-F)*OD(J) + F*OD(J+1)
C             JFSH=1
C            ENDIF
C            if(idiag.gt.0)print*,'h',J,OD(J)
C           ENDDO

C          The following section was found not to be as accurate as
C          desired due to misalignments at boundaries and so needs some 
C          post-processing in gsetrad.f
           DO J=1,NPRO
C            OD(J)=XDEEP*OD(J)/XOD
C            Q(J)=XDEEP*Q(J)/XOD
C            IF(H(J).LT.HKNEE)THEN
C             IF(H(J+1).GE.HKNEE)THEN
C              Q(J)=Q(J)*(1.0 - (HKNEE-H(J))/(H(J+1)-H(J)))
C             ELSE
C              Q(J) = 0.0
C             ENDIF
C            ENDIF
            IF(Q(J).GT.1e10)Q(J)=1e10
            IF(Q(J).LT.1e-36)Q(J)=1e-36
            NTEST=ISNAN(Q(J))
            IF(NTEST)THEN
             if(idiag.gt.0)then
              print*,'Error in subprofretg.f, cloud density is NAN'
	      print*,'Setting Q(J) to 1e-36'
             endif
             Q(J)=1.0e-36
            ENDIF

            X1(J)=Q(J)
            if(idiag.gt.0)print*,'H(J),X1(J),Q(J)'
            if(idiag.gt.0)print*,H(J),X1(J),Q(J)
            JCONT=1
     

           ENDDO

C        25       CONTINUE


         ENDIF
         NP = 3
         
        ELSEIF(VARIDENT(IVAR,1).EQ.333)THEN
C         if(idiag.gt.0)print*,'Surface gravity (log10(g))'
         IPAR = -1
         NP = 1
        ELSEIF(VARIDENT(IVAR,1).EQ.222)THEN
C         if(idiag.gt.0)print*,'Sromovsky cloud layering'
         IPAR = -1
         NP = 8
        ELSEIF(VARIDENT(IVAR,1).EQ.223)THEN
C         if(idiag.gt.0)print*,'Sromovsky cloud layering with methane'
         IPAR = -1
         DO I=1,NVMR
           IF(IDGAS(I).EQ.6)IPAR=I
         ENDDO
         IF(IPAR.LT.0)THEN
           print*,'Error in subprofretg. Model 223 defined, but'
           print*,'no CH4 in .ref file'
           STOP
         ENDIF
         DO I=1,NPRO
           X1(I)=VMR(I,IPAR)
         ENDDO
         PCH4 = EXP(XN(NXTEMP+5))/1.013
         XFAC = EXP(XN(NXTEMP+6))
         IF(XFAC.GT.1.0)THEN
          if(idiag.gt.0)then
           print*,'Error in subprofretg, model 223. XFAC > 1.'
           print*,'Limiting to 1.0'
          endif
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
C         if(idiag.gt.0)print*,'Sromovsky cloud layering with extended UTC'
         IPAR = -1
         NP = 9
        ELSEIF(VARIDENT(IVAR,1).EQ.225)THEN
C         if(idiag.gt.0)print*,'Sromovsky cloud layering with extended UTC, cut-off'
C         if(idiag.gt.0)print*,'and variable methane'
         IPAR = -1
         DO I=1,NVMR
           IF(IDGAS(I).EQ.6)IPAR=I
         ENDDO
         IF(IPAR.LT.0)THEN
           print*,'Error in subprofretg. Model 225 defined, but'
           print*,'no CH4 in .ref file'
           STOP
         ENDIF
         DO I=1,NPRO
           X1(I)=VMR(I,IPAR)
         ENDDO
         PCH4 = EXP(XN(NXTEMP+5))/1.013
         XFAC = EXP(XN(NXTEMP+11))
C         IF(XFAC.GT.1.0)THEN
C          if(idiag.gt.0)print*,'Error in subprofretg, model 223. XFAC > 1.'
C          if(idiag.gt.0)print*,'Limiting to 1.0'
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
C         if(idiag.gt.0)print*,'Two cloud layering'
         IPAR = -1
         NP = 8

        ELSEIF(VARIDENT(IVAR,1).EQ.227)THEN
C         if(idiag.gt.0)print*,'Creme Brulee layering'
         IPAR = -1
         NP = 7

        ELSEIF(VARIDENT(IVAR,1).EQ.228)THEN
C         if(idiag.gt.0)print*,'Creme Brulee layering'
         IPAR = -1
         NP = 7

        ELSEIF(VARIDENT(IVAR,1).EQ.229)THEN
C         if(idiag.gt.0)print*,'Creme Brulee layering'
         IPAR = -1
         NP = 7

        ELSE

         print*,'SUBPROFRETG: VARTYPE NOT RECOGNISED'
         STOP
        ENDIF
 
       ENDIF
       

       IF(IPAR.GT.0)THEN
        IF(IPAR.LE.NVMR)THEN
         DO I=1,NPRO
          VMR(I,IPAR)=X1(I)
         ENDDO

C        Extra section for combined cloud/gas profile - Model 10 and
C        model 41
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
          IF(VARIDENT(IVAR,3).GE.0)THEN
           DO I=1,NPRO
            CONT(JCONT,I)=X1(I)
           ENDDO
C           DO I=1,NPRO
C            CONT(JCONT,I)=X1(I)/XRHO(I)
C           ENDDO
          ENDIF
          IF(VARIDENT(IVAR,3).EQ.40)THEN
           DO I=1,NPRO
            CONT(JCONT,I)=X1(I)
            CONT(JCONT+1,I)=X2(I)
           ENDDO
          ENDIF
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

       if(idiag.gt.0)print*,'Reading in svp flags from : ',IPFILE
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
                 if(idiag.gt.0)then
                  print*,'Subprofretg: gas predicted to condense'
                  print*,'setting to VMR to SVP x VP / PRESS'
                  print*,IDGAS(IGAS),ISOGAS(IGAS)
                 endif
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
                 if(idiag.gt.0)print*,IGAS,IDGAS(IGAS)
C                *** Extra hack code here for Methane on Neptune
                 IF(IDGAS(IGAS).EQ.6)THEN
                  DO J=1,NPRO
C                  methane and above tropopause
                   IF(P(J).LT.0.1.AND.VMR(J,IGAS).GT.1.5E-3)THEN
                    VMR(J,IGAS)=1.5E-3
                   ENDIF
                  ENDDO
                 ENDIF
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

C     Now make sure the resulting VMRs add up to 1.0 for an
C     AMFORM=1 profile
      IF(AMFORM.EQ.1)THEN
        if(idiag.gt.0)print*,'ISCALE : ',(ISCALE(J),J=1,NVMR)
        CALL ADJUSTVMR(NPRO,NVMR,VMR,ISCALE,IERR)
        if(idiag.gt.0)print*,'IERR  = ',IERR

        IF(IERR.EQ.1)THEN
         if(idiag.gt.0)then
          print*,'Warning from subprofretg. VMRS do not add to 1'
          print*,'Resetting to reference'       
         endif
         DO I=1,NPRO
          DO J=1,NVMR
           VMR(I,J)=(1.0-FLAT)*VMRREF(JLAT,I,J)+
     &		FLAT*VMRREF(JLAT+1,I,J)
          ENDDO
         ENDDO
         CALL ADJUSTVMR(NPRO,NVMR,VMR,ISCALE,IERR)
         if(idiag.gt.0)print*,'IERRX  = ',IERR

        ENDIF

      ENDIF


C     ********  Modify profile with hydrostatic equation ********
      IF(JHYDRO.EQ.0)THEN
       if(idiag.gt.0)print*,'Calling xhydrostath - modified'
       if(idiag.gt.0)then
        do i=1,npro
         print*,'i, P(i),H(i) = ',i,P(i),H(i)
        enddo
       endif
C       print*,AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT
       CALL XHYDROSTATH(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1  IDGAS,ISOGAS,H,P,T,VMR,SCALE)
       if(idiag.gt.0)then
        do i=1,npro
         print*,'modA - i, P(i),H(i) = ',i,P(i),H(i)
        enddo
       endif
      ELSE
       CALL XHYDROSTATP(AMFORM,IPLANET,LATITUDE,NPRO,NVMR,MOLWT,
     1  IDGAS,ISOGAS,H,P,T,VMR,HTAN,PTAN,SCALE)
      ENDIF

      
C     ************* Write out modified profiles *********

      CALL FILE(IPFILE,IPFILE,'prf')
      OPEN(UNIT=2,FILE=IPFILE,STATUS='UNKNOWN',ERR=52)
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
C          if(idiag.gt.0)print*,'sub7',I,P(I),T(I)
505   CONTINUE
C
52    CONTINUE
      CLOSE(2)


      OPEN(UNIT=2,FILE='aerosol.prf',STATUS='UNKNOWN')
C      if(idiag.gt.0)print*,'subprofretg. Writing aerosol.prf'
      BUFFER='# aerosol.prf'   
      WRITE(2,10)BUFFER
      WRITE(2,*)NPRO, NCONT
      DO 41 I=1,NPRO
        WRITE(2,*) H(I),(CONT(J,I),J=1,NCONT)
C        if(idiag.gt.0)print*,H(I),(CONT(J,I),J=1,NCONT)
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


      SUBROUTINE XPROJ(XLAT,XLON,V)
      IMPLICIT NONE
      REAL V(3),XLAT,XLON,THETA,PHI,DTR

      DTR=3.1415927/180.0

      PHI = XLON*DTR
      THETA = (90.0-XLAT)*DTR

      V(1)=SIN(THETA)*COS(PHI)
      V(2)=SIN(THETA)*SIN(PHI)
      V(3)=COS(THETA)

      RETURN

      END
