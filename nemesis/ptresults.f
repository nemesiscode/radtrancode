      PROGRAM PTRESULTS
C     $Id: ptresults.f,v 1.8 2018-10-19 $
C--------------------------------------------------------------
C_TITLE:  PTRESULTS: Calculates the pressure and its uncertainty at different heights
C                    from the hydrostatic equilibrium equation, given a retrieved temperature
C                    profile and retrieved pressure at a given tangent height
C
C_KEYS:   RADTRAN,PROG,VMS
C
C_DESCR:  Reads in the .ref file of a Nemesis run. It also reads a .unc file, in which the 
C         retrieved temperatures and pressures with the corresponding uncertainties are stored.
C         The rest of the pressure levels are then computed, as well with the uncertainties
C         calculated using the propagation of errors of both temperature and pressure. Finally,
C         the density profiles are also calculated, with the corresponding uncertainties
C
C_ARGS:
C
C_FILES : unit 1 - input profile file
C         unit 2 - output profile file
C
C_CALLS:  FILE,PROMPT,ASKYN
C
C_BUGS:
C
C_HIST:   28jun92 SBC Original version
C
C_END:
C--------------------------------------------------------------
      IMPLICIT NONE
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/planrad.f'
      include '../radtran/includes/constdef.f'
      include '../radtran/includes/dbcom.f'

C     MAXPRO is the maximum number of vertical points which can be stored.
C     MAXGAS is the maximum number of vertical mixing ratio profiles.
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXGAS),RHUM
      REAL TRET(MAXPRO),TRETERR(MAXPRO),RHORET(MAXPRO),RHORETERR(MAXPRO)
      REAL NUMDENS(MAXPRO),NUMDENSERR(MAXPRO)
      REAL XVMR(MAXGAS),SCALE(MAXPRO),SCALEERR(MAXPRO)
      REAL HTAN,PTAN,PTANERR,PRET(MAXPRO),PRETERR(MAXPRO)
      REAL XSUM,XP1,XP2,WCOL,PX,TX,DELH,SH,SHERR
      REAL XPR,XSCALE(MAXPRO),RADIUS
      INTEGER ISCALE(MAXGAS),IPLANET,IVMR,JW,IERR,JTAN,NLATREF
      REAL XV(MAXPRO),XXMASS(MAXPRO),CALCMOLWT,XMOLWT
      INTEGER NPRO,NVMR,IDGAS(MAXGAS),ISOGAS(MAXGAS),NPRO1
C     H is the height in kilometres above some NOMINAL zero.
C     P is the pressure in atmospheres (not bar).
C     T is the temperature in Kelvin.
C     VMR holds the NVMR volume mixing ratio profiles for each of the gases.
C     There are NPRO points in each profile.
C     ID and ISO hold the local identifier and isotope identifier for
C     each gas. Note that this program does not check that you only include
C     each gas once or that the identifiers are valid.
C----------------------------------------------------------------------------
      REAL GCONST,R,PCONV,VMRCONV,FRAC,XFMIN,KB
      PARAMETER (GCONST=6.668E-11, R=8.3143, KB=1.3806485E-23)
      PARAMETER (XFMIN=0.001)
      INTEGER NSIMP, NCONV, AMFORM
      PARAMETER (NSIMP=101)
      CHARACTER*100 BUFFER,REFFILE,RUNNAME,UNCFILE,PTPFILE
      CHARACTER COMM
      REAL VMRPRO(MAXPRO*6),VMRH(MAXPRO*6)
      INTEGER I,J,K,N,IFREEZE
      LOGICAL VMRFIL,Q,CONMIX,NEWFILE,AMFLAG,REFEXIST,UNCEXIST,GASGIANT
      REAL TMIN,TMAX,LGVMIN,LGVMAX,HMAX,PR,HT,TEMP,DX,DY,X(3),Y(3),X0,Y0
      REAL A,B,C,D,SVP,HMIN,PP,LATITUDE,TCORR,X1,X2
      REAL WT(NSIMP),DH,G,SUM,P0,FNORM,MOLWT,XCORR
      CHARACTER*12 XFORM,YFORM
      CHARACTER*100 TEXT
      CHARACTER*8 PNAME
      CHARACTER*8 SNAME(20)
      CHARACTER*1 ANS
      INTEGER IXFORM,IYFORM
      LOGICAL ASKYN

      DATA (SNAME(J),J=1,12) /'Mercury','Venus','Earth','Mars',
     & 'Jupiter','Saturn', 'Uranus','Neptune','Pluto',
     & 'Sun','Titan','Prplanet'/
C----------------------------------------------------------------------------
C
C
C     READING .REF FILE
C
C
C----------------------------------------------------------------------------

      CALL prompt('Enter run name : ')
      READ(5,1)buffer
1     FORMAT(a)
      runname = buffer(1:36)

      CALL readrefhead(runname,npro,nvmr,gasgiant)
      if(npro.gt.maxpro)then
       print*,'Error in Nemesis. npro > maxpro : ',npro,maxpro
       stop
      endif

      call file(runname,reffile,'ref')
      inquire(file=reffile,exist=refexist)
      if(refexist)then
        OPEN(UNIT=1,FILE=reffile,STATUS='OLD')
        MOLWT=-1.
C       First skip header
54      READ(1,1)BUFFER
        IF(BUFFER(1:1).EQ.'#') GOTO 54
        READ(BUFFER,*)AMFORM
        READ(1,*)NLATREF
        IF(NLATREF.GT.1)THEN
         PRINT*,'error :: number of latitudes in PTresults must be 1'
         STOP
        ENDIF
C       Reading AMFORM
        IF(AMFORM.EQ.0)THEN
         READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
        ELSE
         READ(1,*)IPLANET,LATITUDE,NPRO,NVMR 
        ENDIF

C       Reading gases
        DO 23 I=1,NVMR
         READ(1,*)IDGAS(I),ISOGAS(I)
         ISCALE(I)=1
23      CONTINUE

C       Skip header
        READ(1,*)
C       Reading profiles
        DO 30 I=1,NPRO
          READ(1,*)H(I),P(I),T(I),(VMR(I,J),J=1,NVMR)
30      CONTINUE
        CLOSE(UNIT=1)

      else
        print*,'Error in PTresults :: .ref file does not exist'
        stop
      endif


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



C----------------------------------------------------------------------------
C
C
C     READING .UNC FILE
C
C
C----------------------------------------------------------------------------


      call file(runname,uncfile,'unc')
      inquire(file=uncfile,exist=uncexist)
      if(uncexist)then
        OPEN(UNIT=2,FILE=uncfile,STATUS='OLD')
         READ(2,*)NPRO1
         if(npro.ne.npro1)then
          print*,'Error in PTresults. npro must be the same in' 
          print*,'.ref and .unc files'
          stop
         endif
C       Reading retrieved temperature profile and its uncertainty
        DO I=1,NPRO
          READ(2,*)TRET(I),TRETERR(I)
        ENDDO      
c       Reading retrieved pressure and its uncertainty
        READ(2,*)HTAN,PTAN,PTANERR
        CLOSE(UNIT=2)
      else
        print*,'Error in PTresults :: .unc file does not exist'
        stop
      endif


C----------------------------------------------------------------------------
C
C
C     CALCULATING PRESSURE LEVES VIA HYDROSTATIC EQUILIBRIUM
C
C
C----------------------------------------------------------------------------


C     First find level immediately below the reference altitude
C     Then calculate scale heights at each level
      JTAN=1
      DO I=1,NPRO
       IF(H(I).LT.HTAN)THEN
        JTAN=I
       ENDIF
       DO J=1,NVMR
        XVMR(J)=VMR(I,J)
       ENDDO

       CALL NEWGRAV(IPLANET,LATITUDE,H(I),RADIUS,G,PNAME)

       IF(AMFORM.EQ.0)THEN
        XMOLWT=MOLWT
       ELSEIF(AMFORM.EQ.1)THEN
        XMOLWT=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
       ENDIF

       SCALE(I)=R*TRET(I)/(XMOLWT*G)
       SCALEERR(I)=R*TRETERR(I)/(XMOLWT*G)

      ENDDO

C     Calculate pressure levels just above and below reference height
      SH = 0.5*(SCALE(JTAN)+SCALE(JTAN+1))
      SHERR = 0.5*SQRT((SCALEERR(JTAN))**2.0 + (SCALEERR(JTAN+1))**2.0)
      DELH = H(JTAN+1)-HTAN
      PRET(JTAN+1)=PTAN*EXP(-DELH/SH)
      PRETERR(JTAN+1)=PRET(JTAN+1)*SQRT( (PTANERR**2.0/PTAN**2.0) + 
     1   ( (SHERR**2.0) * (DELH**2.0) / SH**4.0 ) ) 
     
      DELH = H(JTAN)-HTAN
      PRET(JTAN)=PTAN*EXP(-DELH/SH)
      PRETERR(JTAN)=PRET(JTAN)*SQRT((PTANERR**2.0/PTAN**2.0) + 
     1  ( (SHERR**2.0) * (DELH**2.0) / SH**4.0 ) )

      
C     Calculate pressure levels from HTAN to the top of the atmosphere
      DO 301 I=JTAN+2,NPRO
        SH = 0.5*(SCALE(I-1) + SCALE(I))
        SHERR = 0.5*SQRT((SCALEERR(I-1))**2.0 + (SCALEERR(I))**2.0)
        DELH = H(I)-H(I-1)
        PRET(I)=PRET(I-1)*EXP(-DELH/SH)
        PRETERR(I)=PRET(I)*SQRT( (PRETERR(I-1)**2.0/PRET(I-1)**2.0) + 
     1    ( (SHERR**2.0) * (DELH**2.0) / SH**4.0 ) )
301   CONTINUE

c     Calculate pressure levels from HTAN to the bottom of the atmosphere
      PRINT*,JTAN
      DO 311 I=JTAN-1,1,-1
        SH = 0.5*(SCALE(I) + SCALE(I+1))
        SHERR = 0.5*SQRT((SCALEERR(I))**2.0 + (SCALEERR(I+1))**2.0)
        DELH = H(I)-H(I+1)
        PRET(I)=PRET(I+1)*EXP(-DELH/SH)
        PRETERR(I)=PRET(I)*SQRT( (PRETERR(I+1)**2.0/PRET(I+1)**2.0) +
     1    ( (SHERR**2.0) * (DELH**2.0) / SH**4.0 ) )
311   CONTINUE
 


C----------------------------------------------------------------------------
C
C
C     CALCULATING THE DENSITIES AND THE PROPAGATED UNCERTAINTIES
C
C
C----------------------------------------------------------------------------
    

      DO I=1,NPRO

       IF(AMFORM.EQ.0)THEN
        XMOLWT=MOLWT
       ELSEIF(AMFORM.EQ.1)THEN
        XMOLWT=CALCMOLWT(NVMR,XVMR,IDGAS,ISOGAS)
       ENDIF

       RHORET(I) = PRET(I) * 101325.0 * XMOLWT / 1000.0 / (R * TRET(I))
       RHORETERR(I) = RHORET(I) * SQRT( ( PRETERR(I) / PRET(I) )**2.0 + 
     1    ( TRETERR(I) / TRET(I) )**2. )

      ENDDO



C----------------------------------------------------------------------------
C
C
C     CALCULATING NUMBER DENSITIES AND THE PROPAGATED UNCERTAINTIES
C
C
C----------------------------------------------------------------------------


      DO I=1,NPRO

       NUMDENS(I) = PRET(I) * 101325.0 / (KB * TRET(I))
       NUMDENSERR(I) = NUMDENS(I) * SQRT( ( PRETERR(I) / PRET(I) )**2.0+
     1    ( TRETERR(I) / TRET(I) )**2. )

      ENDDO


C----------------------------------------------------------------------------
C
C
C     WRITING .PTP FILE WITH THE RESULTS 
C
C
C----------------------------------------------------------------------------

      
      call file(runname,ptpfile,'ptp')
      open(UNIT=3,file=ptpfile,status='unknown')
      write(3,*)npro
      write(3,*)' H(km) ',' T(K) ',' T_err(K) ',' P(atm) ',
     1  ' P_err(atm) ',' RHO(kg/m3) ',' RHO_ERR(kg/m3) ',
     2  ' Dens(m-3) ',' Dens_ERR(m-3) '
      DO I=1,NPRO
        write(3,1000)H(I),TRET(I),TRETERR(I),PRET(I),PRETERR(I),
     1             RHORET(I),RHORETERR(I),NUMDENS(I),NUMDENSERR(I)
 1000   FORMAT(1X,F13.6,1X,F13.6,1X,F13.6,1X,E13.6,1X,E13.6,1X,
     1         E13.6,1X,E13.6,1X,E13.6,1X,E13.6)

      ENDDO
      close(unit=3)

      END

