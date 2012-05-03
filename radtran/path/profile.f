      PROGRAM PROFILE
C     $Id: profile.f,v 1.8 2011-01-26 10:27:12 irwin Exp $
C--------------------------------------------------------------
C_TITLE:  PROFILE: manipulates and plots vertical atmospheric profiles
C
C_KEYS:   RADTRAN,PROG,VMS
C
C_DESCR:  Reads in an atmospheric profile of temperature, pressure,
C         absorber volume mixing ratios versus height and allows manipulation and plotting
C         and plotting of the profile.
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
      INTEGER MAXPRO,MAXVMR
      PARAMETER(MAXPRO=300,MAXVMR=20)
C     MAXPRO is the maximum number of vertical points which can be stored.
C     MAXVMR is the maximum number of vertical mixing ratio profiles.
      REAL H(MAXPRO),P(MAXPRO),T(MAXPRO),VMR(MAXPRO,MAXVMR),RHUM
      REAL H1(MAXPRO),P1(MAXPRO),T1(MAXPRO),VMR1(MAXPRO,MAXVMR)
      REAL XV(MAXPRO)
      INTEGER NPRO,NVMR,ID(MAXVMR),ISO(MAXVMR),NPRO1
C     H is the height in kilometres above some NOMINAL zero.
C     P is the pressure in atmospheres (not bar).
C     T is the temperature in Kelvin.
C     VMR holds the NVMR volume mixing ratio profiles for each of the gases.
C     There are NPRO points in each profile.
C     ID and ISO hold the local identifier and isotope identifier for
C     each gas. Note that this program does not check that you only include
C     each gas once or that the identifiers are valid.
C
C----------------------------------------------------------------------------
      REAL GCONST,R,PCONV,VMRCONV,FRAC,XFMIN
      PARAMETER (PI=3.1415927, GCONST=6.668E-11, R=8.3143)
      PARAMETER (XFMIN=0.001)
      INTEGER NSIMP, NCONV
      PARAMETER (NSIMP=101)
      CHARACTER*100 IPFILE,OPFILE,BUFFER,IDFILE
      CHARACTER COMM
      REAL VMRPRO(MAXPRO*6),VMRH(MAXPRO*6)
      INTEGER I,J,K,N,IFREEZE
      LOGICAL VMRFIL,Q,CONMIX,NEWFILE
      REAL TMIN,TMAX,LGVMIN,LGVMAX,HMAX,PR,HT,TEMP,DX,DY,X(3),Y(3),X0,Y0
      REAL A,B,C,D,SVP,HMIN,PP,LATITUDE,TCORR,X1,X2
      REAL WT(NSIMP),DH,G,SUM,P0,FNORM,MOLWT,SCALE,XCORR
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
C     First reading in an existing set of profiles which could have been
C     produced using a text editor or a previous use of this program or a
C     project specific program.
C     The nominal format (IFORM=0) is:
C     IFORM
C     IPLANET,LATITUDE,NPRO NVMR MOLWT
C     ID(1) ISO(1)
C     :
C     :
C     ID(NVMR) ISO(NVMR)
C     "H        P           T           gas1         gas2         gas3"
C      H(1)     P(1)        T(1)        VMR(1,1)     VMR(1,2)     VMR(1,3)
C      :        :           :           :            :            :
C      :        :           :           :            :            :
C      H(NPRO   P(NPRO)     T(NPRO)     VMR(NPRO,1)  VMR(NPRO,2)  VMR(NPRO,3)
C     "gas4        gas5.....     gasVMR   "
C      VMR(1,4)    VMR(1,5)......VMR(1,NVMR)
C      :           :             :
C      :           :             :
C      VMR(NPRO,4) VMR(NPRO,5)   VMR(NPRO,NVMR)
C     i.e. profiles are stored in blocks of 6 columns each with a descriptive
C     header. Other value of IFORM allow minor variations to format, mainly
C     for data input.
C     IFORM=1 assumes pressure profile unknown
      NEWFILE=ASKYN('Enter profile from an existing file?')
      IF(NEWFILE)THEN
       CALL PROMPT('Enter existing filename : ')
       READ(*,10)IPFILE
10     FORMAT(A)
       CALL FILE(IPFILE,IPFILE,'prf')
       OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
C     First skip header
54     READ(1,1)BUFFER
       IF(BUFFER(1:1).EQ.'#') GOTO 54
       READ(BUFFER,*)IFORM
1      FORMAT(A)
       IF(AMFORM.NE.2)THEN
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
       ELSE
        READ(1,*)IPLANET,LATITUDE,NPRO,NVMR
       ENDIF
       PNAME = SNAME(IPLANET)
       DO 20 I=1,NVMR
       READ(1,*)ID(I),ISO(I)
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

       IF(IFORM.EQ.2)THEN
        DO 303 I=1,NPRO
         FRAC=0.
         DO 43 J=1,NVMR
          FRAC=FRAC+VMR(I,J)
43       CONTINUE
         FRAC=ABS(FRAC-1.0)
         IF(FRAC.GT.XFMIN)THEN
       PRINT*,'Error in Profile.f. Sum of vmrs does not add up to 1.0'
          STOP
         ENDIF
303      CONTINUE 
       ENDIF
C      all processing below assumes that heights are in ascending order
C      so sorting just in case
       DO 12 J=1,NPRO
       DO 12 I=1,NPRO-1
       IF(ABS(H(I)-H(I+1)).LT.0.01)THEN
	WRITE(*,14)
14      FORMAT(' identical height values found')
C	STOP
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
12     CONTINUE
      ELSE
       IFORM=1
       NVMR=0
       CALL PROMPT('Enter planet ID : ')
       READ*,IPLANET
       PNAME = SNAME(IPLANET)
       CALL PROMPT('Enter planetocentric latitude (deg) : ')
       READ*,LATITUDE
       CALL PROMPT('Enter mean mol. wt of atm (g) : ')
       READ*,MOLWT
       CALL PROMPT('Enter number of levels : ')
       READ*,NPRO
       DO I=1,NPRO
        PRINT*,'Point ',I
        CALL PROMPT('Enter ht (km) and temperature (K) : ')
        READ*,H(I),T(I)
       END DO
      END IF
C
C     now entering main control loop
C     a simple menu list is used to ensure program will work anywhere
40    WRITE(*,41)
      WRITE(*,42)
      WRITE(*,417)
      WRITE(*,44)
      WRITE(*,45)
      WRITE(*,46)
      WRITE(*,47)
      WRITE(*,1002)
      WRITE(*,77)
      WRITE(*,49)
      WRITE(*,433)
      WRITE(*,434)
      WRITE(*,435)
      WRITE(*,436)
      WRITE(*,437)
      WRITE(*,2010)
      WRITE(*,50)
      WRITE(*,41)
41    FORMAT(' ---------------------------------------------------')
42    FORMAT(' A - display summary of profiles')
417   FORMAT(' B - add a VMR profile')
44    FORMAT(' C - remove a VMR profile')
45    FORMAT(' D - force VMR to saturation vapour pressure')
46    FORMAT(' E - compute temperature profile from hydrostatic equil.')
47    FORMAT(' F - compute pressure profile from hydrostatic equil.')
1002  FORMAT(' K - compute height profile from hydrostatic equil.')
49    FORMAT(' H - output profiles and exit')
77    FORMAT(' L - add temperature offset tp T - profile')
50    FORMAT(' Q - quit')
435   FORMAT(' M - multiply vmr profile by constant')
433   FORMAT(' I - multiply pressure by const (e.g. scale bar-atm)')
434   FORMAT(' P - output IDL compatable profile data')
436   FORMAT(' G - Map profile on to new pressure grid')
437   FORMAT(' W - Set water vmr profile to pr-um')
 2010 FORMAT(' V - output single gas VMR profile')
      CALL PROMPT('command?')
      READ(*,10)COMM
      CALL UPCASE(COMM)

C --------------------------------------------------------------------
      IF(COMM.EQ.'A')THEN
	WRITE(*,101)NPRO,PNAME
101     FORMAT(' profiles contain',I4,
     1  ' points,  planet name =',A8)
	WRITE(*,102)NVMR
102     FORMAT(1X,I2,' volume mixing ratio profile(s) included')
	DO 103 I=1,NVMR
	WRITE(*,104)I,ID(I),ISO(I)
104     FORMAT(1X,'gas:',I2,'  identifier:',I2,'  isotope:',I4)
103     CONTINUE

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'B')THEN
C       VMR profiles can be read in or calculated
	VMRFIL=ASKYN('read in profile from file?')
	IF(VMRFIL)THEN
	  CALL PROMPT('Enter name of VMR file : ')
	  READ(*,10)IPFILE
	  CALL FILE(IPFILE,IPFILE,'vmr')
	  OPEN(UNIT=1,FILE=IPFILE,STATUS='OLD')
	  READ(1,*)ID(NVMR+1),ISO(NVMR+1)
	  J=0
105       J=J+1
          READ(1,*,END=106)VMRH(J),VMRPRO(J)
	  GOTO 105
106       CONTINUE
          J=J-1
	  IF(J.LT.2)THEN
	    WRITE(*,108)
108         FORMAT(' error reading profile')
	    GOTO 109
	  END IF
C         now sorting
	  DO 114 K=1,J
	  DO 114 I=1,J-1
	  IF(ABS(VMRH(I)-VMRH(I+1)).LT.0.01)THEN
	    WRITE(*,115)
115         FORMAT(' identical heights found')
	    GOTO 109
	  END IF
	  IF(VMRH(I).GT.VMRH(I+1))THEN
	    TEMP=VMRH(I+1)
	    VMRH(I+1)=VMRH(I)
	    VMRH(I)=TEMP
	    TEMP=VMRPRO(I+1)
	    VMRPRO(I+1)=VMRPRO(I)
	    VMRPRO(I)=TEMP
	    END IF
114       CONTINUE
C         now interpolating the input array to find values at profile
C         heights
	  NVMR=NVMR+1
	  DO 107 I=1,NPRO
	  CALL VERINT(VMRH,VMRPRO,J,VMR(I,NVMR),H(I))
107       CONTINUE
109       CLOSE(UNIT=1)
	 ELSE
	  CALL PROMPT('enter gas and isotope identifier')
	  NVMR=NVMR+1
	  READ(*,*)ID(NVMR),ISO(NVMR)
          CONMIX=ASKYN('Assume constant mixing ratio')
          IF(CONMIX)THEN
 	   CALL PROMPT('enter volume mixing ratio')
	   READ(*,*)VMR(1,NVMR)
	   DO 113 I=2,NPRO
	   VMR(I,NVMR)=VMR(1,NVMR)
113        CONTINUE
          ELSE
          PRINT*,'Enter height and vmr in ascending order'
          PRINT*,'Enter any number of points (>2) and terminate'
          PRINT*,'with -1,-1'
	  J=0
205       J=J+1
	  READ*,VMRH(J),VMRPRO(J)
          IF(VMRH(J).EQ.-1)THEN
           J=J-1
           GOTO 206
          END IF
	  GOTO 205
206       CONTINUE
	  IF(J.LT.2)THEN
	    WRITE(*,208)
208         FORMAT('Profile must have at least 2 points')
	    GOTO 209
	  END IF
C         now sorting
	  DO 214 K=1,J
	  DO 214 I=1,J-1
	  IF(ABS(VMRH(I)-VMRH(I+1)).LT.0.01)THEN
	    WRITE(*,115)
	    GOTO 209
	  END IF
	  IF(VMRH(I).GT.VMRH(I+1))THEN
	    TEMP=VMRH(I+1)
	    VMRH(I+1)=VMRH(I)
	    VMRH(I)=TEMP
	    TEMP=VMRPRO(I+1)
	    VMRPRO(I+1)=VMRPRO(I)
	    VMRPRO(I)=TEMP
	    END IF
214       CONTINUE
C         now interpolating the input array to find values at profile
C         heights
	  DO 207 I=1,NPRO
	  CALL VERINT(VMRH,VMRPRO,J,VMR(I,NVMR),H(I))
207       CONTINUE
209       CONTINUE
          END IF
	 END IF

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'C')THEN
	 CALL PROMPT('which profile?')
	 READ(*,*)K
	 IF(K.GT.NVMR.OR.K.LT.1)THEN
	   WRITE(*,110)
110        FORMAT('no such profile')
	 ELSE
	  WRITE(*,111)K,ID(K),ISO(K)
111       FORMAT(' profile',I3,' is for gas:',I2,' isotope:',I4)
	  Q=ASKYN('is this the one you want to remove?')
	  IF(Q)THEN
	    DO 112 J=K+1,NVMR
	    ID(J-1)=ID(J)
	    ISO(J-1)=ISO(J)
	    DO 117 I=1,NPRO
	    VMR(I,J-1)=VMR(I,J)
117         CONTINUE
112         CONTINUE
	    NVMR=NVMR-1
	    END IF
	  END IF

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'D')THEN
	CALL PROMPT('which profile?')
        READ(*,*)K
        IF(K.GT.NVMR.OR.K.LT.1)THEN
	   WRITE(*,178)
178        FORMAT('no such profile')
	ELSE
	  WRITE(*,111)K,ID(K),ISO(K)
	  Q=ASKYN('is this the one you want to force?')
	  IF(Q)THEN
	   WRITE(*,401)
	   WRITE(*,402)
	   WRITE(*,403)
401     FORMAT(/' the form assumed for calculating vapour pressure is')
402     FORMAT(' svp = exp( a + b/T + cT + dT^2 )')
403     FORMAT(' svp is in atmospheres, the temperature T in kelvin')
	   CALL PROMPT('enter a, b, c, d')
	   READ(*,*,ERR=40)A,B,C,D
           CALL PROMPT('Enter rel. humid. after condens. (0.0-1.0):')
           READ(*,*,ERR=40)RHUM
           CALL PROMPT('Freeze upper levels at cold trap values Y/N? ')
           READ(*,10)ANS
           IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')THEN
            IFREEZE=1
           ELSE
            IFREEZE=0
           ENDIF
           PRINT*,'I,P(I),PP,SVP'
           JSWITCH=0
	   DO 404 I=1,NPRO
	   SVP=RHUM*EXP(A+B/T(I)+C*T(I)+D*T(I)*T(I))
	   PP=VMR(I,K)*P(I)
           PRINT*,I,P(I),PP,SVP
           IF(PP.GT.SVP)THEN
            JSWITCH=1
            VMR(I,K)=SVP/P(I)
           END IF
           I1=MAX(1,I-1)
           IF(JSWITCH.EQ.1.AND.VMR(I,K).GT.VMR(I1,K))THEN
             IF(IFREEZE.EQ.1)VMR(I,K)=VMR(I1,K)
           END IF
404        CONTINUE
	  END IF
        END IF

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'E')THEN

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'F')THEN
C       hydrostatic equilibrium
        CALL PROMPT('Enter level number of reference pressure')
        READ*,JZERO
	CALL PROMPT('enter reference pressure (atm)')
	READ(*,*)P0
        P(JZERO)=P0
	SUM=0.
	DO 405 J=2,NSIMP-1
	IF(J.EQ.(2*(J/2)))THEN
	  WT(J)=4.
	 ELSE
	  WT(J)=2.
	  END IF
	WT(1)=1.
	WT(NSIMP)=1.
405     CONTINUE
	WRITE(*,422)
422     FORMAT(1X,' height ',' temp   ','    g    ',' scale HT',
     1  '   pressure  ')
	DO 410 I=JZERO+1,NPRO
	DH=(H(I)-H(I-1))/FLOAT(NSIMP-1)
	FNORM=DH/3.
	DO 420 J=1,NSIMP
	HT=H(I-1)+FLOAT(J-1)*DH
	CALL VERINT(H,T,NPRO,TEMP,HT)
        CALL NEWGRAV(IPLANET,LATITUDE,HT,RADIUS,G,PNAME)
	SCALE=R*TEMP/(MOLWT*G)
	SUM=SUM+WT(J)*FNORM/SCALE
420     CONTINUE
	P(I)=P0*EXP(-SUM)
	WRITE(*,421)HT,TEMP,G,SCALE,P(I)
421     FORMAT(1X,F8.3,F9.3,F9.3,F9.3,E13.5)
410     CONTINUE

        SUM=0.0

	DO 411 I=JZERO-1,1,-1
	DH=(H(I)-H(I+1))/FLOAT(NSIMP-1)
	FNORM=DH/3.
	DO 423 J=1,NSIMP
	HT=H(I+1)+FLOAT(J-1)*DH
	CALL VERINT(H,T,NPRO,TEMP,HT)
        CALL NEWGRAV(IPLANET,LATITUDE,HT,RADIUS,G,PNAME)
	SCALE=R*TEMP/(MOLWT*G)
	SUM=SUM+WT(J)*FNORM/SCALE
423     CONTINUE
	P(I)=P0*EXP(-SUM)
	WRITE(*,421)HT,TEMP,G,SCALE,P(I)
411     CONTINUE

        IFORM=0

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'K')THEN
C       hydrostatic equilibrium
	CALL PROMPT('enter reference level and height : ')
	READ*,JZERO,HZERO
        CALL PROMPT('Enter fixed gravity (1), or calculate (2) ?')
        READ*,IGRAV 
        IF(IGRAV.EQ.1)THEN
         CALL PROMPT('Enter gravitational acceleration (ms-2) : ')
         READ*,G
        ENDIF
        H(JZERO)=HZERO
	WRITE(*,422)
	DO 412 I=JZERO+1,NPRO
           TEMP = 0.5*(T(I-1) + T(I))
           IF(IGRAV.EQ.2)THEN
	    CALL NEWGRAV(IPLANET,LATITUDE,H(I-1),RADIUS,G,PNAME)
           ENDIF
	   SCALE=R*TEMP/(MOLWT*G)
           H(I) = H(I-1) - SCALE*LOG(P(I)/P(I-1))
  	   WRITE(*,421)H(I),TEMP,G,SCALE,P(I)
412     CONTINUE
        DO 413 I = JZERO-1,1,-1
           TEMP = 0.5*(T(I) + T(I+1))
	   IF(IGRAV.EQ.2)THEN
             CALL NEWGRAV(IPLANET,LATITUDE,H(I+1),RADIUS,G,PNAME)
           ENDIF
	   SCALE=R*TEMP/(MOLWT*G)
           H(I) = H(I+1) -SCALE*LOG(P(I)/P(I+1))
  	   WRITE(*,421)H(I),TEMP,G,SCALE,P(I)
413     CONTINUE
        IFORM=0

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'H')THEN
53      CALL PROMPT('output filename?')
	READ(*,10)OPFILE
	CALL FILE(OPFILE,OPFILE,'prf')
	OPEN(UNIT=2,FILE=OPFILE,STATUS='UNKNOWN',ERR=52)
        WRITE(2,*)IFORM
        IF(IFORM.NE.2)THEN
 	 WRITE(2,501)IPLANET,LATITUDE,NPRO,NVMR,MOLWT
        ELSE
 	 WRITE(2,511)IPLANET,LATITUDE,NPRO,NVMR
        ENDIF
501     FORMAT(1X,I3,1X,F6.2,1X,I3,I3,F8.3)
511     FORMAT(1X,I3,1X,F6.2,1X,I3,I3)
	DO 503 I=1,NVMR
	WRITE(2,502)ID(I),ISO(I)
502     FORMAT(1X,I3,I5)
503     CONTINUE
C       writing the first block of profiles
	N=MIN(NVMR,3)
	WRITE(2,504)(I,I=1,N)
504     FORMAT(1X,' height (km) ',' press (atm) ','  temp (K)   ',
     1  3(' VMR gas',I3,2X))
	DO 505 I=1,NPRO
        IF(IFORM.EQ.0)THEN
	  WRITE(2,506)H(I),P(I),T(I),(VMR(I,J),J=1,N)
506       FORMAT(1X,F13.3,E13.5,F13.4,3(E13.5))
         ELSE IF(IFORM.EQ.1)THEN
	  WRITE(2,5101)H(I),T(I),(VMR(I,J),J=1,N)
5101      FORMAT(1X,F13.3,13X,F13.4,3(E13.5))
          END IF
505     CONTINUE
C       writing additional blocks if any
C       N VMR profiles have been written so far
507     IF(NVMR.GT.N)THEN
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
	STOP
52      CONTINUE
	GOTO 53

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'I')THEN
        CALL PROMPT('Enter constant : ')
        READ*,PCONV
        DO 609 I=1,NPRO
         P(I)=P(I)*PCONV
609     CONTINUE

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'M')THEN
        CALL PROMPT('Enter vmr profile number : ')
        READ*,NCONV
        CALL PROMPT('Enter constant : ')
        READ*,VMRCONV
        DO 610 I=1,NPRO
         VMR(I,NCONV)=VMR(I,NCONV)*VMRCONV
610     CONTINUE

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'L')THEN
        CALL PROMPT('Enter temperature offset : ')
        READ*,TCORR
        CALL PROMPT('Enter height range : ')
        READ*,X1,X2
        DO 611 I=1,NPRO
         IF(H(I).GE.X1.AND.H(I).LE.X2)THEN
          T(I)=T(I)+TCORR
         ENDIF
611     CONTINUE

C --------------------------------------------------------------------
       ELSE IF(COMM.EQ.'P')THEN
        CALL PROMPT('Enter IDL filename : ')
        READ(*,10)IDFILE
        CALL FILE(IDFILE,IDFILE,'idl')
        OPEN(44,FILE=IDFILE,STATUS='UNKNOWN')
        WRITE(44,*)NVMR
        DO 637 I=1,NVMR
         WRITE(44,*)ID(I),ISO(I)
637     CONTINUE

        WRITE(44,*)NPRO
        DO 639 I=1,NPRO
         WRITE(44,*)H(I),P(I),T(I)
639     CONTINUE

        DO 638 I=1,NPRO
         WRITE(44,*)(VMR(I,J),J=1,NVMR)
638     CONTINUE
        CLOSE(44)
C --------------------------------------------------------------------
      ELSE IF(COMM.EQ.'V')THEN
        
 2011    CALL PROMPT('output filename?')
         READ(*,10)OPFILE
         CALL FILE(OPFILE,OPFILE,'vmr')
         OPEN(UNIT=2,FILE=OPFILE,STATUS='UNKNOWN',ERR=2012)
         
	 CALL PROMPT('which profile?')
	 READ(*,*)K
	 IF(K.GT.NVMR.OR.K.LT.1)THEN
            WRITE(*,2013)
 2013       FORMAT('no such profile')
         ELSE
            WRITE(2,2014) ID(K),ISO(K)
 2014       FORMAT(1X,I2,1X,I4)
            
            DO 2015 I=1,NPRO
               WRITE(2,2016) H(I),VMR(I,K)
 2016          FORMAT(1X,F13.3,1X,E13.5)
 2015       CONTINUE
         ENDIF

         GOTO 40
 2012    CONTINUE
         GOTO 2011

C --------------------------------------------------------------------
      ELSE IF(COMM.EQ.'G')THEN
         CALL PROMPT('Enter pmax,pmin (bars) and new NPRO : ')
         READ*,XP1,XP2,NPRO1
	 XP1 = XP1/1.013
         XP2 = XP2/1.013
         DO I=1,NPRO1
          P1(I) = EXP(LOG(XP1) + (LOG(XP2)-LOG(XP1))*(I-1)/(NPRO1-1))
          CALL VERINT(P,T,NPRO,T1(I),P1(I))
          CALL VERINT(P,H,NPRO,H1(I),P1(I))
          DO J=1,NVMR
           DO K=1,NPRO
            XV(K) = VMR(K,J)
           ENDDO
           CALL VERINT(P,XV,NPRO,VMR1(I,J),P1(I))
          ENDDO
         ENDDO

         PRINT*,'Mapping done. Now modifying data'
         NPRO=NPRO1
         DO I=1,NPRO
          P(I)=P1(I)
          T(I)=T1(I)
          H(I)=H1(I)
          DO J=1,NVMR
           VMR(I,J)=VMR1(I,J)
          ENDDO
         ENDDO

         Print*,'Regridding done'

C --------------------------------------------------------------------
      ELSE IF(COMM.EQ.'W')THEN
        JW = 0
        DO I=1,NVMR
         IF(ID(I).EQ.1.AND.ISO(I).EQ.0)JW=I
        ENDDO
        IF(JW.LT.1)THEN
         PRINT*,'No water profile present'
        ELSE
         PRINT*,'Water in column ',JW
C        Make sure water is condensed
         A = 17.073
         B = -6110.6
         C = 6.8881E-4
         D = 0.0
         RHUM=1.0
         JSWITCH=0
	 DO I=1,NPRO
	   SVP=RHUM*EXP(A+B/T(I)+C*T(I)+D*T(I)*T(I))
	   PP=VMR(I,JW)*P(I)
           IF(PP.GT.SVP)THEN
            JSWITCH=1
            VMR(I,JW)=SVP/P(I)
           END IF
           I1=MAX(1,I-1)
           IF(JSWITCH.EQ.1.AND.VMR(I,JW).GT.VMR(I1,J))THEN
             VMR(I,JW)=VMR(I1,JW)
           END IF
         ENDDO
         WCOL = 0.0
         DO I=1,NPRO-1
          print*,I,p(i),t(i),h(i),vmr(i,jw)
          PX = 0.5*(P(I)*VMR(I,JW)+P(I+1)*VMR(I+1,JW))
          TX = 0.5*(T(I)+T(I+1))
          DELH = (H(I+1)-H(I))
          WCOL = WCOL + (DELH*1E5)*PX*1.013E5*6.023E17/(8.31*TX)
         ENDDO        
         WCOL = 1000*WCOL/2.68E19
         PRINT*,'Current water column = ',WCOL
         CALL PROMPT('Enter desired water-column (pr-um) : ')
         READ*,XPR
         PRINT*,'Correction factor = ',XPR/WCOL
         DO I=1,NPRO
          VMR(I,JW)=VMR(I,JW)*XPR/WCOL
         ENDDO
        ENDIF
      ELSE IF(COMM.EQ.'Q')THEN
         STOP
      ELSE
         WRITE(*,51)
 51      FORMAT(' no such command')
         GOTO 40
      END IF
      GOTO 40
      END
