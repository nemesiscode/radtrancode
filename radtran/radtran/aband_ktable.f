      PROGRAM ABAND_KTABLE
C     ***************************************************************
C     Program to calculate a k-table from band data over a regular grid 
C     of wavelengths/wavenumbers.
C
C     Pat Irwin	?/?/??	Original version
C     Pat Irwin	26/4/12	Added comments.
C
C     ***************************************************************

      IMPLICIT NONE

      INTEGER MG,MDATA,MGAS,MBIN,MPAR,ICALC
      PARAMETER (MG=10,MDATA=20,MGAS=20,MBIN=8000,MPAR=20)
      INTEGER IDGAS(MGAS),ISOGAS(MGAS),NGAS,NPOINT,I,J,NPOINT1
      REAL TWAVEN(MBIN,2),TPOUT(MBIN,MGAS,7)
      REAL G_ORD(MG),K_G(MG),DEL_G(MG)
      REAL U(MDATA),TRAN(MDATA),TFIT(MDATA),SIG(MDATA),CHISQ
      REAL P(MPAR),T(MPAR),PRESS,TEMP,LPMIN,LPMAX,TMIN,TMAX
      REAL KNU0,DELAD,Y0,EL,SFB,C1,C2,VSTART,VEND
      REAL VMIN,DELV,FWHM,STEP,TMP,QROT,ALAMBDA,ERR,Q
      REAL DP,DT,X,BDUM(9),TPART(4),XFAC,T2(MPAR,MPAR)
      INTEGER NG,JGAS,NP,NT,NTRAN,IMOD,IAV(MBIN),IFORM
      INTEGER IP,IT,K,ISTEP,IMETHOD,IBIN,NAV,BANDTYP(MGAS)
      LOGICAL IODD
      INTEGER NODD,NEVEN
      INTEGER IWAVE,KWAVE,NBIN
      REAL WMIN,WMAX,WFWHM,WCENTRAL,W1,W2,X1,X2
      REAL WAVEIN(MBIN),FWHMIN(MBIN)
      REAL FBIN(MBIN)
      CHARACTER*100 IPFILE,OPFILE1,OPFILE2
      CHARACTER*20 HEAD
      CHARACTER*10 BUFFER
      CHARACTER*1 ANS

      DATA G_ORD/0.013047, 0.067468, 0.160295, 0.283302, 0.425563,
     1          0.574437, 0.716698, 0.839705, 0.932532, 0.986953/

      DATA DEL_G/0.033336, 0.074726, 0.109543, 0.134633, 0.147762,
     1          0.147762, 0.134633, 0.109543, 0.074726, 0.033336/


      CALL PROMPT('Enter name of input band file : ')
      READ(5,1)IPFILE
1     FORMAT(A)

      CALL FILE(IPFILE,IPFILE,'ban')

      OPEN(12,FILE=IPFILE,STATUS='OLD')
c     First skip header     
11    READ(12,500)BUFFER
      if(BUFFER.NE.'**********')GOTO 11

      READ(12,401)VMIN
      READ(12,402)DELV
      READ(12,400)FWHM
      READ(12,403)NPOINT
      READ(12,404)NGAS
      PRINT*,'NGAS = ',NGAS
      DO I=1,NGAS
       READ(12,405)IDGAS(I),ISOGAS(I)
       PRINT*,I,IDGAS(I),ISOGAS(I)
      ENDDO

      CALL PROMPT('Select gas number : ')
      READ*,JGAS
      READ(12,406)HEAD
      
      DO 105 I=1,NGAS
       READ(12,406)HEAD
       READ(12,406)HEAD
       IF(HEAD(1:3).EQ.'Ex.')THEN
        DO 111 J=1,NPOINT
         READ(12,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         READ(12,933)TPOUT(J,I,6),TPOUT(J,I,7)
         WAVEIN(J)=TWAVEN(J,1)
         FWHMIN(J)=TWAVEN(J,2)
111     CONTINUE
        BANDTYP(I)=1
       ELSEIF(HEAD(1:3).EQ.'2-E')THEN
        DO 109 J=1,NPOINT
         READ(12,934)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5),TPOUT(J,I,6),TPOUT(J,I,7)
         WAVEIN(J)=TWAVEN(J,1)
         FWHMIN(J)=TWAVEN(J,2)
109     CONTINUE
        IF(HEAD(1:4).EQ.'2-E1')THEN
         BANDTYP(I)=2
        ELSE
         BANDTYP(I)=3
        ENDIF
        ELSEIF(HEAD(1:4).EQ.'Kark')THEN
        DO 179 J=1,NPOINT
         READ(12,*)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4)
        WAVEIN(J)=TWAVEN(J,1)
        FWHMIN(J)=FWHM
179     CONTINUE
        BANDTYP(I)=4 
       ELSE
        DO 110 J=1,NPOINT
C         READ(12,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
C     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         READ(12,*)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         TPOUT(J,I,6)=0.0
         TPOUT(J,I,7)=0.0
         WAVEIN(J)=TWAVEN(J,1)
         FWHMIN(J)=TWAVEN(J,2)
110     CONTINUE
       ENDIF
105   CONTINUE

      CLOSE(12)


      PRINT*,'BANDTYP has been set to ',BANDTYP(JGAS)
      PRINT*,' -1 = T_GV - EKS'
      PRINT*,' 0 = T_GV - extended voigt'
      PRINT*,' 1 = T_GV - Kam C1,C2'
      PRINT*,' 2 = T_GV, 2-EL, T**QROT'
      PRINT*,' 3 = T_GV, 2-EL, tpart(4)'
      PRINT*,' 4 = Karkoschka09'
      PRINT*,' 5 = Goody-Lorentz'

      PRINT*,'OK (Y/N) : '
      READ(5,1)ANS
      IF(ANS.EQ.'n'.OR.ANS.EQ.'N')THEN
       PRINT*,'Enter new value : '
       READ*,BANDTYP(JGAS)
      ENDIF


      CALL PROMPT('Enter number of pressure points ( < 20 ) : ')
      READ*,NP  
3     CALL PROMPT('Enter range(1) or individual pressures(2) : ')
      READ*,IP
      IF(IP.LT.1.OR.IP.GT.2)GOTO 3

      IF(IP.EQ.1)THEN
       CALL PROMPT('Enter log(pmin), log(pmax) : ')
       READ*,LPMIN,LPMAX
       DP=(LPMAX-LPMIN)/FLOAT(NP-1)     
       DO 5 J=1,NP
         X=LPMIN + (J-1)*DP
         P(J)=EXP(X)
         PRINT*,J,P(J)
5      CONTINUE

      ELSE
     
       PRINT*,'Enter pressures (in bars)'
       DO 6 J=1,NP
         READ*,X
         P(J)=X/1.013
         PRINT*,J,P(J)
6      CONTINUE
      ENDIF

      PRINT*,'Entering NT < 0 activates new pressure-dependent'
      PRINT*,'temperatures'
      CALL PROMPT('Enter number of temperature points ( < 20 ) : ')
      READ*,NT
      IF(NT.GT.0)THEN
4       CALL PROMPT('Enter range(1) or individual temperatures(2) : ')
        READ*,IT
        IF(IT.LT.1.OR.IT.GT.2)GOTO 4
        IF(IT.EQ.1)THEN 
          CALL PROMPT('Enter Tmin, Tmax : ')
          READ*,TMIN,TMAX
          DT=(TMAX-TMIN)/(NT-1)
          DO J=1,NT
            T(J)=TMIN + (J-1)*DT
            PRINT*,J,T(J)
          ENDDO
     
        ELSE

          PRINT*,'Enter Temperatures (in K)'
          DO J=1,NT
            READ*,T(J)   
            PRINT*,J,T(J) 
          ENDDO

        ENDIF
        
      ELSE

       CALL PROMPT('Enter range(1) or individual temperatures(2) : ')
       READ*,IT

        DO I=1,NP
         PRINT*,'Pressure level = ',I,P(I)
         IF(IT.EQ.1)THEN
           CALL PROMPT('Enter Tmin, Tmax : ')
           READ*,TMIN,TMAX
           DT=(TMAX-TMIN)/(ABS(NT)-1)
           DO J=1,ABS(NT)
             T2(I,J)=TMIN + (J-1)*DT
             PRINT*,J,T2(I,J)
           ENDDO
         ELSE
           PRINT*,'Enter Temperatures (in K)'
           DO J=1,ABS(NT)
             READ*,T2(I,J)
             PRINT*,J,T2(I,J)
           ENDDO
         ENDIF

        ENDDO

      ENDIF

C      CALL PROMPT('Enter number of points in transmission curves : ')
C      READ*,NTRAN 
      NTRAN=20

      DO J=1,4
       TPART(J)=0.0
      ENDDO
      QROT=1.5
       
      IF(BANDTYP(JGAS).LT.3.OR.BANDTYP(JGAS).EQ.5)THEN   
       CALL PROMPT('Enter QROT : ')
       READ*,QROT   
      ELSEIF(BANDTYP(JGAS).EQ.3)THEN   
       CALL PROMPT('Enter TPART(1-4) : ')
       READ*,(TPART(J),J=1,4)
      ENDIF

      CALL PROMPT('Enter Q (0=foreign, 1=self) : ')
      READ*,Q

C      CALL PROMPT('Enter initial ALAMBDA : ')
C      READ*,ALAMBDA

      ALAMBDA=20000.0
      NG = 10

      PRINT*,'Enter wave space of input band file'
      CALL PROMPT('(0=wavenumber, 1=wavelength (um) : ')
      READ*,KWAVE
           
      PRINT*,'Enter wave space of desired ASCII k-table'
      CALL PROMPT('k-table (0=wavenumber, 1=wavelength (um) : ')
      READ*,IWAVE
      print*,KWAVE,IWAVE

      CALL PROMPT('Enter min, max, step and FWHM : ')
      READ*,WMIN,WMAX,STEP,WFWHM
      PRINT*,WMIN,WMAX,STEP,WFWHM

      NBIN = 1 + INT((WMAX-WMIN)/STEP)


      IF(BANDTYP(JGAS).LT.4)THEN
       CALL PROMPT('Enter trans. model. -1=EKS,0=others : ')
       READ*,IMOD
       IF(BANDTYP(JGAS).LT.2.AND.IMOD.LT.0)BANDTYP(JGAS)=-1
      ENDIF

      CALL PROMPT('Enter output ASCII k-table name : ')
      READ(5,1)OPFILE1
      CALL FILE(OPFILE1,OPFILE1,'asc')

      CALL PROMPT('Enter output ASCII transmission file name : ')
      READ(5,1)OPFILE2
      CALL FILE(OPFILE2,OPFILE2,'tra')

      print*,OPFILE1
      print*,OPFILE2

      OPEN(12,FILE=OPFILE1,STATUS='UNKNOWN')
      OPEN(13,FILE=OPFILE2,STATUS='UNKNOWN')

      WRITE(12,1)IPFILE
      WRITE(12,*)IMOD,' ! IMOD (-1=EKS, 0=others)'
      WRITE(12,*)BANDTYP(JGAS),' ! Band type'
      WRITE(12,*)IDGAS(JGAS),ISOGAS(JGAS),' ! ID,ISO'
      WRITE(12,*)WMIN,WMAX,' ! WMIN,WMAX'
      WRITE(12,*)STEP,WFWHM,NBIN,' ! STEP,WFWHM,NBIN'
      WRITE(12,*)NP,' ! Number of pressures'
      WRITE(12,*)(P(I),I=1,NP)
      WRITE(12,*)NT,' ! Number of temperatures'
      IF(NT.GT.0)THEN
       WRITE(12,*)(T(I),I=1,NT)
      ELSE
       DO J=1,NP
        WRITE(12,*)(T2(J,I),I=1,ABS(NT))
       ENDDO
      ENDIF
      WRITE(12,*)NG,' ! Number of g-ordinates'
      WRITE(12,*)(DEL_G(I),I=1,NG)
      WRITE(12,*)NTRAN,' ! Number of transmission points'
      WRITE(12,*)QROT,Q,' ! QROT,Q'
      WRITE(12,*)(TPART(I),I=1,4),' ! Tpart'

      WRITE(13,1)IPFILE
      WRITE(13,*)IMOD,' ! IMOD (-1=EKS, 0=others)'
      WRITE(13,*)BANDTYP(JGAS),' ! Band type'
      WRITE(13,*)WMIN,WMAX,' ! WMIN,WMAX'
      WRITE(13,*)STEP,WFWHM,NBIN,' ! STEP,WFWHM,NBIN'
      WRITE(13,*)NP,' ! Number of pressures'
      WRITE(13,*)(P(I),I=1,NP)
      WRITE(13,*)NT,' ! Number of temperatures'
      IF(NT.GT.0)THEN
       WRITE(13,*)(T(I),I=1,NT)
      ELSE
       DO J=1,NP
        WRITE(13,*)(T2(J,I),I=1,ABS(NT))
       ENDDO
      ENDIF
      WRITE(13,*)NG,' ! Number of g-ordinates'
      WRITE(13,*)(DEL_G(I),I=1,NG)
      WRITE(13,*)NTRAN,' ! Number of transmission points'
      WRITE(13,*)QROT,Q,' ! QROT,Q'
      WRITE(13,*)(TPART(I),I=1,4),' ! Tpart'

      DO I=1,9
       BDUM(I)=0.0
      ENDDO

      DO 100 IBIN = 1,NBIN
       WCENTRAL = WMIN + (IBIN-1)*STEP
       W1 = WCENTRAL-0.5*WFWHM
       W2 = W1 + WFWHM
       print*,'IBIN,W1,W2,WCEN',IBIN,W1,W2,WCENTRAL
       print*,'delv = ',delv
       CALL GETAVBINS(NPOINT,KWAVE,WAVEIN,FWHMIN,DELV,W1,W2,
     1 IWAVE,NAV,IAV,FBIN)

       WRITE(12,*)I,WCENTRAL,'! IPOINT, WAVE'
       WRITE(12,934)(BDUM(J),J=1,9)
       WRITE(13,*)I,WCENTRAL,'! IPOINT, WAVE'
       WRITE(13,934)(BDUM(J),J=1,9)

       KNU0 = 0.0
C       print*,nav
       DO I=1,NAV 
        IF(TPOUT(IAV(I),JGAS,1).GT.KNU0)THEN
         KNU0 = TPOUT(IAV(I),JGAS,1)
        ENDIF
C        print*,i,iav(i),fbin(i)
       ENDDO      

       DO 150 IP=1,NP
         PRESS = P(IP)
         WRITE(6,*)IP,PRESS,' ! IP,PRESS'
         WRITE(12,*)IP,PRESS,' ! IP,PRESS'
         WRITE(13,*)IP,PRESS,' ! IP,PRESS'
         DO 200 IT=1,ABS(NT)
          IF(NT.GT.0)THEN
            TEMP=T(IT)
          ELSE
            TEMP=T2(IP,IT)
          ENDIF
          WRITE(6,*)IT,TEMP,' ! IT,TEMP'
          WRITE(12,*)IT,TEMP,' ! IT,TEMP'
          WRITE(13,*)IT,TEMP,' ! IT,TEMP'

          IF(KNU0.EQ.0.0)THEN
            DO K=1,NG
             K_G(K)=0.0
            ENDDO
            DO K=1,NTRAN
             U(K)=0.0
             TRAN(K)=0.0
             TFIT(K)=0.0
            ENDDO
            ERR = 0.0
          ELSE 

            CALL TRANSET(BANDTYP,TPOUT,JGAS,NAV,IAV,FBIN,QROT,TPART,
     1       PRESS,TEMP,Q,U,TRAN,SIG,NTRAN)

            CALL ROUGHK(NTRAN,U,TRAN,NG,DEL_G,K_G)

            CALL FINEFITK(NTRAN,U,TRAN,SIG,NG,DEL_G,K_G,ICALC)


C           Add in nasty check for weird things happening for K_G(NG)
            XFAC=K_G(NG)/K_G(NG-1)
            IF(XFAC.GT.1E7)THEN
             K_G(NG)=K_G(NG-1)+(K_G(NG-1)-K_G(NG-2))
            ENDIF

            CHISQ=0.0
            DO J=1,NTRAN
             TFIT(J)=0.0
             DO K=1,NG
              TFIT(J)=TFIT(J)+DEL_G(K)*EXP(-K_G(K)*U(J))
             ENDDO
             CHISQ=CHISQ + (TRAN(J)-TFIT(J))**2
            ENDDO
            ERR = SQRT(CHISQ/FLOAT(NTRAN))

          ENDIF

          WRITE(12,*)(K_G(K),K=1,NG)
          WRITE(13,*)(U(K),K=1,NTRAN)
          WRITE(13,*)(TRAN(K),K=1,NTRAN)
          WRITE(13,*)(TFIT(K),K=1,NTRAN)
          WRITE(12,*)ERR,IMOD,ICALC
          WRITE(13,*)ERR,IMOD,ICALC
200     CONTINUE
150    CONTINUE
100   CONTINUE

      CLOSE(12)
      CLOSE(13)

      PRINT*,'Aband_ktable. Run OK'

500   FORMAT(1X,A10)

401   FORMAT(8X,F8.2)
402   FORMAT(8X,F8.2)
400   FORMAT(8X,F8.2)
403   FORMAT(10X,I5)
404   FORMAT(8X,I3)      
405   FORMAT(7X,I3,I3)


C401   FORMAT(1X,'VMIN = ',F8.2)
C402   FORMAT(1X,'DELV = ',F8.2)
C400   FORMAT(1X,'FWHM = ',F8.2)
C403   FORMAT(1X,'NPOINT = ',I5)
C404   FORMAT(1X,'NGAS = ',I3)      
C405   FORMAT(1X,'Gas : ',I3,I3)

406   FORMAT(1X,A)
932   FORMAT(1X,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5)
933   FORMAT(1x,2(e12.5))
934   FORMAT(1X,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5)

      END


      SUBROUTINE TRANSLATE(IWAVE,W1,W2,KWAVE,X1,X2)
      INTEGER IWAVE,KWAVE
      REAL W1,W2,X1,X2
       
      IF(IWAVE.EQ.0)THEN 
C      Required table is in wavenumbers
       IF(KWAVE.EQ.0)THEN
        X1=W1
        X2=W2
       ELSEIF(KWAVE.EQ.1)THEN
        X1 = 1E4/W2
        X2 = 1E4/W1
       ELSE
        PRINT*,'TRANSLATE: KWAVE NOT DEFINED',KWAVE
        STOP
       ENDIF

      ELSEIF(IWAVE.EQ.1)THEN
C      Required table is in wavelength (microns)
       IF(KWAVE.EQ.0)THEN
        X1=1E4/W2  
        X2=1E4/W1  
       ELSEIF(KWAVE.EQ.1)THEN
        X1 = W1
        X2 = W2
       ELSE 
        PRINT*,'TRANSLATE: KWAVE NOT DEFINED',KWAVE
        STOP
       ENDIF

      ELSE
       PRINT*,'TRANSLATE: IWAVE NOT DEFINED',IWAVE
       STOP
      ENDIF
        
      RETURN
      END

