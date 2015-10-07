      PROGRAM ABAND_KTABLEC
C     ***************************************************************
C     Program to calculate a k-table from band data at a set of output 
C     wavelengths/wavenumbers read in from an input file.
C
C     Pat Irwin ?/?/??   Original version
C     Pat Irwin 26/4/12  Added comments.
C     SPH       15/08/13 Read from filter file include
C
C     ***************************************************************

      IMPLICIT NONE

      INTEGER MG,MDATA,MGAS,MBIN,MPAR,ICALC
      PARAMETER (MG=10,MDATA=20,MGAS=20,MBIN=8000,MPAR=20)
      INTEGER IDGAS(MGAS),ISOGAS(MGAS),NGAS,NPOINT,I,J,NPOINT1
      REAL TWAVEN(MBIN,2),TPOUT(MBIN,MGAS,7),WCEN(MBIN)
      REAL G_ORD(MG),K_G(MG),DEL_G(MG)
      REAL U(MDATA),TRAN(MDATA),TFIT(MDATA),SIG(MDATA),CHISQ
      REAL P(MPAR),T(MPAR),PRESS,TEMP,LPMIN,LPMAX,TMIN,TMAX
      REAL KNU0,DELAD,Y0,EL,SFB,C1,C2,VSTART,VEND
      REAL VMIN,DELV,FWHM,STEP,TMP,QROT,ALAMBDA,ERR,Q
      REAL DP,DT,X,BDUM(9),TPART(4)
      INTEGER NG,JGAS,NP,NT,NTRAN,IMOD,IAV(MBIN),IFORM
      INTEGER IP,IT,K,ISTEP,IMETHOD,IBIN,NAV,BANDTYP(MGAS)
      LOGICAL IODD
      INTEGER NODD,NEVEN
      INTEGER IWAVE,KWAVE,NBIN,ISHAPE,ISHAPE1
      REAL WMIN,WMAX,WFWHM,WCENTRAL
      REAL WAVEIN(MBIN),FWHMIN(MBIN)
      REAL FBIN(MBIN)
      CHARACTER*100 IPFILE,OPFILE1,OPFILE2,GFILE
      CHARACTER*20 HEAD
      CHARACTER*10 BUFFER

c ------------------------------------------------------
c Added by SPH. New variables, please put in a more 
c sensible place after debugging.
c ------------------------------------------------------
      INTEGER NFILBIN
      REAL WFIL(MBIN), TFIL(MBIN)
c ------------------------------------------------------

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
         READ(12,932)TWAVEN(J,1),TWAVEN(J,2),TPOUT(J,I,1),TPOUT(J,I,2),
     1 TPOUT(J,I,3),TPOUT(J,I,4),TPOUT(J,I,5)
         TPOUT(J,I,6)=0.0
         TPOUT(J,I,7)=0.0
         WAVEIN(J)=TWAVEN(J,1)
         FWHMIN(J)=TWAVEN(J,2)
110     CONTINUE
       ENDIF
       print*,'Gas ',I,'Bandtype: ',BANDTYP(I)
105   CONTINUE

      CLOSE(12)

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


      CALL PROMPT('Enter number of temperature points ( < 20 ) : ')
      READ*,NT
4     CALL PROMPT('Enter range(1) or individual temperatures(2) : ')
      READ*,IT
      IF(IT.LT.1.OR.IT.GT.2)GOTO 4
      IF(IT.EQ.1)THEN 
       CALL PROMPT('Enter Tmin, Tmax : ')
       READ*,TMIN,TMAX
       DT=(TMAX-TMIN)/(NT-1)
       DO 7 J=1,NT
         T(J)=TMIN + (J-1)*DT
         PRINT*,J,T(J)
7      CONTINUE
     
      ELSE

       PRINT*,'Enter Temperatures (in K)'
       DO 8 J=1,NT
         READ*,T(J)   
         PRINT*,J,T(J) 
8      CONTINUE

      ENDIF
        
C      CALL PROMPT('Enter number of points in transmission curves : ')
C      READ*,NTRAN 
      NTRAN=20

      DO J=1,4
       TPART(J)=0.0
      ENDDO
      QROT=1.5
       
      IF(BANDTYP(JGAS).LT.3)THEN   
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

      PRINT*,'Enter instrument line shape code'
      PRINT*,'0 = square'
      PRINT*,'1 = triangular'
      PRINT*,'2 = gaussian'
      PRINT*,'3 = individual filters'
      CALL PROMPT('Enter code : ')
      READ*,ISHAPE
      CALL PROMPT('Enter name of wavelengths or filter file ')
      READ(5,1)GFILE
      OPEN(12,FILE=GFILE,STATUS='OLD')
      IF(ISHAPE.LT.3)THEN
       READ(12,*)NBIN,WFWHM,ISHAPE1
       IF(ISHAPE.NE.ISHAPE1)THEN
        PRINT*,'ISHAPE <> ISHAPE1'
        PRINT*,ISHAPE,ISHAPE1
        STOP
       ENDIF
       DO I=1,NBIN
        READ(12,*)WCEN(I)
       ENDDO
       NFILBIN=0
       CLOSE(12)

      ELSE
       WFWHM=-1.
       READ(12,*)NBIN
C      Read in the central wavelengths. Need to read through file, close 
C      and then open again.
       DO I=1,NBIN
        READ(12,*)WCEN(I)
        READ(12,*)NFILBIN
        DO J=1,NFILBIN
         READ(12,1)BUFFER
        ENDDO
       ENDDO
       CLOSE(12)
       OPEN(12,FILE=GFILE,STATUS='OLD')
       READ(12,*)NBIN

      ENDIF
c ---------------------------------------------------------------------

      IMOD=0 
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

      OPEN(14,FILE=OPFILE1,STATUS='UNKNOWN')
      OPEN(13,FILE=OPFILE2,STATUS='UNKNOWN')

      WRITE(14,1)IPFILE
      WRITE(14,*)IMOD,' ! IMOD (-1=EKS, 0=others)'
      WRITE(14,*)BANDTYP(JGAS),' ! Band type'
      WRITE(14,*)IDGAS(JGAS),ISOGAS(JGAS),' ! ID,ISO'
      WRITE(14,*)NP,' ! Number of pressures'
      WRITE(14,*)(P(I),I=1,NP)
      WRITE(14,*)NT,' ! Number of temperatures'
      WRITE(14,*)(T(I),I=1,NT)
      WRITE(14,*)NG,' ! Number of g-ordinates'
      WRITE(14,*)(DEL_G(I),I=1,NG)
      WRITE(14,*)NTRAN,' ! Number of transmission points'
      WRITE(14,*)QROT,Q,' ! QROT,Q'
      WRITE(14,*)(TPART(I),I=1,4),' ! Tpart'
      WRITE(14,*)NBIN,WFWHM,ISHAPE,' ! NBIN, WFWHM, ISHAPE'
      DO I=1,NBIN
        WRITE(14,*)WCEN(I)
      ENDDO

      WRITE(13,1)IPFILE
      WRITE(13,*)IMOD,' ! IMOD (-1=EKS, 0=others)'
      WRITE(13,*)BANDTYP(JGAS),' ! Band type'
      WRITE(13,*)NP,' ! Number of pressures'
      WRITE(13,*)(P(I),I=1,NP)
      WRITE(13,*)NT,' ! Number of temperatures'
      WRITE(13,*)(T(I),I=1,NT)
      WRITE(13,*)NG,' ! Number of g-ordinates'
      WRITE(13,*)(DEL_G(I),I=1,NG)
      WRITE(13,*)NTRAN,' ! Number of transmission points'
      WRITE(13,*)QROT,Q,' ! QROT,Q'
      WRITE(13,*)(TPART(I),I=1,4),' ! Tpart'
      WRITE(13,*)NBIN,WFWHM,ISHAPE,' ! NBIN, WFWHM, ISHAPE'
      DO I=1,NBIN
       WRITE(13,*)WCEN(I)
      ENDDO

      DO I=1,9
       BDUM(I)=0.0
      ENDDO

      DO 100 IBIN = 1,NBIN
c ---------------------------------------------------------------------
c Added by SPH. Reading the filter file has been moved here for ishape=3.
C Remember to then close unit 12 if ishape=3 (otherwise it is already closed).
c ---------------------------------------------------------------------     
       print*,'IBIN',IBIN
       IF (ISHAPE.EQ. 3) THEN
         READ(12,*)WCEN(IBIN)
         READ(12,*)NFILBIN
         print*,WCEN(IBIN),NFILBIN
         IF (NFILBIN.GT.MBIN) THEN
            PRINT*,'***Error: NFILBIN > MBIN',NFILBIN, MBIN
            STOP
         ENDIF
         DO J=1,NFILBIN
            READ(12,*)WFIL(J),TFIL(J)
         ENDDO
         IF (IBIN.EQ.NBIN) CLOSE(12)
       ENDIF
 
       WCENTRAL = WCEN(IBIN)
       print*,'IBIN,WCEN',IBIN,WCENTRAL

       CALL GETAVBINSC(NPOINT,KWAVE,WAVEIN,FWHMIN,DELV,WCENTRAL,
     1     WFWHM,ISHAPE,NFILBIN,WFIL,TFIL,IWAVE,NAV,IAV,FBIN)

       WRITE(14,*)I,WCENTRAL,'! IPOINT, WAVE'
       WRITE(14,934)(BDUM(J),J=1,9)
       WRITE(13,*)I,WCENTRAL,'! IPOINT, WAVE'
       WRITE(13,934)(BDUM(J),J=1,9)

       print*,'NFILBIN,Range : ',NFILBIN,WFIL(1),WFIL(NFILBIN)
       print*,'NAV = ',NAV
          
       KNU0 = 0.0
       tmp=0.
       DO I=1,NAV 
        IF(TPOUT(IAV(I),JGAS,1).GT.KNU0)THEN
         KNU0 = TPOUT(IAV(I),JGAS,1)
        ENDIF
C        print*,I,IAV(I),WAVEIN(IAV(I)),FBIN(I)
        tmp=tmp+fbin(i)
       ENDDO      
C       print*,'Sum of weights = ',tmp

       DO 150 IP=1,NP
         PRESS = P(IP)
C         WRITE(6,*)IP,PRESS,' ! IP,PRESS'
         WRITE(14,*)IP,PRESS,' ! IP,PRESS'
         WRITE(13,*)IP,PRESS,' ! IP,PRESS'
         DO 200 IT=1,NT
          TEMP=T(IT)
C          WRITE(6,*)IT,TEMP,' ! IT,TEMP'
          WRITE(14,*)IT,TEMP,' ! IT,TEMP'
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

          WRITE(14,*)(K_G(K),K=1,NG)
          WRITE(13,*)(U(K),K=1,NTRAN)
          WRITE(13,*)(TRAN(K),K=1,NTRAN)
          WRITE(13,*)(TFIT(K),K=1,NTRAN)
          WRITE(14,*)ERR,IMOD,ICALC
          WRITE(13,*)ERR,IMOD,ICALC
200     CONTINUE
150    CONTINUE
100   CONTINUE

      CLOSE(14)
      CLOSE(13)

      PRINT*,'Aband_ktablec. Run OK'

500   FORMAT(1X,A10)

401   FORMAT(8X,F8.2)
402   FORMAT(8X,F8.2)
400   FORMAT(8X,F8.2)
403   FORMAT(10X,I5)
404   FORMAT(8X,I3)      
405   FORMAT(7X,I3,I3)


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

