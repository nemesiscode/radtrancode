      PROGRAM PAR_KTABLE
C     **************************************************************
C     Program to convert a band_ktable ASCII output file into a
C     more formatted .par file
C
C     Pat Irwin		29/1/04
C
C     **************************************************************

      IMPLICIT NONE

      INTEGER MG,MDATA,MBIN,MPAR,ICALC
      PARAMETER (MG=10,MDATA=20,MBIN=8000,MPAR=20)
      INTEGER IDGAS,ISOGAS,NPOINT,I,J,NPOINT1
      INTEGER LUN0,LUN1,ITAB,IS1,IS2
      REAL X1,X2
      REAL TWAVEN(MBIN,2),TPOUT(MBIN,7),VV,TE1,PE1,V1,V2
      REAL G_ORD(MG),K_G(MG),DEL_G(MG),KSTORE(MPAR,MPAR,MG)
      REAL ERRSTORE(MPAR,MPAR),DELV,FWHM,QROT,Q,TPART(4)
      REAL U(MDATA),TRAN(MDATA),TFIT(MDATA),SIG(MDATA),CHISQ
      REAL P(MPAR),T(MPAR),PRESS,TEMP,LPMIN,LPMAX,TMIN,TMAX
      REAL KNU0,DELAD,Y0,EL,SFB,C1,C2,VSTART,VEND
      INTEGER TMOD(8000,20,20),TCALC(8000,20,20),IFORM
      INTEGER NG,NP,NT,NTRAN,IMOD,IMOD1,IBAND,BANDTYP
      INTEGER IP,IT,K,IMETHOD,I1,IP1,IT1
      CHARACTER*100 IPFILE,OPFILE1,OPFILE2
      CHARACTER*100 TNAME
      CHARACTER*20 HEAD
      CHARACTER*10 BUFFER
      CHARACTER*1 ANS
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      DATA G_ORD/0.013047, 0.067468, 0.160295, 0.283302, 0.425563,
     1          0.574437, 0.716698, 0.839705, 0.932532, 0.986953/

      DATA DEL_G/0.033336, 0.074726, 0.109543, 0.134633, 0.147762,
     1          0.147762, 0.134633, 0.109543, 0.074726, 0.033336/


      idiag=1

      CALL PROMPT('Enter input ASCII k-table name : ')
      READ(5,1)OPFILE1
1     FORMAT(A)
      CALL FILE(OPFILE1,OPFILE1,'asc')

      CALL PROMPT('Old(0) or New(1) format? ')
      READ*,IFORM

      CALL PROMPT('Enter output ASCII .par name : ')
      READ(5,1)OPFILE2
      CALL FILE(OPFILE2,OPFILE2,'par')


      LUN0 = 12
      LUN1 = 13
      OPEN(LUN0,FILE=OPFILE1,STATUS='OLD')
      OPEN(LUN1,FILE=OPFILE2,STATUS='UNKNOWN')


      READ(LUN0,1)IPFILE
      READ(LUN0,*)IMOD
      IF(IFORM.EQ.1)READ(LUN0,*)BANDTYP
      READ(LUN0,*)IDGAS,ISOGAS
      READ(LUN0,*)VSTART,VEND
      READ(LUN0,*)DELV,FWHM,NPOINT
      READ(LUN0,*)NP
      READ(LUN0,*)(P(I),I=1,NP)
      READ(LUN0,*)NT
      READ(LUN0,*)(T(I),I=1,NT)
      READ(LUN0,*)NG
      READ(LUN0,*)(DEL_G(I),I=1,NG)
      READ(LUN0,*)NTRAN
      READ(LUN0,*)QROT,Q
      IF(IFORM.EQ.1)READ(LUN0,*)(TPART(I),I=1,4)       

      PRINT*,'Input range is ',VSTART,VEND

      CALL PROMPT('Enter required output range : ')
      READ*,X1,X2

      CALL PROMPT('Output band parameters (Y/N)?')
      READ(5,1)ANS
      IBAND=0
      IF(ANS.EQ.'y'.OR.ANS.EQ.'Y')IBAND=1

      IS1 = 1+INT((X1-VSTART)/DELV)
      IS2 = 1+INT((X2-VSTART)/DELV)

      V1 = VSTART + (IS1-1)*DELV
      V2 = VSTART + (IS2-1)*DELV

      print*,IS1,IS2,V1,V2
      NPOINT1 = 1 + IS2-IS1
      CALL PROMPT('Enter title of data : ')
      READ(5,1)TNAME
      WRITE(LUN1,1)TNAME
      CALL PROMPT('Enter date : ')
      READ(5,1)TNAME
      WRITE(LUN1,*)'Patrick Irwin, AOPP, University of Oxford ',TNAME
      WRITE(LUN1,*)'VMIN,DELV,FWHM = ',V1,DELV,FWHM
      WRITE(LUN1,*)'NPOINT,VMAX = ',NPOINT1,V2
      WRITE(LUN1,*)'Gas ID,ISO : ',IDGAS,ISOGAS

      WRITE(LUN1,*)'Number of g-ordinates: ',NG
      WRITE(LUN1,*)'G_ORD, DEL_G'
      DO I=1,NG
       WRITE(LUN1,*)G_ORD(I),DEL_G(I)
      ENDDO

      WRITE(LUN1,*)'Number of pressures : ',NP 
      WRITE(LUN1,*)'Pressures (atm)'
      DO I=1,NP
       WRITE(LUN1,*)P(I)
      ENDDO

      WRITE(LUN1,*)'Number of temperatures : ',NT
      WRITE(LUN1,*)'Temperatures (K)'
      DO I=1,NT
       WRITE(LUN1,*)T(I)
      ENDDO
      WRITE(LUN1,*)
      WRITE(LUN1,*)'Absorption coefficients are in (km-amagat)^-1'
      WRITE(LUN1,*)
      CALL PROMPT('Data tab. by wavenumber (0) or wavelength (1)? : ')
      READ*,ITAB
      WRITE(LUN1,*)'Band parameter format:'
      IF(IFORM.EQ.0)THEN
       WRITE(LUN1,*)'First line contains V,DELV,KNU0,DELTA/AD,Y0,
     1EL,SFB'
       WRITE(LUN1,*)'Second line contains two far-wing continuum
     1coefficents C1,C2'
      ELSE
       WRITE(LUN1,*)'Line contains V,DELV,KNU0,DELTA/AD,Y0,E1,SFB,
     1E2,A'
      ENDIF
      WRITE(LUN1,*)

      WRITE(LUN1,*)'K-table format:'
      WRITE(LUN1,*)'For each pressure/temperature there are two lines:'
      WRITE(LUN1,*)'  The first line contains PRESS,(K_G(J),J=1,5)'
      WRITE(LUN1,*)'  The second line contains ERR,(K_G(J),J=6,10)'
      WRITE(LUN1,*)'N.B. ERR is the standard deviation of the'
      WRITE(LUN1,*)'     transmission fitting error' 

      DO 100 I=1,NPOINT
       READ(LUN0,*)I1,TWAVEN(I,1)
       VV = TWAVEN(I,1)
C       PRINT*,'I1',I1,TWAVEN(I,1)
       IF(IFORM.EQ.0)THEN
        READ(LUN0,932)TWAVEN(I,1),TWAVEN(I,2),TPOUT(I,1),
     1   TPOUT(I,2),TPOUT(I,3),TPOUT(I,4),
     2   TPOUT(I,5)
        READ(LUN0,933)TPOUT(I,6),TPOUT(I,7)
       ELSE
        READ(LUN0,934)TWAVEN(I,1),TWAVEN(I,2),TPOUT(I,1),
     1   TPOUT(I,2),TPOUT(I,3),TPOUT(I,4),
     2   TPOUT(I,5),TPOUT(I,6),TPOUT(I,7)
       ENDIF

       DO 150 IP=1,NP
         READ(LUN0,*)IP1,PRESS
C         print*,'IP1',IP1,PRESS
         DO 200 IT=1,NT
          READ(LUN0,*)IT1,TEMP
C          PRINT*,'IT1',IT1,TEMP
          READ(LUN0,*)(KSTORE(IP,IT,K),K=1,NG)
          READ(LUN0,*)ERRSTORE(IP,IT),IMOD1,ICALC
200     CONTINUE
150    CONTINUE

       IF(I.GE.IS1.AND.I.LE.IS2)THEN
        DO 301 IT=1,NT
         TE1 = T(IT)
         WRITE(LUN1,*)' '
         IF(ITAB.EQ.0)THEN
           WRITE(LUN1,996)VV,TE1
         ELSE
           WRITE(LUN1,997)VV,TE1
         ENDIF
996      FORMAT('          RESULTS FOR',F7.1,' CM-1,',F6.1,
     &' DEGREES K')
997      FORMAT('          RESULTS FOR',F7.4,' MICRON,',F6.1,
     &' DEGREES K')
         WRITE(LUN1,*)' '
         IF(IBAND.EQ.1)THEN
          WRITE(LUN1,*)'Band parameters'
          WRITE(LUN1,932)TWAVEN(I,1),TWAVEN(I,2),TPOUT(I,1),
     1    TPOUT(I,2),TPOUT(I,3),TPOUT(I,4),
     2    TPOUT(I,5)
          WRITE(LUN1,933)TPOUT(I,6),TPOUT(I,7)
         ENDIF

         WRITE(LUN1,*)'Pressure, fitted K-coeffs and transmission error'

         DO 401 IP=1,NP
          PE1 = P(IP)
          WRITE(LUN1,707) PE1,(26850.0*KSTORE(IP,IT,K),K=1,5)
          WRITE(LUN1,708)ERRSTORE(IP,IT),
     1      (26850.0*KSTORE(IP,IT,K),K=6,10)
401      CONTINUE

301     CONTINUE
 
       ENDIF

100   CONTINUE

707   FORMAT(1x,E11.4,5(E13.6))
708   FORMAT(3x,F9.7,5(E13.6))

      CLOSE(LUN0)
      CLOSE(LUN1)

932   FORMAT(1x,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5)
933   FORMAT(1x,2(e12.5))
934   FORMAT(1x,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5)

      END
