      PROGRAM DIRECT_KTABLE
C     ****************************************************************
C     Converts ASCII k-tables outputted by Band_ktable to
C     direct-access .kta file.
C
C     Pat Irwin	 Documented and debugged	15/1/04
C     ****************************************************************
      IMPLICIT NONE

      INTEGER MG,MBIN,MPAR,NLMAX
      PARAMETER (MG=21,MBIN=8000,MPAR=20,NLMAX=2000)
      INTEGER IDGAS,ISOGAS,NPOINT,I,LUNOUT
      REAL TWAVEN(MBIN,2),TPOUT(MBIN,7),WAVEIN(MBIN)
      REAL G_ORD(MG),K_G(MG),DEL_G(MG),FWHMIN(MBIN)
      REAL P(MPAR),T(MPAR),PRESS,TEMP,T2(MPAR,MPAR)
      REAL VSTART,VEND,XFAC
      REAL DELV,FWHM,STEP,QROT,ERR,Q
      REAL KSTORE(MBIN,MPAR,MPAR,MG)
      REAL ERRM(MBIN,MPAR,MPAR)
      REAL CONT(NLMAX),FAC(NLMAX),FRAC(MBIN)
      REAL WCENTRAL,W1,W2,V1,V2,X1,X2,SUM,WMIN,WMAX,WFWHM,TPART(4)
      INTEGER IMODM(MBIN,MPAR,MPAR),ILOOP,IBIN,NBIN,IREC
      INTEGER NG,NP,NT,NTRAN,IMOD,IMOD1,ICALC,BANDTYP
      INTEGER IP,IT,K,IS1,IS2,IWAVE,KWAVE,IFORM
      INTEGER I1,IP1,IT1,IAV
      CHARACTER*100 IPFILE
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      DATA G_ORD/0.013047, 0.067468, 0.160295, 0.283302, 0.425563,
     1          0.574437, 0.716698, 0.839705, 0.932532, 0.986953,
     2		0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./


      idiag=1

      LUNOUT=30

      CALL PROMPT('Enter name of input ASCII k-table file : ')
      READ(5,1)IPFILE
1     FORMAT(A)

      CALL FILE(IPFILE,IPFILE,'asc')

      OPEN(12,FILE=IPFILE,STATUS='OLD')

      CALL PROMPT('Old(0) or New(1)format? ')
      READ*,IFORM

      READ(12,1)IPFILE
      READ(12,*)IMOD
      IF(IFORM.EQ.1)READ(12,*)BANDTYP
      READ(12,*)IDGAS,ISOGAS
      READ(12,*)VSTART,VEND
      READ(12,*)DELV,FWHM,NPOINT
      READ(12,*)NP
      READ(12,*)(P(I),I=1,NP)
      READ(12,*)NT
      READ(12,*)(T(I),I=1,NT)
      READ(12,*)NG
      READ(12,*)(DEL_G(I),I=1,NG)
      READ(12,*)NTRAN
      READ(12,*)QROT,Q
      IF(IFORM.EQ.1)READ(12,*)(TPART(I),I=1,4)
C      print*,'OK to here'

      DO 100 I=1,NPOINT
       READ(12,*)I1,TWAVEN(I,1)
C       print*,I1,TWAVEN(I,1)

       IF(IFORM.EQ.0)THEN
        READ(12,932)TWAVEN(I,1),TWAVEN(I,2),TPOUT(I,1),
     1   TPOUT(I,2),TPOUT(I,3),TPOUT(I,4),TPOUT(I,5)
        READ(12,933)TPOUT(I,6),TPOUT(I,7)
       ELSE
        READ(12,934)TWAVEN(I,1),TWAVEN(I,2),TPOUT(I,1),
     1   TPOUT(I,2),TPOUT(I,3),TPOUT(I,4),TPOUT(I,5),
     2   TPOUT(I,6),TPOUT(I,7)
       ENDIF
       WAVEIN(I)=TWAVEN(I,1)
       FWHMIN(I)=TWAVEN(I,2)
       
       DO 150 IP=1,NP
         READ(12,*)IP1,PRESS
         DO 200 IT=1,NT
          READ(12,*)IT1,TEMP

          READ(12,*)(KSTORE(I,IP,IT,K),K=1,NG)
          READ(12,*)ERR,IMOD1,ICALC
          ERRM(I,IP,IT)=ERR
          IMODM(I,IP,IT)=IMOD1

200     CONTINUE
150    CONTINUE
100   CONTINUE

      CLOSE(12)

      CALL OPEN_K(LUNOUT,VSTART,NPOINT,DELV,FWHM,IDGAS,ISOGAS,
     1  NG,NP,NT,P,T,T2,DEL_G,G_ORD,IREC)

      DO 10 I=1,NPOINT

       DO 20 IP=1,NP
        DO 30 IT=1,ABS(NT)
         DO K=1,NG
          K_G(K)=KSTORE(I,IP,IT,K)
         ENDDO

         CALL WRITE_K(LUNOUT,NG,K_G,IREC)

30      CONTINUE
20     CONTINUE

10    CONTINUE 

      CLOSE(LUNOUT)


932   FORMAT(1X,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5)
933   FORMAT(1x,2(e12.5))
934   FORMAT(1X,F8.2,F8.2,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5,E12.5)

      END

