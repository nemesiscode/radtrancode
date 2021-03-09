      PROGRAM TABLE_TRANS
C     ***************************************************************
C     Program to calculate the transmission for a set of pressures,
C     temperatures and path lengths. Path conditions (and laboratory
C     spectra are read in from Neil Bowles original measurements and
C     then synthetic transmission spectra generated for the same
C     conditions and compared.
C
C     Pat Irwin 26/6/17  Original version
C
C     ***************************************************************

      IMPLICIT NONE
      INCLUDE '../includes/constdef.f'
      INTEGER MG,MDATA,MGAS,MBIN,MPAR,ICALC
      PARAMETER (MG=10,MDATA=20,MGAS=20,MBIN=8000,MPAR=20)
      INTEGER IDGAS(MGAS),ISOGAS(MGAS),NGAS,NPOINT,I,J,NPOINT1
      REAL TWAVEN(MBIN,2),TPOUT(MBIN,MGAS,7),WCEN(MBIN)
      REAL G_ORD(MG),K_G(MG),DEL_G(MG),T1
      REAL U, CHISQ,SUM
      REAL P(MPAR),T(MPAR),PRESS,TEMP,LPMIN,LPMAX,TMIN,TMAX
      REAL KNU0,DELAD,Y0,EL,SFB,C1,C2,VSTART,VEND
      REAL VMIN,DELV,FWHM,STEP,TMP,QROT,ALAMBDA,ERR,Q
      REAL DP,DT,X,BDUM(9),TPART(4),T2(MPAR,MPAR)
      INTEGER NG,JGAS,NP,NT,NTRAN,IMOD,IAV(MBIN),IFORM
      INTEGER IP,IT,K,ISTEP,IMETHOD,IBIN,NAV,BANDTYP(MGAS)
      LOGICAL IODD
      INTEGER NODD,NEVEN,ifile,nfile,LE,LEN,L
      INTEGER IWAVE,KWAVE,NBIN,ISHAPE,ISHAPE1
      REAL WMIN,WMAX,WFWHM,WCENTRAL
      REAL WAVEIN(MBIN),FWHMIN(MBIN)
      REAL FBIN(MBIN),xpath,Terr,w1,trans,etrans
      CHARACTER*100 IPFILE,OPFILE1,OPFILE2,GFILE,FLIST
      CHARACTER*100 FTMP,XTMP
      CHARACTER*20 HEAD
      CHARACTER*10 BUFFER
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

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

      idiag=1

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
        IF(HEAD(1:2).EQ.'GL')THEN
         BANDTYP(I)=5
        ELSE
         BANDTYP(I)=0
        ENDIF
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
       print*,'Allowed bandtypes are '
       print*,'-1 = Goody-Voigt - Kim Strong Type'
       print*,' 0 = Goody-Voigt - extended'
       print*,' 1 = Goody-Voigt - Kam Sihra type (C1,C2)'
       print*,' 2 = Goody-Voigt - 2-EL, T**QROT'
       print*,' 3 = Goody-Voigt - 2-EL, tpart(4)'
       print*,' 4 = Karkoschka and Tomasko (2009)'
       print*,' 5 = Goody-Lorentz (Neil Bowles)'
105   CONTINUE

      CLOSE(12)


      DO J=1,4
       TPART(J)=0.0
      ENDDO
      QROT=1.5
       
      IF(BANDTYP(JGAS).LT.3.OR.BANDTYP(JGAS).GT.4)THEN   
       CALL PROMPT('Enter QROT : ')
       READ*,QROT   
      ELSEIF(BANDTYP(JGAS).EQ.3)THEN
       CALL PROMPT('Enter TPART(1-4) : ')
       READ*,(TPART(J),J=1,4)
      ENDIF

      CALL PROMPT('Enter Q (0=foreign, 1=self) : ')
      READ*,Q

      nav=1
      fbin(1)=1.

C     Read list of laboratory transmission spectra filenames.
      CALL PROMPT('Enter name of files list : ')
      READ(5,1)FLIST

C     Need number if names in list
      call prompt('Enter nfile : ')
      read(5,*)nfile

c      open(12,file=flist,status='old',readonly)
      open(12,file=flist,status='old',action='read')

      SUM=0.
      do ifile=1,nfile
       read(12,1)ipfile
       call file(ipfile,opfile1,'tct')
       print*,ipfile
C       print*,opfile1
c       open(13,file=ipfile,status='old',readonly)
       open(13,file=ipfile,status='old',action='read')
       open(14,file=opfile1,status='unknown')

       read(13,*)xpath
       write(14,*)xpath
       chisq=0.
       do i=1,npoint      
        read(13,*)TEMP,Terr
C        print*,TEMP,Terr
        read(13,*)PRESS
C        print*,PRESS
        read(13,1)ftmp
        if(i.eq.1)then
         print*,'xpath, press, temp = ',xpath,press,temp
        endif
C        print*,ftmp
        call remsp(ftmp)
        le=len(ftmp)
C        print*,'len = ',le
C       Need this bit of code as this line in Neil's file has an
C       annoying semi-colon at the end of the line that needs to be
C       ignored.
        DO 40 L=LE,1,-1
         IF(FTMP(L:L).eq.';')THEN
          LE=L-1
          GOTO 50
         END IF
40      CONTINUE
50      CONTINUE
        xtmp=ftmp(1:le)
C        print*,xtmp
        read(xtmp,*)w1,trans,etrans     
        iav(1)=i
C	convert mb to atm
        PRESS=PRESS/1013.0
C       Calculate number of molecules/cm2 per km
        U=MODBOLTZ*PRESS/TEMP
C       Calculate number of molecules/cm2 in path (1e-20 factor because
C       strengths already normalised) 
        U = 1e-20*U*xpath/1000.
        call calc_tau(bandtyp,tpout,jgas,nav,iav,fbin,U,
     1   PRESS,TEMP,q,qrot,
     1   tpart,T1)
C        print*,w1,twaven(i,1),trans,T1,etrans
        chisq = chisq + ((T1-trans)/etrans)**2
C        print*,U,PRESS,TEMP,T1
        write(14,*)TEMP,Terr
        write(14,*)PRESS
        etrans=0.
        write(14,*)w1,T1,etrans     
       enddo
       close(13)
       close(14)
       print*,'CHISQ, CHISQ/ny = ',chisq,chisq/float(npoint)
       SUM=SUM+CHISQ
      enddo
      print*,'Total : '
      print*,'SUM, SUM/(ny*nfile) = ',SUM,SUM/float(npoint*nfile)

      
      close(12)

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


