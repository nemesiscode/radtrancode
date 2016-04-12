      PROGRAM CALC_FNKTABLEC
C     $Id:
C---------------------------------------------------------------------------
C_TITLE:  CALC_FNKTABLEC
C
C_ARGS:
C
C_KEYS:
C
C_DESCR:  calculates the absorption coefficient look-up table for a gas of
C	  the user's choice.
C
C_FILES:  UNIT LUN: output file
C
C_CALLS:  PROMPT	Similar to WTEXT, or PRINT in that the argument is
C			written to the screen.
C	  FILE		Forces correct VMS style file extension for a
C			filename. i.e. assumes a <4 character extension
C			after a dot separator.
C	  RDKEY
C	  RDGAS
C	  RDISO
C	  LBL_CONT
C	  LBL_FKNEW	calculates the cumulative K-Distribution for a
C			spectral interval for a mixture of gases at a
C			number of pressures and temperatures. This is done
C			by first generating the lbl absorption coefficient
C			spectra and then analysing these according to the 
C			equation given by Lacis and Oinas (1991):
C				f(k) = (1/(V2-V1))*SUM(ABS(dV/DK))
C			f(k) is then summed to give the cumulative k
C			distribution.
C	  WTEXT		Similar to PROMPT, or PRINT in that the argument
C			is written to the screen.
C	  REMSP		Removes leading spaces from text string.
C
C_BUGS:   
C
C_HIST:	30.3.00		PGJI	ORIGINAL VERSION
C	26/4/12		PGJI	Updated for new radtrans/nemesis and added 
C				option to calculate k-table for spectral 
C				channels.
C   
C---------------------------------------------------------------------------
C     note dbcom defines the linedata base variables. it is not normally stored
C     in the same directory as the rest of the code
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
      INCLUDE '../includes/lincom.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/bincom.f'
C-----------------------------------------------------------------------------
      INTEGER LUN,LOOP,LUN0,LUN1,MDATA,MFIL,NGMAX
      PARAMETER (LUN=2,LUN0=30,LUN1=31,MFIL=1000)
      PARAMETER (NGMAX=51)
      REAL X,PMIN,PMAX,TMIN,TMAX,DT,DP
      REAL VMINC,VMAXC
      INTEGER IREC,IREC0,I,IWAVE,NFIL,IMULTI
      INTEGER IEXO,IPTF
      real tot_time
      real rate
      integer c1,c2,cr,time1,time2,cm
      REAL TFIL(MFIL),VFIL(MFIL),V1(MFIL),V2(MFIL),T1(MFIL)
      CHARACTER*100 KTAFIL,FILFILE,OPFILE1
      CHARACTER*1 ANS
      REAL KNU0,DELAD,Y0,EL
      REAL QROT,FRAC1,MAXDV

      REAL P1,TE1,TEMP1(MAXK),PRESS1(MAXK),VCEN(MPOINT)
      REAL TABLE(MAXK,MAXK,MAXG)
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),ERRK(MAXG)
C     G_ORD: Gauss-Legendre ordinates for calculating the k-distribution.
C     K_G: Calculated k-distribution.
C     DEL_G: Gauss-Legendre weights for integration.
C     **** all these are now defined by zgauleg.f ***

      CALL system_clock(count_rate=cr)
      CALL system_clock(count_max=cm)
      rate = REAL(cr)
      call system_clock(time1)

      PRINT*,'Enter IEXO, IPTF'
      PRINT*,'(IEXO=1 uses temperature dependent line databases'
      PRINT*,'relevant for exoplanet k-tables spanning a very large'
      PRINT*,'temperature range. Set to 0 otherwise.)'
      PRINT*,'(IPTF=1 uses partition functions for gases listed in'
      PRINT*,'partfextra.dat file. IPTF=0 uses default partition'
      PRINT*,'functions listed in gasinfo file'
      CALL PROMPT('Enter IEXO, IPTF : ')
      READ*,IEXO,IPTF

      CALL PROMPT('Enter NG : ')
      READ*,NG

c  ** calc g_ord and del_g **
      call zgauleg(g_ord,del_g,ng,ngmax)


      CALL PROMPT('Use Wavelengths(0) or Wavenumbers(1): ')
      READ*,IWAVE

      PRINT*,'For bin centres do you want a regular wave step (0)'
      CALL PROMPT('or irregular (1) : ')
      READ*,IREG

      IF(IREG.EQ.0)THEN
       CALL PROMPT('Enter minumum wavenumber/wavelength X: ')
       READ*,VMIN
       CALL PROMPT('Enter DELX and NPOINT : ')
       READ*,DELV,NPOINT
       VMAX=VMIN + (NPOINT-1)*DELV
       PRINT*,' XMIN -> XMAX by DELX: ',VMIN,VMAX,DELV
       DO I=1,NPOINT
        VCEN(I)=VMIN+(I-1)*DELV
       ENDDO
      ELSE
       CALL PROMPT('Enter number of points : ')
       READ*,NPOINT
       CALL PROMPT('Enter central wavenumbers/wavelengths : ')
       READ*,(VCEN(I),I=1,NPOINT)
       VMIN=VCEN(1)
       VMAX=VCEN(NPOINT)
       DELV=-1.0
      ENDIF

      CALL PROMPT('Square bins (0), or non-square(1) : ')
      READ*,SQBIN

      IF(SQBIN.EQ.0)THEN
       CALL PROMPT('Enter FWHM : ')
       READ*,FWHM

C      Calculate min/max wavelength/wavenumbers for continuum calculation
       IF(IWAVE.EQ.0)THEN
          VMIN1 = 1E4/(VMAX+0.5*FWHM)
          VMAX1 = 1E4/(VMIN-0.5*FWHM)
       ELSE
          VMIN1=VMIN-0.5*FWHM
          VMAX1=VMAX+0.5*FWHM
       ENDIF

       print*,'test: VMIN1,VMAX1 = ',VMIN1,VMAX1
       
      ELSE

       FWHM=0.0

C      Read min/max wavelength/wavenumbers for continuum calculation
       PRINT*,'Enter min,max wavenumber/wavelengths for'
       PRINT*,'Continuum calculation. These should cover the range'
       PRINT*,'VMIN-(extreme lower spread of filter width) to'
       PRINT*,'VMAX+(extreme upper spread of filter width).' 
       PRINT*,'(N.B. it is better to overestimate rather than'
       PRINT*,'underestimate these values)'
       CALL PROMPT('Enter values : ')
       READ*,VMINC,VMAXC
C      Convert wavelength range to wavenumber range if IWAVE=0
       IF(IWAVE.EQ.0)THEN
          VMIN1 = 1E4/VMAXC
          VMAX1 = 1E4/VMINC
       ELSE
          VMIN1 = VMINC
          VMAX1 = VMAXC
       ENDIF

       print*,'test: VMIN1,VMAX1 = ',VMIN1,VMAX1

       PRINT*,'Same bin shape for all output points(0) or'
       CALL PROMPT('new filter profile each time(1) : ')
       READ*,IMULTI

       CALL PROMPT('Enter name of filter function file : ')
       READ(5,23)FILFILE
       IFIL=13
       OPEN(IFIL,FILE=FILFILE,STATUS='OLD')
       READ(IFIL,*)NCHAN
       IF(IMULTI.EQ.1)THEN
        IF(NCHAN.NE.NPOINT)THEN
         PRINT*,'Not same number of filter functions as wavecentres'
         PRINT*,'NCHAN,NPOINT',NCHAN,NPOINT
         STOP
        ENDIF
       ENDIF

      ENDIF 
 
      CALL PROMPT('Enter gas ID,ISO,IPROC : ')
      READ*,IDGAS(1),ISOGAS(1),IPROC(1)
      NGAS=1

      CALL PROMPT('Enter WING,VREL : ')
      READ*,WING,VREL

      CALL PROMPT('Enter number of pressure points ( < 20 ) : ')
      READ*,NP
      CALL PROMPT('Enter log(pmin), log(pmax) : ')
      READ*,PMIN,PMAX
      DP=(PMAX-PMIN)/(NP-1)

      DO 5 J=1,NP
        X=PMIN + (J-1)*DP
        PRESS1(J)=EXP(X)
        PRINT*,J,press1(J)
5     CONTINUE

      CALL PROMPT('Enter number of temperature points ( < 20 ) : ')
      READ*,NT
      CALL PROMPT('Enter Tmin, Tmax : ')
      READ*,TMIN,TMAX
      DT=(TMAX-TMIN)/(NT-1)

      DO 6 J=1,NT
        TEMP1(J)=TMIN + (J-1)*DT
        PRINT*,J,TEMP1(J)
6     CONTINUE

      WRITE(*,*)'Enter fractional abundance of absorber'
      WRITE(*,*)'0.0 will set the line width to be completely foreign'
      WRITE(*,*)'broadened. 1.0 will set the line width to be'
      CALL PROMPT('completely self-broadened : ')
      READ*,FRAC1
      Q=FRAC1

      CALL PROMPT('Enter line wing cut-off (cm^-1) : ')
      READ*,MAXDV

      IF(IEXO.EQ.0)THEN
       CALL PROMPT('Enter database name: ')
       READ(5,23)OPFILE
23     FORMAT(A)
       CALL FILE(OPFILE,KEYFIL,'key')
       CALL RDKEY(LUN)
       CALL RDGAS
       CALL RDISO
      ELSE
C      If IEXO<>0, then we need to read in temperature-dependent database
C      (for exoplanet k-tables)
       CALL PROMPT('Enter database root name: ')
       READ(5,23)OPFILE1
      ENDIF
      

      CALL PROMPT('Enter output filename : ')
      READ(5,23)OPFILE
      CALL FILE(OPFILE,KTAFIL,'kta')

      IRECL=ISYS()
      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=IRECL)
      IREC0=11 + 2*NG + 2 + NP + NT + 2
C     Add in extra buffer to list wavelengths if a non-uniform grid is
C     specified
      IF(DELV.LE.0.0)THEN 
        IREC0=IREC0+NPOINT
      ENDIF
c     PRINT*,'IREC0 = ',irec0
      WRITE(LUN0,REC=1)IREC0
      WRITE(LUN0,REC=2)NPOINT
      WRITE(LUN0,REC=3)VMIN
      WRITE(LUN0,REC=4)DELV
      WRITE(LUN0,REC=5)FWHM
      WRITE(LUN0,REC=6)NP
      WRITE(LUN0,REC=7)NT
      WRITE(LUN0,REC=8)NG
      WRITE(LUN0,REC=9)IDGAS(1)
      WRITE(LUN0,REC=10)ISOGAS(1)
      IREC=11
      DO 299 J=1,NG
        WRITE(LUN0,REC=IREC)G_ORD(J)
        IREC=IREC+1
299   CONTINUE
      DO 399 J=1,NG
        WRITE(LUN0,REC=IREC)DEL_G(J)
        IREC=IREC+1
399   CONTINUE
      IREC=11 + 2*NG + 2
      DO 301 J=1,NP
        WRITE(LUN0,REC=IREC)PRESS1(J)
        IREC=IREC+1
301   CONTINUE
      DO 302 J=1,NT
        WRITE(LUN0,REC=IREC)TEMP1(J)
        IREC=IREC+1
302   CONTINUE
C     Write out central wavelengths if non-uniform grid
      IF(DELV.LE.0.0)THEN
       DO 303 J=1,NPOINT
        WRITE(LUN0,REC=IREC)VCEN(J)
C        PRINT*,J,VCEN(J)
        IREC=IREC+1
303    CONTINUE
      ENDIF

C     Calculate continuum absorption for all pressures/temperatures

      PRINT*,'Calculating Continuum' 
      PRINT*,'VMIN1,VMAX1,WING,VREL',VMIN1,VMAX1,WING,VREL

      DO 102 K=1,NT

         IF(IEXO.NE.0)THEN
C          Read in temperature specific linedata file.
           OPFILE=OPFILE1
           DO I=1,LEN(OPFILE)
            IF(OPFILE(I:I).NE.' ')KK=I
           ENDDO
           I1 = INT(K/10)
           I2 = K-10*I1 
           OPFILE(KK+1:KK+1)=CHAR(I1+48)
           OPFILE(KK+2:KK+2)=CHAR(I2+48)
    
           CALL FILE(OPFILE,KEYFIL,'key')
           CALL RDKEY(LUN)
           CALL RDGAS
           CALL RDISO
         ENDIF


         DO 105 J=1,NP

            PRINT*,'Continuum. J,K = ',J,K,' P = ',PRESS1(J),' T = ',
     & 	    TEMP1(K)
   
cc            WRITE(*,*)'CALLING lbl_kcont'
            CALL LBL_KCONT(VMIN1,VMAX1,WING,VREL,PRESS1(J),TEMP1(K),
     1      IDGAS(1),ISOGAS(1),FRAC1,IPROC(1),J,K,MAXDV,IPTF,IEXO)
cc            WRITE(*,*)'lbl_kcont COMPLETE'
cc            WRITE(*,*)' '
105      CONTINUE
102   CONTINUE
 
      PRINT*,'Continuum OK'

      IREC=IREC0

      OPEN(UNIT=LUN1,FILE='tempk.dat',STATUS='UNKNOWN')
      WRITE(LUN1,*)VMIN,VMAX,DELV
      WRITE(LUN1,*)NP,(PRESS1(J),J=1,NP)
      WRITE(LUN1,*)NT,(TEMP1(J),J=1,NT)

      PRINT*,'Calculating LBL spectra and k-distribution'

      DO 10 IPOINT=1,NPOINT
        WRITE(*,*)'Current Wave: ',VCEN(IPOINT)
        WRITE(LUN1,*)'Current Wave: ',VCEN(IPOINT)

        IF(SQBIN.EQ.0)THEN
            VSTART = VCEN(IPOINT) - 0.5*FWHM
            VEND = VSTART+FWHM

            IF(IWAVE.EQ.0)THEN
                VV = VSTART
                VSTART = 1E4/VEND
                VEND = 1E4/VV
                NFIL=25
                DO I=1,NFIL
                 VFIL(I)=VSTART+FLOAT(I-1)*(VEND-VSTART)/FLOAT(NFIL-1)
                 TFIL(I)=(VFIL(I)/VFIL(1))**2
                ENDDO
            ELSE
                NFIL=2
                VFIL(1)=VSTART
                VFIL(2)=VEND
                TFIL(1)=1.0
                TFIL(2)=1.0
            ENDIF
        ELSE
            IF(IPOINT.EQ.1.OR.IMULTI.EQ.1)THEN
             READ(IFIL,*)VFILC
             READ(IFIL,*)NFIL
             IF(VFILC.NE.VCEN(IPOINT).AND.IMULTI.EQ.1)THEN
              PRINT*,'*** WARNING *** Calc_ktablec'
              PRINT*,'VFILC,VCEN not the same'
              PRINT*,VFILC,VCEN(IPOINT)
             ENDIF
             IF(NFIL.GT.MFIL)THEN
              PRINT*,'*** WARNING *** Calc_ktablec'
              PRINT*,'Instrument Function has NFIL>MFIL'
              PRINT*,'Aborting'
              STOP
             ENDIF
             DO I=1,NFIL
              READ(IFIL,*)V1(I),T1(I)
              V2(I)=V1(I)-VFILC
              IF(T1(I).LT.0.0)THEN
               PRINT*,'*** WARNING *** Calc_ktablec'
C               PRINT*,'Instrument Function has negative elements'
C               PRINT*,'Aborting'
C               STOP
              ENDIF 
             ENDDO
            ENDIF 

C           Centre filter
            DO I=1,NFIL
              V1(I) = V2(I)+VCEN(IPOINT)
            ENDDO

            IF(IWAVE.EQ.1)THEN
             DO I=1,NFIL
              VFIL(I)=V1(I)
              TFIL(I)=T1(I)
             ENDDO
            ELSE
C            Modify filter function if in wavelength space
             DO I=1,NFIL
               J=NFIL+1-I
               VFIL(I)=1E4/V1(J)
               TFIL(I)=T1(J)*(VFIL(I)/VFIL(1))**2
             ENDDO
            ENDIF

            VSTART=VFIL(1)
            VEND=VFIL(NFIL)

        ENDIF
        PRINT*,'VSTART,VEND = ',VSTART,VEND
C        PRINT*,'NFIL = ',NFIL
C        DO I=1,NFIL
C         PRINT*,VFIL(I),TFIL(I)
C        ENDDO
           
        DO 30 K=1,NT

          IF(IEXO.NE.0)THEN
C            Read in temperature specific line data file
             OPFILE=OPFILE1
             DO I=1,LEN(OPFILE)
              IF(OPFILE(I:I).NE.' ')KK=I
             ENDDO
             I1 = INT(K/10)
             I2 = K-10*I1 
             OPFILE(KK+1:KK+1)=CHAR(I1+48)
             OPFILE(KK+2:KK+2)=CHAR(I2+48)
    
             CALL FILE(OPFILE,KEYFIL,'key')
             CALL RDKEY(LUN)
             CALL RDGAS
             CALL RDISO
C            Read in the lines into the bins (bins have already been
C            defined by lbl_kcont.f) 
             CALL LOADBINS(WING,NGAS,IDGAS,ISOGAS)
             PRINT*,'Reading in lines for temperature : ',TEMP1(K)
          ENDIF

          DO 20 J=1,NP
            P1=PRESS1(J)
            TE1=TEMP1(K)

            WRITE(*,*)'Pressure, temperature: ',P1,TE1
            WRITE(LUN1,*)'Pressure, temperature: ',P1,TE1
            CALL LBL_FKNEW(IWAVE,VSTART,VEND,P1,TE1,
     1          IDGAS(1),ISOGAS(1),IPROC(1),J,K,FRAC1,MAXDV,IPTF,
     2		NPOINT)

            DELV=(VEND-VSTART)/FLOAT(NPOINT)

            CALL CALC_FKDIST_WAVEC(IWAVE,VSTART,DELV,NPOINT,
     1    NFIL,VFIL,TFIL,G_ORD,DEL_G,K_G,NGMAX,NG)

            WRITE(LUN1,*)(K_G(LOOP),LOOP=1,NG)

            DO 40 LOOP=1,NG
             TABLE(J,K,LOOP)=K_G(LOOP)          
40          CONTINUE

20        CONTINUE
30      CONTINUE


C       K-tables assume the loop is pressure followed by temperature
C       unlike the order we have followed here, which is temperature
C       followed by pressure. Hence, we need to reverse the order for 
C       output

        DO J=1,NP
         DO K=1,NT

          DO LOOP=1,NG
           K_G(LOOP)=TABLE(J,K,LOOP)
           WRITE(LUN0,REC=IREC)K_G(LOOP)
           IREC=IREC+1
          ENDDO

         ENDDO
        ENDDO

10    CONTINUE


C-------------------------------------------------------------------------
C
C	Close files and shut down the program.
C
C-------------------------------------------------------------------------
C If the code succeeds in completion, delete the "tempk.dat" file since 
C its main usefulness is when the code crashes prior to completion.
      CLOSE(UNIT= LUN1,STATUS= 'DELETE')

      CALL WTEXT('%CALC_KTABLEC.f :: calculation complete')
      call system_clock(time2)
      tot_time=(time2-time1)/rate
      WRITE(*,244)tot_time
244   FORMAT('%Elapsed time (s) = ',F8.1)

      CLOSE(IFIL)
      STOP

      END
************************************************************************
************************************************************************
