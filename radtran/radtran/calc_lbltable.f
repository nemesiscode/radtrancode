      PROGRAM CALC_LBLTABLE
C     $Id:
C---------------------------------------------------------------------------
C_TITLE:  CALC_LBLTABLE
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
C     IMPLICIT NONE
C     note dbcom defines the linedata base variables. it is not normally stored
C     in the same directory as the rest of the code
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/dbcom.f'
C--------------------------------------------------------------
      INCLUDE '../includes/lincom.f'
      INCLUDE '../includes/pathcom.f'
      INCLUDE '../includes/bincom.f'
      INCLUDE '../includes/lcocom.f'
C-----------------------------------------------------------------------------
      INTEGER LUN,LOOP,LUN0,MDATA,MFIL,NGMAX,N1
      PARAMETER (LUN=2,LUN0=30,MFIL=1000)
      PARAMETER (NGMAX=51)
      REAL X,PMIN,PMAX,TMIN,TMAX,DT,DP
      REAL VMINC,VMAXC
      INTEGER IREC,IREC0,I,IWAVE,NFIL,IMULTI
      INTEGER IEXO,IPTF
      real tot_time
      real rate
      integer c1,c2,cr,time1,time2,cm
      CHARACTER*100 KTAFIL,OPFILE1,LCOFIL
      CHARACTER*1 ANS
      LOGICAL FEXIST
      REAL KNU0,DELAD,Y0,EL
      REAL QROT,FRAC1,MAXDV

      REAL P1,TE1,TEMP1(MAXK),PRESS1(MAXK),VCEN(MPOINT)
      REAL TEMP2(MAXK,MAXK),OUTPUT(MPOINT)
      REAL TABLE(MAXK,MAXK,MAXG)
      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),ERRK(MAXG)
C     G_ORD: Gauss-Legendre ordinates for calculating the k-distribution.
C     K_G: Calculated k-distribution.
C     DEL_G: Gauss-Legendre weights for integration.
C     **** all these are now defined by zgauleg.f ***
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      CALL system_clock(count_rate=cr)
      CALL system_clock(count_max=cm)
      rate = REAL(cr)
      call system_clock(time1)

      IJLCO=0

      idiag=1

      PRINT*,'Enter IEXO, IPTF'
      PRINT*,'(IEXO=1 uses temperature dependent line databases'
      PRINT*,'relevant for exoplanet k-tables spanning a very large'
      PRINT*,'temperature range. Set to 0 otherwise.)'
      PRINT*,'(IPTF=1 uses partition functions for gases listed in'
      PRINT*,'partfextra.dat file. IPTF=0 uses default partition'
      PRINT*,'functions listed in gasinfo file'
      CALL PROMPT('Enter IEXO, IPTF : ')
      READ*,IEXO,IPTF


      CALL PROMPT('Use Wavelengths(0) or Wavenumbers(1): ')
      READ*,IWAVE

      CALL PROMPT('Enter minumum wavenumber/wavelength X: ')
      READ*,VMIN
      CALL PROMPT('Enter DELX and NPOINT : ')
      READ*,DELV,NPOINT
      VMAX=VMIN + (NPOINT-1)*DELV
      PRINT*,' XMIN -> XMAX by DELX: ',VMIN,VMAX,DELV

      IF(NPOINT.GT.MPOINT)THEN 
       print*,'NPOINT exceeds MPOINT',NPOINT,MPOINT
       print*,'Increase MPOINT and recompile'
       stop
      ENDIF

      VMIN1=VMIN
      VMAX1=VMAX
      IF(IWAVE.EQ.0)THEN
       VMIN1=1E4/VMAX
       VMAX1=1E4/VMIN
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
      IF(NP.GT.1)THEN
       DP=(PMAX-PMIN)/(NP-1)
      ELSE
       DP=0.
      ENDIF
      DO 5 J=1,NP
        X=PMIN + (J-1)*DP
        PRESS1(J)=EXP(X)
        PRINT*,J,press1(J)
5     CONTINUE

      PRINT*,'Entering NT < 0 activates new pressure-dependent'
      PRINT*,'temperatures'
      CALL PROMPT('Enter number of temperature points ( < 20 ) : ')
      READ*,NT
      IF(NT.GT.0.AND.NT.LE.MAXK)THEN
       CALL PROMPT('Enter Tmin, Tmax : ')
       READ*,TMIN,TMAX
       IF(NT.GT.1)THEN
        DT=(TMAX-TMIN)/(NT-1)
       ELSE
        DT=0.
       ENDIF
       DO 6 J=1,NT
         TEMP1(J)=TMIN + (J-1)*DT
         PRINT*,J,TEMP1(J)
6      CONTINUE
      ELSE
       IF(NT.LT.0)THEN
        DO 7 I=1,NP
         PRINT*,'Pressure level = ',I,PRESS1(I)
         CALL PROMPT('Enter Tmin, Tmax : ')
         READ*,TMIN,TMAX
         DT=(TMAX-TMIN)/(ABS(NT)-1)
         DO 8 J=1,ABS(NT)
          TEMP2(I,J)=TMIN+(J-1)*DT
8        CONTINUE
7       CONTINUE
       ENDIF
      ENDIF

      WRITE(*,*)'Enter fractional abundance of absorber'
      WRITE(*,*)'0.0 will set the line width to be completely foreign'
      WRITE(*,*)'broadened. 1.0 will set the line width to be'
      CALL PROMPT('completely self-broadened : ')
      READ*,FRAC1
      Q=FRAC1

      CALL PROMPT('Enter line wing cut-off (cm^-1) : ')
      READ*,MAXDV

      IF(MAXDV.NE.VREL)THEN
       PRINT*,'*Warning from calc_lbltable.f: Cut-off and VREL'
       PRINT*,'are different. It is usual to set these to be the same.'
       PRINT*,'V_CUTOFF = ',MAXDV
       PRINT*,'VREL = ',VREL
      ENDIF

      IF(IEXO.EQ.0)THEN
       CALL PROMPT('Enter database name: ')
       READ(5,23)OPFILE
23     FORMAT(A)
       CALL FILE(OPFILE,KEYFIL,'key')
       CALL RDKEY(LUN)
       CALL RDGAS
       CALL RDISO

       CALL FILE(OPFILE,LCOFIL,'lco')
       INQUIRE(FILE=LCOFIL,EXIST=FEXIST)
       IF(FEXIST)THEN
           CALL INIT_LCO(LCOFIL)
       ENDIF

      ELSE
C      If IEXO<>0, then we need to read in temperature-dependent database
C      (for exoplanet k-tables)
       CALL PROMPT('Enter database root name: ')
       READ(5,23)OPFILE1
      ENDIF
      

      CALL PROMPT('Enter output filename : ')
      READ(5,23)OPFILE
      CALL FILE(OPFILE,KTAFIL,'lta')

      IRECL=ISYS()
      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=IRECL)
      print*,'NT = ',NT
      IF(NT.LT.0)THEN
       IREC0=10 + NP + NP*ABS(NT) + 2
      ELSE
       IREC0=10 + NP + NT + 2
      ENDIF
      PRINT*,'IREC0 = ',irec0
      print*,IREC0,NPOINT,VMIN,DELV
      print*,NP,NT,IDGAS(1),ISOGAS(1)

      WRITE(LUN0,REC=1)IREC0
      WRITE(LUN0,REC=2)NPOINT
      WRITE(LUN0,REC=3)VMIN
      WRITE(LUN0,REC=4)DELV
      WRITE(LUN0,REC=5)NP
      WRITE(LUN0,REC=6)NT
      WRITE(LUN0,REC=7)IDGAS(1)
      WRITE(LUN0,REC=8)ISOGAS(1)
      IREC=9
      DO 301 J=1,NP
        WRITE(LUN0,REC=IREC)PRESS1(J)
        IREC=IREC+1
301   CONTINUE
      IF(NT.GT.0)THEN
        DO 302 J=1,NT
          WRITE(LUN0,REC=IREC)TEMP1(J)
          IREC=IREC+1
302     CONTINUE
      ELSE
       DO 307 I=1,NP
        DO 304 J=1,ABS(NT)
          WRITE(LUN0,REC=IREC)TEMP2(I,J)
          IREC=IREC+1         
304     CONTINUE
307    CONTINUE
      ENDIF

C     Calculate continuum absorption for all pressures/temperatures

      PRINT*,'Calculating Continuum' 
      PRINT*,'VMIN1,VMAX1,WING,VREL',VMIN1,VMAX1,WING,VREL

      N1=ABS(NT)
      DO 102 K=1,N1

         IF(IEXO.NE.0)THEN
C          Read in temperature specific linedata file.
           IF(NT.LT.0)THEN
             PRINT*,'IEXO <> 0 AND NT<0 not yet implemented'
             STOP
           ENDIF

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

            P1=PRESS1(J)
            IF(NT.LT.0)THEN
             TE1=TEMP2(J,K)
            ELSE
             TE1=TEMP1(K)
            ENDIF

            PRINT*,'Continuum. J,K = ',J,K,' P = ',P1,' T = ',
     & 	    TE1
   
C            WRITE(*,*)'CALLING lbl_kcont'
            CALL LBL_KCONT(VMIN1,VMAX1,WING,VREL,P1,TE1,
     1      IDGAS(1),ISOGAS(1),FRAC1,IPROC(1),J,K,MAXDV,IPTF,IEXO)
C            WRITE(*,*)'lbl_kcont COMPLETE'
C            WRITE(*,*)' '

105      CONTINUE
102   CONTINUE

      PRINT*,'Continuum OK'

      IREC=IREC0

      PRINT*,'Calculating LBL spectra and k-distribution'

      N1=ABS(NT)
      DO 30 K=1,N1

          IF(IEXO.NE.0)THEN
C            Read in temperature specific line data file
             IF(NT.LT.0)THEN
               STOP
             ENDIF
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
            IF(NT.LT.0)THEN
             TE1=TEMP2(J,K)
            ELSE
             TE1=TEMP1(K)
            ENDIF

            WRITE(*,*)'Pressure, temperature: ',P1,TE1
            CALL LBL_KNOW(IWAVE,VMIN,DELV,NPOINT,P1,TE1,
     1          IDGAS(1),ISOGAS(1),IPROC(1),J,K,FRAC1,MAXDV,IPTF,
     2		OUTPUT)

             DO I=1,NPOINT
              IREC=IREC0 + (I-1)*NP*N1+(J-1)*N1+K-1
C              print*,I,VMIN+(I-1)*DELV,IREC0,IREC,OUTPUT(I)
              WRITE(LUN0,REC=IREC)OUTPUT(I)
             ENDDO

20        CONTINUE
30    CONTINUE


C-------------------------------------------------------------------------
C
C	Close files and shut down the program.
C
C-------------------------------------------------------------------------

      CALL WTEXT('%CALC_LBLTABLE.f :: calculation complete')
      call system_clock(time2)
      tot_time=(time2-time1)/rate
      WRITE(*,244)tot_time
244   FORMAT('%Elapsed time (s) = ',F8.1)

      CLOSE(IFIL)
      STOP

      END
************************************************************************
************************************************************************
