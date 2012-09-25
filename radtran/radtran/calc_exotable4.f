      PROGRAM CALC_EXOTABLE4
C     $Id:
C***********************************************************************
C_TITL:	CALC_EXOTABLE4.f
C
C_DESC: Calculates the absorption coefficient look-up table for a gas of
C       the user's choice. Code is similar to Calc_fnktable, but is modified
C       to use a range of line databases calculated at different temperatures 
C       by Makedbloop. This was necessary as some linedatabase, such as BT2
C       have far too many lines to deal with efficiently. It turns out that 
C       at different temperatures, different sets of the lines account for
C       most of the absorption. Hence, the non-negligible lines at each 
C       temperature are output to a separate database and the k-table for that
C       temperature calculated from that database.
C
C	Calc_exotable is ALMOST same as calc_exotable3. However, to
C       use Calc_exotable4, we need to create a blank k-table first 
C       and then run the code and fill it up with k-coefficients for 
C       each grid point. However,this process is totally manual, not 
C       automated, so is probably not a usable code for all users. 
C
C
C_ARGS:	See the definitions below.
C
C_FILE:	unit=LUN0	output file "<ktafil>.kta"
C
C_CALL:	FILE		Forces correct VMS style file extension for a
C			filename. i.e. assumes a <4 character extension
C			after a dot separator.
C	RDKEY		Reads in line data key file.
C	RDGAS		Reads in details of gases in data base.
C	RDISO		Reads in isotope details for line data base.
C	RDKEY_BAND	Reads in line data key file for band data.
C	LBL_KCONT_EXO3	Subroutine to calculate the line continuum in the
C			range VMIN - VMAX in bins of width WING.
C	LBL_FKNEW_EXO3	Calculates the cumulative K-Distribution for a
C			spectral interval for a mixture of gases at a
C			number of pressures and temperatures. This is done
C			by first generating the lbl absorption coefficient
C			spectra and then analysing these according to the 
C			equation given by Lacis and Oinas (1991):
C				f(k) = (1/(V2-V1))*SUM(ABS(dV/DK))
C			f(k) is then summed to give the cumulative k
C			distribution.
C	REMSP		Removes leading spaces from text string.
C
C_HIST: 30.3.00	PGJI	ORIGINAL VERSION
C	09jun10 PGJI    Converted from calc_fnktable to deal with
C 			exoplanetary line databases which have so many
C			lines that the important lines for each table
C			temperature must be selected first with
C			SELECTTEMPLOOP
C	00nnn00 JML	Adapted from calc_exotable3.f
C	26apr12 PGJI	Updated and commented.
C
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/bincom.f'
C ../includes/bincom.f stores the line bin variables (including NLINES,
C FSTLIN, LSTLIN) and band parameters.
      INCLUDE '../includes/dbcom.f'
C ../includes/dbcom.f stores the line database variables.
      INCLUDE '../includes/lincom.f'
C ../includes/lincom.f stores the linedata variables (including MAXLIN,
C VLIN, SLIN, ALIN, ELIN, SBLIN, TDW, TDWS and that lot).
      INCLUDE '../includes/pathcom.f'
C ../includes/parcom.f stores the parameter values such as MAXLAY,

      INTEGER I,II,IV,IV0,J,K,LI,LJ,IWAVE
      INTEGER IREC,IREC0,IRECL,ISYS
      INTEGER IDGAS1,ISOGAS1,IPROC1
C IDGAS1: The local gas identifier.
C ISOGAS1: The local gas-isotopic identifier; if zero all isotopes of the
C gas are included.
C IPROC1: Line wing identifier.

      INTEGER NP,NT,LOOP
C NP: Number of pressures.
C NT: Number of temperatures.
      INTEGER LUN,LUN0,LUN1,MFIL,NFIL
      PARAMETER (LUN=2,LUN0=30,LUN1=31,MFIL=1000)


      REAL K_BOLTZ,C_LIGHT,R_GAS,PI,VFIL(MFIL),TFIL(MFIL)
      PARAMETER (K_BOLTZ=1.38E-23,C_LIGHT=2.998E8,R_GAS=8.3143)        
      PARAMETER (PI=3.1415927)

C K_BOLTZ: Boltzman constant [J K^-1 molecule^-1].
C C_LIGHT: Speed of light [m sec^-1].
C R_GAS: Universal gas constant [J K^-1 mol^-1].
      REAL QROT
      INTEGER MDATA,NG,NGMAX,BOTREC,STAREC,NXTREC
      PARAMETER (MDATA=20,QROT=1.5,NG=20)
C NG: Number of ordinates in k-distribution.

      REAL TOT_TIME
C TOT_TIME: The total system time it took to complete program execution.
      DOUBLE PRECISION TIME,TIME1,TIME2
C TIME: Temporary variable returned by GETTIME containing the system time.
C TIME1: System time at the beginning of program execution.
C TIME2: System time at the end of program execution.

      REAL VSTART,VEND,DELVSF,V1
C VSTART: Beginning of spectral range wavenumber [cm-1].
C VEND: End of spectral range wavenumber [cm-1].
C DELVSF: DELV scale factor. Used to avoid non-integer DO loops.

      REAL VMAX,PMIN,PMAX,TMIN,TMAX,VMIN1,VMAX1,VMINX
C VMAX: Wavenumber [cm-1] maximum (VMIN is already declared elsewhere --
C likely in some COMMON block in one of those ../include files).
C PMIN: Pressure [atm] minimum.
C PMAX: Pressure [atm] maximum.
C TMIN: Temperature [Kelvin] minimum.
C TMAX: Temperature [Kelbin] maximum.

      REAL P1,TE1,DT,DP,XP,TEMP1(MAXK),PRESS1(MAXK)
C P1: Pressure [atm] at level J.
C TE1: Temperature [Kelvin] at level K.
C PRESS1: Pressure [atm].
C TEMP1: Temperature [Kelvin].

      REAL FRAC,MAXDV
C FRAC: Required fraction (0=air broadened,1=self).
C MAXDV: Line wing cut-off parameter [cm-1]: The maximum line width away
C within which to consider the contribution of the line wings.

      REAL U,XE(MDATA),YE(MDATA),SIGE(MDATA)
      REAL SDES,SUM1,ALAMDA1,LNABSCO
      REAL KNU0,DELAD,Y0,EL,SFB,CB1,CB2

      REAL G_ORD(MAXG),K_G(MAXG),DEL_G(MAXG),ERRK(MAXG)
C G_ORD: Gauss-Legendre ordinates for calculating the k-distribution.
C K_G: Calculated k-distribution.
C DEL_G: Gauss-Legendre weights for integration.


      CHARACTER*100 KTAFIL,LINFIL,LINFIL1

      REAL XMASS,GETMASS,VBOT,VTOP,RANGE,TCORS1,TCORS2,TCORDW,TRATIO
      INTEGER ITEMP,IPRESS,INORM,FIRST,LAST,NLIN,LINE,LINE1
      INTEGER OFSTBIN,OLSTBIN,FSTBIN1,IBIN,JBIN,IPTF
      REAL XPW,MXWID,V2,V0,XPD,PARTF,DPEXP
      INTEGER I1,I2,IBIN1,IBIN2,IFULL,IBINMAX,ILIN
      REAL RANGE1,VS1,VE1,VFIRST,VLAST
      INTEGER NVMAX,NPOINTX
C      PARAMETER (NPOINT1=801)

C******************************** CODE *********************************

c  ** calc g_ord and del_g **
      ngmax=maxg
      call zgauleg(g_ord,del_g,ng,ngmax)

C Obtain the system time for use in determining the total time used by the
C program since execution.
      CALL GETTIME(TIME)
      TIME1 = TIME

C     Set NFIL to -1 to stop calc_fkdist_wavec trying to average over a
C     filter function

      NFIL=-1

      IPTF=0

      CALL PROMPT('Use Wavelengths(0) or Wavenumbers(1): ')
      READ*,IWAVE

      WRITE(*,*)'Enter minumum : '
      READ*,VMIN

      WRITE(*,*)'Enter FWHM, DELV and NPOINT : '
      READ*,FWHM,DELV,NPOINT
      VMAX = VMIN + (NPOINT - 1)*DELV

      WRITE(*,*)' VMIN --> VMAX by DELV: ',VMIN,VMAX,DELV

      WRITE(*,*)'Enter gas ID,ISO,IPROC : '
      READ*,IDGAS1,ISOGAS1,IPROC1


      NGAS=1
      IDGAS(1)=IDGAS1
      ISOGAS(1)=ISOGAS1

      WRITE(*,*)'Enter WING,VREL : '
      READ*,WING, VREL

      WRITE(*,*)'Enter number of pressure points ( <= 20 ) : '
      READ*,NP
      WRITE(*,*)'Enter log(pmin), log(pmax) : '
      READ*,PMIN,PMAX
      DP = (PMAX - PMIN)/(NP - 1)
      DO 5 J=1,NP
        XP = PMIN + (J - 1)*DP
        PRESS1(J) = EXP(XP)
        WRITE(*,*)J,PRESS1(J)
5     CONTINUE

      WRITE(*,*)'Enter number of temperature points ( <= 20 ) : '
      READ*,NT
      WRITE(*,*)'Enter Tmin, Tmax : '
      READ*,TMIN,TMAX
      DT = (TMAX - TMIN)/(NT - 1)
      DO 6 J=1,NT
        TEMP1(J) = TMIN + (J - 1)*DT
        WRITE(*,*)J,TEMP1(J)
6     CONTINUE

      WRITE(*,*)'---------- Hydrogen CIA Absorption details ----------'
      WRITE(*,*)'Enter fractional abundance of absorber'
      WRITE(*,*)'0.0 will set the line width to be completely foreign'
      WRITE(*,*)'broadened. 1.0 will set the line width to be'
      WRITE(*,*)'completely self-broadened : '
      READ*,FRAC

      WRITE(*,*)'Enter line wing cut-off (cm^-1) : '
      READ*,MAXDV

      WRITE(*,*)'Enter root database name : '
      READ(5,23)LINFIL1
23    FORMAT(A)

      WRITE(*,*)'Enter output filename : '
      READ(5,23)OPFILE
      CALL FILE(OPFILE,KTAFIL,'kta')
      IRECL = ISYS()
      OPEN(UNIT=LUN0,FILE=KTAFIL,STATUS='UNKNOWN',ACCESS='DIRECT',
     1 RECL=IRECL)
      
C JM  wavenumber space : 300-30000 cm-1, DELV=5 cm-1 (NPX=5941)
C                                        DELV=1 cm-1 (NPX=27001)
C JM  wavelength space : 0.3-30 um, DELV=0.005 um (NPX=5941)
	NPOINTX=5941
	VMINX=300.0
      IF(IWAVE.EQ.0)THEN
	NPOINTX=5941
	VMINX=0.3
      ENDIF

      IREC0 = 11 + 2*NG + 2 + NP*NT + 2
cc      WRITE(*,*)' CALC_EXOTABLE3.f :: IREC0 = ',irec0
      WRITE(LUN0,REC=1)IREC0
      WRITE(LUN0,REC=2)NPOINTX
      WRITE(LUN0,REC=3)VMINX
      WRITE(LUN0,REC=4)DELV
      WRITE(LUN0,REC=5)FWHM
      WRITE(LUN0,REC=6)NP
      WRITE(LUN0,REC=7)NT
      WRITE(LUN0,REC=8)NG
      WRITE(LUN0,REC=9)IDGAS1
      WRITE(LUN0,REC=10)ISOGAS1

      IREC = 11
      DO 299 J=1,NG
        WRITE(LUN0,REC=IREC)G_ORD(J)
        IREC = IREC + 1
299   CONTINUE
      DO 399 J=1,NG
        WRITE(LUN0,REC=IREC)DEL_G(J)
        IREC = IREC + 1
399   CONTINUE
      IREC = 11 + 2*NG + 2
      DO 301 J=1,NP
        WRITE(LUN0,REC=IREC)PRESS1(J)
        IREC = IREC + 1
301   CONTINUE
      DO 302 J=1,NT
        WRITE(LUN0,REC=IREC)TEMP1(J)
        IREC = IREC + 1
302   CONTINUE


      IREC = IREC0
      TE1 = 0.0
C     First generate table and pack with zeros so we can insert new
C     values anywhere
c      DO IV=1,NPOINT
c       DO J=1,NP
c        DO K=1,NT
c        DO LOOP=1,NG
c         WRITE(LUN0,REC=IREC)TE1
c        IREC=IREC+1
c        ENDDO
c        ENDDO
c       ENDDO
c      ENDDO
c	print*,irec
c	close(lun0)
c	stop
C Calculate continuum absorptions


       IF(IWAVE.EQ.0)THEN
          V1 = VMIN
          VMIN1 = 1E4/(VMAX+0.5*FWHM)
          VMAX1 = 1E4/(V1-0.5*FWHM)
       ELSE
          VMIN1=VMIN-0.5*FWHM
          VMAX1=VMAX+0.5*FWHM
       ENDIF
c    	print*, VMIN1,VMAX1,V1,'test1'
       VBOT = VMIN1 - VREL
       VTOP = VMAX1 + VREL
       IF(VBOT.LT.0.0)VBOT = 0.0
       RANGE = VTOP - VBOT
       NBIN = INT(RANGE/WING) + 1
       print*, VMIN,VMAX,VBOT,VTOP,'v-cork'
       IF(NBIN.GT.MAXBIN)THEN
         WRITE(*,*)'CALC_EXOTABLE3.f :: *ERROR* NBIN > MAXBIN'
         WRITE(*,*)'Stopping program.'
         WRITE(*,*)' '
         WRITE(*,*)'NBIN, MAXBIN = ',NBIN,MAXBIN
         STOP
       ENDIF
       DO 127 I=1,NBIN
        VBIN(I) = VBOT + (I - 1)*WING
127    CONTINUE

        IREC = IREC0

       DO 30 ITEMP=1,NT

        TE1 = TEMP1(ITEMP)

        LINFIL=LINFIL1
        DO I=1,LEN(LINFIL)
         IF(LINFIL(I:I).NE.' ')J=I
        ENDDO

        I1 = INT(ITEMP/10)
        I2 = ITEMP-10*I1

        LINFIL(J+1:J+1)=CHAR(I1+48)
        LINFIL(J+2:J+2)=CHAR(I2+48)

        PRINT*,'ITEMP, TEMP(ITEMP) = ',ITEMP, TE1

        CALL FILE(LINFIL,KEYFIL,'key')
        CALL RDKEY(LUN)
        CALL RDGAS
        CALL RDISO
        OPEN(UNIT=DBLUN,FILE=DBFILE,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='DIRECT',RECL=DBRECL)
c	print*,dbrecl,dbfile,1.5


C Checking isotope included in model and setting mass for doppler width
C calculation

        XMASS=GETMASS(IDGAS1,ISOGAS1)
 
        PRINT*,'XMASS = ',XMASS

 

C       NOTE: TCORS1 includes a factor of 1.e-27 for scaling of stored 
C       lines
        IF(ISOGAS1.EQ.0)THEN
           K = 1                   
        ELSE
           K = ISOGAS1
        ENDIF
        TCORS1 = PARTF(IDGAS1,K,TE1,IPTF)*1.E-27   
        TCORS2 = 1.439*(TE1 - 296.)/(296.*TE1)
        TCORDW = 4.301E-7*SQRT(TE1/XMASS)
        TRATIO = 296./TE1


C       Now cycle through bins
C       Find starting record in line data file
c	print*,dbrec,1.11
        CALL FNDWAV(VBOT)
        NXTREC= DBREC
c	print*,dbrec,1.12
        V1=VBOT
        INORM=1

145     CONTINUE
C       Read in lines until buffer full
        FIRST=1
        RANGE1=VTOP-V1
        print*,'NXTREC,V1,RANGE1,DBREC',NXTREC,V1,RANGE1,DBREC
        CALL LOADLINEBUFFER(NXTREC,V1,RANGE1,MAXLIN,MAXLIN,VLIN,
     1   SLIN,ALIN,ELIN,IDLIN,SBLIN,PSHIFT,DOUBV,TDW,TDWS,LLQ,NLIN,
     2   FIRST,LAST,NGAS,IDGAS,ISOGAS,IFULL,V2)
c        print*,'nlin = ',nlin
c	print*,slin,elin,13
        DO 101 LINE1=1,NLIN
             LINE=FIRST+LINE1-1
             IF(LINE.GT.MAXLIN)LINE=LINE-MAXLIN
             tstim_arr(line) = (1.0 - dpexp(-1.439*vlin(line)/TE1))/
     >          (1.0 - dpexp(-1.439*vlin(line)/296.0))
             lnabsco = log(slin(line))+log(tcors1)+tcors2*elin(line)+
     >         log(tstim_arr(line))
             absco_arr(line)=exp(lnabsco)
             ad_arr(line) = tcordw*vlin(line)
c	print*,10,vlin(line),tstim_arr(line),absco_arr(line)

c   	 print*,SLIN(line),tcors1,tcors2,elin(line),
c     >     tstim_arr(line),9.9
c         print*,LINE,VLIN(LINE),IDLIN(LINE),TE1,P1,
c     >     ABSCO_arr(line),9.9
c         stop

101     CONTINUE

        DO 20 IPRESS=1,NP
          P1 = PRESS1(IPRESS)
       
          DO LINE1=1,NLIN
           LINE=FIRST+LINE1-1
           IF(LINE.GT.MAXLIN)LINE=LINE-MAXLIN
           y_arr(line) = (alin(line)*(1.-frac)*tratio**tdw(line)+
     >         (alin(line) - sblin(line))*frac*tratio**tdws(line))*
     >         p1/ad_arr(line)

          ENDDO
          CALL LBL_KCONT_EXO3(INORM,FIRST,NLIN,VBOT,VMIN1,VMAX1,
     1       WING,IPRESS,P1,TE1,IDGAS1,ISOGAS1,FRAC,IPROC1,MAXDV)
             INORM=0

20      CONTINUE
c------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        IF(NLIN.EQ.MAXLIN)THEN
         V1=V2
         GOTO 145
        ENDIF

       
c        IREC = IREC0

C The use of a non-integer DO loop variable or expression is obsolescent
C in Fortran 90 and deleted in Fortran 95. As DELV is often LT 1.0, a
C scale factor is used below to allow IV to be declared as an integer and
C to avoid a non-integer DO loop.

        DELVSF = 1/DELV

        IBINMAX=0

        IV0=INT((VMINX*DELVSF))

        DO 10 IV=(VMIN*DELVSF),(VMAX*DELVSF),(DELV*DELVSF)

         WRITE(*,*)'CALC_EXOTABLE3.f :: Current Wavenumber = ',IV/DELVSF

         VS1 = (IV/DELVSF) - 0.5*FWHM
         VE1 = VS1 + FWHM

         IF(IWAVE.EQ.0)THEN
               V1=VS1
               VS1=1E4/VE1
               VE1=1E4/V1
         ENDIF

c         print*,VS1,VE1,V1,'test2'

C        Figure out which bins are needed here
         FSTBIN = INT((VS1-VBOT)/WING)
         LSTBIN = INT((VE1-VBOT)/WING) + 2

         V1 = VBOT+(FSTBIN-1)*WING
         V2 = VBOT+LSTBIN*WING

c         print*,FSTBIN,V1,LSTBIN,V2,IBINMAX,VBOT,WING

c         IF(LSTBIN.GT.IBINMAX)THEN
C         Need to read in new lines
          print*,'Reading in new lines'
          VFIRST=V1
          RANGE1=V2-VFIRST 
          CALL FNDWAV(VFIRST)
          NXTREC= DBREC
          FIRST=1 
c          print*,'vlin',vlin
c          print*,'NXTREC,VFIRST,RANGE1',NXTREC,VFIRST,RANGE1     
          CALL LOADLINEBUFFER(NXTREC,VFIRST,RANGE1,MAXLIN,MAXLIN,
     1     VLIN,SLIN,ALIN,ELIN,IDLIN,SBLIN,PSHIFT,DOUBV,TDW,TDWS,
     2     LLQ,NLIN,FIRST,LAST,NGAS,IDGAS,ISOGAS,IFULL,VLAST)
          print*,'NLIN = ',NLIN
          DO IBIN=1,NBIN  
           NLINES(IBIN)=0
           FSTLIN(IBIN)=0
           LSTLIN(IBIN)=-1
          ENDDO

          DO ILIN=1,NLIN
           J = FIRST+ILIN-1
           IF(J.GT.MAXLIN) J=J-MAXLIN
c       print*,j,nlin,vlin(j),elin(j),'test5'
c       if (vlin(j).lt.2000) stop
           tstim_arr(J) = (1.0 - dpexp(-1.439*vlin(J)/te1))/
     >          (1.0 - dpexp(-1.439*vlin(J)/296.0))
           lnabsco = log(slin(J))+log(tcors1)+tcors2*elin(J)
     >        +log(tstim_arr(J))
           absco_arr(J)=exp(lnabsco)
C           absco_arr(J) = slin(J)*tcors1*
C     >          dpexp(tcors2*elin(J))*tstim_arr(J)
           ad_arr(J) = tcordw*vlin(J)

           I= 1 + INT((VLIN(J) - VBIN(1))/WING)
           NLINES(I)= NLINES(I) + 1
c           PRINT*,vlin(j),I
           IF(FSTLIN(I).EQ.0)FSTLIN(I)= J
           LSTLIN(I)= J
c           PRINT*,fstlin(i),lstlin(i)
C          Last bin which is definitely full
           IBINMAX=I-1
          ENDDO

c          PRINT*,'Lines Loaded : '
          DO IBIN=1,NBIN 
            IF(NLINES(IBIN).GT.0)THEN 
c             PRINT*,IBIN,VBIN(IBIN),NLINES(IBIN),FSTLIN(IBIN),
c     1        LSTLIN(IBIN)
            ENDIF
          ENDDO



c        ENDIF


         DO 21 IPRESS=1,NP
          P1 = PRESS1(IPRESS)

        print*,'temp,press',te1,p1
          DO LINE1=1,NLIN
           LINE=FIRST+LINE1-1
           IF(LINE.GT.MAXLIN)LINE=LINE-MAXLIN
           y_arr(line) = (alin(line)*(1.-frac)*tratio**tdw(line)+
     >         (alin(line) - sblin(line))*frac*tratio**tdws(line))*
     >         p1/ad_arr(line)
c     	print*,line,vlin(line),tstim_arr(line),y_arr(line),'test4'
          ENDDO
          print*,'Calling LBL_FKNEW_EXO3'
	print*, 11,VS1,VE1,'test6'
          CALL LBL_FKNEW_EXO3(IWAVE,VS1,VE1,P1,IPRESS,TE1,IDGAS1,
     1       ISOGAS1,IPROC1,XMASS,FRAC,MAXDV,NPOINT)

          DELV=(VE1-VS1)/FLOAT(NPOINT)

          CALL CALC_FKDIST_WAVEC(IWAVE,VS1,DELV,NPOINT,
     1    NFIL,VFIL,TFIL,G_ORD,DEL_G,K_G,NGMAX,NG)

          IREC=IREC0 + (ITEMP-1+(IPRESS-1+(IV-IV0)*NP)*NT)*NG
          DO 40 LOOP=1,NG
              WRITE(LUN0,REC=IREC)K_G(LOOP)
              IREC = IREC + 1
	print*,IREC,itemp,iv,iv0,ipress,loop
40        CONTINUE

21       CONTINUE
10      CONTINUE


        CLOSE(DBLUN)

30     CONTINUE

C-----------------------------------------------------------------------
C
C	Close files and shut down the program.
C
C-----------------------------------------------------------------------

c      CLOSE(LUN0)

      WRITE(*,*)' CALC_EXOTABLE3.f :: calculation complete.'
      CALL GETTIME(TIME)
      TIME2 = TIME
      TOT_TIME = SNGL(TIME2 - TIME1)
      WRITE(*,244)TOT_TIME
244   FORMAT(/' Elapsed time (s) = ',F8.1)

      END
************************************************************************
************************************************************************

