      SUBROUTINE EDLINE(BUFFER,QIDENT,ACO2,NCO2,AH2O,NH2O,YACO2,
     1  YNCO2,YAN2,YNN2)
C     $Id: edline.f,v 1.12 2011-06-23 09:09:02 irwin Exp $
C***********************************************************************
C_TITL:	EDLINE.f
C
C_DESC:	Change line parameters as copied by SELECT. Before select.f writes
C	a line data record to the new data base file it calls this routine
C	to allow editing of the record. The action of this routine will
C	vary a lot and multiple versions are expected to exist. edline.f
C	is a generic version which is linked by default. other versions
C	are named edline_xxxxx.f or similar and must be linked explicitly.
C	The file also contains the routine EDSET which prompts the user as
C	necessary for decisions regarding the line editing.
C
C_ARGS:	Input variables:
C	BUFFER		CHARA*(*)	Holds ascii format line data
C					record. CHARACTER*(*) declares an
C					incoming character variable whose
C					length is unknown.
C	QIDENT(190)	CHARA*6
C	ACO2(190)	REAL
C	NCO2(190)	REAL
C	AH2O(190)	REAL
C	NH2O(190)	REAL
C	YACO2(200)	REAL
C	YNCO2(200)	REAL
C	YAN2(200)	REAL
C	YNN2(200)	REAL
C	QGEISA80	CHARA*35	Used to store the upper/lower 
C                                       vibrational state and upper/lower
C                                       rotational state for 80-character
C					format GEISA.
C	QGEISA120	CHARA*36	Used to store the upper/lower
C					vibrational state and upper/lower
C					rotational state for
C					120-character format GEISA.
C
C_CALL:	No calls made.
C
C_HIST:  15oct92  SBC	ORIGINAL VERSION just sets self to air broadened
C			widths.
C        04feb93  PGJI	add option to multiply H2O air width by 1.3
C        10feb93  PGJI	add option to calculate CO2 broadening of H2O from
C                 	delaye data and replace air broadened width and
C			temperature dependence. Also calculate temperature
C			dependence of self broadening and place in
C			probability column (which we will assume is zero.
C			Convert BUFFER to numbers at beginning and
C			convert back to ASCII at end to simplify
C	5/1998	CAN	Removed "Hack to cure line strength bug with CH3D
C			in GEISA" ...
C			        IF(TNID.EQ.23)THEN
C			          TNSTR= TNSTR*RELABU(3,6)
C			        ENDIF
C			becuase happens AGAIN later in rdline.f.
C	28feb00	CAN	altered read/write format statements to handle
C			the problem of the line strength (TNSTR) sometines
C			being smaller than smallest machine-storable
C			number. The method used here is to read TNSTR in
C			and out as a 10-character string, so if you need
C			to numerically manipulate it you will need to work
C			around this, remembering that the resulting number
C			may not be machine storable. This can also cause
C			problems when you write back to BUFFER.
C			Note 1: this method is different from the method
C			  in rdline, which reads in mantissa and exponent
C                         separately, and then scales by 1e+47 to get 
C                         a usable number. If you want to do that here,
C                         you will need to re-scale before output.
C			Note 2: the number TNPROB has not caused problems
C                         with smallness so far, so here is read in
C                         in the regular manner. Again this is different
C                         from the method in RDLINE which handles mant./exp.
C	20jan04	NT	Added option to broaden PH3 for Saturn
C				h2 and he abundances.
C	12jan05	NT	added output format options to allow TNWIDS
C			to be >1 (but still less than 10). this is necessary
C			for the Hitran04 HCN data which has widths of
C			over 1.0 in some cases.
C***************************** VARIABLES *******************************

C ../includes/dbcom.f stores the line database variables .
      INCLUDE '../includes/dbcom.f' 

      INTEGER I,J
      INTEGER TNID,TNISO,TNUGQI,TNLGQI,TNACC(3),TNREF(3)
      INTEGER TNACC04(6),TNREF04(6)

      REAL ACO2(190),NCO2(190),AH2O(190),NH2O(190)
      REAL YACO2(200),YNCO2(200),YAN2(200),YNN2(200)
      REAL WIDC,WIDC1,EXPC,EXPC1,EXPH,EXPH1,LNCONS
      REAL WIDH,WIDH1,WIDCT0,WIDCT1,WIDHT0,WIDHT1
      REAL TNPROB,TNWIDA,TNWIDS,TNLSE,TNTDEP,TNSHIF,TNEINA
      REAL GAMMA_H2,GAMMA_HE,TNUWGHT,TNLWGHT
      DOUBLE PRECISION TNWAVE
      CHARACTER*(*) BUFFER
      CHARACTER*6 QIDENT(190)
      CHARACTER*9 TNULQ,TNLLQ
      CHARACTER*15 TNULQ04,TNLLQ04,TNLGQ04,TNUGQ04
      CHARACTER*10 TNSTR
      CHARACTER*35 QGEISA80
      CHARACTER*36 QGEISA120
      CHARACTER*1 TNFLAG


      LOGICAL WSBTAB,CO2TAB,H2OTAB,SCOFLD,WEDAD,BEZARD,BERGH,PMIRR
      LOGICAL H2HePH3J,H2HePH3S,BERGC,H2HeCH4
      COMMON /EDLOG/WSBTAB,CO2TAB,H2OTAB,SCOFLD,WEDAD,BEZARD,BERGH,
     $   BERGC,PMIRR,H2HePH3J,H2HePH3S,H2HeCH4

C******************************** CODE *********************************

1     FORMAT(A)
      IF(DBFORM.EQ.0)THEN
C=======================================================================
C
C	HITRAN (either of 100-, 112-, or 160-character formats)
C
C=======================================================================
        IF(DBRECL.EQ.100)THEN
          READ(BUFFER,100)TNID,TNISO,TNWAVE,TNSTR,TNPROB,TNWIDA,TNWIDS,
     1    TNLSE,TNTDEP,TNUGQI,TNLGQI,TNULQ,TNLLQ,TNACC,TNREF
100       FORMAT(I2,I1,F12.6,A10,E10.3,F5.4,F5.4,F10.4,F4.2,8X,I3,I3,
     1    A9,A9,3I1,3I2)
          TDOUBV= 0.0
        ELSEIF(DBRECL.EQ.112)THEN
          READ(BUFFER,105)TNID,TNISO,TNWAVE,TNSTR,TNPROB,TNWIDA,TNWIDS,
     1    TNLSE,TNTDEP,TNUGQI,TNLGQI,TNULQ,TNLLQ,TNACC,TNREF,TDOUBV
105       FORMAT(I2,I1,F12.6,A10,E10.3,F5.4,F5.4,F10.4,F4.2,8X,I3,I3,
     1    A9,A9,3I1,3I2,F12.7)
        ELSEIF(DBRECL.EQ.52)THEN
           READ(BUFFER,207)TNID,TNISO,TNWAVE,TNSTR,
     1     TNLSE,TNWIDA,TNTDEP,TNWIDS,TNTDEPS
207        FORMAT(I2,I1,F12.6,A10,F10.4,F5.4,F3.2,F6.4,F3.2)
        ELSEIF(DBRECL.EQ.160)THEN
          READ(BUFFER,107)TNID,TNISO,TNWAVE,TNSTR,TNEINA,TNWIDA,TNWIDS,
     1     TNLSE,TNTDEP,TNSHIF,TNUGQ04,TNLGQ04,TNULQ04,TNLLQ04,TNACC04,
     2     TNREF04,TNFLAG,TNUWGHT,TNLWGHT
107       FORMAT(I2,I1,F12.6,A10,E10.3,F5.4,F5.4,F10.4,F4.2,F8.6,
     1     A15,A15,A15,A15,6I1,6I2,A1,F7.1,F7.1)        
        ELSE
          WRITE(*,*)' EDLINE.f :: HITRAN format not recognised.'
          WRITE(*,*)' Stopping program.'
          WRITE(*,*)' '
          WRITE(*,*)' Accepted DBRECL for DBFORM = 0:100,112,160.'
          WRITE(*,*)' DBFORM, DBRECL = ',DBFORM,DBRECL
          STOP
        ENDIF 
      ELSEIF(DBFORM.EQ.1)THEN
C=======================================================================
C
C	GEISA (either of 80-, 120- or 211-character formats)
C
C=======================================================================
        IF(DBRECL.EQ.82)THEN
          READ(BUFFER,200)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA120,
     1    TNTDEP,TNISO,TNID
200       FORMAT(F10.3,A10,F5.3,F10.3,A36,F4.2,I4,I3)
C Note: the above assumes that GEISA version with temp dependence is used.
C Also ignores quantum numbers in this format.
          TNPROB= 0.0
          TNWIDS= 0.0
          TNPSH= 0.0
          TNUGQI= 0
          TNLGQI= 0
          TNULQ= ' '
          TNLLQ= ' '
          TNACC(1)= 0
          TNACC(2)= 0
          TNACC(3)= 0
          TNREF(1)= 0
          TNREF(2)= 0
          TNREF(3)= 0
        ELSEIF(DBRECL.EQ.80)THEN
          READ(BUFFER,210)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA80,
     1    TNTDEP,TNISO,TNID
210       FORMAT(F10.3,A10,F5.3,F10.3,A35,A3,I4,I3)
C Note: the above assumes that GEISA version with temp dependence is used.
C Also ignores quantum numbers in this format.
          TNPROB= 0.0
          TNWIDS= 0.0
          TNPSH= 0.0
          TNUGQI= 0
          TNLGQI= 0
          TNULQ= ' '
          TNLLQ= ' '
          TNACC(1)= 0
          TNACC(2)= 0
          TNACC(3)= 0
          TNREF(1)= 0
          TNREF(2)= 0
          TNREF(3)= 0
        ELSEIF(DBRECL.EQ.120)THEN
          READ(BUFFER,205)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA120,
     1    TNTDEP,TNISO,TNID,TNPROB,TNWIDS,TNPSH,TNACC,TNREF
205       FORMAT(F10.3,A10,F5.3,F10.3,A36,F4.2,I4,I3,6X,
     1    E10.3,F5.4,F8.6,3I1,3I2)
          TNUGQI= 0
          TNLGQI= 0
          TNULQ= ' '
          TNLLQ= ' '
        ELSEIF(DBRECL.EQ.211)THEN
          READ(BUFFER,206)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA120,
     1    TNTDEP,TNISO,TNID,TNPROB,TNWIDS,TNPSH,TNACC,TNREF
206       FORMAT(F12.6,1X,A10,F6.4,F10.4,A36,F4.2,I3,I3,6X,
     1    E10.3,F5.4,F8.6,3I1,3I2)
          TNUGQI= 0
          TNLGQI= 0
          TNULQ= ' '
          TNLLQ= ' '
        ELSE
         WRITE(*,*)' EDLINE.f :: GEISA format not recognised.'
         WRITE(*,*)' Stopping program.'
         WRITE(*,*)' '
         WRITE(*,*)' Accepted DBRECL for DBFORM = 1, either 80,120.'
         WRITE(*,*)' DBFORM, DBRECL = ',DBFORM,DBRECL
         STOP
        ENDIF
      ENDIF

      TNSHIF= 0.0
C-----------------------------------------------------------------------
C
C     SET NEW DEFAULT THAT TNPROB = TEMPERATURE DEPENDENCE OF SB WIDTH
C     
C     IF NOT SET THEN = TNTDEP
C     (took next line out: why was it ever in? CAN 5/1998)
cc     TNPROB=TNTDEP
C     -----------------------

      LNCONS=0.6797241
C     LOG(296/150)
      
C=======================================================================
C
C	WSBTAB=TRUE means fudge self broadened widths to be same as air 
C	broadened. Useful to compare with old calculations.
C
C=======================================================================
      IF(WSBTAB)THEN
        IF(DBFORM.EQ.0)THEN
C HITRAN format
          TNWIDS= TNWIDA
        ELSEIF(DBFORM.EQ.1)THEN
C GEISA format; GEISA does not include self broadened widths
        ELSE
          WRITE(*,*)' EDLINE.f :: invalid line format.'
          WRITE(*,*)' EDLINE.f :: Stopping program.'
          STOP
        ENDIF
      ENDIF
      
C=======================================================================
C
C	CO2TAB=TRUE means multiply air broadened widths of H2O by a
C	factor of 1.3 to simulate CO2 broadening
C
C=======================================================================
      IF(CO2TAB)THEN
        IF(DBFORM.EQ.0)THEN
C Only do for HITRAN format for time being
          IF(TNID.EQ.1) THEN
            TNWIDA= TNWIDA*1.3
          ENDIF
        ENDIF
      ENDIF

C=======================================================================
C
C	H2OTAB=TRUE means replace air broadened width and temperature
C	dependence of water with Delaye calculated value. Also replace
C	transition probability column by temperature dependence of self
C	broadening calculated from delaye. If self broadening is set to
C	zero, put in Delaye calculated value.
C
C=======================================================================
      IF(H2OTAB)THEN
        IF(DBFORM.EQ.0)THEN
C Only do for HITRAN format for time being
          IF(TNID.EQ.1)THEN
C Check that we have water vapour ...
            WIDC= 0.0
            WIDH= 0.0
            EXPC= 0.0
            EXPH= 0.0
            WIDC1= WIDC
            WIDH1= WIDH
            EXPC1= EXPC
            EXPH1= EXPH
            IF(TNULQ(1:2).EQ.TNLLQ(1:2))THEN
C Q branch lines *******************************************************
C JL=JU. Use TNLLQ to find widths.
              DO 10 I=1,187
                IF(TNLLQ.EQ.QIDENT(I))THEN
                  WIDC= ACO2(I)
                  WIDH= AH2O(I)
                  EXPC= NCO2(I)
                  EXPH= NH2O(I)
                  WIDC= 0.001*WIDC*(300.0/296.0)**EXPC
                  WIDH= 0.001*WIDH*(300.0/296.0)**EXPH
                ENDIF
 10           CONTINUE
            ELSE
C P,R branch lines *****************************************************
              DO 20 I=1,187
                IF(TNULQ.EQ.QIDENT(I))THEN
                  WIDC= ACO2(I)
                  WIDH= AH2O(I)
                  EXPC= NCO2(I)
                  EXPH= NH2O(I)
                  WIDC= WIDC*(300.0/296.0)**EXPC
                  WIDH= WIDH*(300.0/296.0)**EXPH
                ENDIF
                IF(TNLLQ.EQ.QIDENT(I))THEN
                  WIDC1= ACO2(I)
                  WIDH1= AH2O(I)
                  EXPC1= NCO2(I)
                  EXPH1= NH2O(I)
                  WIDC1= WIDC1*(300.0/296.0)**EXPC1
                  WIDH1= WIDH1*(300.0/296.0)**EXPH1
                ENDIF
 20           CONTINUE
C Calculation of P,R lines from Delaye: For CO2, n=4. For H2O, n=3.
              WIDCT0= ((WIDC**3+WIDC1**3)*0.5)**0.33333333
              WIDHT0= ((WIDH**2+WIDH1**2)*0.5)**0.5
C WIDT0 = Width at 296K
              WIDC= WIDC*(296.0/150.0)**EXPC
              WIDC1= WIDC1*(296.0/150.0)**EXPC1
              WIDH= WIDH*(296.0/150.0)**EXPH
              WIDH1= WIDH1*(296.0/150.0)**EXPH1
                  
              WIDCT1= ((WIDC**3+WIDC1**3)*0.5)**0.33333333
              WIDHT1= ((WIDH**2+WIDH1**2)*0.5)**0.5
C WIDT1 = Width at 150K
                  
                  
C Calculate temperature exponent n from estimates of width at 296,150K
              IF(WIDCT1.NE.0.AND.WIDCT0.NE.0)THEN
                EXPC= (LOG(WIDCT1) - LOG(WIDCT0))/LNCONS
              ELSE
                WIDCT0= 0
              ENDIF
              IF(WIDHT1.NE.0.AND.WIDHT0.NE.0)THEN
                EXPH= (LOG(WIDHT1) - LOG(WIDHT0))/LNCONS
              ENDIF
              WIDC= 0.001*WIDCT0
              WIDH= 0.001*WIDHT0
            ENDIF

C Put results back into database line. Put T dependence of sb into
C transition probability column. If no match has been found, do nothing
            IF(WIDC.NE.0)THEN
              TNWIDA= WIDC
              TNTDEP= EXPC
              TNPROB= EXPH
              IF(TNWIDS.LT.0.0001)THEN
                TNWIDS= WIDH
              ENDIF
            ELSE
              WRITE(*,*)' EDLINE.f :: Delaye: No Match. No correction'
              WRITE(*,*)' EDLINE.f :: made. TNULQ,TNLLQ: ',TNULQ,TNLLQ
            ENDIF
          ENDIF
        ENDIF
      ENDIF

C=======================================================================
C
C	SCOFLD=TRUE means ?
C
C=======================================================================
      IF(SCOFLD)THEN
        IF(TNID.EQ.2)THEN
          IF(TNLLQ(7:7).EQ.' ')THEN
            I= ICHAR(TNLLQ(8:8))-48
          ELSE
            I= ICHAR(TNLLQ(8:8))-48 + 10*(ICHAR(TNLLQ(7:7))-48)
          ENDIF
          IF(TNLLQ(5:5).EQ.'R')THEN
            I= 2*I + 2
          ELSEIF(TNLLQ(5:5).EQ.'P')THEN
            I= 2*I
          ELSEIF(TNLLQ(5:5).EQ.'Q')THEN
            I= 2*I + 1
          ENDIF
          TNWIDS= YACO2(I)*1.1/1.013246
          TNWIDA= YAN2(I)*1.1/1.013246
          TNTDEP= YNN2(I)
          TNPROB= YNCO2(I)
        ELSEIF(TNID.EQ.1)THEN
          TNTDEP= 0.5
          TNWIDA= TNWIDA/1.013246
          TNWIDS= TNWIDA*6.5
          TNPROB= TNTDEP    
        ENDIF
      ENDIF

C=======================================================================
C
C	WEDAD=TRUE means ?
C
C=======================================================================
      IF(WEDAD)THEN
        IF(TNID.EQ.2)THEN
          IF(TNLLQ(7:7).EQ.' ')THEN
            I= ICHAR(TNLLQ(8:8))-48
          ELSE
            I= ICHAR(TNLLQ(8:8))-48 + 10*(ICHAR(TNLLQ(7:7))-48)
          ENDIF
          IF(TNLLQ(5:5).EQ.'R')THEN
            I= 2*I + 2
          ELSEIF(TNLLQ(5:5).EQ.'P')THEN
            I= 2*I
          ELSEIF(TNLLQ(5:5).EQ.'Q')THEN
            I= 2*I + 1
          ENDIF
          IF(TNWIDS.EQ.0.0)TNWIDS= YACO2(I)/1.013246
          IF(TNTDEP.EQ.0.0)TNTDEP= YNN2(I)
          TNPROB= YNCO2(I)
        ELSE IF(TNID.EQ.1)THEN
          TNTDEP= 0.5
          TNWIDA= TNWIDA/1.013246
          TNWIDS= TNWIDA*6.5
          TNPROB= TNTDEP    
        ENDIF
      ENDIF

C=======================================================================
C
C	BEZARD=TRUE means to set TNWIDA and TNTDEP variables of certain
C	species to values that were mutually agreed upon for the purposes
C	of RT model inter-comparisons.
C
C	TNWIDA (also known as "LNWIDA") is the air-broadened halfwidth
C	HWHM [cm-1/atm] at 296K.
C	TNTDEP (also known as "LNTDEP") is the coefficient of
C	temperature dependence of HWHM.
C=======================================================================
      IF(BEZARD)THEN
         IF(LOCID(TNID).EQ.6)THEN
C Modify CH4 foreign broadening parameters
           TNWIDA= 0.061
           TNTDEP= 0.45
         ELSE IF(LOCID(TNID).EQ.11)THEN
C Modify NH3 foreign broadening parameters
          TNWIDA= 0.072
          TNTDEP= 0.73
         ELSE IF(LOCID(TNID).EQ.28)THEN
C Modify PH3 foreign broadening parameters
           TNWIDA= 0.093
           TNTDEP= 0.68
         ELSE IF(LOCID(TNID).EQ.26)THEN
C Modify C2H2 foreign broadening
           TNWIDA= 0.085
           TNTDEP= 0.75
         ELSE IF(LOCID(TNID).EQ.27)THEN
C Modify C2H6 foreign broadening
           TNWIDA= 0.105
           TNTDEP= 0.94
         ELSE IF(LOCID(TNID).EQ.32)THEN
C Modify C2H4 foreign broadening
           TNWIDA= 0.101
           TNTDEP= 0.63
         ENDIF
      ENDIF





C=======================================================================
C
C	BERGH=TRUE means to set linewidths and temperature exponents as 
C       suggested for hydrogen-broadening of methane at near-IR wavelengths
C
C	TNWIDA (also known as "LNWIDA") is the air-broadened halfwidth
C	HWHM [cm-1/atm] at 296K.
C	TNTDEP (also known as "LNTDEP") is the coefficient of
C	temperature dependence of HWHM.
C=======================================================================
      IF(BERGH)THEN
         IF(LOCID(TNID).EQ.6)THEN
C Modify CH4 foreign broadening parameters
           TNWIDA= 0.06
           TNTDEP= 0.44
           IF(TNISO.EQ.212)THEN
            print*,'CH3D'
            TNWIDA=0.07
            TNTDEP= 0.6
           ENDIF
         ENDIF
      ENDIF


C=======================================================================
C
C	H2HeCH4=TRUE means to set linewidths and temperature exponents as 
C       suggested for hydrogen-broadening of methane from ExoMOL
C
C	TNWIDA (also known as "LNWIDA") is the air-broadened halfwidth
C	HWHM [cm-1/atm] at 296K.
C	TNTDEP (also known as "LNTDEP") is the coefficient of
C	temperature dependence of HWHM.
C=======================================================================
      IF(H2HeCH4)THEN
         IF(LOCID(TNID).EQ.6)THEN
C Modify CH4 foreign broadening parameters. Assumes H2/He = 0.845/0.155
           TNWIDA = 0.062*0.845 + 0.035*0.155
           TNTDEP = 0.5*0.845 + 0.3*0.155
         ENDIF
      ENDIF
      

C=======================================================================
C
C	BERGC=TRUE means to set linewidths and temperature exponents as 
C       suggested for hydrogen-broadening of methane at near-IR wavelengths
C
C       Also modifies the CH3D fraction from 0.005 assumed in the 
C       Campargue line data to
C
C	TNWIDA (also known as "LNWIDA") is the air-broadened halfwidth
C	HWHM [cm-1/atm] at 296K.
C	TNTDEP (also known as "LNTDEP") is the coefficient of
C	temperature dependence of HWHM.
C=======================================================================
      IF(BERGC)THEN
         IF(LOCID(TNID).EQ.6)THEN
C Modify CH4 foreign broadening parameters
           TNWIDA= 0.06
           TNTDEP= 0.44
           IF(TNISO.EQ.212)THEN
            print*,'CH3D'
            TNWIDA=0.07
            TNTDEP= 0.6
C           Need to modify CH3D linestength here...
           ENDIF
         ENDIF
      ENDIF
      

C=======================================================================
C
C	PMIRR=TRUE means that for water vapour: replace air broadened
C	width and temperature dependence of water with Delaye calculated
C	CO2 value. Also replace transition probability column by
C	temperature dependence of self broadening calculated from
C	delaye. If self broadening is set to zero, put in Delaye
C	calculated value.
C     
C	PMIRR=TRUE means that for CO2: leave air broadened CO2
C	unchanged. If temperature dependence of sb = 0, put gamma_air =
C	gamma_N2 from Yamamoto. If self broadened width = 0, substitute
C	Yamamoto value. Alter transition probability column to hold
C	Yamamoto temperature dependence of self broadening. 
C
C=======================================================================
      IF(PMIRR) THEN
        IF(DBFORM.EQ.0)THEN
C Only do for HITRAN format for time being
          IF(TNID.EQ.2)THEN
            IF(TNLLQ(7:7).EQ.' ')THEN
              I= ICHAR(TNLLQ(8:8))-48
            ELSE
              I= ICHAR(TNLLQ(8:8))-48 + 10*(ICHAR(TNLLQ(7:7))-48)
            ENDIF
            IF(TNLLQ(5:5).EQ.'R')THEN
              I= 2*I + 2
            ELSE IF(TNLLQ(5:5).EQ.'P')THEN
              I= 2*I
            ELSE IF(TNLLQ(5:5).EQ.'Q')THEN
              I= 2*I + 1
            ENDIF
            IF(TNWIDS.EQ.0.0)TNWIDS= YACO2(I)/1.013246
            IF(TNTDEP.EQ.0.0)TNTDEP= YNN2(I)
            TNPROB= YNCO2(I)
          ELSE IF(TNID.EQ.1)THEN
C Check that we have water vapour ...
            WIDC= 0.0
            WIDH= 0.0
            EXPC= 0.0
            EXPH= 0.0
            WIDC1= WIDC
            WIDH1= WIDH
            EXPC1= EXPC
            EXPH1= EXPH
            IF(TNULQ(1:2).EQ.TNLLQ(1:2))THEN
C Q branch lines *******************************************************
C JL=JU. Use TNLLQ to find widths.
              DO 11 I=1,187
                IF(TNLLQ.EQ.QIDENT(I))THEN
                  WIDC= ACO2(I)
                  WIDH= AH2O(I)
                  EXPC= NCO2(I)
                  EXPH= NH2O(I)
                  WIDC= 0.001*WIDC*(300.0/296.0)**EXPC
                  WIDH= 0.001*WIDH*(300.0/296.0)**EXPH
                ENDIF
 11           CONTINUE
            ELSE
C P,R branch lines *****************************************************
              DO 21 I=1,187
                IF(TNULQ.EQ.QIDENT(I))THEN
                  WIDC= ACO2(I)
                  WIDH= AH2O(I)
                  EXPC= NCO2(I)
                  EXPH= NH2O(I)
                  WIDC= WIDC*(300.0/296.0)**EXPC
                  WIDH= WIDH*(300.0/296.0)**EXPH
                ENDIF
                IF(TNLLQ.EQ.QIDENT(I))THEN
                  WIDC1= ACO2(I)
                  WIDH1= AH2O(I)
                  EXPC1= NCO2(I)
                  EXPH1= NH2O(I)
                  WIDC1= WIDC1*(300.0/296.0)**EXPC1
                  WIDH1= WIDH1*(300.0/296.0)**EXPH1
                ENDIF
 21           CONTINUE
C Calculation of P,R lines from Delaye: For CO2, n=4. For H2O, n=3,
              WIDCT0= ((WIDC**3 + WIDC1**3)*0.5)**0.33333333
              WIDHT0= ((WIDH**2 + WIDH1**2)*0.5)**0.5
C WIDT0 = Width at 296K
              WIDC= WIDC*(296.0/150.0)**EXPC
              WIDC1= WIDC1*(296.0/150.0)**EXPC1
              WIDH= WIDH*(296.0/150.0)**EXPH
              WIDH1= WIDH1*(296.0/150.0)**EXPH1
              WIDCT1= ((WIDC**3 + WIDC1**3)*0.5)**0.33333333
              WIDHT1= ((WIDH**2 + WIDH1**2)*0.5)**0.5
C WIDT1 = Width at 150K


C Calculate temperature exponent n from estimates of width at 296,150K
              IF(WIDCT1.NE.0.AND.WIDCT0.NE.0)THEN
                EXPC= (LOG(WIDCT1)-LOG(WIDCT0))/LNCONS
              ELSE
                WIDCT0= 0
              ENDIF
              IF(WIDHT1.NE.0.AND.WIDHT0.NE.0)THEN
                EXPH= (LOG(WIDHT1)-LOG(WIDHT0))/LNCONS
              ENDIF
              WIDC= 0.001*WIDCT0
              WIDH= 0.001*WIDHT0
            ENDIF
C***********************************************************************
C Put results back into database line. Put T dependence of sb into 
C transition probability column. If no match has been found, do nothing.
            IF(WIDC.NE.0)THEN
              TNWIDA= WIDC
              TNTDEP= EXPC
              TNPROB= EXPH
              IF(TNWIDS.LT.0.0001)THEN
                TNWIDS= WIDH
              ENDIF
            ELSE
              WRITE(*,*)' EDLINE.f :: Delaye: No Match. No correction'
              WRITE(*,*)' made. TNULQ,TNLLQ: ',TNULQ,TNLLQ
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      
C=======================================================================
C
C	H2HePH3=TRUE means to H2-HE broaden PH3 linedata according to the
C	formulation in Levy etal 1993 (Journal of Molecular Spectroscopy).
C	H2HePH3J=TRUE => used h2/he abundance for jupiter
C	H2HePH3S=TRUE => used h2/he abundance for saturn
C
C	The formulation is based on the lower-state
C	rotational/local-quanta index information (in other words, the
C	J-value) stored in the respective linedatabases. Jovian abundances
C	of H2 (0.863) and He (0.134) are currently used in the
C	calculation, but can be changed in the future.
C
C	The coefficient of temperature-dependence of TNWIDA ("TNTDEP")
C	was left as it was in Bruno Bezard's assumptions above.
C
C=======================================================================
c  ** jupiter case **
      IF(H2HePH3J)THEN
        READ(QGEISA120(24:27),*)J
	  IF(LOCID(TNID).EQ.28)THEN
          GAMMA_H2= 0.863*(0.1078 - (0.0014*J))
          GAMMA_HE= 0.134*(0.0618 - (0.0012*J))
          TNWIDA= GAMMA_H2 + GAMMA_HE
C PH3 coefficient of HWHM temperature dependance ...
	    TNTDEP= 0.68
        ELSE IF(LOCID(TNID).NE.28)THEN
           WRITE(*,*)' EDLINE.f :: Error: Can only Levy, A etal (1993)'
           WRITE(*,*)' H2-He broaden PH3. Stopping program.'
           WRITE(*,*)' '
           WRITE(*,*)' LOCID(TNID) = ',LOCID(TNID)
           STOP
        ENDIF
      ENDIF
c  ** saturn case **
      IF(H2HePH3S)THEN
        READ(QGEISA120(24:27),*)J
	  IF(LOCID(TNID).EQ.28)THEN
c	**  conrath and gautier 2000 values for vmr H2 and He **
          GAMMA_H2= 0.881*(0.1078 - (0.0014*J))
          GAMMA_HE= 0.119*(0.0618 - (0.0012*J))
          TNWIDA= GAMMA_H2 + GAMMA_HE
C PH3 coefficient of HWHM temperature dependance ...
	    TNTDEP= 0.68
        ELSE IF(LOCID(TNID).NE.28)THEN
           WRITE(*,*)' EDLINE.f :: Error: Can only Levy, A etal (1993)'
           WRITE(*,*)' H2-He broaden PH3. Stopping program.'
           WRITE(*,*)' '
           WRITE(*,*)' LOCID(TNID) = ',LOCID(TNID)
           STOP
        ENDIF
      ENDIF

C=======================================================================
C
C	Write modified data back to new database
C
C=======================================================================

      IF(DBFORM.EQ.0)THEN
C=======================================================================
C
C	HITRAN (either of 100-, 112- or 160-character formats)
C
C=======================================================================
        IF(DBRECL.EQ.100)THEN
	   IF (TNWIDS.LT.1.0) THEN
          WRITE(BUFFER,106)TNID,TNISO,TNWAVE,TNSTR,TNPROB,TNWIDA,TNWIDS,
     1    TNLSE,TNTDEP,TNSHIF,TNUGQI,TNLGQI,TNULQ,TNLLQ,TNACC,TNREF
         ELSE
         WRITE(BUFFER,1061)TNID,TNISO,TNWAVE,TNSTR,TNPROB,TNWIDA,TNWIDS,
     1    TNLSE,TNTDEP,TNSHIF,TNUGQI,TNLGQI,TNULQ,TNLLQ,TNACC,TNREF
         ENDIF
106       FORMAT(I2,I1,F12.6,A10,E10.3,F5.4,F5.4,F10.4,F4.2,F8.5,I3,I3,
     1    A9,A9,3I1,3I2)
1061      FORMAT(I2,I1,F12.6,A10,E10.3,F5.4,F5.3,F10.4,F4.2,F8.5,I3,I3,
     1    A9,A9,3I1,3I2)
        ELSE IF(DBRECL.EQ.52)THEN
           WRITE(BUFFER,207)TNID,TNISO,TNWAVE,TNSTR,
     1     TNLSE,TNWIDA,TNTDEP,TNWIDS,TNTDEPS
        ELSE IF(DBRECL.EQ.112)THEN
	   IF (TNWIDS.LT.1.0) THEN
          WRITE(BUFFER,111)TNID,TNISO,TNWAVE,TNSTR,TNPROB,TNWIDA,TNWIDS,
     1    TNLSE,TNTDEP,TNSHIF,TNUGQI,TNLGQI,TNULQ,TNLLQ,TNACC,TNREF,
     2    TDOUBV
         ELSE
         WRITE(BUFFER,1111)TNID,TNISO,TNWAVE,TNSTR,TNPROB,TNWIDA,TNWIDS,
     1    TNLSE,TNTDEP,TNSHIF,TNUGQI,TNLGQI,TNULQ,TNLLQ,TNACC,TNREF,
     2    TDOUBV
         ENDIF
111       FORMAT(I2,I1,F12.6,A10,E10.3,F5.4,F5.4,F10.4,F4.2,F8.5,I3,I3,
     1    A9,A9,3I1,3I2,F12.7)
1111      FORMAT(I2,I1,F12.6,A10,E10.3,F5.4,F5.3,F10.4,F4.2,F8.5,I3,I3,
     1    A9,A9,3I1,3I2,F12.7)
        ELSE
	  IF (TNWIDS.LT.1.0) THEN
          WRITE(BUFFER,108)TNID,TNISO,TNWAVE,TNSTR,TNEINA,TNWIDA,TNWIDS,
     1    TNLSE,TNTDEP,TNSHIF,TNUGQ04,TNLGQ04,TNULQ04,TNLLQ04,
     2    TNACC04,TNREF04,TNFLAG,TNUWGHT,TNLWGHT
          ELSE
         WRITE(BUFFER,1081)TNID,TNISO,TNWAVE,TNSTR,TNEINA,TNWIDA,TNWIDS,
     1    TNLSE,TNTDEP,TNSHIF,TNUGQI04,TNLGQI04,TNULQ04,TNLLQ04,
     2    TNACC04,TNREF04,TNFLAG,TNUWGHT,TNLWGHT
          ENDIF
108       FORMAT(I2,I1,F12.6,A10,E10.3,F5.4,F5.4,F10.4,F4.2,F8.6,A15,   
     1    A15,A15,A15,6I1,6I2,A1,F7.1,F7.1)
1081      FORMAT(I2,I1,F12.6,A10,E10.3,F5.4,F5.3,F10.4,F4.2,F8.6,A15,   
     1    A15,A15,A15,6I1,6I2,A1,F7.1,F7.1)
        ENDIF
      ELSE
C=======================================================================
C
C	GEISA (either of 80-, 82- or 120-character formats)
C
C=======================================================================
        IF(DBRECL.EQ.80)THEN
          WRITE(BUFFER,210)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA80,
     1    TNTDEP,TNISO,TNID
        ELSEIF(DBRECL.EQ.82)THEN
          WRITE(BUFFER,200)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA120,
     1    TNTDEP,TNISO,TNID
        ELSEIF(DBRECL.EQ.120)THEN
	   IF (TNWIDS.LT.1.0) THEN
          WRITE(BUFFER,205)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA120,
     1    TNTDEP,TNISO,TNID,TNPROB,TNWIDS,TNPSH,TNACC,TNREF
         ELSE
          WRITE(BUFFER,2051)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA120,
     1    TNTDEP,TNISO,TNID,TNPROB,TNWIDS,TNPSH,TNACC,TNREF
         ENDIF
        ELSE
	   IF (TNWIDS.LT.1.0) THEN
          WRITE(BUFFER,206)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA120,
     1    TNTDEP,TNISO,TNID,TNPROB,TNWIDS,TNPSH,TNACC,TNREF
         ELSE
          WRITE(BUFFER,2061)TNWAVE,TNSTR,TNWIDA,TNLSE,QGEISA120,
     1    TNTDEP,TNISO,TNID,TNPROB,TNWIDS,TNPSH,TNACC,TNREF
         ENDIF
        ENDIF
      ENDIF

2051  FORMAT(F10.3,A10,F5.3,F10.3,A36,F4.2,I4,I3,6X,
     1    E10.3,F5.3,F8.6,3I1,3I2)
2061  FORMAT(F12.6,1X,A10,F6.4,F10.4,A36,F4.2,I3,I3,6X,
     1    E10.3,F5.3,F8.6,3I1,3I2)

      RETURN

      END
