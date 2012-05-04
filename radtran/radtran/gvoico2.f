C***********************************************************************
C
C  PROGRAM        GVOICO2     SUBROUTINE          
C
C     $Id: gvoico2.f,v 1.2 2011-06-17 15:40:26 irwin Exp $
C
C  PURPOSE        COMPUTE CO2 VOIGT LINESHAPE MODIFIED FOR 
C                 SUB-LORENTZIAN LINE WINGS AND LINE MIXING 
C
C  VERSION        3.0   D.P. EDWARDS   27/02/91
C                 (C) COPYRIGHT 1991 UCAR/NCAR
C                 ALL RIGHTS RESERVED
C                 GENLN2 SOFTWARE AND RELATED MATERIALS MAY BE 
C                 USED ONLY UNDER AN EXECUTED VALID LICENCE AGREEMENT
C 
C                 MODIFIED 16/2/93 FOR INCORPORATION INTO GENLBL
C                 P.G.J.IRWIN
C                 NOTE - GENLBL MUST BE MODIFIED TO SEND TEMPERATURE
C                 AND DOPPLER WIDTH TO THE FUNCTION SHAPE
C
C  DESCRIPTION    THIS ROUTINE CALCULATES THE CHI FACTOR MODIFIED
C                 VOIGT FUNCTION OVER A FREQUENCY MESH DUE TO
C                 A LINE CENTRED AT SWN CM-1 USING 
C
C                 TEMPERATURE DEPENDENT CHI FACTOR FOR THE CO2 
C                 4.3 micron BAND IS USED
C                 REF: C.Cousin,R.Le Doucen,C.Boulet and A.Henry, (1985)
C                 Temperature dependence of the absorption in the
C                 region beyond the 4.3 mic band head of CO2.
C                 2: N2 and O2 broadening
C                 Appl.Opt. 24 (22) 3899-3907
C
C                 OR USES LINE MIXING FOR APPROPRIATE LINES
C
C  NOTE:          MUST BE CAREFUL ABOUT INCLUDING LINE MIXING A LONG 
C                 WAY FROM THE MIXED LINES. THIS CAN LEAD TO NUMERICAL 
C                 PROBLEMS. LINE MIXING IS INCLUDED HERE IF THE Y COEF
C                 IS SET AND IF THE LINE IS LESS THAN 20cm-1 (SOMEWHAT
C                 ARBITRARY CUTOFF) FROM THE FREQUENCY POINT OF THE 
C                 CALCULATION. IF IT IS GREATER THAN 20cm-1 THEN THE
C                 LINE WING IS ASSUMED TO BE SUB-LORENTZIAN.
C               
C  ARGUMENTS      NUM    I*4 I/P NUMBER OF FREQUENCY INTERVALS
C                 FF     R*4 I/P FINE WAVENUMBER GRID [cm-1]
C                 SWN    R*4 I/P WAVENUMBER OF LINE CENTRE [cm-1] 
C                 STRPAR R*4 I/P PATH ADJUSTED LINE STRENGTH [cm-1]
C                 WIDPAR R*4 I/P PATH ADJUSTED LORENTZ HALFWIDTH [cm-1] 
C                 DOPWID R*4 I/P DOPPLER LINE WIDTH [cm-1]
C                 YMIX   R*4 I/P LINE MIXING Y COEFFICIENT
C                 T      R*4 I/P PATH TEMPERATURE [K]
C                 XABS   R*4 O/P FREQUENCY MESH POINT ABSORPTIONS
C
C  SUBROUTINES    VOIVEC
C***********************************************************************
C

       REAL FUNCTION GVOICO2(X,Y,T,AD)
       IMPLICIT NONE
       INTEGER MXTCH,MXREG,MXDIV
       PARAMETER (MXTCH=3, MXREG=6, MXDIV=4001)
C  
       INTEGER NREG(MXTCH)
       REAL XABS(1),FF(1)
       REAL TCH(MXTCH), WNL(MXREG,MXTCH), WNH(MXREG,MXTCH),
     +      ACO(MXREG,MXTCH), BCO(MXREG,MXTCH),X,Y,T,AD       
       COMPLEX VT(MXDIV)

       INTEGER NUM,NT,IR,IC1,IC2,I,J,IP
       REAL WIDPAR,SWN,STRPAR,YMIX,LCONS,DOPWID,REPWID,H0
       REAL FDIF,CHI,YM,CHI1,CHI2
C-----------------------------------------------------------------------
       DATA NT/3/ TCH/193.,238.,296./ NREG/6,4,4/
C
       DATA (WNL(IR,1),IR=1,6)/0.,  9., 23., 28., 50., 135./
       DATA (WNL(IR,2),IR=1,4)/0.,  5., 22., 50./  
       DATA (WNL(IR,3),IR=1,4)/0.,  .5, 20., 50./
C
       DATA (WNH(IR,1),IR=1,6)/9., 23., 28., 50., 135., 1.E5/
       DATA (WNH(IR,2),IR=1,4)/5., 22., 50., 1.E5/
       DATA (WNH(IR,3),IR=1,4)/.5, 20., 50., 1.E5/
C
       DATA (ACO(IR,1),IR=1,6)/1.0, 3.908, 0.207, 0.219, 0.146, 1.164/
       DATA (ACO(IR,2),IR=1,4)/1.0, 1.968, 0.160, 0.162/
       DATA (ACO(IR,3),IR=1,4)/1.0, 1.064, 0.125, 0.146/
C
       DATA (BCO(IR,1),IR=1,6)/0.0, 0.1514, 3.778E-3, 0.0276, 
     +   0.0196, 0.0350/
       DATA (BCO(IR,2),IR=1,4)/0.0, 0.1354, 0.0214,   0.0216/
       DATA (BCO(IR,3),IR=1,4)/0.0, 0.1235, 0.0164,   0.0196/ 
C
C  CONSTANTS: SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.56418958
C
C      FUDGE GENLN2 DEFINED VARIABLES FROM GENLBL ONES
C  ---------------*

       NUM=0
       LCONS=0.83255461
       SWN=4000.0
       STRPAR=1.0
       DOPWID=AD*LCONS

       FF(1)=SWN+X*DOPWID/LCONS
       YMIX=0.0
       WIDPAR=Y*DOPWID/LCONS
C ----------------*

       REPWID = LCONS/DOPWID
       H0 = REPWID*0.56418958*STRPAR
       Y = WIDPAR*REPWID
C
C  CALCULATE THE COMPLEX PROBABILITY FUNCTION
C
       CALL VOIVEC(NUM,FF,SWN,REPWID,Y,VT)
C
C  DETERMINE TEMPERATURE REGION FOR CHI-FACTOR INTERPOLATION
C
       IF (T .GE. TCH(NT)) THEN
         IC1 = NT - 1
         IC2 = NT
       ELSEIF (T .LT. TCH(1)) THEN
         IC1 = 1
         IC2 = 2
       ELSE
         DO 5 I=1,NT-1
           IF (T .GE. TCH(I) .AND. T .LT. TCH(I+1)) THEN
             IC1 = I
             IC2 = I+1
             GOTO 6
           ENDIF
    5    CONTINUE
       ENDIF
C
C  LOOP OVER FREQUENCY MESH TO CALCULATE CHI FACTOR
C
    6  DO 10 IP=1,NUM+1
         FDIF = ABS(FF(IP) - SWN)
C
C  LINE MIXING OR CHI FACTOR REGIME
C
         IF (YMIX .NE. 0.0 .AND. FDIF .LT. 10.0) THEN
C
C  COMPUTE ABSORPTION DUE TO VOIGT LINE SHAPE WITH 1st ORDER LINE MIXING
C
           CHI = 1.0
           YM = YMIX
         ELSE
C
C  CALCULATE CHI FACTORS AT ADJACENT TEMP POINTS FOR THIS FDIF
C
           DO 20 J=1,NREG(IC1)
             IF (FDIF .GE. WNL(J,IC1) .AND. FDIF .LT. WNH(J,IC1)) THEN
               IF (IC1 .EQ. 1 .AND. J .EQ. 3) THEN
                 CHI1 = ACO(J,IC1) - BCO(J,IC1)*FDIF
                 GOTO 25
               ELSE
                 CHI1 = ACO(J,IC1)*EXP(-BCO(J,IC1)*FDIF)
                 GOTO 25
               ENDIF
             ENDIF
   20      CONTINUE
C
   25      DO 30 J=1,NREG(IC2)
             IF (FDIF .GE. WNL(J,IC2) .AND. FDIF .LT. WNH(J,IC2)) THEN
               IF (IC2 .EQ. 1 .AND. J .EQ. 3) THEN
                 CHI2 = ACO(J,IC2) - BCO(J,IC2)*FDIF
                 GOTO 35
               ELSE
                 CHI2 = ACO(J,IC2)*EXP(-BCO(J,IC2)*FDIF)
                 GOTO 35
               ENDIF
             ENDIF
   30      CONTINUE
C
   35      CHI = CHI1 + (CHI2-CHI1)*(T-TCH(IC1))/(TCH(IC2)-TCH(IC1))
           IF (CHI .GT. 1.0) CHI = 1.0
           YM = 0.0
         ENDIF
C
C  COMPUTE ABSORPTION DUE TO VOIGT LINE SHAPE
C
         XABS(IP) = H0*CHI*(REAL(VT(IP)) + YM*AIMAG(VT(IP)))
C
   10  CONTINUE       
C
       GVOICO2=XABS(1)/(LCONS/DOPWID)
       RETURN                                              
       END                                      

C***********************************************************************
C
C  PROGRAM        VOIVEC     SUBROUTINE          
C
C  PURPOSE        COMPUTE COMPLEX PROBABILITY FUNCTION
C
C  VERSION        3.0   D.P. EDWARDS   27/02/91
C                 (C) COPYRIGHT 1991 UCAR/NCAR
C                 ALL RIGHTS RESERVED
C                 GENLN2 SOFTWARE AND RELATED MATERIALS MAY BE 
C                 USED ONLY UNDER AN EXECUTED VALID LICENCE AGREEMENT
C
C                 ALTERED TO GET RID OF PARRAY FILE 
C                 P.G.J.IRWIN    16/2/93
C 
C  DESCRIPTION    THIS ROUTINE CALCULATES THE COMPLEX PROBABILITY
C                 FUNCTION USING A VECTORIZED VERSION OF THE
C                 HUMLICEK JQSRT V27 437 1982 PAPER.
C                 THE CALCULATION IS PERFORMED FOR THE ARRAY OF X,Y
C                 PAIRS FOR A GIVEN LINE OVER THE FINE MESH POINTS
C                 OF THE CURRENT WIDE MESH. 
C
C  ARGUMENTS      NUM    I*4 I/P NUMBER OF FREQUENCY INTERVALS
C                 FF     R*4 I/P FINE WAVENUMBER GRID [cm-1]
C                 SWN    R*4 I/P WAVENUMBER OF LINE CENTRE [cm-1] 
C                 REPWID R*4 I/P SQRT(ln2)/(DOPPLER WIDTH [cm-1])
C                 Y      R*4 I/P VOIGT Y PARAMETER OF LINE
C                 VT     COM O/P COMPLEX PROBABILITY FUNCTION
C
C***********************************************************************
C
       SUBROUTINE VOIVEC(NUM,FF,SWN,REPWID,Y,VT)
C
       PARAMETER (MXDIV=4001)

       REAL      FF(NUM+1)
       COMPLEX   U,T,VT(MXDIV)
C-----------------------------------------------------------------------
C
C  SORT THE (X,Y) PAIRS INTO THE 4 REGIONS OF THE HUMLICEK 
C  EXPRESSIONS OF THE VOIGT LINE PROFILE.
C
       DO 20 I=1,NUM+1
         X =  (FF(I) - SWN)*REPWID
         S1V = ABS(X) + Y
         S2V = (0.195*ABS(X)) - 0.176
         T = CMPLX(Y,-X)
C
C  FOR REGION 1 OF HUMLICEK
C
         IF (S1V .GE. 15.0) THEN
           VT(I) = T*0.5641896/(0.5+(T*T))
C
C  REGION 2 OF HUMLICEK
C
         ELSEIF (S1V .GE. 5.5) THEN
           U = T*T
           VT(I) = T*(1.410474 + U*.5641896)/(.75 + U*(3.+U))
C
C  REGION 3 OF HUMLICEK
C
         ELSEIF (Y .GE. S2V) THEN
           VT(I) =
     1     (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+
     2     T*.5642236))))/
     3     (16.4955+T*(38.82363+T*(39.27121+
     4     T*(21.69274+T*(6.699398+T)))))
C
C  REGION 4 OF HUMLICEK
C
         ELSE
           U = T*T
           VT(I)=CEXP(U)-T*(36183.31-U*(3321.9905-
     1     U*(1540.787-U*(219.0313-U*
     2     (35.76683-U*(1.320522-U*.56419))))))/
     3     (32066.6-U*(24322.84-U*
     4     (9022.228-U*(2186.181-U*(364.2191-
     5     U*(61.57037-U*(1.841439-U)))))))
         ENDIF
C
   20  CONTINUE
C
       RETURN                                              
       END                                                            
