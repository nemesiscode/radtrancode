      BLOCK DATA INPUT                                                  
c
c	Clean F77 version of an old CYBER program
c	Runs on VAX ( tested 17-NOV-1988) OK
C                                                                       
c	==================================================
c	If you use this program, please cite original work:
c     A. Borysow and L. Frommhold,
c     "Collision induced rototranslational absorption spectra of
c     CH4-CH4 pairs at temperatures from 50 to 300K",
c	Astrophysical Journal, vol.318, p.940-943, (1987). 
c	==================================================

C     THESE DATA SERVE AS THE INPUT TO THE MAIN PROGRAM ADDEM.          
C     TEMP = TEMPERATURE IN KELVIN, SHOULD BE BETWEEN 40. AND 300.      
C     FNUMIN = LOWEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU)      
C     FNUMAX = HIGHEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU)     
C     DNU = FREQUENCY INCREMENT IN CM-1. DNU SHOULD BE CHOSEN SO        
C           THAT NOT MORE THAN 600 STEPS ARE NEEDED TO GO FROM          
C           FNUMIN TO FNUMAX (ELSE ARRAY DIMENSIONS OF FREQ,ABSCOEF     
C           MUST BE ADJUSTED IN ADDEM).                                 
C                                                                       
      implicit double precision (a-h,o-z)
      CHARACTER*10 LGAS                                                 
C      COMMON /BLOCKIN/ TEMP,FNUMIN,FNUMAX,DNU                           
      COMMON /LIKEA/ LIKE,LGAS                                           

C these are now parameters
C      DATA TEMP/163./                                                   
C      DATA FNUMIN/0./                                                   
C      DATA FNUMAX/600./                                                 
C      DATA DNU/10./                                                     
      DATA LIKE/1/                                                      
      DATA LGAS/' CH4 - CH4' /                                          
      END                                                               

      subroutine ch4ch4_s(temp,fnumin,fnumax,dnu,nf1,freq,alfatot)

C       ****************************************************************
C     subroutine adapted from Borysow code to calculate CIA absorption
C     of CH4-CH4 collisions (for Titan atmosphere).
C
C     Input variables:
C       temp            double  Temperature (K)
C       fnumin          double  Lowest wavenumber in spectrum (cm-1)
C       fnumax          double  Highest wavenumber in spectrum (cm-1)
C       dnu             double  Wavenumber step (cm-1)
C
C     Output variables
C       nf              integer Number of points in spectrum
C       freq(601)       double  Frequency grid (cm-1)
C       alfatot(601)    double  Absorption coefficient (cm-1 amagat-2)
C
C     C. Nixon 	30-1-97	Original
C     Pat Irwin	2/3/12	Checked for Radtrans2.0
C
C       ****************************************************************

C       PROGRAM GENERATES THE CH4-CH4 CIA SPECTRA at temperatures 50-300
C
C       WRITTEN BY ALEKSANDRA BORYSOW, University of Texas  at AUSTIN 
c	(MARCH 25, 1987); now at Joint Institute of Laboratory Astrophysics, 
c	University of Colorado, Boulder.
C       Program compatible with work described in the paper by
c	A. Borysow & L. Frommhold;
c	ASTROPHYSICAL JOURNAL, VOL.318, pp. 940-943, 1987.
C                                                                       
C       File CORRCH4. contains the STATISTICAL (NUCLEAR) CORRECTIONS    
C                                                                       
	implicit double precision (a-h,o-z)
	character*10 LGAS
      COMMON /ROT/ AMUL0,DELTA                                          
      COMMON /STATT/ QW3(51,3),QW4(51,4)                                 
      COMMON /LIKEA/ LIKE,LGAS                                           
      DIMENSION FREQ(601), ABSCOEF(601), ALFATOT(601)                   
      DATA OMEGA0,DELOM/2.599,0.0071/                                   
      DATA PHI,DELPHI/-8.19,0.033/                                      
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      Y(X,A,B,C)=A*DEXP((C*X+B)*X)                                      

      NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1                                 
C save the value of nf
      nf1 = nf
      IF (NF.GT.601) NF=601                                             
      FNUMAX=FNUMIN+FLOAT(NF-1)*DNU                                     

      if (temp .lt. 50 .or. temp .gt. 300) then
         if(idiag.gt.0)print*, 'ch4ch4_s: Warning'
         if(idiag.gt.0)print*, 'Temperature should be 50 < T < 300 K'
      end if

      CALL PARTSUM8 (TEMP)                                               

C     THE CH4 - CH4 SPECTRA   FOR 50-300K                               
C     =================                                                 
C                                                                       
      X=DLOG(TEMP)                                                      
      DO 10 I=1,NF                                                      
         FREQ(I)=FNUMIN+FLOAT(I-1)*DNU                                  
         ALFATOT(I)=0.0                                                 
   10 ABSCOEF(I)=0.                                                     


      DATA (QW3(J,1),J=1,51) /
     & -1.00000, -1.00000, -1.00000,  4.09091, -1.00000,
     & -1.00000,  .79679,  .27141, -.85139,  .59819,
     & -.00412, -.31878,  .18283,  .15133, -.32477, 
     & .21392, .02191, -.14717,  .07524,  .08569,
     & -.16823,  .10771,  .01794, -.08386,  .04014,
     & .05426, -.10240,  .06447,  .01345, -.05399,
     & .02474,  .03726, -.06876,  .04281,  .01019, 
     & -.03761,  .01670,  .02711, -.04932,  .03046,
     & .00792, -.02769,  .01201,  .02059, -.03708,  
     & .02276,  .00630, -.02122,  .00904,  .01616,
     & -.02889/

      DATA (QW3(J,2),J=1,51) /
     & -1.00000,-1.00000,-1.00000,-1.00000, 1.46154,
     & -1.00000,  .52036, -.38246,  .48307, -.58061,
     &  .54775, -.42488,  .32187, -.30004,  .33377,
     & -.35746,  .33367, -.27900,  .23619, -.22949,
     & .24628, -.25572,  .24034, -.20959, .18643, 
     & -.18371,  .19376, -.19855,  .18803, -.16835,
     & .15390, -.15262,  .15930, -.16211,  .15451,
     & -.14085,  .13100, -.13033,  .13510, -.13691,
     & .13118, -.12115,  .11400, -.11364,  .11722,
     & -.11846,  .11400, -.10631,  .10090, -.10070,
     & .10348/

      DATA (QW3(J,3),J=1,51) /
     & 11.00000,-1.00000,-1.00000, 1.01399,  .09091,
     & -1.00000,  .94475, -.16780, -.51321,  .52338,
     & .01421, -.50745,  .49964, -.06185, -.34591,
     & .34888,  .00564, -.34073,  .33830, -.03176,
     & -.26007,  .26130,  .00304, -.25650,  .25545,
     & -.01926, -.20817,  .20879,  .00191, -.20566,
     & .20511, -.01291, -.17348,  .17383,  .00131,
     & -.17164,  .17132, -.00925, -.14867,  .14889,
     & .00095, -.14727,  .14707, -.00695, -.13006,
     & .13021,  .00073, -.12896,  .12882, -.00541,
     & -.11559/


      DATA (QW4(J,1),J=1,51) /
     & -1.00000,-1.00000,-1.00000,-1.00000,-1.00000,
     & -1.00000, 1.82353,-1.00000, -.00478,  .22586,
     & .19088, -.58488,  .55212, -.21168, -.06443, 
     & .05729,  .14181, -.28953,  .24581, -.07845,
     & -.04207,  .01948,  .09314, -.16876,  .13649,
     & -.03831, -.02777,  .00772,  .06418, -.10976,
     & .08625, -.02197, -.01938,  .00323,  .04653,
     & -.07688,  .05925, -.01398, -.01420,  .00128,
     & .03515, -.05678,  .04315, -.00956, -.01082,
     & .00037,  .02744, -.04362,  .03279, -.00690,
     & -.00850/
      DATA (QW4(J,2),J=1,51) /
     & -1.00000,-1.00000,-1.00000,-1.00000, 2.35664,
     & -1.00000, -.43987,  .63467,  .00794, -.53268,
     & .47071, -.04638, -.25892,  .20786,  .05459,
     & -.23119,  .17454,  .01650, -.14054,  .09758,
     & .03999, -.12597,  .08757,  .02020, -.08639,  
     & .05576,  .02883, -.07869,  .05185,  .01719,
     & -.05811,  .03585,  .02073, -.05366,  .03404,
     & .01393, -.04166,  .02491,  .01572, -.03888,
     & .02397,  .01128, -.03129,  .01828,  .01229,
     & -.02945,  .01775,  .00924, -.02434,  .01397,
     & .00985/
      DATA (QW4(J,3),J=1,51) /
     & -1.00000,-1.00000,-1.00000, 2.16484,-1.00000,
     & -1.00000, 1.60491, -.46218,-1.00000, 1.31017,
     & -.34627, -.78497,  .98223, -.18287, -.72146,
     & .84929, -.16187, -.59662,  .69125, -.09565,
     & -.55198,  .62115, -.09262, -.47526,  .53051,
     & -.05845, -.44470,  .48796, -.05974, -.39348,
     & .42963, -.03932, -.37166,  .40124, -.04166,
     & -.33521,  .36069, -.02824, -.31896,  .34046,
     & -.03069, -.29177,  .31068, -.02125, -.27924,
     & .29557, -.02353, -.25819,  .27279, -.01656,
     & -.024825/
      DATA (QW4(J,4),J=1,51) /
     & 11.00000,-1.00000,-1.00000,  .40260,  .62896,
     & -1.00000,  .67713, -.17620, -.04947, -.07182,
     & .30704, -.40281,  .29894, -.12495,  .04080, 
     & -.09160,  .19422, -.23764,  .18677, -.09887,
     & .05507, -.08289,  .14023, -.16496,  .13482, 
     & -.08186,  .05504, -.07257,  .10916, -.12512,
     & .10520, -.06982,  .05172, -.06377,  .08914,
     & -.10029,  .08615, -.06085,  .04781, -.05661,
     & .07522, -.08345,  .07290, -.05391,  .04406,
     & -.05077,  .06501, -.07133,  .06315, -.04838,
     & .04069/
C                                                                       
C     THESE STATISTICAL FACTORS ARE VALID ONLY FOR DELTA J>0            
C     FOR DELTA J<0 THESE WILL BE EQUAL TO 1.                           
C     METHANE INDUCED COMPONENTS                                        
C     ===============================                                   
C     OCTOPOLE- INDUCED TERM (43)                                       
C                                                                       
      AMUL0=OMEGA0                                                      
      DELTA=DELOM                                                       
      S=Y(X,0.31933D-51/7.,-8.24624D0,0.73213D0)                        
      E=Y(X,-0.61214D0,0.18496D0,-0.04047D0)
      T1=Y(X,0.36548D-11,-0.63873D0,0.00818D0) 
      T2=Y(X,0.16277D-11,-0.52965D0,0.00404D0)                          
      T3=Y(X,0.46747D-11,-0.73107D0,0.01892D0)                          
      T4=Y(X,0.20588D-11,-0.04288D0,-0.04548D0)                         
      CALL ADSPEC1 (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,0,3,3,4)
      DO 60 I=1,NF                                                      
   60 ALFATOT(I)=ABSCOEF(I)+ALFATOT(I)                                  
C                                                                       
C     HEXADECAPOLE-INDUCED TERM                                         
C     =============================                                     
C                                                                       
      AMUL0=PHI                                                         
      DELTA=DELPHI                                                      
      S=Y(X,0.12504D-53/9.,-7.59574D0,0.66118D0)                        
      E=Y(X,-0.5D0,0.0D0,0.0D0)                                         
      T1=Y(X,0.40711D-11,-1.20681D0,0.099D0)                            
      T2=Y(X,0.1743D-12,0.37304D0,-0.08709D0)                           
      T3=Y(X,0.42016D-12,0.09897D0,-0.05468D0)                          
      T4=Y(X,0.15957D-7,-5.10481D0,0.55655D0)                           
      CALL ADSPEC1 (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,0,4,4,5)
      DO 70 I=1,NF                                                      
   70 ALFATOT(I)=ABSCOEF(I)+ALFATOT(I)                                  

   80 FORMAT ( 23H1ABSORPTION SPECTRA OF ,A10,  3H AT,F8.1,  2H K,/1X,43
     1(1H=),/, 11H MIN.FREQ.=,F8.1,  5H CM-1,10X, 10HMAX.FREQ.=,F8.1,  5
     2H CM-1,10X, 15HFREQ.INCREMENT=,F8.2,  5H CM-1,5X,  2HIN,I5,  6H ST
     3EPS,//)                                                           
   90 FORMAT (F8.5)                                                     
  100 FORMAT (10F8.5)                                                   
  110 FORMAT (/,' ABSORPTION COEFFICIENT ALPHA(FNU), FROM', F10.2,      
     1' CM-1 TO', F10.2,' CM-1',/,' AT', F6.2,' CM-1 INCREMENTS',
     1 ' TEMPERATURE =',F7.2, ' K,  IN UNITS OF [CM-1 AMAGAT-2]',/)     
  120 FORMAT (  101 (6E13.5,/))                                         
      return
      end

C ********************* START OF SUBS *****************************

      SUBROUTINE ADSPEC1 (G0,EP,TAU1,TAU2,TAU5,TAU6,TEMP,NF,FREQ,ABSCOEF
     1,MP,LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE)                           
C                                                                       
C     FOR INDUCTION BY  A  T E T R A H E D R A L M O L E C U L E        
C     THIS PROGRAM GENERATES LISTING OF CIA TR ALFA(OMEGA)              
C                                                                       
      implicit double precision (a-h,o-z)
      COMMON /CHPART/ Q1,WCH(2),B01,D01,JRANGE2                         
      COMMON /STATT/ QQW3(51,3),QQW4(51,4)                               
      COMMON /ROT/ AMUL0,DELTA                                          
      DIMENSION ABSCOEF(NF), FREQ(NF)                                   
      DATA CLOSCHM,BOLTZWN/2.68675484D19,.6950304/                      
      DATA HBAR,PI,CLIGHT/1.054588757D-27,3.1415926535898,2.9979250D10/ 
      ECH4(I)=B01*FLOAT(I)                                              
      PCH4(J,T)=DEXP(-1.4387859/T*ECH4(J*(J+1)))                        
      TWOPIC=2.*PI*CLIGHT                                               
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2             
      CALIB=CALIB/FLOAT(1+LIKE)                                         
      BETA=1./(BOLTZWN*TEMP)                                            
      LIST=NF                                                           
      DO 10 I=1,LIST                                                    
   10 ABSCOEF(I)=0.0                                                    

      DO 60 I=1,JRANGE2                                                 
         J=I-1                                                          
         P=FLOAT(2*J+1)*PCH4(J,TEMP)/Q1                                 
         AMUL2=(AMUL0+DELTA*FLOAT(J*(J+1)))**2                          
         DO 40 IP=1,LAMBDA                                              
C                                                                       
C         POSITIVE DELTA J                                              
C                                                                       
            JP=J+IP                                                     
            OMEGA1=ECH4(JP*(JP+1))-ECH4(J*I)                            
            IF (LAMBDA.EQ.3) CC=(1.+QQW3(I,IP)/4.)                      
            IF (LAMBDA.EQ.4) CC=(1.+QQW4(I,IP)/4.)                      
            FAC=CALIB*P*CC                                              
            DO 20 IQ=1,LIST                                             
               FRQ=FREQ(IQ)-OMEGA1                                      
               WKF=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))*FAC*
     1         (2.*DFLOAT(JP)+1.) 
               XBG=G0*BGAMA8(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)*AMUL2      
               ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF                          
   20       CONTINUE                                                    
C                                                                       
C     NEGATIVE DELTA J                                                  
C                                                                       
            JP=J-IP                                                     
            IF (JP.LT.0) GO TO 40                                       
            OMEGA1=ECH4(JP*(JP+1))-ECH4(J*I)                            
            FAC=CALIB*P                                                 
            DO 30 IQ=1,LIST                                             
               FRQ=FREQ(IQ)-OMEGA1                                      
               WKF=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))*FAC*
     1         (2.*DFLOAT(JP)+1.)                                       
               XBG=G0*BGAMA8(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)*AMUL2      
               ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF                          
   30       CONTINUE                                                    
   40    CONTINUE                                                       
C                                                                       
C     DELTA J=0                                                         
C                                                                       
         JP=J                                                           
         FAC=CALIB*P                                                    
         DO 50 IQ=1,LIST                                                
            FRQ=FREQ(IQ)                                                
            WKF=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))*FAC*
     1      (2.*FLOAT(JP)+1.)                                           
            XBG=G0*BGAMA8(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)*AMUL2         
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF                             
   50    CONTINUE                                                       
   60 CONTINUE
      RETURN                                                            
C                                                                       
   70 FORMAT (/, 32H LAMBDA1,LAMBDA2, LAMBDA,LVALUE=,2I3,2X,2I3, 12H COM
     1PONENT .,/15X, 22HLINE SHAPE PARAMETERS:,6E12.3,5X,  5HG(0)=,E12.3
     2/)                                                                
   80 FORMAT ((1X,10E12.4,/))                                           
      END                                                               

      FUNCTION BGAMA8 (FNU,T1,T2,EPS,T3,T4,TEMP)                         
C                                                                       
	implicit double precision (a-h,o-z)
      double precision K0                                                         
      DATA PI,CLIGHT/3.1415926535898,2.99792458D10/                     
      DATA HBAR,BOLTZ/1.0545887D-27,1.380662D-16/                       
      P1(X)=((((((.0045813*X+.0360768)*X+.2659732)*X+1.2067492)*X+3.0899
     1424)*X+3.5156229)*X+1.)                                           
      P2(X)=((((((.00000740*X+.00010750)*X+.00262698)*X+.03488590)*X+.23
     1069756)*X+.42278420)*X-.57721566)                                 
      P3(X)=((((((.00032411*X+.00301532)*X+.02658733)*X+.15084934)*X+.51
     1498869)*X+.87890594)*X+.5)                                        
      P4(X)=((((((-.00004686*X-.00110404)*X-.01919402)*X-.18156897)*X-.6
     17278579)*X+.15443144)*X+1.)                                       
      P5(X)=((((((.00053208*X-.00251540)*X+.00587872)*X-.01062446)*X+.02
     1189568)*X-.07832358)*X+1.25331414)                                
      P6(X)=((((((-.00068245*X+.00325614)*X-.00780353)*X+.01504268)*X-.0
     13655620)*X+.23498619)*X+1.25331414)                               
C                                                                       
      OMEGA=2.*PI*CLIGHT*FNU                                            
      T0=HBAR/(2.*BOLTZ*TEMP)                                           
      Z=DSQRT((1.+(OMEGA*T1)**2)*(T2*T2+T0*T0))/T1                      
      IF (Z-2.) 10,10,20                                                
   10 XK1=Z*Z*DLOG(Z/2.)*P3((Z/3.75)**2)+P4((Z/2.)**2)                  
      GO TO 30                                                          
   20 XK1=DSQRT(Z)*DEXP(-Z)*P6(2./Z)                                    
   30 ZP=DSQRT((1.+(OMEGA*T4)**2)*(T3*T3+T0*T0))/T4                     
      IF (ZP-2.) 40,40,50                                               
   40 K0=-DLOG(ZP/2.)*P1((ZP/3.75)**2)+P2((ZP/2.)**2)                   
      GO TO 60                                                          
   50 K0=DEXP(-ZP)*P5(2./ZP)/DSQRT(ZP)                                  
   60 BGAMA8=((T1/PI)*DEXP(T2/T1+T0*OMEGA)*XK1/
     1  (1.+(T1*OMEGA)**2)+EPS*(T3/
     1PI)*DEXP(T3/T4+T0*OMEGA)*K0)/(1.+EPS)                             
      RETURN                                                            
      END                                                               

      SUBROUTINE PARTSUM8 (TEMP)                                         
	implicit double precision (a-h,o-z)
      COMMON /CHPART/ Q1,WCH(2),B01,D01,JRANGE2                         
      DATA B01,D01,WCH(1),WCH(2)/5.24,0.,1.,1./                         
      ECH4(I)=B01*FLOAT(I)                                              
      PCH4(J,T)=DEXP(-1.4387859/T*ECH4(J*(J+1)))                        
C                                                                       
C     Q1,B01,D01,WCH - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR CH4   
C     *** PARTITION FUNCTION FOR CH4 *******************                
C                                                                       
      Q1=0.                                                             
      J=0                                                               
   10 DQ=PCH4(J,TEMP)*FLOAT((2*J+1)**2)                                 
      Q1=Q1+DQ                                                          
      J=J+1                                                             
      IF (DQ.GT.Q1/1000.) GO TO 10                                      
      JRANGE2=J+3                                                       
      RETURN                                                            
      end
      
