      subroutine h2ch4_sub(temp,fnumin,fnumax,dnu,
     >					nf1,freq,alfatot,slit1)
C       ****************************************************************
C     subroutine adapted from Borysow code to calculate CIA absorption
C     of H2-CH4 collisions.
c
c	_9 added to all function/sunroutines/common block names to avoid confusion
c	with similar routines when compiled together into a library.
c
cIMPLICIT DOUBLE PRECISION replaced with IMPLICIT double precision (to be same as n2n2_s.f etc)
C
C     Input variables:
C	temp			double	Temperature (K)
C	fnumin		double	Lowest wavenumber in spectrum (cm-1)
C	fnumax		double	Highest wavenumber in spectrum (cm-1)
C	dnu			double	Wavenumber step (cm-1)
C     slit1            	double  	Slit size
C
C     Output variables
C	nf			integer	Number of points in spectrum
C	freq(601)		double	Frequency grid (cm-1)
C	alfatot(601)	double	Absorption coefficient (cm-1 amagat-2)
C
C     Nick Teanby	23-06-04
C     Pat Iriwn		2/3/12	Updated for Radtrans2.0
c
C
c	=============================================
c	Copyright, Aleksandra Borysow, 1988
c	=============================================
c	If you use this program, please cite original work:
c       A. Borysow and L. Frommhold,
c       "Theoretical collision induced rototranslational absorption
c       spectra for the outer planets: H2-CH4 pairs",
c	Astrophysical Journal, vol. 304, p.849-865, (1986)
c	=============================================

C     THIS PROGRAM MUST BE ACCOMPANIED BY FILE corrch4
C     THESE DATA SERVE AS THE INPUT TO THE MAIN PROGRAM ADDEM.          
C     TEMP = TEMPERATURE IN KELVIN, SHOULD BE BETWEEN 40. AND 300.      
C     FNUMIN = LOWEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU)      
C     FNUMAX = HIGHEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU)     
C     DNU = FREQUENCY INCREMENT IN CM-1. DNU SHOULD BE CHOSEN SO        
C           THAT NOT MORE THAN 600 STEPS ARE NEEDED TO GO FROM          
C           FNUMIN TO FNUMAX (ELSE ARRAY DIMENSIONS OF FREQ,ABSCOEF     
C           MUST BE ADJUSTED IN ADDEM).                                 
C                                                                       
C                                                                       
C       ****************************************************************
C       PROGRAM PREPARED BY ALEKSANDRA BORYSOW, UNIVERSITY OF TEXAS @AUS
C       AND JOINT INSTITUTE FOR LABORATORY ASTROPHYSICS, UNIV. COLORADO 
C       LAST CORRECTION DATE: 14 NOVEMBER 1988                          
C       THIS PROGRAM IS COMPATIBLE WITH PAPER: A. BORYSOW & L. FROMMHOLD
C       ASTROPHYSICAL JOURNAL; VOL. 304, PP.849-865; 1986.              
C       ****************************************************************
C      PROGRAM GENERATES H2-CH4 FREE-FREE, BOUND-FREE & BOUND-BOUND
C      COLLISION INDUCED ABSORPTION SPECTRA  
C                                                                       
C     OUTPUT: file output 						ADD00036
C     File corrch4 contains STATISTICAL (NUCLEAR) CORRECTIONS           
C     FOR METHANE INDUCTION                                             
C                                                                       
       implicit double precision (A-H,O-Z)
      Character*10 LGAS
      COMMON /APP3_9/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT                
      COMMON /RSILO_9/ RSILO(201)                                         
      COMMON /BOU43_9/ INITB                                              
      COMMON /BB_9/ OMEG,RSI,RSIGG,ALFA,SCAL,NSOL                        
      COMMON /BF_9/ G0BF,DELBF,OM0                                        
      COMMON /STAT_9/ QW3(51,3),QW4(51,4)                                 
      COMMON /LIKE_9/ LIKE,LGAS                                           
      DIMENSION FREQ(601), ABSCOEF(601), ALFATOT(601)                   
      DIMENSION RSI(201), RSIGG(201), TT(2), SS(1), OMEG(201), AIG(201) 
      Y(X,A,B,C)=A*DEXP((C*X+B)*X)                                      
C                                                                       
      LIKE = 0                                                      
      LGAS = ' H2 - CH4'                                           

      if (slit1.lt.4.29999999D0) then
        pause 'ERROR: h2ch4_sub: DONT PUT SLIT < 4.3!'
      endif
      slit=slit1
      NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1                                 
      IF (NF.GT.601) NF=601 
c     save value of nf
      NF1=NF                                            
      FNUMAX=FNUMIN+DFLOAT(NF-1)*DNU                                    
      CALL PARTSUM_9 (TEMP)                                               
C                                                                       
C     THE H2 - CH4 SPECTRA   FOR 50-300K                                
C     =================                                                 
C                                                                       
C     ONLY B-F + F-F TERMS INCLUDED IN MODELLING                        
C     TEMPERATURE INTERPOLATION GOOD FOR 45-300K                        
C                                                                       
      X=DLOG(TEMP)                                                      
      DO 10 I=1,NF                                                      
         FREQ(I)=FNUMIN+DFLOAT(I-1)*DNU                                 
         ALFATOT(I)=0.0                                                 
   10 ABSCOEF(I)=0.                                                     

      call setqw_9(qw3,qw4)

                                                                      
C     THESE STATISTICAL FACTORS ARE VALID ONLY FOR DELTA J>0            
C     FOR DELTA J<0 THESE WILL BE EQUAL TO 1.                           
C     HYDROGEN'S QUADRUPOLE                                             
C                                                                       
      EPS=1.D-5                                                         
      TT(1)=10.D0                                                       
      CALL BOUND32_9 (TEMP,RSI,NSOL)                                      
      DO 60 I=1,NSOL                                                    
         RSILO(I)=DLOG(RSI(I)*1.D80)                                    
   60 OMEG(I)=DFLOAT(I-1)*DX                                            
      CALL SPLINE_9 (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)           
      SCAL=1.D0                                                         
      S=Y(X,.22916D-59,-1.27134D0, .12423D0)
      E=Y(X,-.2019D1,-.02259d0,-.05891D0)                               
      T1=Y(X,.48254D-13,.6421d0,-.10109D0)                              
      T2=Y(X,.97826D-12,-.48654d0,-.0146D0)                             
      T3=Y(X,.35247D-12,.10164d0,-.07879D0)                             
      T4=Y(X,.13961D-13,1.11146d0,-.09561D0)                            
C                                                                       
C     THIS PART FOR MODELING A LOW FREQUENCY PART OF BOUND-FREE         
C     TRANSLATIONAL SPECTRAL FUNCTION, BY A DESYMMETRIZED GAUSSIAN      
C     PROFILE                                                           
C                                                                       
      G0BF=Y(X,0.73536D-72,-.79815D0,-.0585D0)                          
      DELBF=Y(X,3.123D0,-.00178D0,0.00021D0)                            
      OM0=Y(X,8.6922D0,0.00785D0,-0.00054D0)                            
      CALL ADDSPEC_9 (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,2,0,2,
     >3)
      DO 70 I=1,NF                                                      
   70 ALFATOT(I)=ABSCOEF(I)+ALFATOT(I)                                  
C                                                                       
C     PARAMETERS FOR 4045  (HYDROGEN HEXADECAPOLE) COMPONENTS           
C                                                                       
      CALL BOUN54C_9 (TEMP,RSI,NSOL)                                      
      DO 80 I=1,NSOL                                                    
         RSILO(I)=DLOG(RSI(I)*1.D80)                                    
   80 OMEG(I)=DFLOAT(I-1)*DX                                            
      CALL SPLINE_9 (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)           
      S=Y(X,(2.7321D-60/9.)*0.038,-2.32012D0,.23082D0)                  
      E=Y(X,-1.8198D0,-.00665D0,-.05626D0)                              
      T1=Y(X,1.3492D-13,.14472D0,-.06506D0)                             
      T2=Y(X,1.6596D-12,-.77679D0,.01401D0)                             
      T3=Y(X,5.9914D-13,-.16208D0,-.05895D0)                            
      T4=Y(X,1.9405D-14,.95965D0,-.10246D0)                             
      SCAL=0.038D0                                                      
      CALL ADDSPEC_9 (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,4,0,4,
     >5)
      DO 90 I=1,NF                                                      
   90 ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)                                  
C                                                                       
C     METHANE INDUCED COMPONENTS                                        
C     ===============================                                   
C     OCTOPOLE- INDUCED TERM (43)                                       
C                                                                       
      EPS=1.D-5                                                         
      TT(1)=10.D0                                                       
      CALL BOUND43_9 (TEMP,RSI,NSOL)                                      
      DO 100 I=1,NSOL                                                   
         RSILO(I)=DLOG(RSI(I)*1.D80)                                    
  100 OMEG(I)=FLOAT(I+INITB-2)*DX                                       
      CALL SPLINE_9 (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)           
      SCAL=1.D0                                                         
      S=Y(X,.91196D-60/7.,-1.56529D0,.15284D0)                          
      E=Y(X,-.5D0,0.D0,0.D0)                                            
      T1=Y(X,.51456D-13,.60523D0,-.10856D0)                             
      T2=Y(X,.62514D-12,-.51384D0,.00754D0)                             
      T3=Y(X,.55346D-12,-.40381D0,.00208D0)                             
      T4=Y(X,.13804D-13,1.9307D0,-.27921D0)                             
      CALL ADSPEC1_9 (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,0,3,3,
     >4)
      DO 110 I=1,NF                                                     
  110 ALFATOT(I)=ABSCOEF(I)+ALFATOT(I)                                  
C                                                                       
C     HEXADECAPOLE-INDUCED TERM                                         
C     =============================                                     
C                                                                       
      EPS=1.D-5                                                         
      TT(1)=10.D0                                                       
      CALL BOUN54C_9 (TEMP,RSI,NSOL)                                      
      DO 120 I=1,NSOL                                                   
         RSILO(I)=DLOG(RSI(I)*1.D80)                                    
  120 OMEG(I)=DFLOAT(I-1)*DX                                            
      CALL SPLINE_9 (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)           
      SCAL=1.D0                                                         
      S=Y(X,2.7321D-60/9.,-2.32012D0,.23082D0)                          
      E=Y(X,-1.8198D0,-.00665D0,-.05626D0)                              
      T1=Y(X,1.3492D-13,.14472D0,-.06506D0)                             
      T2=Y(X,1.6596D-12,-.77679D0,.01401D0)                             
      T3=Y(X,5.9914D-13,-.16208D0,-.05895D0)                            
      T4=Y(X,1.9405D-14,.95965D0,-.10246D0)                             
      CALL ADSPEC1_9 (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,0,4,4,
     >5)
      DO 130 I=1,NF                                                     
  130 ALFATOT(I)=ABSCOEF(I)+ALFATOT(I)                                  

  140 FORMAT (' ABSORPTION SPECTRA OF ', A10,' AT', F8.1,
     1  'K',/1X,
     1 ' MIN.FREQ.=', F8.1, ' CM-1', 10X, ' MAX.FREQ.=', F8.1,  
     2 ' CM-1',/, ' FREQ.INCREMENT=', F8.2, ' CM-1', 5X, 'IN', I5,  
     1  ' STEPS',/)                                                     
  150 FORMAT (F8.5)                                                     
  160 FORMAT (10F8.5)                                                   
  170 FORMAT (/, ' ABSORPTION COEFFICIENT ALPHA(FNU), FROM ',F5.1,      
     1' CM-1 TO', F7.1,  9H CM-1, AT, F6.2, /, 
     1  23H CM-1 INCREMENTS, AT T=,F7.2                                 
     2, 29H K, IN UNITS OF CM-1 AMAGAT-2,/)                             
  180 FORMAT (  1H ,10E13.5)                                            
  190 FORMAT (4F10.5,I5)                                                
C        
	return                                                               
      END
      
                                                                     
      SUBROUTINE ADDSPEC_9 (G0,EP,TAU1,TAU2,TAU5,TAU6,TEMP,NF,FREQ,
     1ABSCOEF,MP,LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE)                           
C                                                                       
C     THIS ENTRY FOR L I N E A R   M O L E C U L E!!!!!                 
C     SET LAMBDA1 NONZERO                                               
C     THIS PROGRAM GENERATES LISTING OF CIA TR ALFA(OMEGA)              
C     IF EITHER LAMBDA1 OR LAMBDA2 EQUAL TO ZERO - SINGLE TRANSITIONS;  
C     LAMBDA1 CORRESPONDS TO H2, LAMBDA2 TO CH4                         
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.              
C                                                                       
       implicit double precision (A-H,O-Z)
      COMMON /BB_9/ OMEG(201),RSI(201),RSIGG(201),BETA,SCAL,NSOL          
      COMMON /RSILO_9/ RSILO(201)                                         
      COMMON /APP3_9/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT                
      COMMON /H2PART_9/ Q,WH2(2),B0,D0,JRANGE1                            
      COMMON /BF_9/ G0BF,DELBF,OM0                                        
      COMMON /CHPART_9/ Q1,WCH(2),B01,D01,JRANGE2                         
      DIMENSION ABSCOEF(NF), FREQ(NF)                                   
      DATA CLOSCHM,BOLTZWN/2.68675484D19,.6950304/                      
      DATA HBAR,PI,CLIGHT/1.054588757D-27,3.1415926535898,2.9979250D10/ 
      EH2(I)=(B0-DFLOAT(I)*D0)*DFLOAT(I)                                
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*DEXP(-1.4387859/T*
     1 EH2(J*(J+1)))                                                    
      TWOPIC=2.*PI*CLIGHT                                               
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2             
      CALIB=CALIB/DFLOAT(1+LIKE)                                        
      BETA=1./(BOLTZWN*TEMP)                                            
      LIST=NF                                                           
      DO 10 I=1,LIST                                                    
   10 ABSCOEF(I)=0.0                                                    
C                                                                       
C     ROTATIONAL SPECTRUM FOR THE DETAILED LISTING   *******************
C                                                                       
c      WRITE (6,50) LAMBDA1,LAMBDA2,LAMBDA,LVALUE,G0,EP,TAU1,TAU2,
c     1 TAU5,TAU6,G0*BGAMA_9(0.,TAU1,TAU2,EP,TAU5,TAU6,TEMP)               
									
      IF (LAMBDA1.EQ.0) STOP 6666                                       
C                                                                       
C     FOR THIS ENTRY LAMBDA1 HAS TO BE NONZERO                          
C     SINGLE TRANSITIONS ON HYDROGEN'S ROTATIONAL FREQUENCIES.          
C     LAMBDA IS EQUAL TO LAMBDA1 FOR SINGLE TRANSITIONS                 
C                                                                       
      JPLUSL=JRANGE1+LAMBDA                                             
      DO 40 I=1,JRANGE1                                                 
         J=I-1                                                          
      DO 40 IP=1,JPLUSL                                                 
         JP=IP-1                                                        
         CGS=CLEBSQR_9(J,LAMBDA,JP)                                       
         IF (CGS) 40,40,20                                              
   20    P=PH2(J,TEMP)/Q                                                
         OMEGA1=EH2(JP*IP)-EH2(J*I)                                     
         FAC=CALIB*P*CGS                                                

         DO 30 IQ=1,LIST                                                
            FRQ=FREQ(IQ)-OMEGA1                                         
            WKI=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))                      
            WKF=WKI*FAC                                                 
            XBG=G0*BGAMA_9(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)               


            IF ( DABS(FRQ).LE.WNRMAX3 ) 
     1       XBG=XBG+SCAL*SPECFCT_9(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)

            IF (LVALUE.EQ.3 .AND. G0BF.NE.0.D0) 
     1      XBG = XBG + BGAUS_9(FRQ, G0BF, DELBF, OM0, TEMP)

C                                                                       
C     THIS MODELS THE PART OF THE BOUND-FREE SPECTRAL FUNCTION BY      
C     MEANS OF THE DESYMMETRIZED GAUSSIAN PROFILE                       
C                                                                       
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF                             
   30    CONTINUE                                                       
   40 CONTINUE                                                          
      RETURN                                                            
                                                                       
   50 FORMAT (/, 32H LAMBDA1,LAMBDA2, LAMBDA,LVALUE=,2I3,2X,2I3, 12H COM
     1PONENT .,/15X, 22HLINE SHAPE PARAMETERS:,6E12.3,5X,  5HG(0)=,E12.3
     2/)                                                                
   60 FORMAT ((1X,10E12.4,/))                                           
                                                                       
      END                                                               

      SUBROUTINE ADSPEC1_9 (G0,EP,TAU1,TAU2,TAU5,TAU6,TEMP,NF,FREQ,
     1ABSCOEF,MP,LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE)                           
C                                                                       
C     FOR INDUCTION BY A  T E T R A H E D R A L   M O L E C U L E!!     
C     THIS PROGRAM GENERATES LISTING OF CIA TR ALFA(OMEGA)              
C     LAMBDA1 CORRESPONDS TO H2, LAMBDA2 TO CH4                         
C     FOR THIS ENTRY SET LAMBDA2 NONZERO                                
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.              
C                                                                       
       implicit double precision (A-H,O-Z)
      COMMON /BB_9/ OMEG(201),RSI(201),RSIGG(201),BETA,SCAL,NSOL          
      COMMON /RSILO_9/ RSILO(201)                                         
      COMMON /APP3_9/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT                
      COMMON /H2PART_9/ Q,WH2(2),B0,D0,JRANGE1                            
      COMMON /CHPART_9/ Q1,WCH(2),B01,D01,JRANGE2                         
      COMMON /STAT_9/ QQW3(51,3),QQW4(51,4)                               
      DIMENSION ABSCOEF(NF), FREQ(NF)                                   
      DATA CLOSCHM,BOLTZWN/2.68675484D19,.6950304/                      
      DATA HBAR,PI,CLIGHT/1.054588757D-27,3.1415926535898,2.9979250D10/ 
      ECH4(I)=B01*DFLOAT(I)                                             
      PCH4(J,T)=DEXP(-1.4387859/T*ECH4(J*(J+1)))                        
      TWOPIC=2.*PI*CLIGHT                                               
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2             
      CALIB=CALIB/DFLOAT(1+LIKE)                                        
      BETA=1./(BOLTZWN*TEMP)                                            
      LIST=NF                                                           
      DO 10 I=1,LIST                                                    
   10 ABSCOEF(I)=0.0                                                    
C                                                                       
C     ROTATIONAL SPECTRUM: DETAILED LISTING   *******************       
C                                                                       
c      WRITE (6,70) LAMBDA1,LAMBDA2,LAMBDA,LVALUE,G0,EP,TAU1,TAU2,
c     1 TAU5,TAU6, G0*BGAMA_9(0.,TAU1,TAU2,EP,TAU5,TAU6,TEMP)              
C                                                                       
      DO 60 I=1,JRANGE2                                                 
         J=I-1                                                          
         P=DFLOAT(2*J+1)*PCH4(J,TEMP)/Q1                                
         DO 40 IP=1,LAMBDA                                              
C                                                                       
C     POSITIVE DELTA J                                                  
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
               XBG=G0*BGAMA_9(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)            
               
               IF ( DABS(FRQ).LE.WNRMAX3 ) 
     1         XBG= XBG + SPECFCT_9(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)
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
               XBG=G0*BGAMA_9(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)            
               
               IF ( DABS(FRQ).LE.WNRMAX3 ) 
     1     XBG=XBG+SPECFCT_9(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)
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
     1      (2.*DFLOAT(JP)+1.) 
            XBG=G0*BGAMA_9(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)               
            
            IF (DABS(FRQ).LE.WNRMAX3) 
     1   XBG=XBG+SPECFCT_9(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)   
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF                             
   50    CONTINUE                                                       
   60 CONTINUE                                                          
c      WRITE(6,80) (ABSCOEF(I),I=1,LIST)                                 
      RETURN                                                            
C                                                                       
   70 FORMAT (/, 32H LAMBDA1,LAMBDA2, LAMBDA,LVALUE=,2I3,2X,2I3, 12H COM
     1PONENT .,/15X, 22HLINE SHAPE PARAMETERS:,6E12.3,5X,  5HG(0)=,E12.3
     2/)                                                                
   80 FORMAT ((1X,10E12.4,/))                                           
C                                                                       
      END                                                               
      FUNCTION CLEBSQR_9 (L,LAMBDA,LP)                                    
C                                                                       
C     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)          
C     FOR INTEGER ARGUMENTS ONLY                                        
C     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.                
C                                                                       
       implicit double precision (A-H,O-Z)
      FC=DFLOAT(2*LP+1)                                                 
      GO TO 10                                                          
C                                                                       
      ENTRY THREEJ2_9                                                    
C                                                                       
C     THIS ENTRY RETURNS THE SQUARED 3-J SYMBOL   L LAMBDA LP           
C                                                 0    0    0           
C     INSTEAD OF THE CLEBSCH-GORDAN COEFFICIENT                         
C     (LIMITATION TO INTEGER ARGUMENTS ONLY)                            
C                                                                       
C     NOTE THAT THE THREE-J SYMBOLS ARE COMPLETELY SYMMETRIC IN THE     
C     ARGUMENTS. IT WOULD BE ADVANTAGEOUS TO REORDER THE INPUT ARGUMENT 
C     LIST SO THAT LAMBDA BECOMES THE SMALLEST OF THE 3 ARGUMENTS.      
C                                                                       
      FC=1.                                                             
   10 CLEBSQR_9=0.                                                        
      IF (((L+LAMBDA).LT.LP).OR.((LAMBDA+LP).LT.L).OR.((L+LP).LT.LAMBDA)
     1) RETURN                                                          
      IF (MOD(L+LP+LAMBDA,2).NE.0) RETURN                               
      IF ((L.LT.0).OR.(LP.LT.0).OR.(LAMBDA.LT.0)) RETURN                
      F=1./DFLOAT(L+LP+1-LAMBDA)                                        
      IF (LAMBDA.EQ.0) GO TO 30                                         
      I1=(L+LP+LAMBDA)/2                                                
      I0=(L+LP-LAMBDA)/2+1                                              
      DO 20 I=I0,I1                                                     
   20 F=F*DFLOAT(I)/FLOAT(2*(2*I+1))                                    
   30 P=FC*F*FCTL_9(LAMBDA+L-LP)*FCTL_9(LAMBDA+LP-L)                        
      CLEBSQR_9=P/(FCTL_9((LAMBDA+L-LP)/2)*FCTL_9((LAMBDA+LP-L)/2))**2        
      RETURN                                                            
C                                                                       
      END                                                               
      FUNCTION FCTL_9 (N)                                                 
       implicit double precision (A-H,O-Z)
      P(Z)=((((-2.294720936D-4)/Z-(2.681327160D-3))/Z+(3.472222222D-3))/
     1Z+(8.333333333D-2))/Z+1.                                          
      FCTL_9=1.                                                           
      IF (N.LE.1) RETURN                                                
      IF (N.GT.15) GO TO 20                                             
      J=1                                                               
      DO 10 I=2,N                                                       
   10 J=J*I                                                             
      FCTL_9=DFLOAT(J)                                                    
      RETURN                                                            
   20 Z=DFLOAT(N+1)                                                     
      FCTL_9=(DEXP(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)                 
      RETURN                                                            
C                                                                       
      END                                                               
      FUNCTION BGAMA_9 (FNU,T1,T2,EPS,T3,T4,TEMP)                         
C                                                                       
C     EQUATION 13, SO-CALLED EBC MODEL, OF BORYSOW,TRAFTON,FROMMHOLD,   
C     AND BIRNBAUM, ASTROPHYS.J., (1985)                                
C                                                                       
       implicit double precision (A-H,O-Z)
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
   60 BGAMA_9=((T1/PI)*DEXP(T2/T1+T0*OMEGA)*XK1/
     1   (1.+(T1*OMEGA)**2)+EPS*(T3/PI)*DEXP(T3/T4+T0*OMEGA)*K0)/
     1   (1.+EPS)                              			
      RETURN                                                            
      END                                                               

      SUBROUTINE PARTSUM_9 (TEMP)                                         
C                                                                       
C     H2 ROTATIONAL PARTITION SUM Q = Q(T).                             
C                                                                       
       implicit double precision (A-H,O-Z)
      COMMON /H2PART_9/ Q,WH2(2),B0,D0,JRANGE1                            
      COMMON /CHPART_9/ Q1,WCH(2),B01,D01,JRANGE2                         
      DATA B0,D0,WH2(1),WH2(2)/59.3392,0.04599,1.,3./                   
      DATA B01,D01,WCH(1),WCH(2)/5.24,0.,1.,1./                         
      EH2(I)=(B0-DFLOAT(I)*D0)*DFLOAT(I)                                
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*DEXP(-1.4387859*
     1 EH2(J*(J+1))/T)                                                  
      ECH4(I)=B01*DFLOAT(I)                                             
      PCH4(J,T)=DEXP(-1.4387859/T*ECH4(J*(J+1)))                        
C                                                                       
C     Q,B01,D01,WCH - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR CH4    
C     Q,B0,D0,WH2 - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR H2       
C                                                                       
      Q=0.                                                              
      J=0                                                               
   10 DQ=PH2(J,TEMP)                                                    
      Q=Q+DQ                                                            
      J=J+1                                                             
      IF (DQ.GT.Q/900.) GO TO 10                                        
      JRANGE1=J                                                         
C                                                                       
C     *** PARTITION FUNCTION FOR CH4 *******************                
C                                                                       
      Q1=0.                                                             
      J=0                                                               
   20 DQ=PCH4(J,TEMP)*DFLOAT((2*J+1)**2)                                
      Q1=Q1+DQ                                                          
      J=J+1                                                             
      IF (DQ.GT.Q1/1000.) GO TO 20                                      
      JRANGE2=J+3                                                       
      RETURN                                                            
C                                                                       
   30 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF H2: Q=,F8.2,10X,  
     17HJ MAX =,I3/)                                                    
C                                                                       
      END                                                               
      SUBROUTINE PROFILE_9 (X,Y)                                          
       implicit double precision (A-H,O-Z)
      COMMON /BL3_9/ RSI(401)                                             
C                                                                       
C     ATRIANGULAR SLIT FUNCTION IS USED.                                
C                                                                       
      COMMON /APP3_9/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT                
      IF (Y) 50,60,10                                                   
   10 CONTINUE                                                          
      X0=(NSRI+1.)+X/DX                                                 
      NC=X0                                                             
      N1=NC+1                                                           
      SLOPE=Y/SLIT                                                      
      NU=X0-NS                                                          
      IF (NU.LT.1) NU=1                                                 
      IF (NU.GT.NSRIUP) RETURN                                          
      NO=X0+NS                                                          
      IF (NO.GT.NSRIUP) NO=NSRIUP                                       
      IF (NO.LT.1) RETURN                                               
      IF (NC.GT.NSRIUP) NC=NSRIUP                                       
      IF (NC.LE.1) GO TO 30                                             
      DO 20 I=NU,NC                                                     
         XI=(I-1.)*DX-WNRMAX3                                           
         DR=SLOPE*(XI-(X-SLIT))                                         
         IF (DR.LE.0.) GO TO 20                                         
         RSI(I)=RSI(I)+DR                                               
   20 CONTINUE                                                          
   30 IF (NC.GE.NSRIUP) RETURN                                          
      IF (N1.LT.1) N1=1                                                 
      DO 40 I=N1,NO                                                     
         XI=(I-1.)*DX-WNRMAX3                                           
         DR=Y-SLOPE*(XI-X)                                              
         IF (DR.LE.0.) GO TO 40                                         
         RSI(I)=RSI(I)+DR                                               
   40 CONTINUE                                                          
      RETURN                                                            
c   50 WRITE(6,70) SLIT                                                  
   50 continue                                                  
   60 CONTINUE
      RETURN                                                            
C                                                                       
   70 FORMAT (/, 30H A TRIANGULAR SLIT FUNCTION OF,F6.3, 23H CM-1 HALFWI
     1DTH IS USED,/)                                                    
C                                                                       
      END                                                               
      FUNCTION SPECFCT_9 (FREQ,OMEGA,PHI,PHI2,N,RTEMP)                    
C                                                                       
C     THIS INTERPOLATES THE SPECTRAL FUNCTION PHI(FREQ) DEFINED AT      
C     OMEGA(N) AS PHI(N). PHI2 IS THE SECOND DERIVATIVE AT OMEGA        
C     WHICH MUST BE OBTAINED FIRST (USE SPLINE_9 FOR THAT PURPOSE).       
C     RTEMP IS THE RECIPROCAL TEMPERATURE IN CM-1 UNITS.                
C     NOTE THAT WE INTERPOLATE 1.E80 TIMES THE LOGARITHM OF PHI(OMEGA)  
C                                                                       
	implicit double precision (A-h,o-z)
      DIMENSION PHI(N), PHI2(N), OMEGA(N)                               
	dimension f(1), gp(1)
c	correction above  1-MAR-1989 18:10:02 in all subroutine
      TFAC=0.                                                           
      F(1)=FREQ                                                         
      IF (F(1) ) 10,20,20                                               
   10 F(1)=dABS(F(1))                                                   
      TFAC=(-RTEMP*F(1))                                                
   20 IF (F(1).LE.OMEGA(N)) GO TO 30                                    
      SPECFCT_9=DEXP(-(PHI(N-1)-PHI(N))*(F(1)-OMEGA(N))/
     1  (OMEGA(N)-OMEGA(N-1))+                                          
     1PHI(N)+TFAC)*(1.D-80)                                             
      RETURN                                                            
c	exchanged ixpolat to spline   1-MAR-1989 18:08:27 
   30 CALL SPLINE_9  (N,1,0,1.D-6,OMEGA,PHI,F,GP,SI,NR,PHI2)              
c	f(1), gp(1)
      SPECFCT_9=DEXP(TFAC+GP(1))*(1.D-80)                                 
      RETURN                                                            
C                                                                       
      END                                                               

      SUBROUTINE SPLINE_9 (L,M,K,EPS,X,Y,T,SS,SI,NR,S2)                   
C                                                                       
C     SPLINE INTERPOLATION AND QUADRATURE, THIRD ORDER AFTER GREVILLE.  
C     INPUT ARGUMENTS L...Y, OUTPUT SS...NR.                            
C     L DATA POINTS X(1), Y(1) ... X(L),Y(L)                            
C     EPS=ERROR CRITERION, TYPICALLY EPS=1.E-5 FOR 5 DECI. PLACES ACCURA
C     M ARGUMENTS T(1)..T(M) FOR WHICH FUNCTION VALUES SS(1)..SS(M), FOR
C     K=0; OR FIRST OR SECOND DERIVATIVE FOR K=1 OR -1, RESPECTIVELY.   
C     NOTE THAT M HAS TO BE AT LEAST EQUAL TO 1.                        
C     SI=INTEGRAL (OVER WHOLE INTERVAL) FOR K=2 ONLY.                   
C     FOR 'NATURAL' SPLINE FUNCTIONS, S2(1)=S2(L)=0. MUST BE INPUT*NOTE*
C     N0 INDICATES THE NUMBER OF OUT-OF-RANGE CALLS. X(1)@T(I)@X(L)     
C     EXTRAPOLATE WITH CAUTION. (ASSUMPTION D2Y/DX2 = 0.)               
C     S2(I) IS THE 2ND DERIVATIVE AT X=X(I) AND IS COMPUTED INTERNALLY. 
C                                                                       
       implicit double precision (A-H,O-Z)
      DIMENSION X(L), Y(L), T(M), SS(M), S2(L)                          
      DELY(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))                               
      B(I)=(X(I)-X(I-1))*0.5/(X(I+1)-X(I-1))                            
      C(I)=3.*(DELY(I)-DELY(I-1))/(X(I+1)-X(I-1))                       
      N=L                                                               
      N1=N-1                                                            
      NR=0                                                              
      DO 10 I=2,N1                                                      
   10 S2(I)=C(I)/1.5                                                    
      OMEGA=1.0717968                                                   
      IC=0                                                              
C                                                                       
C     'NATURAL' SPLINE FUNCTIONS OF THIRD ORDER.                        
C                                                                       
      S2(N)=0.                                                          
      S2(1)=S2(N)                                                       
   20 ETA=0.                                                            
      IC=IC+1                                                           
      SM=DABS(S2(1))                                                    
      DO 30 I=2,N                                                       
         IF (DABS(S2(I)).GT.SM) SM=DABS(S2(I))                          
   30 CONTINUE                                                          
      EPSI=EPS*SM                                                       
      DO 50 I=2,N1                                                      
         W=(C(I)-B(I)*S2(I-1)-(0.5-B(I))*S2(I+1)-S2(I))*OMEGA           
         IF (DABS(W)-ETA) 50,50,40                                      
   40    ETA=DABS(W)                                                    
   50 S2(I)=S2(I)+W                                                     
      IF (ETA-EPSI) 60,20,20                                            
      ENTRY IXPOLAT_9                                                     
C                                                                       
C     THIS ENTRY USEFUL WHEN ITERATION PREVIOUSLY COMPLETED             
C                                                                       
      N=L                                                               
      N1=N-1                                                            
      NR=0                                                              
      IC=-1                                                             
   60 IF (K.EQ.2) GO TO 260                                             
      GO TO 70                                                          
   70 DO 250 J=1,M                                                      
         I=1                                                            
         IF (T(J)-X(1)) 110,210,80                                      
   80    IF (T(J)-X(N)) 100,190,150                                     
   90    IF (T(J)-X(I)) 200,210,100                                     
  100    I=I+1                                                          
         GO TO 90                                                       
  110    NR=NR+1                                                        
         HT1=T(J)-X(1)                                                  
         HT2=T(J)-X(2)                                                  
         YP1=DELY(1)+(X(1)-X(2))*(2.*S2(1)+S2(2))/6.                    
         IF (K) 140,130,120                                             
  120    SS(J)=YP1+HT1*S2(1)                                            
         GO TO 250                                                      
  130    SS(J)=Y(1)+YP1*HT1+S2(1)*HT1*HT1/2.                            
         GO TO 250                                                      
  140    SS(J)=S2(I)                                                    
         GO TO 250                                                      
  150    HT2=T(J)-X(N)                                                  
         HT1=T(J)-X(N1)                                                 
         NR=NR+1                                                        
         YPN=DELY(N1)+(X(N)-X(N1))*(S2(N1)+2.*S2(N))/6.                 
         IF (K) 180,170,160                                             
  160    SS(J)=YPN+HT2*S2(N)                                            
         GO TO 250                                                      
  170    SS(J)=Y(N)+YPN*HT2+S2(N)*HT2*HT2/2.                            
         GO TO 250                                                      
  180    SS(J)=S2(N)                                                    
         GO TO 250                                                      
  190    I=N                                                            
  200    I=I-1                                                          
  210    HT1=T(J)-X(I)                                                  
         HT2=T(J)-X(I+1)                                                
         PROD=HT1*HT2                                                   
         S3=(S2(I+1)-S2(I))/(X(I+1)-X(I))                               
         SS2=S2(I)+HT1*S3                                               
         DELSQS=(S2(I)+S2(I+1)+SS2)/6.                                  
         IF (K) 240,220,230                                             
  220    SS(J)=Y(I)+HT1*DELY(I)+PROD*DELSQS                             
         GO TO 250                                                      
  230    SS(J)=DELY(I)+(HT1+HT2)*DELSQS+PROD*S3/6.                      
         GO TO 250                                                      
  240    SS(J)=SS2                                                      
  250 CONTINUE                                                          
  260 SI=0.                                                             
      DO 270 I=1,N1                                                     
         H=X(I+1)-X(I)                                                  
  270 SI=SI+0.5*H*(Y(I)+Y(I+1))-H**3*(S2(I)+S2(I+1))/24.                
      IF (K.EQ.2) NR=IC                                                 
      RETURN                                                            
      END                                                               

      SUBROUTINE BOUND32_9 (TEMP,RSI,NSOL)                                
       implicit double precision (A-H,O-Z)
      DIMENSION EB(2,9), A(9,9), RSI(201)                               
C                                                                       
C     EB(I,K) - BOUND ENERGIES FOR VIB.NR.V=I-1, ROT.NR J=K-1           
C     A(I,J)=( C(I,L,J)**2) * (MTX.EL. (I,BETA,J)**2 )                  
C     I,J-DENOTE ROTATIONAL STATES FOR VIB. LEVEL V=0                   
C                                                                       
      COMMON /APP3_9/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT                
      COMMON /BL3_9/ RSIBB(401)                                           
      DATA (EB(1,J),J=1,9)/-26.1308,-24.9681,-22.6532,-19.2080,
     1  -14.669,-9.0915,-2.564,4.7286,12.4947/
      DATA (EB(2,J),J=1,9)/-0.516,-0.1246,0.,0.,0.,0.,0.,0.,0./         
      DATA TWOPIC/1.88365183D11/                                        
      NSRI=190                                                          
      WNRMAX3=25.D0                                                     
C                                                                       
C     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN      
C     WNRMAX3 - THE FREQUENCY RANGE OF B-B CONTRIBUTION                 
C     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH       
C     B-B SPECTRUM, IN [CM-1].                                          
C                                                                       
      NSRIUP=2*NSRI+1                                                   
      DX=WNRMAX3/DFLOAT(NSRI)                                           
      DO 10 I=1,401                                                     
   10 RSIBB(I)=0.0                                                      
      NS=INT(SLIT/DX)                                                   
      DO 20 I=1,9                                                       
      DO 20 J=1,9                                                       
   20 A(I,J)=0.0                                                        
C                                                                       
C     A(I,J) FOR L=32 CONTRIBUTION H2-CH4 SYSTEM                        
C                                                                       
      A(1,4)=.14127D-39                                                 
      A(2,3)=.18414D-39                                                 
      A(2,5)=.23414D-39                                                 
      A(3,4)=.18534D-39                                                 
      A(3,6)=.30894D-39                                                 
      A(4,5)=.21745D-39                                                 
      A(4,7)=.36564D-39                                                 
      A(5,6)=.24659D-39                                                 
      A(5,8)=.40032D-39                                                 
      A(6,7)=.26656D-39                                                 
      A(7,8)=.27376D-39                                                 
      A(8,9)=.26805D-39                                                 
      A(6,9)=.41233D-39                                                 
      ALFA=1./(0.69519*TEMP)                                            
      RM=1.7768*1.672649D-24                                            
      PI=3.141592654                                                    
      PF=((RM*(1.380662D-16)*TEMP*2.*PI/(6.626176D-27)**2)**1.5)        
C                                                                       
C     DELTA V = 0 BELOW                                                 
C                                                                       
      DO 40 L=1,8                                                       
         STOKE1=EB(1,L+1)-EB(1,L)                                       
C                                                                       
C     FREQUENCY SHIFT FOR DELTA J=+1,DELTA V=0                          
C                                                                       
         IF (L.GT.6) GO TO 30                                           
         STOKE3=EB(1,L+3)-EB(1,L)                                       
C                                                                       
C     STOKE3- FREQUENCY SHIFT FOR DELTA J=+3, DELTA V=0                 
C                                                                       
         STOKI=A(L,L+3)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE_9 (STOKE3,STOKI)                                    
         STOKIP=A(L,L+3)*DEXP(-ALFA*EB(1,L+3))/PF                       
         CALL PROFILE_9 (-STOKE3,STOKIP)                                  
   30    STOKI=A(L,L+1)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE_9 (STOKE1,STOKI)                                    
         STOKIP=A(L,L+1)*DEXP(-ALFA*EB(1,L+1))/PF                       
         CALL PROFILE_9 (-STOKE1,STOKIP)                                  
   40 CONTINUE                                                          
C                                                                       
C     DELTA V=1 BELOW                                                   
C                                                                       
      AV1=.39933D-41                                                    
      STOKE=22.5287                                                     
      STOKI=AV1*DEXP(-ALFA*EB(1,3))/PF                                  
      CALL PROFILE_9 (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF                                 
      CALL PROFILE_9 (-STOKE,STOKIP)                                      
      AV1=0.3079D-41                                                    
      STOKE=18.6921                                                     
      STOKI=AV1*DEXP(-ALFA*EB(1,4))/PF                                  
      CALL PROFILE_9 (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,1))/PF                                 
      CALL PROFILE_9 (-STOKE,STOKIP)                                      
      AV1=0.4281D-41                                                    
      STOKE=14.544                                                      
      STOKI=AV1*DEXP(-ALFA*EB(1,5))/PF                                  
      CALL PROFILE_9 (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF                                 
      CALL PROFILE_9 (-STOKE,STOKIP)                                      
      DO 50 N=1,NSRIUP                                                  
   50 RSIBB(N)=RSIBB(N)/TWOPIC/SLIT                                     
      NSOL=NSRI+1                                                       
      DO 60 I=1,NSOL                                                    
      K=(NSRI+1)+(I-1)                                                  
   60 RSI(I)=RSIBB(K)                                                   
C                                                                       
C     RSI - BOUND-BOUND CONTRIBUTION FOR POSITIVE FREQUENCY SHIFTS      
C                                                                       
      RETURN                                                            
C                                                                       
      END                                                               
      SUBROUTINE BOUND43_9 (TEMP,RSI,NSOL)                                
       implicit double precision (A-H,O-Z)
      DIMENSION EB(2,9), A(9,9), RSI(201)                               
C                                                                       
C     EB(I,K) - BOUND ENERGIES FOR VIB.NR.V=I-1, ROT.NR J=K-1           
C     A(I,J)=( C(I,L,J)**2) * (MTX.EL. (I,BETA,J)**2 )                  
C     I,J-DENOTE ROTATIONAL STATES FOR VIB. LEVEL V=0                   
C                                                                       
      COMMON /APP3_9/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT                
      COMMON /BL3_9/ RSIBB(401)                                           
      COMMON /BOU43_9/ INITB                                              
      DATA (EB(1,J),J=1,9)/-26.1308,-24.9681,-22.6532,-19.2080,-14.669,-
     19.0915,-2.564,4.7286,12.4947/                                     
      DATA (EB(2,J),J=1,9)/-0.516,-0.1246,0.,0.,0.,0.,0.,0.,0./         
      DATA TWOPIC/1.88365183D11/                                        
      NSRI=190                                                          
      WNRMAX3=30.                                                       
      INITB=1                                                           
C                                                                       
C     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN      
C     WNRMAX3 - THE FREQUENCY RANGE OF B-B CONTRIBUTION                 
C     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH       
C     B-B SPECTRUM, IN [CM-1].                                          
C                                                                       
      NSRIUP=2*NSRI+1                                                   
      DX=WNRMAX3/DFLOAT(NSRI)                                           
      DO 10 I=1,401                                                     
   10 RSIBB(I)=0.0                                                      
      NS=INT(SLIT/DX)                                                   
      DO 20 I=1,9                                                       
      DO 20 J=1,9                                                       
   20 A(I,J)=0.0                                                        
C                                                                       
C     A(I,J) FOR L,Lambda={3,4} CONTRIBUTION H2-CH4 SYSTEM              
C     MATRIX ELEMENTS AS IF M(ll')/7                                    
C                                                                       
      A(1,5)=.33038D-41                                                 
      A(2,4)=.45246D-41                                                 
      A(2,6)=.52213D-41                                                 
      A(3,5)=.42125D-41                                                 
      A(3,7)=.65899D-41                                                 
      A(4,6)=.46727D-41                                                 
      A(4,8)=.74203D-41                                                 
      A(5,7)=.50413D-41                                                 
      A(5,9)=.77286D-41                                                 
      A(6,8)=.51657D-41                                                 
      A(7,9)=.50297D-41                                                 
      ALFA=1./(0.69519*TEMP)                                            
      RM=1.7768*1.672649D-24                                            
      PI=3.141592654                                                    
      PF=((RM*(1.380662D-16)*TEMP*2.*PI/(6.626176D-27)**2)**1.5)        
C                                                                       
C     DELTA V = 0 BELOW                                                 
C                                                                       
      DO 40 L=1,7                                                       
         STOKE1=EB(1,L+2)-EB(1,L)                                       
C                                                                       
C     FREQUENCY SHIFT FOR DELTA J=+2,DELTA V=0                          
C                                                                       
         IF (L.GT.5) GO TO 30                                           
         STOKE3=EB(1,L+4)-EB(1,L)                                       
C                                                                       
C     STOKE3- FREQUENCY SHIFT FOR DELTA J=+4, DELTA V=0                 
C                                                                       
         STOKI=A(L,L+4)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE_9 (STOKE3,STOKI)                                    
         STOKIP=A(L,L+4)*DEXP(-ALFA*EB(1,L+4))/PF                       
         CALL PROFILE_9 (-STOKE3,STOKIP)                                  
   30    STOKI=A(L,L+2)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE_9 (STOKE1,STOKI)                                    
         STOKIP=A(L,L+2)*DEXP(-ALFA*EB(1,L+2))/PF                       
         CALL PROFILE_9 (-STOKE1,STOKIP)                                  
   40 CONTINUE                                                          
C                                                                       
C     DELTA V=1 BELOW                                                   
C                                                                       
      AV1=.12153D-42                                                    
      STOKE=19.0835                                                     
      STOKI=AV1*DEXP(-ALFA*EB(1,4))/PF                                  
      CALL PROFILE_9 (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF                                 
      CALL PROFILE_9 (-STOKE,STOKIP)                                      
      AV1=0.89592D-43                                                   
      STOKE=14.1528                                                     
      STOKI=AV1*DEXP(-ALFA*EB(1,5))/PF                                  
      CALL PROFILE_9 (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,1))/PF                                 
      CALL PROFILE_9 (-STOKE,STOKIP)                                      
      AV1=0.1182D-42                                                    
      STOKE=8.9669                                                      
      STOKI=AV1*DEXP(-ALFA*EB(1,6))/PF                                  
      CALL PROFILE_9 (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF                                 
      CALL PROFILE_9 (-STOKE,STOKIP)                                      
      DO 50 N=1,NSRIUP                                                  
   50 RSIBB(N)=RSIBB(N)/TWOPIC/SLIT                                     
      NSOL=NSRI+1                                                       
      DO 60 I=1,NSOL                                                    
      K=(NSRI+1)+(I-1)                                                  
   60 RSI(I)=RSIBB(K)                                                   
C                                                                       
C     RSI - BOUND-BOUND CONTRIBUTION FOR POSITIVE FREQUENCY SHIFTS      
C                                                                       
      I=0                                                               
   70 I=I+1                                                             
      IF (RSI(I).EQ.0.) GO TO 70                                        
      IF (I.EQ.1) GO TO 90                                              
      NSOL=NSOL-I+1                                                     
      DO 80 J=1,NSOL                                                    
   80 RSI(J)=RSI(J+I-1)                                                 
      INITB=I                                                           
   90 CONTINUE                                                          
      RETURN                                                            
      END                                                               

      SUBROUTINE BOUN54C_9 (TEMP,RSI,NSOL)                                
       implicit double precision (A-H,O-Z)
      DIMENSION EB(2,9), A(9,9), RSI(201)                               
C                                                                       
C     EB(I,K) - BOUND ENERGIES FOR VIB.NR.V=I-1, ROT.NR J=K-1           
C     A(I,J)=( C(I,L,J)**2) * (MTX.EL. (I,BETA,J)**2 )                  
C     I,J-DENOTE ROTATIONAL STATES FOR VIB. LEVEL V=0                   
C                                                                       
      COMMON /APP3_9/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT                
      COMMON /BL3_9/ RSIBB(401)                                           
      DATA (EB(1,J),J=1,9)/-26.1308,-24.9681,-22.6532,-19.2080,-14.669,-
     19.0915,-2.564,4.7286,12.4947/                                     
      DATA (EB(2,J),J=1,9)/-0.516,-0.1246,0.,0.,0.,0.,0.,0.,0./         
      DATA TWOPIC/1.88365183D11/                                        
      NSRI=190                                                          
      WNRMAX3=34.                                                       
C                                                                       
C     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN      
C     WNRMAX3 - THE FREQUENCY RANGE OF B-B CONTRIBUTION                 
C     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH       
C     B-B SPECTRUM, IN [CM-1].                                          
C                                                                       
      NSRIUP=2*NSRI+1                                                   
      DX=WNRMAX3/DFLOAT(NSRI)                                           
      DO 10 I=1,401                                                     
   10 RSIBB(I)=0.0                                                      
      NS=INT(SLIT/DX)                                                   
      DO 20 I=1,9                                                       
      DO 20 J=1,9                                                       
   20 A(I,J)=0.0                                                        
C                                                                       
C     A(I,J) FOR L=54(CH4) CONTRIBUTION H2-CH4 SYSTEM                   
C     MATRIX ELEMENTS AS IF M(LL')/9                                    
C                                                                       
      A(1,6)=.99184D-42                                                 
      A(2,5)=.14168D-41                                                 
      A(2,7)=.14969D-41                                                 
      A(3,4)=.16093D-41                                                 
      A(3,6)=.12455D-41                                                 
      A(4,5)=.12587D-41                                                 
      A(4,7)=.13120D-41                                                 
      A(5,6)=.12714D-41                                                 
      A(5,8)=.13379D-41                                                 
      A(6,7)=.12818D-41                                                 
      A(7,8)=.12378D-41                                                 
      A(6,9)=.12963D-41                                                 
      A(8,9)=.11352D-41                                                 
      A(3,8)=.17927D-41                                                 
      A(4,9)=.19144D-41                                                 
      ALFA=1./(0.69519*TEMP)                                            
      RM=1.7768*1.672649D-24                                            
      PI=3.141592654                                                    
      PF=((RM*(1.380662D-16)*TEMP*2.*PI/(6.626176D-27)**2)**1.5)        
C                                                                       
C     DELTA V = 0 BELOW                                                 
C                                                                       
      DO 50 L=1,8                                                       
         STOKE1=EB(1,L+1)-EB(1,L)                                       
C                                                                       
C     FREQUENCY SHIFT FOR DELTA J=+1,DELTA V=0                          
C                                                                       
         IF (L.GT.6) GO TO 40                                           
         STOKE3=EB(1,L+3)-EB(1,L)                                       
C                                                                       
C     STOKE3- FREQUENCY SHIFT FOR DELTA J=+3, DELTA V=0                 
C                                                                       
         IF (L.GT.4) GO TO 30                                           
         STOKE5=EB(1,L+5)-EB(1,L)                                       
         STOKI=A(L,L+5)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE_9 (STOKE5,STOKI)                                    
         STOKIP=A(L,L+5)*DEXP(-ALFA*EB(1,L+5))/PF                       
         CALL PROFILE_9 (-STOKE5,STOKIP)                                  
   30    STOKI=A(L,L+3)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE_9 (STOKE3,STOKI)                                    
         STOKIP=A(L,L+3)*DEXP(-ALFA*EB(1,L+3))/PF                       
         CALL PROFILE_9 (-STOKE3,STOKIP)                                  
   40    STOKI=A(L,L+1)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE_9 (STOKE1,STOKI)                                    
         STOKIP=A(L,L+1)*DEXP(-ALFA*EB(1,L+1))/PF                       
         CALL PROFILE_9 (-STOKE1,STOKIP)                                  
   50 CONTINUE                                                          
C                                                                       
C     DELTA V=1 BELOW                                                   
C                                                                       
      AV1=0.45814D-43                                                   
      STOKE=14.5442                                                     
      STOKI=AV1*DEXP(-ALFA*EB(1,5))/PF                                  
      CALL PROFILE_9 (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF                                 
      CALL PROFILE_9 (-STOKE,STOKIP)                                      
      AV1=.32581D-43                                                    
      STOKE=8.5755                                                      
      STOKI=AV1*DEXP(-ALFA*EB(1,6))/PF                                  
      CALL PROFILE_9 (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,1))/PF                                 
      CALL PROFILE_9 (-STOKE,STOKIP)                                      
      AV1=.41005D-43                                                    
      STOKE=2.4395                                                      
      STOKI=AV1*DEXP(-ALFA*EB(1,7))/PF                                  
      CALL PROFILE_9 (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF                                 
      CALL PROFILE_9 (-STOKE,STOKIP)                                      
      DO 60 N=1,NSRIUP                                                  
   60 RSIBB(N)=RSIBB(N)/TWOPIC/SLIT                                     
      NSOL=NSRI+1                                                       
      DO 70 I=1,NSOL                                                    
      K=(NSRI+1)+(I-1)                                                  
   70 RSI(I)=RSIBB(K)                                                   
C                                                                       
C     RSI - BOUND-BOUND CONTRIBUTION FOR POSITIVE FREQUENCY SHIFTS      
C                                                                       
      RETURN                                                            
      END                                                               

      FUNCTION BGAUS_9 (FNU,G0,DELTA,OMEGA0,TEMP)                         
C                                                                       
C     THIS IS DESYMMETRIZED GAUSSIAN PROFILE, G0 IS A ZEROTH MOMENT
C     FNU,DELTA AND OMEGA0 ARE IN CM-1                                  
C     G0 IN CGS                                                         
C                                                                       
      implicit double precision (A-H,O-Z)
      DATA BKW,TWOPI/0.6950304256D0,6.2831853D0/                        
      D=G0/(DELTA*DSQRT(TWOPI))                                         
      DESYM=2./(1. + DEXP(-FNU/(BKW*TEMP)))                             
      FEXP=(FNU-OMEGA0)**2/(2.*DELTA**2)                                
      FEXP=DMIN1(FEXP,300.D0)                                           
      BGAUS_9=D*DESYM*DEXP(-FEXP)                                         
      RETURN                                                            
      END                                                               

c-----------------------------------------------------------------------
	subroutine setqw_9(qw3,qw4)
c-----------------------------------------------------------------------
c	set quantum corrections
c
c	numbers are from A. Borysow's file 'corrch4' dowloaded from the 
c	website on 23-06-04
c
c	a data statement cannot be used for commom block variables so
c	must be set individually via qw3=a, qw4=b
c
c	n.teanby	original version	23-06-04
c
c-----------------------------------------------------------------------
      implicit none
      double precision QW3(51,3),QW4(51,4)
      double precision a1(51),a2(51),a3(51),b1(51),b2(51),b3(51),b4(51)
      integer i
      data a1 /-1.00000D0,-1.00000D0,-1.00000D0,4.09091D0,-1.00000D0,
     >-1.00000D0,.79679D0,.27141D0,-.85139D0,.59819D0,
     >-.00412D0,-.31878D0,.18283D0,.15133D0,-.32477D0,
     >.21392D0,.02191D0,-.14717D0,.07524D0,.08569D0,
     >-.16823D0,.10771D0,.01794D0,-.08386D0,.04014D0,
     >.05426D0,-.10240D0,.06447D0,.01345D0,-.05399D0,
     >.02474D0,.03726D0,-.06876D0,.04281D0,.01019D0,
     >-.03761D0,.01670D0,.02711D0,-.04932D0,.03046D0,
     >.00792D0,-.02769D0,.01201D0,.02059D0,-.03708D0,
     >.02276D0,.00630D0,-.02122D0,.00904D0,.01616D0,
     >-.02889D0/
	data a2 /-1.00000D0,-1.00000D0,-1.00000D0,-1.00000D0,1.46154D0,
     >-1.00000D0,.52036D0,-.38246D0,.48307D0,-.58061D0,
     >.54775D0,-.42488D0,.32187D0,-.30004D0,.33377D0,
     >-.35746D0,.33367D0,-.27900D0,.23619D0,-.22949D0,
     >.24628D0,-.25572D0,.24034D0,-.20959D0,.18643D0,
     >-.18371D0,.19376D0,-.19855D0,.18803D0,-.16835D0,
     >.15390D0,-.15262D0,.15930D0,-.16211D0,.15451D0,
     >-.14085D0,.13100D0,-.13033D0,.13510D0,-.13691D0,
     >.13118D0,-.12115D0,.11400D0,-.11364D0,.11722D0,
     >-.11846D0,.11400D0,-.10631D0,.10090D0,-.10070D0,
     >.10348D0/
      data a3 /11.00000D0,-1.00000D0,-1.00000D0,1.01399D0,.09091D0,
     >-1.00000D0,.94475D0,-.16780D0,-.51321D0,.52338D0,
     >.01421D0,-.50745D0,.49964D0,-.06185D0,-.34591D0,
     >.34888D0,.00564D0,-.34073D0,.33830D0,-.03176D0,
     >-.26007D0,.26130D0,.00304D0,-.25650D0,.25545D0,
     >-.01926D0,-.20817D0,.20879D0,.00191D0,-.20566D0,
     >.20511D0,-.01291D0,-.17348D0,.17383D0,.00131D0,
     >-.17164D0,.17132D0,-.00925D0,-.14867D0,.14889D0,
     >.00095D0,-.14727D0,.14707D0,-.00695D0,-.13006D0,
     >.13021D0,.00073D0,-.12896D0,.12882D0,-.00541D0,
     >-.11559D0/
      data b1 /-1.00000D0,-1.00000D0,-1.00000D0,-1.00000D0,-1.00000D0,
     >-1.00000D0,1.82353D0,-1.00000D0,-.00478D0,.22586D0,
     >.19088D0,-.58488D0,.55212D0,-.21168D0,-.06443D0,
     >.05729D0,.14181D0,-.28953D0,.24581D0,-.07845D0,
     >-.04207D0,.01948D0,.09314D0,-.16876D0,.13649D0,
     >-.03831D0,-.02777D0,.00772D0,.06418D0,-.10976D0,
     >.08625D0,-.02197D0,-.01938D0,.00323D0,.04653D0,
     >-.07688D0,.05925D0,-.01398D0,-.01420D0,.00128D0,
     >.03515D0,-.05678D0,.04315D0,-.00956D0,-.01082D0,
     >.00037D0,.02744D0,-.04362D0,.03279D0,-.00690D0,
     >-.00850D0/
      data b2 /-1.00000D0,-1.00000D0,-1.00000D0,-1.00000D0,2.35664D0,
     >-1.00000D0,-.43987D0,.63467D0,.00794D0,-.53268D0,
     >.47071D0,-.04638D0,-.25892D0,.20786D0,.05459D0,
     >-.23119D0,.17454D0,.01650D0,-.14054D0,.09758D0,
     >.03999D0,-.12597D0,.08757D0,.02020D0,-.08639D0,
     >.05576D0,.02883D0,-.07869D0,.05185D0,.01719D0,
     >-.05811D0,.03585D0,.02073D0,-.05366D0,.03404D0,
     >.01393D0,-.04166D0,.02491D0,.01572D0,-.03888D0,
     >.02397D0,.01128D0,-.03129D0,.01828D0,.01229D0,
     >-.02945D0,.01775D0,.00924D0,-.02434D0,.01397D0,
     >.00985D0/
      data b3 /-1.00000D0,-1.00000D0,-1.00000D0,2.16484D0,-1.00000D0,
     >-1.00000D0,1.60491D0,-.46218D0,-1.00000D0,1.31017D0,
     >-.34627D0,-.78497D0,.98223D0,-.18287D0,-.72146D0,
     >.84929D0,-.16187D0,-.59662D0,.69125D0,-.09565D0,
     >-.55198D0,.62115D0,-.09262D0,-.47526D0,.53051D0,
     >-.05845D0,-.44470D0,.48796D0,-.05974D0,-.39348D0,
     >.42963D0,-.03932D0,-.37166D0,.40124D0,-.04166D0,
     >-.33521D0,.36069D0,-.02824D0,-.31896D0,.34046D0,
     >-.03069D0,-.29177D0,.31068D0,-.02125D0,-.27924D0,
     >.29557D0,-.02353D0,-.25819D0,.27279D0,-.01656D0,
     >-.024825D0/
      data b4 /11.00000D0,-1.00000D0,-1.00000D0,.40260D0,.62896D0,
     >-1.00000D0,.67713D0,-.17620D0,-.04947D0,-.07182D0,
     >.30704D0,-.40281D0,.29894D0,-.12495D0,.04080D0,
     >-.09160D0,.19422D0,-.23764D0,.18677D0,-.09887D0,
     >.05507D0,-.08289D0,.14023D0,-.16496D0,.13482D0,
     >-.08186D0,.05504D0,-.07257D0,.10916D0,-.12512D0,
     >.10520D0,-.06982D0,.05172D0,-.06377D0,.08914D0,
     >-.10029D0,.08615D0,-.06085D0,.04781D0,-.05661D0,
     >.07522D0,-.08345D0,.07290D0,-.05391D0,.04406D0,
     >-.05077D0,.06501D0,-.07133D0,.06315D0,-.04838D0,
     >.04069D0/
      
c	assign qw3 and qw4
	do i=1,51
        qw3(i,1)=a1(i)
        qw3(i,2)=a2(i)
        qw3(i,3)=a3(i)
        qw4(i,1)=b1(i)
        qw4(i,2)=b2(i)
        qw4(i,3)=b3(i)
        qw4(i,4)=b4(i)
      enddo
      
      return
      end
