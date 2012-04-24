      subroutine n2h2_s(temp,fnumin,fnumax,dnu,nf1,freq,alfatot,slit1)
C       ****************************************************************
C     subroutine adapted from Borysow code to calculate CIA absorption
C     of N2-H2 collisions (for Titan atmosphere).
C
C     Input variables:
C	temp		double	Temperature (K)
C	fnumin		double	Lowest wavenumber in spectrum (cm-1)
C	fnumax		double	Highest wavenumber in spectrum (cm-1)
C	dnu		double	Wavenumber step (cm-1)
C       slit            double  Slit size
C
C     Output variables
C	nf		integer	Number of points in spectrum
C	freq(601)	double	Frequency grid (cm-1)
C	alfatot(601)	double	Absorption coefficient (cm-1 amagat-2)
C
C     C. Nixon 	22-1-97
C     Pat Irwin	2/3/12	Updated for Radtrans2.0
C
C       ****************************************************************
C       PROGRAM PREPARED BY ALEKSANDRA BORYSOW; UNIVERSITY OF TEXAS AT A
C       AND JOINT INSTITUTE FOR LABORATORY ASTROPHYSICS, UNIVERSITY OF C
C       LAST revision DATE: 14 NOVEMBER 1988                            
C       PROGRAM COMPATIBLE WITH PAPER: A. BORYSOW & L. FROMMHOLD;       
C       ASTROPHYSICAL JOURNAL; VOL. 303, PP. 495-510, (1986).           
C       ****************************************************************
C     PROGRAM GENERATES THE H2-N2, FREE-FREE, BOUND-FREE & B-B CONTRIBUT
C     TO THE CIA SPECTRA                                                

      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      COMMON /APP3/ SLIT,DX,WNRMAX,NSRI,NS,NSRIUP,JSLIT  
      COMMON /RSILO/ RSILO(201)                                         
      COMMON /BB/ OMEG,RSI,RSIGG,ALFA,SCAL,NSOL                         
      DIMENSION FREQ(601), ABSCOEF(601), ALFATOT(601)                   
      DIMENSION RSI(201), RSIGG(201), TT(2), SS(1), OMEG(201), AIG(201) 
      Y(X,A,B,C)=A*DEXP((C*X+B)*X)                                      

      slit = slit1
C      SLIT=4.3
c	The parameters taken in subroutines BOUNDxx are specially chosen
c	for SLIT=4.3 cm-1; if a much smaller slit is chosen, parameters 
c	WNRMAX, DX and NRSI have to be adjusted
C     ( the slit width is now a parameter. The WNRMAX etc is now 
C      calculated for each SLIT value. 
C      see BOUND32 and BOUND54 subroutines - C.A.N. )

      if (slit .lt. 0.50) then
         print*, 'n2h2_s: Warning'
         print*, 'Slit width approaches minimum grid scale of'
         print*, 'integration: 0.12 cm-1 for bound54 transitions.'
      end if

c	NOTE: Take care that all b-b intensities are larger than zero
c	in case of separated lines, (extremelty small slit-width) put 
c	RSI(i)=1.d-99 or so, to prevent a failure in LOG(RSI).

C     now allows for variable slit width, and fixed the LOG(RSI) to
C     check for log(0). Slit width is passed as a parameter. C.A.N.

      NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1   
C save the value of nf
      nf1 = nf
      IF (NF.GT.601) NF=601      
      FNUMAX=FNUMIN+DFLOAT(NF-1)*DNU

      if (temp .lt. 45 .or. temp .gt. 300) then
         print*, 'n2h2_s: Warning'
         print*, 'Temperature should be 45 < T < 300 K'
      end if

      CALL PARTSUM5 (TEMP)                                               
C                                                                       
C     THE H2-N2 SPECTRA   FOR 50-300K                                   
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
      EPS=1.D-5                                                         
      TT(1)=10.                                                         

      CALL BOUND32 (TEMP,RSI,NSOL)                                      

      DO 20 I=1,NSOL

C this if-clause added by C.Nixon 28-1-97 to account for the fact that
C since I included a variable slit width (instead of fixed 4.3), rsi
C can now take zero values, and the log of a 0 is an error!
C         print*, rsi(i)*1.d0
        if (rsi(i) .eq. 0.0) rsi(i) = 1.e-100

        RSILO(I)=DLOG(RSI(I)*1.D80)

        OMEG(I)=DFLOAT(I-1)*DX
 20   continue
      CALL SPLINE5 (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,RSIGG,NR)           
C                                                                       
C     THE LAMBDA1,LAMBDA2,LAMBDA,L = 2023 AND 0223 COMPONENTS:          
C     2023 - CORRESPONDS TO HYDROGEN'S SINGLE TRANSITIONS               
C     0223 - CORRESPONDS TO NITROGEN'S TRANSITIONS                      
C     HYDROGEN'S QUADRUPOLE                                             
C                                                                       
      SCAL=33.74d0                                                      
      IF (TEMP-145.d0) 30,30,40                                         
   30 S=Y(X,33.74*5.9340D-62,-1.5429D0,.15677D0)                        
      E=-Y(X,4.239838D-3,2.76328D0,-.38222D0)                           
      T1=Y(X,8.66209D-14,.59129D0,-.11062D0)                            
      T2=Y(X,4.576639D-13,-.20236D0,-.03646D0)                          
      T3=Y(X,1.564106D-13,.42334D0,-.10565D0)                           
      T4=Y(X,1.37791D-11,-1.61496D0,.19558D0)                           
      GO TO 50                                                          
   40 S=Y(X,33.74*5.5440D-63,-0.58016D0,.05844D0)                       
      T1=Y(X,4.4024D-12,-.92834D0,.03489D0)                             
      T2=Y(X,2.0536D-13,-.0549D0,-.03709D0)                             
      E=-Y(X,3.79637D-6,4.41619D0,-.43833D0)                            
      T3=Y(X,1.1773D-13,.22391D0,-.06086D0)                             
      T4=Y(X,1.28155D-4,-6.30489D0,.47603D0)                            
   50 CONTINUE                                                          

      CALL ADDSP5C (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,0,2,0,2,3)   
      DO 60 I=1,NF                                                      
   60 ALFATOT(I)=ABSCOEF(I)                                             
C                                                                       
C     NITROGEN'S QUADRUPOLE   (0,4,4,5)                                 
C                                                                       
      SCAL=34.9472D0                                                    
      IF (TEMP-145.) 70,70,80                                           
   70 S=Y(X,34.9472*5.9340D-62,-1.5429D0,.15677D0)                      
      GO TO 90                                                          
   80 S=Y(X,34.9472*5.5440D-63,-0.58016D0,.05844D0)                     
   90 CONTINUE                                                          
      CALL ADDSP5C (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,0,0,2,2,3)   
      DO 100 I=1,NF                                                     
  100 ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)                                  
C                                                                       
C     DOUBLE  TRANSITIONS (2,2,3,3)                                     
C                                                                       
      SCAL=3.3810D0                                                     
      IF (TEMP-145.D0) 110,110,120                                      
  110 S=Y(X,3.3810*5.9340D-62,-1.5429D0,.15677D0)                       
      GO TO 130                                                         
  120 S=Y(X,3.3810*5.5440D-63,-0.58016D0,.05844D0)                      
  130 CONTINUE                                                          
      CALL ADDSP5C (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,0,2,2,3,3)   
      DO 140 I=1,NF                                                     
  140 ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)                                  
      CALL BOUND54 (TEMP,RSI,NSOL)                                      
      DO 150 I=1,NSOL                                                   

C once again, as for the bound32 transition, have to allow for the 
C fact that rsi can be 0. If-clause added by C.Nixon 28-1-97

         if (rsi(i) .eq. 0.0) rsi(i) = 1.e-100
         RSILO(I)=DLOG(RSI(I)*1.D80)
         
         OMEG(I)=DFLOAT(I-1)*DX                                            
 150  continue

      CALL SPLINE5 (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,RSIGG,NR)           
C                                                                       
C     PARAMETERS FOR 4045  (HYDROGEN HEXADECAPOLE) COMPONENTS           
C                                                                       
      S=Y(X,17.7052*9.7163D-65,-1.8937D0,2.01078D-1)                    
      E=-Y(X,9.8311D-7,4.9861D0,-4.61612D-1)                            
      T1=Y(X,2.5235D-15,1.6788D0,-1.86225D-1)                           
      T2=Y(X,3.27189D-13,-0.2089D0,-3.3745D-2)                          
      T3=Y(X,8.5546D-9,-3.8389D0,3.3309D-1)                             
      T4=Y(X,1.33844D-13,0.3935D0,-7.49592D-2)                          
      SCAL=17.7052d0                                                    
      CALL ADDSP5C (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,0,4,0,4,5)   
      DO 160 I=1,NF                                                     
  160 ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)                                  
C                                                                       
C     PARAMETERS FOR 0445 - NITROGEN'S HEXADECAPOLE                     
C     ALL AS FOR 4045 EXCEPT PARAMETER S                                
C                                                                       
      S=Y(X,3181.456*9.7163D-65,-1.8937D0,2.01078D-1)                   
      SCAL=3181.456D0                                                   
      CALL ADDSP5C (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,0,0,4,4,5)   
      DO 170 I=1,NF                                                     
  170 ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)                                  
          
       return
       end

C ************************ END OF MAIN PROGRAM **************************
C ***************** beginning of subroutines and functions **************



      SUBROUTINE ADDSP5C (G0,EP,TAU1,TAU2,TAU5,TAU6,TEMP,NF,FREQ,ABSCOEF
     1,MP,LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE)                           
C                                                                       
C     THIS PROGRAM GENERATES LISTING OF CIA TR ALFA(OMEGA)              
C     IF EITHER LAMBDA1 OR LAMBDA2 EQUAL TO ZERO - SINGLE TRANSITIONS;  
C     LAMBDA1 CORRESPONDS TO H2, LAMBDA2 TO N2                          
C     DOUBLE TRANSITIONS ARE ASSUMED OTHERWISE.                         
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.              
C                                                                       
      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      LOGICAL COND                                                      
      COMMON /BB/ OMEG(201),RSI(201),RSIGG(201),BETA,SCAL,NSOL          
      COMMON /RSILO/ RSILO(201)                                         
      COMMON /APP3/ SLIT,DX,WNRMAX,NSRI,NS,NSRIUP,JSLIT                 
      COMMON /H2PARTC/ Q,WH2(2),B0,D0,JRANGE1                            
      COMMON /N2PARTB/ Q1,WN2(2),B01,D01,JRANGE2                         
      DIMENSION ABSCOEF(601), FREQ(601)                                 
      DATA CLOSCHM,BOLTZWN/2.68675484D19,.6950304/                      
      DATA HBAR,PI,CLIGHT/1.054588757D-27,3.1415926535898,2.9979250D10/ 
      EH2(I)=(B0-DFLOAT(I)*D0)*DFLOAT(I)                                
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*DEXP(-1.4387859/T*
     1  EH2(J*(J+1)))                                                         
      EN2(I)=(B01-DFLOAT(I)*D01)*DFLOAT(I)                              
      PN2(J,T)=DFLOAT(2*J+1)*WN2(1+MOD(J,2))*DEXP(-1.4387859/T*
     1  EN2(J*(J+1)))                                                         


      TWOPIC=2.*PI*CLIGHT                                               
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2             
      CALIB=CALIB/DFLOAT(1+LIKE)                                        
      BETA=1./(BOLTZWN*TEMP)                                            
      LIST=NF                                                           
      DO 10 I=1,LIST                                                    
   10 ABSCOEF(I)=0.0                                                    
C                                                                       
C     ROTATIONAL SPECTRUM; DETAILED LISTING   *******************       
C                       

      IF ((LAMBDA1.EQ.0).OR.(LAMBDA2.EQ.0)) GO TO 70                    
      JPLUSL=JRANGE1+MAX0(LAMBDA1,LAMBDA2)                              
      JPLU2=JRANGE2+MAX0(LAMBDA1,LAMBDA2)                               

      DO 60 I1=1,JRANGE1                                                
         J1=I1-1                                                        
      DO 60 IP1=1,JPLUSL                                                
         JP1=IP1-1                                                      
         CG1S=CLEBSQR5(J1,LAMBDA1,JP1)                                   
         IF (CG1S) 60,60,20                                             
   20    P1=PH2(J1,TEMP)/Q                                              
         OMEGA1=EH2(JP1*IP1)-EH2(J1*I1)                                 
         DO 50 I2=1,JRANGE2                                             
            J2=I2-1                                                     
         DO 50 IP2=1,JPLU2                                              
            JP2=IP2-1                                                   
            CG2S=CLEBSQR5(J2,LAMBDA2,JP2)                                
            IF (CG2S) 50,50,30                                          
   30       P2=PN2(J2,TEMP)/Q1                                          
            OMEGA2=EN2(JP2*IP2)-EN2(J2*I2)                              
            FAC=CALIB*P1*P2*CG1S*CG2S                                   

            DO 40 I=1,LIST                                              
               FRQ=FREQ(I)-OMEGA1-OMEGA2                                
               WKI=FREQ(I)*(1.-DEXP(-BETA*FREQ(I)))                     
               WKF=WKI*FAC                                              
               XBG=G0*BGAMA5(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)            

               COND=ABS(FRQ).LE.WNRMAX                                  
               IF (COND) XBG=XBG+SCAL*SPECFCT(FRQ,OMEG,RSILO,RSIGG,NSOL,
     1         BETA)                                                    

               ABSCOEF(I)=ABSCOEF(I)+XBG*WKF                            
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60 CONTINUE                                                          
      GO TO 150                                                         
   70 IF (LAMBDA1.EQ.0) GO TO 110                                       
C                                                                       
C     SINGLE TRANSITIONS ON HYDROGEN'S ROTATIONAL FREQUENCIES.          
C                                                                       
      JPLUSL=JRANGE1+LAMBDA                                             
      DO 100 I=1,JRANGE1                                                
         J=I-1                                                          
      DO 100 IP=1,JPLUSL                                                
         JP=IP-1                                                        
         CGS=CLEBSQR5(J,LAMBDA,JP)                                       


         IF (CGS) 100,100,80                                            
   80    P=PH2(J,TEMP)/Q                                                
         OMEGA1=EH2(JP*IP)-EH2(J*I)                                     
         FAC=CALIB*P*CGS                                                
         DO 90 IQ=1,LIST                                                
            FRQ=FREQ(IQ)-OMEGA1                                         
            WKI=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))                      
            WKF=WKI*FAC                                                 
            XBG=G0*BGAMA5(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)               
            COND=ABS(FRQ).LE.WNRMAX                                     
            IF (COND) XBG=XBG+SCAL*SPECFCT(FRQ,OMEG,RSILO,RSIGG,NSOL,BET
     1      A)                                                          

            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF                             
   90    CONTINUE                                                       
  100 CONTINUE                                                          
      GO TO 150                                                         
C                                                                       
C     SINGLE TRANSITIONS ON NITROGEN'S ROTATIONAL FREQUENCIES.          
C                                                                       
  110 JPLU2=JRANGE2+LAMBDA                                              
      DO 140 I=1,JRANGE2                                                
         J=I-1                                                          
      DO 140 IP=1,JPLU2                                                 
         JP=IP-1                                                        
         CGS=CLEBSQR5(J,LAMBDA,JP)                                       
         IF (CGS) 140,140,120                                           
  120    P=PN2(J,TEMP)/Q1                                               
         OMEGA1=EN2(JP*IP)-EN2(J*I)                                     
         FAC=CALIB*P*CGS                                                
         DO 130 IQ=1,LIST                                               
            FRQ=FREQ(IQ)-OMEGA1                                         
            WKI=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))                      
            WKF=WKI*FAC                                                 
            XBG=G0*BGAMA5(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)               
            COND=ABS(FRQ).LE.WNRMAX                                     
            IF (COND) XBG=XBG+SCAL*SPECFCT(FRQ,OMEG,RSILO,RSIGG,NSOL,BET
     1      A)                                                          
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF                             
  130    CONTINUE                                                       
  140 CONTINUE                                                          
  150 CONTINUE                  


      RETURN                                                            

C                                                                       
  163 FORMAT (/, 32H LAMBDA1,LAMBDA2, LAMBDA,LVALUE=,2I3,2X,2I3,
     &12H COMPONENT .,/15X, 22HLINE SHAPE PARAMETERS:,6E12.3,5X,
     &  5HG(0)=,E12.3/)
  170 FORMAT ((1X,10E12.4,/))                                           
C                                                                       
      END                                                               

      FUNCTION CLEBSQR5 (L,LAMBDA,LP)                                    
C                                                                       
C     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)          
C     FOR INTEGER ARGUMENTS ONLY                                        
C     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.                
C                                                                       
      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      FC=DFLOAT(2*LP+1)                                                 
      GO TO 10                                                          
C                                                                       
      ENTRY THREEJ25                                                     
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
   10 CLEBSQR5=0.                                                        
      IF (((L+LAMBDA).LT.LP).OR.((LAMBDA+LP).LT.L).OR.((L+LP).LT.LAMBDA)
     1) RETURN                                                          
      IF (MOD(L+LP+LAMBDA,2).NE.0) RETURN                               
      IF ((L.LT.0).OR.(LP.LT.0).OR.(LAMBDA.LT.0)) RETURN                
      F=1./DFLOAT(L+LP+1-LAMBDA)                                        
      IF (LAMBDA.EQ.0) GO TO 30                                         
      I1=(L+LP+LAMBDA)/2                                                
      I0=(L+LP-LAMBDA)/2+1                                              
      DO 20 I=I0,I1                                                     
   20 F=F*DFLOAT(I)/DFLOAT(2*(2*I+1))                                   
   30 P=FC*F*FCTL5(LAMBDA+L-LP)*FCTL5(LAMBDA+LP-L)                        
      CLEBSQR5=P/(FCTL5((LAMBDA+L-LP)/2)*FCTL5((LAMBDA+LP-L)/2))**2        
      RETURN                                                            
C                                                                       
      END                                                               
      FUNCTION FCTL5 (N)                                                 
      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      P(Z)=((((-2.294720936D-4)/Z-(2.681327160D-3))/Z+(3.472222222D-3))/
     1Z+(8.333333333D-2))/Z+1.                                          
      FCTL5=1.                                                           
      IF (N.LE.1) RETURN                                                
      IF (N.GT.15) GO TO 20                                             
      J=1                                                               
      DO 10 I=2,N                                                       
   10 J=J*I                                                             
      FCTL5=DFLOAT(J)                                                    
      RETURN                                                            
   20 Z=DFLOAT(N+1)                                                     
      FCTL5=(DEXP(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)                 
      RETURN                                                            
C                                                                       
      END                                                               
      FUNCTION BGAMA5 (FNU,T1,T2,EPS,T3,T4,TEMP)                         
C                                                                       
C     EQUATION 13, SO-CALLED EBC MODEL, OF BORYSOW,TRAFTON,FROMMHOLD,   
C     AND BIRNBAUM, ASTROPHYS.J., TO BE PUBLISHED (1985)                
C                                                                       
      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
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
   20 XK1=dSQRT(Z)*DEXP(-Z)*P6(2./Z)                                    
   30 ZP=DSQRT((1.+(OMEGA*T4)**2)*(T3*T3+T0*T0))/T4                     
      IF (ZP-2.) 40,40,50                                               
   40 K0=-DLOG(ZP/2.)*P1((ZP/3.75)**2)+P2((ZP/2.)**2)                   
      GO TO 60                                                          
   50 K0=DEXP(-ZP)*P5(2./ZP)/DSQRT(ZP)                                  
   60 BGAMA5=((T1/PI)*DEXP(T2/T1+T0*OMEGA)*XK1/(1.+(T1*OMEGA)**2)+       
     1  EPS*(T3/PI)*DEXP(T3/T4+T0*OMEGA)*K0)/(1.+EPS)                   
      RETURN                                                            
C                                                                       
      END                                                               
      SUBROUTINE PARTSUM5 (TEMP)                                         
C                                                                       
C     H2 ROTATIONAL PARTITION SUM Q = Q(T).                             
C                                                                       
      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      COMMON /H2PARTC/ Q,WH2(2),B0,D0,JRANGE1                            
      COMMON /N2PARTB/ Q1,WN2(2),B01,D01,JRANGE2                         
      DATA B0,D0,WH2(1),WH2(2)/59.3392,0.04599,1.,3./                   
      DATA B01,D01,WN2(1),WN2(2)/1.98957,.58D-5,2.,1./                  
      EH2(I)=(B0-DFLOAT(I)*D0)*FLOAT(I)                                 
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*DEXP(-1.4387859*
     1  EH2(J*(J+1))/ T)                                                
      EN2(I)=(B01-DFLOAT(I)*D01)*FLOAT(I)                               
      PN2(J,T)=DFLOAT(2*J+1)*WN2(1+MOD(J,2))*DEXP(-1.4387859*
     1  EN2(J*(J+1))/T)                                                 
                                                                        
C                                                                       
C     Q,B01,D01,WN2 - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR N2     
C     Q,B0,D0,WH2 - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR H2       
C                                                                       
      Q=0.                                                              
      J=0                                                               
   10 DQ=PH2(J,TEMP)                                                    
      Q=Q+DQ                                                            
      J=J+1                                                             
      IF (DQ.GT.Q/900.) GO TO 10                                        
      JRANGE1=J                                                         
      Q1=0.                                                             
      J=0                                                               
   20 DQ1=PN2(J,TEMP)                                                   
      Q1=Q1+DQ1                                                         
      J=J+1                                                             
      IF (DQ1.GT.Q1/900.) GO TO 20                                      
      JRANGE2=J                                                         
      RETURN                                                            
C                                                                       
   30 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF H2: Q=,F8.2,10X,  
     17HJ MAX =,I3/)                                                    
   40 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF N2: Q=,F8.2,10X,  
     17HJ MAX =,I3/)                                                    
C                                                                       
      END                                                               
      SUBROUTINE PROFILE (X,Y)                                          
      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      COMMON /BL3/ RSI(401)                                             
C                                                                       
C     ATRIANGULAR SLIT FUNCTION IS USED.                                
C                                                                       
      COMMON /APP3/ SLIT,DX,WNRMAX,NSRI,NS,NSRIUP,JSLIT                 
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
         XI=(I-1.)*DX-WNRMAX                                            
         DR=SLOPE*(XI-(X-SLIT))                                         
         IF (DR.LE.0.) GO TO 20                                         
         RSI(I)=RSI(I)+DR                                               
   20 CONTINUE                                                          
   30 IF (NC.GE.NSRIUP) RETURN                                          
      IF (N1.LT.1) N1=1                                                 
      DO 40 I=N1,NO                                                     
         XI=(I-1.)*DX-WNRMAX                                            
         DR=Y-SLOPE*(XI-X)                                              
         IF (DR.LE.0.) GO TO 40                                         
         RSI(I)=RSI(I)+DR                                               
   40 CONTINUE                                                          
      RETURN                                                            
   50 continue
C      WRITE(12,70) SLIT                                                 
   60 CONTINUE                                                          
      RETURN                                                            
C                                                                       
   70 FORMAT (/, 30H A TRIANGULAR SLIT FUNCTION OF,F6.3, 23H CM-1 HALFWI
     1DTH IS USED,/)                                                    
C                                                                       
      END                                                               

      FUNCTION SPECFCT (FREQ,OMEGA,PHI,PHI2,N,RTEMP)                    
C                                                                       
C     THIS INTERPOLATES THE SPECTRAL FUNCTION PHI(FREQ) DEFINED AT      
C     OMEGA(N) AS PHI(N). PHI2 IS THE SECOND DERIVATIVE AT OMEGA        
C     WHICH MUST BE OBTAINED FIRST (USE SPLINE5 FOR THAT PURPOSE).       
C     RTEMP IS THE RECIPROCAL TEMPERATURE IN CM-1 UNITS.                
C     NOTE THAT WE INTERPOLATE 1.E80 TIMES THE LOGARITHM OF PHI(OMEGA)  
C                                                                       
      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      DIMENSION PHI(N), PHI2(N), OMEGA(N)                               
      dimension f(1), gp(1)
      DATA SCALEF/1.D-80/
      DELY(I)=(PHI(I+1)-PHI(I))/(OMEGA(I+1)-OMEGA(I))
c
      TFAC=0.                                                           
      F(1)=FREQ                                                            
      IF (F(1) ) 10,20,20                                                   
   10 F(1)=dABS(F(1))                                                         
      TFAC=(-RTEMP*F(1))                                                   
   20 IF (F(1).LE.OMEGA(N)) GO TO 30                                       
      SPECFCT=DEXP(-(PHI(N-1)-PHI(N))*(F(1)-OMEGA(N))/
     1  (OMEGA(N)-OMEGA(N-1))+                                          
     1PHI(N)+TFAC)*(1.D-80)                                             
      RETURN                                                            
C			INTERPOLATION IS BEING MADE BELOW:
   30 I = 2
   90 IF (f(1)-OMEGA(I))200,200,100
  100 I = I+1
      GOTO 90
  200 I = I-1
      HT1=f(1)-OMEGA(I)
      HT2=f(1)-OMEGA(I+1)
      PROD=HT1*HT2
      S3=(PHI2(I+1)-PHI2(I))/(OMEGA(I+1)-OMEGA(I))
      SS2=PHI2(I)+HT1*S3
      DELSQS=(PHI2(I)+PHI2(I+1)+SS2)/6.
      gp(1)=PHI(I)+HT1*DELY(I)+PROD*DELSQS
      SPECFCT=DEXP( tfac + gp(1) ) * SCALEF
      RETURN                                                            
      END                                                               

      SUBROUTINE SPLINE5 (L,M,K,EPS,X,Y,T,SS,SI,S2,NR)                   

C                                                                       
C     SPLINE5 INTERPOLATION AND QUADRATURE, THIRD ORDER AFTER GREVILLE.  
C     INPUT ARGUMENTS L...Y, OUTPUT SS...NR.                            
C     L DATA POINTS X(1), Y(1) ... X(L),Y(L)                            
C     EPS=ERROR CRITERION, TYPICALLY EPS=1.E-5 FOR 5 DECI. PLACES ACCURA
C     M ARGUMENTS T(1)..T(M) FOR WHICH FUNCTION VALUES SS(1)..SS(M), FOR
C     K=0; OR FIRST OR SECOND DERIVATIVE FOR K=1 OR -1, RESPECTIVELY.   
C     NOTE THAT M HAS TO BE AT LEAST EQUAL TO 1.                        
C     SI=INTEGRAL (OVER WHOLE INTERVAL) FOR K=2 ONLY.                   
C     FOR 'NATURAL' SPLINE5 FUNCTIONS, S2(1)=S2(L)=0. MUST BE INPUT*NOTE*
C     N0 INDICATES THE NUMBER OF OUT-OF-RANGE CALLS. X(1)@T(I)@X(L)     
C     EXTRAPOLATE WITH CAUTION. (ASSUMPTION D2Y/DX2 = 0.)               
C     S2(I) IS THE 2ND DERIVATIVE AT X=X(I) AND IS COMPUTED INTERNALLY. 
C                                                                       
      IMPLICIT  DOUBLE  PRECISION (A-H,O-Z)
      INTEGER L,M,K,NR
      DIMENSION X(201),Y(201),T(2),SS(1),S2(201)
      DOUBLE PRECISION EPS,SI
      DIMENSION DELY(201),B(201),C(201)

      N=L
      N1=N-1
      DO I=1,N1
       DELY(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))                              
      END DO
      DO I=2,N1
       B(I)=(X(I)-X(I-1))*0.5/(X(I+1)-X(I-1))                           
       C(I)=3.*(DELY(I)-DELY(I-1))/(X(I+1)-X(I-1))                      
      END DO

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
      SM=ABS(S2(1))                                                     
      DO 30 I=2,N                                                       
         IF (ABS(S2(I)).GT.SM) SM=ABS(S2(I))                            
   30 CONTINUE                                                          
      EPSI=EPS*SM                                                       
      DO 50 I=2,N1                                                      
         W=(C(I)-B(I)*S2(I-1)-(0.5-B(I))*S2(I+1)-S2(I))*OMEGA           
         IF (ABS(W)-ETA) 50,50,40                                       
   40    ETA=ABS(W)                                                     
   50 S2(I)=S2(I)+W                                                     
      IF (ETA-EPSI) 60,20,20                                            
      ENTRY IXPOLAT                                                     
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
C                                                                       
      END                                                               
      SUBROUTINE BOUND32 (TEMP,RSI,NSOL)                                
      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      DIMENSION EB(2,8), A(8,8), RSI(201)                               
C                                                                       
C     EB(I,K) - BOUND ENERGIES FOR VIB.NR.V=I-1, ROT.NR J=K-1           
C     A(I,J)=( C(I,L,J)**2) * (MTX.EL. (I,BETA,J)**2 )                  
C     I,J-DENOTE ROTATIONAL STATES FOR VIB. LEVEL V=0                   
C       SEE TABLE 10 OF AP. J. VOL.303, P.508 (1986)                    
C                                                                       
      COMMON /APP3/ SLIT,DX,WNRMAX ,NSRI,NS,NSRIUP,JSLIT                

      COMMON /BL3/ RSIBB(401)                                           
      DATA (EB(1,J),J=1,8)/-17.0467,-15.9605,-13.8056,-10.6201,-6.4709,-
     11.4752,4.0588,10.1923/                                            
      DATA (EB(2,J),J=1,8)/-0.2735,0.,0.,0.,0.,0.,0.,0./                
      DATA TWOPIC/1.88365183D11/                                        
C                                                                       
C     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN      
C     WNRMAX  - THE FREQUENCY RANGE OF B-B CONTRIBUTION                 
C     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH       
C     B-B SPECTRUM, IN [CM-1].                                          
C                                                                       
      WNRMAX = 16.66 + SLIT 
C      WNRMAX = 20.64           (old line - replaced to allow general SLIT)
c	This is an estimated range of frequencies (from 0 to WNRMAX)
c	from the line center, in cm-1, where bound-bound intensities matter
c	at the given slit width. The number 16.66 is the maximum delta-E
C       from Table 10 of the paper from which this code derives.
C                                                          C.A.N. 27-1-97

      NSRI = 200
C      NSRI=129                               (old version)
C       This is the number of integration points. The maximum is 200,
C       may as well take the max. However, this means that DX now
C       varies (used to be fixed)

      DX = WNRMAX / NSRI
C      DX=WNRMAX /DFLOAT(NSRI)         (old version)       
C       This is the resolution of the integration grid (see subroutine
C       PROFILE). C.A.N.

      NSRIUP=2*NSRI+1

      DO 10 I=1,401                                                     
   10 RSIBB(I)=0.0                                                      
      NS=INT(SLIT/DX)                                                   
      DO 20 I=1,8                                                       
      DO 20 J=1,8                                                       
C                                                                       
C     A(I,J) FOR L=3 CONTRIBUTION (REDUCED BETA FUNCTION)               
C                                                                       
   20 A(I,J)=0.0                                                        
      A(1,4)=.1722                                                      
      A(2,3)=.2272                                                      
      A(2,5)=.2767                                                      
      A(3,4)=.2234                                                      
      A(3,6)=.3472                                                      
      A(4,5)=.2519                                                      
      A(4,7)=.3784                                                      
      A(5,6)=.2686                                                      
      A(5,8)=.4007                                                      
      A(6,7)=.2639                                                      
      A(7,8)=.2550                                                      
      DO 30 I=1,8                                                       
      DO 30 J=1,8                                                       
   30 A(I,J)=A(I,J)*1.D-41                                              
      ALFA=1./(0.69519*TEMP)                                            
      RM=1.8657*1.672649D-24                                            
      PI=3.141592654                                                    
      PF=((RM*(1.380662D-16)*TEMP*2.*PI/(6.626176D-27)**2)**1.5)        
      DO 50 L=1,7                                                       
         STOKE1=EB(1,L+1)-EB(1,L)                                       
C                                                                       
C     FREQUENCY SHIFT FOR DELTA J=+1,DELTA V=0                          
C                                                                       
         IF (L.GT.5) GO TO 40                                           
         STOKE3=EB(1,L+3)-EB(1,L)                                       
C                                                                       
C     STOKE3- FREQUENCY SHIFT FOR DELTA J=+3, DELTA V=0                 
C                                                                       
         STOKI=A(L,L+3)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE (STOKE3,STOKI)                                    
         STOKIP=A(L,L+3)*DEXP(-ALFA*EB(1,L+3))/PF                       
         CALL PROFILE (-STOKE3,STOKIP)                                  
   40    STOKI=A(L,L+1)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE (STOKE1,STOKI)                                    
         STOKIP=A(L,L+1)*DEXP(-ALFA*EB(1,L+1))/PF                       
         CALL PROFILE (-STOKE1,STOKIP)                                  
   50 CONTINUE                                                          
      AV1=0.3255D-43                                                    
      STOKE=EB(2,1)-EB(1,4)                                             
      STOKI=AV1*DEXP(-ALFA*EB(1,4))/PF                                  
      CALL PROFILE (STOKE,STOKI)                                        
      STOKIP=AV1*DEXP(-ALFA*EB(2,1))/PF                                 
      CALL PROFILE (-STOKE,STOKIP)                                      
      DO 60 N=1,NSRIUP                                                  
   60 RSIBB(N)=RSIBB(N)/TWOPIC/SLIT                                     
      NSOL=NSRI+1                                                       
      DO 70 I=1,NSOL                                                    
         K=(NSRI+1)+(I-1)                                            
   70 RSI(I)=RSIBB(K)                                                   
C                                                                       
C     RSI - BOUND-BOUND CONTRIBUTION FOR POSITIVE FREQUENCY SHIFTS     
C                                                                      
C      WRITE(6,80) SLIT                                                 
   80 FORMAT ( 10H SLITWIDTH,F10.4, 13H CM-1 ASSUMED,/)                 
      RETURN                                                            
      END                                                               
      SUBROUTINE BOUND54 (TEMP,RSI,NSOL)                                
      IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      DIMENSION EB(2,8), A(8,8), RSI(201)                               
C                                                                       
C     EB(I,K) - BOUND ENERGIES FOR VIB.NR.V=I-1, ROT.NR J=K-1           
C     A(I,J)=( C(I,L,J)**2) * (MTX.EL. (I,BETA,J)**2 )                  
C     I,J-DENOTE ROTATIONAL STATES FOR VIB. LEVEL V=0                   
C       SEE TABLE 10 OF AP. J. VOL.303, P.508 (1986)                    
C                                                                       
      COMMON /APP3/ SLIT,DX,WNRMAX,NSRI,NS,NSRIUP,JSLIT                 
      COMMON /BL3/ RSIBB(401)                                           
      DATA (EB(1,J),J=1,8)/-17.0467,-15.9605,-13.8056,-10.6201,-6.4709,-
     11.4752,4.0588,10.1923/                                            
      DATA (EB(2,J),J=1,8)/-0.2735,0.,0.,0.,0.,0.,0.,0./                
      DATA TWOPIC/1.88365183D11/                                        
C                                                                       
C     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN      
C     WNRMAX - THE FREQUENCY RANGE OF B-B CONTRIBUTION                  
C     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH       
C     B-B SPECTRUM, IN [CM-1].                                          
C                                                                       

      WNRMAX = 24.00 + SLIT   
C      WNRMAX = 28.16          (old line - replaced to allow general SLIT)
c	This is an estimated range of frequencies (from 0 to WNRMAX)
c	from the line center, in cm-1, where bound-bound intensities matter
c	at the given slit width. The number 24.00 is the maximum delta-E
C       from Table 10 (G5) of the paper from which this code derives.
C                                                          C.A.N. 27-1-97

      NSRI = 200
C      NSRI=175                         (old version)
C       This is the number of integration points. The maximum is 200,
C       may as well take it to the limit. However, this means that DX now
C       varies (used to be fixed)

      DX = WNRMAX / NSRI
C      DX=WNRMAX /DFLOAT(NSRI)          (old version)       
C       This is the resolution of the integration grid (see subroutine
C       PROFILE). C.A.N.

C      print*, 'WNRMAX=', WNRMAX
C      print*, 'NSRI  =', NSRI
C      print*, 'DX    =', DX

      NSRIUP=2*NSRI+1                                                   

      DO 10 I=1,401                                                     
   10 RSIBB(I)=0.0                                                      
      NS=INT(SLIT/DX)                                                   
      DO 20 I=1,8                                                       
      DO 20 J=1,8                                                       
C                                                                       
C     A(I,J) FOR L=5 CONTRIBUTION (REDUCED BETA FUNCTION)               
C                                                                       
   20 A(I,J)=0.0                                                        
      A(1,6)=0.08605                                                    
      A(2,5)=0.1294                                                     
      A(2,7)=0.1178                                                     
      A(3,4)=0.1506                                                     
      A(3,6)=0.1065                                                     
      A(3,8)=0.1360                                                     
      A(4,5)=0.1116                                                     
      A(4,7)=0.1008                                                     
      A(5,6)=0.1037                                                     
      A(5,8)=0.09753                                                    
      A(6,7)=0.09175                                                    
      A(7,8)=0.08088                                                    
      DO 30 I=1,8                                                       
      DO 30 J=1,8                                                       
   30 A(I,J)=A(I,J)*1.D-44                                              
      ALFA=1./(0.69519*TEMP)                                            
      RM=1.8657*1.672649D-24                                            
      PI=3.141592654                                                    
      PF=((RM*(1.380662D-16)*TEMP*2.*PI/(6.626176D-27)**2)**1.5)        
      DO 60 L=1,7                                                       
         STOKE1=EB(1,L+1)-EB(1,L)                                       
C                                                                       
C     FREQUENCY SHIFT FOR DELTA J=+1,DELTA V=0                          
C                                                                       
         IF (L.GT.5) GO TO 50                                           
         STOKE3=EB(1,L+3)-EB(1,L)                                       
C                                                                       
C     STOKE3- FREQUENCY SHIFT FOR DELTA J=+3, DELTA V=0                 
C                                                                       
         IF (L.GT.3) GO TO 40                                           
         STOKE5=EB(1,L+5)-EB(1,L)                                              
         STOKI=A(L,L+5)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE (STOKE5,STOKI)                                    
         STOKIP=A(L,L+5)*DEXP(-ALFA*EB(1,L+5))/PF                       
         CALL PROFILE (-STOKE5,STOKIP)                                  
   40    STOKI=A(L,L+3)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE (STOKE3,STOKI)                                    
         STOKIP=A(L,L+3)*DEXP(-ALFA*EB(1,L+3))/PF                       
         CALL PROFILE (-STOKE3,STOKIP)                                  
   50    STOKI=A(L,L+1)*DEXP(-ALFA*EB(1,L))/PF                          
         CALL PROFILE (STOKE1,STOKI)                                    
         STOKIP=A(L,L+1)*DEXP(-ALFA*EB(1,L+1))/PF                       
         CALL PROFILE (-STOKE1,STOKIP)                                  
   60 CONTINUE                                                          
      DO 70 N=1,NSRIUP                                                  
   70 RSIBB(N)=RSIBB(N)/TWOPIC/SLIT                                     
      NSOL=NSRI+1                                                       
      DO 80 I=1,NSOL                                                    
         K=(NSRI+1)+(I-1)                                               
   80 RSI(I)=RSIBB(K)                                                   
C                                                                       
      RETURN                                                            
      END                                                               
