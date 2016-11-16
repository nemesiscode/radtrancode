      SUBROUTINE GET_KG(NLAYER,PRESS,TEMP,NGAS,VWAVE)
C***********************************************************************
C_TITL:	GET_KG
C
C_DESC:	Interpolate to get k value at selected pressure and 
C	temperature. K-coefficients on the P/T grid have already
C	been read in. Interpolated coefficients are passed back to 
C	calling routine via a common block.
C
C_ARGS:	Input variables:
C	NLAYER		INTEGER	Number of layers.
C	PRESS(NLAYER)	REAL	Layer pressures [atm].
C	TEMP(NLAYER)	REAL	Layer temperatures [Kelvin].
C	NGAS		INTEGER	Number of gases.
C	VWAVE		REAL	Desired 'calculation' wavenumber.
C
C	Output variables (via common block):
C	KOUT(NLAYER,NGAS,NG)	REAL	K-coefficients for each layer and
C					each gas.	
C	DKOUTDT(NLAYER,NGAS,NG)	REAL	Rate of change of the above
C					k-coefficients with temperature.
C
C_FILE:	No files openned.
C
C_CALL:	LOCATE
C	RANKK		Subroutine to sort randomised k-coefficients.
C
C_HIST:	30/7/01	PGJI	Original version - conversion from get_k
C	7aug03	NT call rank changed to call rankk to avoid address
C			error caused by multiple occurances of rank.f in various
C			libraries. rankk.f in cirsradg/ is a copy of rank.f in 
C			cirsrad/
C	29/2/12	PGJI	Updated for Radtrans2.0
C
C***************************** VARIABLES *******************************

      IMPLICIT NONE

C The include file ...
      INCLUDE '../includes/arrdef.f'
C Defines the maximum values for a series of variables (layers, bins,
C paths, etc.)


C The input and output variables ...
      INTEGER nlayer,ngas,npoint
      REAL press(nlayer),temp(nlayer),vwave,frac
      REAL kout(maxlay,maxgas),dkoutdt(maxlay,maxgas)


C General variables ...
      INTEGER NP,NT,NG,CP,CT,I,I1,I2,I3,I4,J,N1,NX
C NP: Number of k-table pressures.
C NT: Number of k-table temperatures.
C NG: Number of ordinates in k-distribution.
      INTEGER igas,ilayer,irec

      INTEGER ntab
      INTEGER mtab
      PARAMETER (MTAB=maxk*maxk)
      REAL TABLE(mtab),TABLE2(mtab),frack

      REAL P1,T1,tmp,eps,KTEST
C T1: Profile temperature at each atmospheric layer.
C P1: LOG profile pressure at each atmospheric layer.
C EPS: Max allowed wavenumber difference between wavenumber and nearest 
C tabulated wavenumber to be considered 'aligned'.
      REAL Y1,Y2,Y3,Y4,U,V,X
      REAL pmax,pmin,tmax,tmin,sum,DXDT,DUDT
C PMIN: K-table pressure minimum [atm].
C PMAX: K-table pressure maximum [atm].
C TMIN: K-table temperature minimum [Kelvin].
C TMAX: K-table temperature maximum [Kelvin].

      INTEGER lun0,lun(maxgas),irec0
      REAL vmin,delv,DKDT,XT,vwavex

      REAL P(MAXK),T(MAXK),T2(MAXK,MAXK),TN(MAXK),X1,X2,U2
C P: K-table pressures [atm].
C T: K-table temperatures [Kelvin].

      INTEGER IOFF(MAXLAY,4),CT1,CT2,DUDT2
      REAL UT(MAXLAY),VT(MAXLAY),TDUDT(MAXLAY),FWHMK,DELVK
      REAL UT2(MAXLAY),TDUDT2(MAXLAY)

      LOGICAL KLOG,NTEST,ISNAN

      COMMON /INTERPKLBL/ LUN,IREC0,VMIN,DELV,NPOINT,P,NP,T,
     1 T2,NT,KOUT,DKOUTDT

C******************************** CODE *********************************

      pmax = p(np)
      pmin = p(1)
      if(NT.GT.0)THEN
       tmax = t(nt)
       tmin = t(1)
      else
       TMAX = T2(NP,ABS(NT))
       TMIN = T2(1,1)
      endif

C Work out where P,T for each layer lies in the tables
      DO 51 ilayer=1,NLAYER
        P1 = LOG(PRESS(ilayer))
C First check range is fine
        IF(P1.LT.PMIN)P1 = PMIN
        IF(P1.GT.PMAX)P1 = PMAX
        CALL locate(p,np,p1,cp)
        IF(cp.LT.1)cp = 1
        IF(cp.GE.np)cp = np - 1
        VT(ilayer) = (p1 - p(cp))/(p(cp + 1) - p(cp))

        IF(NT.GT.0)THEN
         T1 = TEMP(ilayer)
         IF(T1.LT.TMIN)T1 = TMIN
         IF(T1.GT.TMAX)T1 = TMAX
         CALL locate(t,nt,t1,ct)
         IF(ct.LT.1)ct = 1
         IF(ct.GE.nt)ct = nt - 1
         UT(ilayer) = (t1 - t(ct))/(t(ct + 1) - t(ct))
         TDUDT(ilayer) = 1.0/(t(ct + 1) - t(ct))

         IOFF(ilayer,1) = nt*(cp - 1) + (ct - 1)
         IOFF(ilayer,2) = nt*cp       + (ct - 1)
         IOFF(ilayer,3) = nt*cp       + ct
         IOFF(ilayer,4) = nt*(cp - 1) + ct

        ELSE

         DO I=1,ABS(NT)
          TN(I)=T2(CP,I)
         ENDDO
         T1=TEMP(ilayer)
         IF(T1.LT.TN(1))T1=TN(1)
         IF(T1.GT.TN(ABS(NT)))T1=TN(ABS(NT))
         CALL LOCATE(TN,ABS(NT),T1,CT1)
         IF(ct1.lt.1)ct1 = 1
         IF(ct1.ge.ABS(nt))ct1=ABS(nt)-1
         UT(iLAYER)=(T1-TN(CT1))/(TN(CT1+1)-TN(CT1))
         TDUDT(ilayer) = 1.0/(tn(ct1 + 1) - tn(ct1))

         DO I=1,ABS(NT)
          TN(I)=T2(CP+1,I)
         ENDDO
         T1=TEMP(ilayer)
         IF(T1.LT.TN(1))T1=TN(1)
         IF(T1.GT.TN(ABS(NT)))T1=TN(ABS(NT))
         CALL LOCATE(TN,ABS(NT),T1,CT2)
         IF(ct2.lt.1)ct2 = 1
         IF(ct2.ge.ABS(nt))ct2=ABS(nt)-1
         UT2(iLAYER)=(T1-TN(CT2))/(TN(CT2+1)-TN(CT2))
         TDUDT2(ilayer) = 1.0/(tn(ct2 + 1) - tn(ct2))
         nx=abs(nt)
         IOFF(ilayer,1) = nx*(cp - 1) + (ct1 - 1)
         IOFF(ilayer,2) = nx*cp       + (ct2 - 1)
         IOFF(ilayer,3) = nx*cp       + ct2
         IOFF(ilayer,4) = nx*(cp - 1) + ct1

        ENDIF
51    CONTINUE

      DO 1000 IGAS=1,NGAS
        LUN0 = LUN(IGAS)
        IF(LUN0.LE.0)THEN
cc          WRITE(*,*)'GET_KGLBL :: No data defined for gas : ',IGAS
          DO ilayer=1,NLAYER
              KOUT(ilayer,IGAS) = 0.0
              DKOUTDT(ilayer,IGAS) = 0.0
          ENDDO
          GOTO 999
        ENDIF

C       Find nearest point in table below current wavelength. Parameter eps is
C       there to prevent small numerical errors in VWAVE screwing things up
C       between platforms.

C       Set minimum closeness to tabulated wavenumbers to be 1/50 of
C       the separation to be considered aligned.
        eps = 0.02*delv

        n1 = INT((vwave + eps - vmin)/delv)
        vwavex = VMIN + (N1-1)*DELV
        XT = ABS(vwave - vwavex)
        IF(XT.GT.EPS)THEN
          print*,'Wavenumber does not match'
          print*,vwave,vwavex
          stop
        ENDIF

        irec = irec0 + np*abs(nt)*(n1-1)

        NTAB = ABS(NT)*NG
        KTEST=0.0   

        DO I=1,NTAB
         READ(LUN0,REC=IREC)TABLE(I) 
         IF(TABLE(I).GT.0.0)THEN
          KTEST=TABLE(I)
          GOTO 202
         ENDIF
         IREC=IREC+1
        ENDDO

202     CONTINUE

        IF(KTEST.EQ.0.0)THEN
          WRITE(*,*)'GET_KBLG :: Zero k-data for GAS: ',IGAS
          WRITE(*,*)'Wavelength/Wavenumber = ',VWAVE
          DO ilayer=1,NLAYER
              KOUT(ilayer,IGAS) = 0.0
              DKOUTDT(ilayer,IGAS) = 0.0
          ENDDO
          GOTO 999
        ELSE
          IREC = IREC0+NTAB*(N1-1)
          DO I=1,NTAB
           READ(LUN0,REC=IREC)TABLE(I)
           IREC=IREC+1
          ENDDO
        ENDIF

C       Now interpolate k-coefficients for conditions in each layer
        DO 1050 ilayer=1,NLAYER
          U = UT(ilayer)
          U2 = UT2(ilayer)
          V = VT(ilayer)
          DUDT = TDUDT(ilayer)
          DUDT2 = TDUDT2(ilayer)

          I1 = IOFF(ilayer,1)
          I2 = IOFF(ilayer,2)
          I3 = IOFF(ilayer,3)
          I4 = IOFF(ilayer,4)

          Y1 = TABLE(I1+1)
          Y2 = TABLE(I2+1)
          Y3 = TABLE(I3+1)
          Y4 = TABLE(I4+1)

          KLOG = .TRUE.
          IF((Y1.LE.0.0).OR.(Y2.LE.0.0).OR.(Y3.LE.0.0).OR.
     1      (Y4.LE.0.0))THEN
              KLOG = .FALSE.
          ENDIF

          IF(KLOG)THEN
                Y1 = LOG(Y1)
                Y2 = LOG(Y2)
                Y3 = LOG(Y3)
                Y4 = LOG(Y4)
          ENDIF

          IF(NT.GT.0)THEN
               X = (1.0 - V)*(1.0 - U)*Y1 + V*(1.0 - U)*Y2 + V*U*Y3 +
     1          (1.0 - V)*U*Y4             
               DXDT = (-(1.0 - V)*Y1 - V*Y2 + V*Y3 + (1 - V)*Y4)*DUDT
          ELSE              
               X1=(1.0-U)*Y1 + U*Y4
               X2=(1.0-U2)*Y2 + U2*Y3
               X = (1.0-V)*X1 + V*X2
               DXDT = (1.0-V)*(-Y1+Y4)*DUDT+V*(-Y2+Y3)*DUDT2 
          ENDIF

          IF(KLOG)THEN
                KOUT(ILAYER,IGAS) = EXP(X)
                DKOUTDT(ILAYER,IGAS) = EXP(X)*DXDT
          ELSE
                KOUT(ILAYER,IGAS) = X
                DKOUTDT(ILAYER,IGAS) = DXDT
          ENDIF
          NTEST=ISNAN(KOUT(ILAYER,IGAS))
          IF(NTEST)THEN
             KOUT(ILAYER,IGAS)=1e-37  
             PRINT*,'Warning, NAN returned by get_klbl.f for gas',igas
             print*,'         VWAVE = ',VWAVE
             print*,'         LAYER,PRESS,TEMP = ',ILAYER,
     1        PRESS(ILAYER),TEMP(ILAYER)
            STOP
          ENDIF
          NTEST=ISNAN(DKOUTDT(ILAYER,IGAS))
          IF(NTEST)THEN
             DKOUTDT(ILAYER,IGAS)=1e-37  
             PRINT*,'Warning, NAN returned by get_klbl.f for gas',igas
             print*,'         VWAVE = ',VWAVE
             print*,'         LAYER,PRESS,TEMP = ',ILAYER,
     1        PRESS(ILAYER),TEMP(ILAYER)
            STOP
          ENDIF

1050    CONTINUE

 999    CONTINUE

1000  CONTINUE

      RETURN

      END
