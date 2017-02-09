************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE GET_KLBL	
C
C	Interpolate to get k value at selected pressure and 
C	temperature. K coefficients on the P/T grid have already
C	been read in.
C
C-----------------------------------------------------------------------
      	SUBROUTINE GET_KLBL(NLAYER,PRESS,TEMP,NGAS,VWAVE)
	IMPLICIT NONE

        INCLUDE '../includes/arrdef.f'

	INTEGER	NP, NT, CP, CT, I, I1, I2, I3, I4, N1,IREC
	INTEGER NGAS,NLAYER,IGAS,LAYER,LUN0,IREC0
        INTEGER NTAB,MTAB,J,CT1,CT2, NX
        LOGICAL NTEST,ISNAN
        parameter (MTAB=maxk*maxk)
        REAL TABLE(MTAB),TABLE2(MTAB),X1,X2,U2,VWAVEX,XT

	REAL	P1, T1,DELV,VMIN, tmp, eps, KTEST,
     1		Y1, Y2, Y3, Y4, U, V, pmax, pmin, tmax, tmin, X, Z
        REAL PRESS(NLAYER),TEMP(NLAYER),VWAVE,FRAC

C       Defines the maximum values for a series of variables (layers, 
C       bins, paths, etc.)
        integer lun(maxgas)
	REAL KOUT(MAXLAY,MAXGAS)
        REAL dkoutdt(maxlay,maxgas)
        REAL P(MAXK),T(MAXK),T2(MAXK,MAXK),TN(MAXK)
	REAL UT(MAXLAY),VT(MAXLAY),UT2(MAXLAY)
        INTEGER IOFF(MAXLAY,4),NPOINT

	COMMON /INTERPKLBL/LUN,IREC0,VMIN,DELV,NPOINT,P,NP,T,
     1   T2,NT,KOUT,DKOUTDT
C-----------------------------------------------------------------------


        eps=0.01
        PMAX = P(NP)
        PMIN = P(1)
        IF(NT.GT.0)THEN
         TMAX = T(NT)
         TMIN = T(1)
        ELSE
         TMAX = T2(NP,ABS(NT))
         TMIN = T2(1,1)
        ENDIF


        DO 1000 IGAS=1,NGAS

         LUN0 = LUN(IGAS)
         IF(LUN0.LE.0) THEN
C	  print*,'No data defined for gas : ',IGAS
          DO LAYER=1,NLAYER
            KOUT(LAYER,IGAS)=0.0
          ENDDO
          GOTO 999
         ENDIF


C        Set minimum closeness to tabulated wavenumbers to be 1/50 of 
C        the separation
         eps = 0.02*delv

C        find nearest point in table below current wavenumber
         N1 = 1 + int((VWAVE + EPS - VMIN)/DELV)
C						  below current wavelength.
C						  Parameter eps is there 
C						  to prevent small numerical
C						  errors in VWAVE screwing
C					          things up between platforms
         if(n1.lt.1.or.n1.gt.npoint)then
          print*,'Wavelength/wavenumber is out of range'
          print*,vwave,vmin,vmin+(npoint-1)*delv
         endif

C         vwavex = VMIN + (N1-1)*DELV
C         XT = ABS(vwave - vwavex)
C         IF(XT.GT.EPS)THEN
C          print*,'Wavenumber does not match'
C          print*,vwave,vwavex
C          stop
C         ENDIF

         IREC = IREC0+NP*ABS(NT)*(N1-1)

         NTAB = ABS(NT)*NP
         if(NTAB.gt.MTAB)then
           print*,'Error in get_kg, NTAB>MTAB'
           print*,NTAB,MTAB
           stop
         endif

         KTEST=0

         DO I=1,NTAB
          READ(LUN0,REC=IREC)TABLE(I)
          IF(TABLE(I).GT.0.0)THEN
           KTEST=TABLE(I)
           GOTO 202
          ENDIF
          IREC=IREC+1
         ENDDO

202      CONTINUE

         IF(KTEST.EQ.0.0)THEN
          print*,'GET_K: Zero k-data for GAS: ',IGAS
          print*,'Wavelength/waveumber = ',VWAVE
          DO LAYER=1,NLAYER
            KOUT(LAYER,IGAS)=0.0
          ENDDO
          GOTO 999
         ELSE
          IREC = IREC0+NTAB*(N1-1)
          DO I=1,NTAB
           READ(LUN0,REC=IREC)TABLE(I)
           IREC=IREC+1
          ENDDO
         ENDIF


C	 Now interpolate k-coefficients for conditions in each layer
         DO 1050 LAYER=1,NLAYER
C          print*,'LAYER = ',LAYER,PRESS(LAYER),
C     1    TEMP(LAYER)

C          IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C           Work out where P,T for each layer lies in the tables
C        ------------------------------------
C	   First check range is fine.
           P1 = LOG(PRESS(LAYER))

C	   Find position of temp and pressure values in k table, then
C	   interpolate
C-----------------------------------------------------------------------
           IF(P1.LT.PMIN)P1=PMIN
           IF(P1.GT.PMAX)P1=PMAX
 	   CALL locate(p, np, p1, cp)
 	   IF(cp.lt.1)cp = 1
	   IF(cp.ge.np)cp = np-1
           VT(LAYER)=(P1-P(CP))/(P(CP+1)-P(CP))

           IF(NT.GT.0)THEN
              T1 = TEMP(LAYER)
              IF(T1.LT.TMIN)T1=TMIN
              IF(T1.GT.TMAX)T1=TMAX
              CALL locate(t, nt, t1, ct)
 	      IF(ct.lt.1)ct = 1
              IF(ct.ge.nt)ct=nt-1
              UT(LAYER)=(T1-T(CT))/(T(CT+1)-T(CT))

              IOFF(LAYER,1) = nt*(cp-1) + ct-1
              IOFF(LAYER,2) = nt*cp     + ct-1
              IOFF(LAYER,3) = nt*cp     + ct
              IOFF(LAYER,4) = nt*(cp-1) + ct

           ELSE
              DO I=1,ABS(NT)
               TN(I)=T2(CP,I)
              ENDDO
              T1=TEMP(LAYER)
              IF(T1.LT.TN(1))T1=TN(1)
              IF(T1.GT.TN(ABS(NT)))T1=TN(ABS(NT))
              CALL LOCATE(TN,ABS(NT),T1,CT1)
 	      IF(ct1.lt.1)ct1 = 1
              IF(ct1.ge.ABS(nt))ct1=ABS(nt)-1
              UT(LAYER)=(T1-TN(CT1))/(TN(CT1+1)-TN(CT1))

              DO I=1,ABS(NT)
               TN(I)=T2(CP+1,I)
              ENDDO
              T1=TEMP(LAYER)
              IF(T1.LT.TN(1))T1=TN(1)
              IF(T1.GT.TN(ABS(NT)))T1=TN(ABS(NT))
              CALL LOCATE(TN,ABS(NT),T1,CT2)
 	      IF(ct2.lt.1)ct2 = 1
              IF(ct2.ge.ABS(nt))ct2=ABS(nt)-1
              UT2(LAYER)=(T1-TN(CT2))/(TN(CT2+1)-TN(CT2))
              nx=abs(nt)
              IOFF(LAYER,1) = nx*(cp-1) + ct1-1
              IOFF(LAYER,2) = nx*cp     + ct2-1
              IOFF(LAYER,3) = nx*cp     + ct2
              IOFF(LAYER,4) = nx*(cp-1) + ct1

           ENDIF            


C          IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	   U=UT(LAYER)
	   U2=UT2(LAYER)
	   V=VT(LAYER)
	   I1 = IOFF(LAYER,1)
	   I2 = IOFF(LAYER,2)
	   I3 = IOFF(LAYER,3)
	   I4 = IOFF(LAYER,4)


C          Check if this should all be +1 (I think OK)
           Y1=TABLE(I1+1)
           Y2=TABLE(I2+1)
           Y3=TABLE(I3+1)
           Y4=TABLE(I4+1)


           IF((Y1.LE.0.0).OR.(Y2.LE.0.0).OR.
     &	     (Y3.LE.0.0).OR.(Y4.LE.0.0))THEN
		        Z = -1.0
           ELSE
                        Y1=LOG(Y1)
                        Y2=LOG(Y2)
                        Y3=LOG(Y3)
                        Y4=LOG(Y4)
			Z = 1.0
           ENDIF

           IF(NT.GT.0)THEN
		 X=(1.0-V)*(1.0-U)*Y1 + V*(1.0-U)*Y2 + 
     1			V*U*Y3 + (1.0-V)*U*Y4
           ELSE
                 X1=(1.0-U)*Y1 + U*Y4
                 X2=(1.0-U2)*Y2 + U2*Y3
                 X = (1.0-V)*X1 + V*X2
           ENDIF

	   IF(Z.LE.0.0)THEN
		KOUT(LAYER,IGAS)=X
	   ELSE
		KOUT(LAYER,IGAS)=EXP(X)
	   ENDIF
           NTEST=ISNAN(KOUT(LAYER,IGAS))
           IF(NTEST)THEN
            KOUT(LAYER,IGAS)=1e-37
            PRINT*,'Warning, NAN returned by get_klbl.f for gas',igas
            stop
           ENDIF

1050     CONTINUE


999      CONTINUE

1000    CONTINUE

C        print*,'kout',kout(1,1,1)
	RETURN

	END

************************************************************************
************************************************************************
