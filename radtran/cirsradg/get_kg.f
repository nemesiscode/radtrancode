      SUBROUTINE GET_KG(NLAYER,PRESS,TEMP,NGAS,IWAVE,VWAVE)
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
C	IWAVE		INTEGER	Ordinate of desired wavenumber (k-tables
C				read in earlier).
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
      INTEGER nlayer,ngas,iwave
      REAL press(nlayer),temp(nlayer),vwave,frac
      REAL kout(maxlay,maxgas,maxg),dkoutdt(maxlay,maxgas,maxg)


C General variables ...
      INTEGER NP,NT,NG,CP,CT,I,I1,I2,I3,I4,J,N1,NX
C NP: Number of k-table pressures.
C NT: Number of k-table temperatures.
C NG: Number of ordinates in k-distribution.
      INTEGER igas,ilayer,irec,IRECX

      INTEGER ntab,loop,count
      INTEGER maxc,mtab
      PARAMETER (maxc=2*maxg,MTAB=maxk*maxk*maxg)
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

      REAL weight(2),fac(maxc),cont(maxc),dcont(maxc)

      INTEGER lun0,lun(maxbin,maxgas),irec0,ireck(maxbin,maxgas)
      REAL vmin,xmin(maxbin,maxgas),delv,delx(maxbin,maxgas)
      REAL fracx(maxbin,maxgas)

      REAL K_G(MAXG),G_ORD(MAXG),DELG(MAXG),DKDT(MAXG)
C K_G: Calculated k-distribution.
C G_ORD: Gauss-Legendre ordinates for calculating the k-distribution.
C DEL_G: Gauss-Legendre weights for integration.

      REAL P(MAXK),T(MAXK),T2(MAXK,MAXK),TN(MAXK),X1,X2,U2
C P: K-table pressures [atm].
C T: K-table temperatures [Kelvin].

      INTEGER IOFF(MAXLAY,4),CT1,CT2,DUDT2
      REAL UT(MAXLAY),VT(MAXLAY),TDUDT(MAXLAY),FWHMK,DELVK
      REAL UT2(MAXLAY),TDUDT2(MAXLAY)

      LOGICAL COINC,KLOG,NTEST,ISNAN

      COMMON /INTERPK/ LUN,IRECK,XMIN,DELX,FRACX,P,NP,T,T2,NT,NG,
     1 DELVK,FWHMK,G_ORD,DELG,KOUT,DKOUTDT

C******************************** CODE *********************************

      pmax = p(np)
      pmin = p(1)
      if(NT.GT.0)THEN
       tmax = t(nt)
       tmin = t(1)
      endif

      print*,'get_kg: ng =',ng
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

         IOFF(ilayer,1) = (nt*ng*(cp - 1)) + (ng*(ct - 1))
         IOFF(ilayer,2) = (nt*ng*cp) + (ng*(ct - 1))
         IOFF(ilayer,3) = (nt*ng*cp) + (ng*ct)
         IOFF(ilayer,4) = (nt*ng*(cp - 1)) + (ng*ct)

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
         IOFF(ilayer,1) = (nx*ng*(cp - 1)) + (ng*(ct1 - 1))
         IOFF(ilayer,2) = (nx*ng*cp) + (ng*(ct2 - 1))
         IOFF(ilayer,3) = (nx*ng*cp) + (ng*ct2)
         IOFF(ilayer,4) = (nx*ng*(cp - 1)) + (ng*ct1)

        ENDIF
51    CONTINUE

      DO 1000 IGAS=1,NGAS
        LUN0 = LUN(IWAVE,IGAS)
        VMIN = XMIN(IWAVE,IGAS)
        DELV = DELX(IWAVE,IGAS)
        FRACK = FRACX(IWAVE,IGAS)
        IREC0 = IRECK(IWAVE,IGAS)
C        print*,'iwave,igas,vmin',IWAVE,IGAS,VMIN
        IF(LUN0.LE.0)THEN
cc          WRITE(*,*)'GET_KG :: No data defined for gas : ',IGAS
          DO ilayer=1,NLAYER
            DO I=1,NG
              KOUT(ilayer,IGAS,I) = 0.0
              DKOUTDT(ilayer,IGAS,I) = 0.0
            ENDDO
          ENDDO
          GOTO 999
        ENDIF

        IF(delv.gt.0)then
C        Find nearest point in table below current wavelength. Parameter eps is
C        there to prevent small numerical errors in VWAVE screwing things up
C        between platforms.

C        Set minimum closeness to tabulated wavenumbers to be 1/50 of
C        the separation to be considered aligned.
         eps = 0.02*delv

          n1 = INT((vwave + eps - vmin)/delv)
          irec = irec0 + np*abs(nt)*ng*n1

        ELSE

C         For irregularly gridded tables IREC0 is assumed to hold the
C         nearest record number in the table to the requested wavelength, not
C         that of the start of the table. This has already been allocated
C         by read_klist.f

          irec = irec0

        ENDIF

        NTAB = ABS(NT)*NP*NG
        IRECX=IREC
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
        IREC=IRECX

        IF(KTEST.EQ.0.0)THEN
          WRITE(*,*)'GET_KG :: Zero k-data for GAS: ',IGAS
          WRITE(*,*)'Wavelength/Wavenumber = ',VWAVE
          DO ilayer=1,NLAYER
            DO I=1,NG
              KOUT(ilayer,IGAS,I) = 0.0
              DKOUTDT(ilayer,IGAS,I) = 0.0
            ENDDO
          ENDDO
          GOTO 999
        ENDIF


        IF(delv.gt.0)THEN
C        Calculate wavelength in table below current wavelength
         tmp = vmin + delv*n1
         frac = (vwave-tmp)/delv
         if(frac.lt.eps)then
C         delv > 0 and requested wavenumber close enough to tabulated to
C         be considered coincident (and thus need to interrogate k-table once
C         no twice.
          COINC=.TRUE.
         else
          COINC=.FALSE.
         endif
C         print*,'GET_KG: vwave, COINC = ',vwave,COINC
C         print*,'Nearest tabulated+ frac: ',tmp,tmp+delv,frac
        ELSE
          if(delv.lt.0)then
C          If delv < 0 then it is assumed that the calculation wavelengths
C          coincide with the central wavelengths and thus IREC0 is assumed to
C          hold the current record number of the nearest wavelength in the
C          table, already set up in read_klist.f.
           COINC=.TRUE.
          else
C          If delv = 0 then it is not assumed that the calculation wavelengths
C          coincide with the central wavelengths and thus IREC0 is assumed to
C          hold the current record number of the nearest wavelength in the
C          table BELOW that requested, already set up in read_klist.f.
           COINC=.FALSE.
          endif
        ENDIF

        ntab = np*abs(nt)*ng
        if(ntab.gt.MTAB)then
         print*,'Error in get_kg, NTAB>MTAB'
         print*,NTAB,MTAB
         stop
        endif
C       Read in table at or below current wavelength
        DO I=1,ntab
          READ(LUN0,REC=IREC)TABLE(I)
          irec = irec + 1
        ENDDO

C       If not in-line, read in k-coeffs for next bin also
        IF(.NOT.coinc)THEN
          DO I=1,ntab
            READ(LUN0,REC=IREC)TABLE2(I)
            irec = irec + 1
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

          IF(coinc)THEN
C=======================================================================
C
C	If wavenumber lines up with table, just read in the k-data for
C	that wavenumber.
C
C	Check to see if any tabulated k-coefficients are zero or 
C	negative. Provided that this is not the case, interpolate
C	log(k-coefficients). This tends to be more accurate(?).
C
C=======================================================================
            Y1 = TABLE(I1+1)
            Y2 = TABLE(I2+1)
            Y3 = TABLE(I3+1)
            Y4 = TABLE(I4+1)

            KLOG = .TRUE.
            IF((Y1.LE.0.0).OR.(Y2.LE.0.0).OR.(Y3.LE.0.0).OR.
     1      (Y4.LE.0.0))THEN
              KLOG = .FALSE.
            ENDIF

            DO I=1,NG
              Y1 = TABLE(I1+I)
              Y2 = TABLE(I2+I)
              Y3 = TABLE(I3+I)
              Y4 = TABLE(I4+I)
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
                K_G(I) = EXP(X)
                DKDT(I) = EXP(X)*DXDT
              ELSE
                K_G(I) = X
                DKDT(I) = DXDT
              ENDIF
              NTEST=ISNAN(K_G(I))
              IF(NTEST)THEN
               print*,'NAN in get_kg'
               print*,'Y1,Y2,Y3,Y4',Y1,Y2,Y3,Y4
               print*,'u,v,x',u,v,x
              ENDIF
            ENDDO
          ELSE
            if(delv.gt.0)then 
             weight(2) = ((vwave - tmp)/delv)
             IF(weight(2).GT.1.0)weight(2) = 1.0
             IF(weight(2).LT.0.0)weight(2) = 0.0
             weight(1) = 1.0 - weight(2)
            else
             weight(1)=frack
             weight(2)=1.0-frack
            endif
C            print*,'Vwave,weights : ',vwave,weight(1),weight(2)
c Fletcher: Initialised cont and dcont arrays, and the parameter X.
	    DO I=1,MAXC
		cont(i)=0.0
		dcont(i)=0.0
            ENDDO
	    X=0.0

            NTAB = np*abs(nt)*ng
            count = 0
            DO loop=1,2
              DO I=1,NG
                IF(loop.EQ.1)THEN
                  Y1 = TABLE(I1+I)
                  Y2 = TABLE(I2+I)
                  Y3 = TABLE(I3+I)
                  Y4 = TABLE(I4+I)
                ELSE
                  Y1 = TABLE2(I1+I)
                  Y2 = TABLE2(I2+I)
                  Y3 = TABLE2(I3+I)
                  Y4 = TABLE2(I4+I)
                ENDIF

                IF(NT.GT.0)THEN
                 X = (1.0 - V)*(1.0 - U)*Y1 + V*(1.0 - U)*Y2 +
     1            V*U*Y3 + (1.0 - V)*U*Y4
                 DXDT = (-(1.0 - V)*Y1 - V*Y2 + V*Y3 + 
     1            (1 - V)*Y4)*DUDT
                ELSE
                 X1=(1.0-U)*Y1 + U*Y4
                 X2=(1.0-U2)*Y2 + U2*Y3
                 X = (1.0-V)*X1 + V*X2
                 DXDT = (1.0-V)*(-Y1+Y4)*DUDT+V*(-Y2+Y3)*DUDT2 
                ENDIF
                count = count + 1
                cont(count) = X
                dcont(count) = DXDT
                fac(count) = delg(I) * weight(loop)
              ENDDO
            ENDDO

            sum = 0.0
            DO I=1,count
              sum = sum + fac(I)
            ENDDO
c	    write(*,*)'Before rankk (count):',count
c	    write(*,*)'Before rankk:',cont
            CALL rankk(delg,ng,cont,dcont,fac,count,k_g,dkdt)
C	    write(*,*)'After rankk:',(k_g(i),i=1,ng)

          ENDIF


          DO I=1,NG
            NTEST=ISNAN(K_G(I))
            IF(NTEST)THEN
             KOUT(ILAYER,IGAS,I)=1e-37  
             PRINT*,'Warning, NAN returned by get_k.f for gas',igas
             print*,'         IWAVE,VWAVE = ',IWAVE,VWAVE
             print*,'         LAYER,PRESS,TEMP = ',ILAYER,
     1        PRESS(ILAYER),TEMP(ILAYER)
            ELSE
             KOUT(ILAYER,IGAS,I)=K_G(I)
            ENDIF

            NTEST=ISNAN(DKDT(I))
            IF(NTEST)THEN
             DKOUTDT(ILAYER,IGAS,I)=1e-37  
             PRINT*,'Warning, Grad NAN returned by get_k.f for gas',igas
             print*,'         IWAVE,VWAVE = ',IWAVE,VWAVE
             print*,'         LAYER,PRESS,TEMP = ',ILAYER,
     1        PRESS(ILAYER),TEMP(ILAYER)
            ELSE
             DKOUTDT(ILAYER,IGAS,I)=DKDT(I)
            ENDIF            
          ENDDO

c 	  write(*,*)'End of get_kg:',k_g

1050    CONTINUE

 999    CONTINUE

1000  CONTINUE

      RETURN

      END
