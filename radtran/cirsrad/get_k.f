************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE GET_K	
C
C	Interpolate to get k value at selected pressure and 
C	temperature. K coefficients on the P/T grid have already
C	been read in.
C
C-----------------------------------------------------------------------
      	SUBROUTINE GET_K(NLAYER,PRESS,TEMP,NGAS,IWAVE,VWAVE)

	IMPLICIT NONE

	INTEGER	NP, NT, NG, CP, CT, I, I1, I2, I3, I4, N1,IREC
	INTEGER IWAVE,NGAS,NLAYER,IGAS,LAYER,LUN0,IREC0,maxc
        INTEGER NTAB,loop,count,MTAB,J
	LOGICAL COINC
        parameter (maxc=40,MTAB=20*20*20)
        REAL TABLE(MTAB),TABLE2(MTAB)

	REAL	P1, T1,DELV,VMIN, tmp, eps, KTEST,
     1		Y1, Y2, Y3, Y4, U, V, pmax, pmin, tmax, tmin, X, Z
        real weight(2),fac(maxc),cont(maxc),sum

        parameter (eps=1e-6)		! Max allowed wavenumber difference
C					  between wavenumber and nearest
C					  tabulated wavenumber to be considered
C					  'aligned'.

        REAL PRESS(NLAYER),TEMP(NLAYER),VWAVE

C       Defines the maximum values for a series of variables (layers, 
C       bins, paths, etc.)
        INCLUDE '../includes/arrdef.f'

        integer lun(maxbin,maxgas), ireck(maxbin,maxgas)
        real xmin(maxbin,maxgas), delx(maxbin,maxgas)
	REAL KOUT(MAXLAY,MAXGAS,MAXG),K_G(MAXG),G_ORD(MAXG)
        REAL DELG(MAXG)
        REAL dkoutdt(maxlay,maxgas,maxg)
        REAL P(MAXK),T(MAXK)
	REAL UT(MAXLAY),VT(MAXLAY),FWHMK,DELVK
        INTEGER IOFF(MAXLAY,4)
        LOGICAL INTERP

	COMMON /INTERPK/LUN,IRECK,XMIN,DELX,P,NP,T,NT,
     1			NG,DELVK,FWHMK,G_ORD,DELG,KOUT,DKOUTDT
C-----------------------------------------------------------------------

C	print*,'GET_K calling parameters'
C        print*,'NLAYER=',NLAYER
C        DO I=1,NLAYER
C         print*,press(I),temp(I)
C        ENDDO
C        print*,'ngas,iwave,vwave = ',ngas,iwave,vwave
C        print*,'igas,lun,ireck,xmin,delx'
C        do igas=1,ngas
C         print*,igas,lun(iwave,igas),ireck(iwave,igas),
C     1    xmin(iwave,igas),delx(iwave,igas)
C        enddo
C        print*,'k-pressures : ',np,(P(I),I=1,NP)
C        print*,'k-temps : ',nt,(T(I),I=1,NT)
C        print*,'ng,delg : ',ng,(delg(i),i=1,ng)
C        print*,'kout(1,1,1)=',kout(1,1,1)

        PMAX = P(NP)
        PMIN = P(1)
        TMAX = T(NT)
        TMIN = T(1)

        INTERP = .TRUE.

        DO 1000 IGAS=1,NGAS

         LUN0 = LUN(IWAVE,IGAS)
         VMIN = XMIN(IWAVE,IGAS)
         DELV = DELX(IWAVE,IGAS)
         IREC0 = IRECK(IWAVE,IGAS)

C         print*,'IGAS,DELV',IGAS,DELV
         IF(LUN0.LE.0) THEN
C	  print*,'No data defined for gas : ',IGAS
          DO LAYER=1,NLAYER
           DO I=1,NG
            KOUT(LAYER,IGAS,I)=0.0
           ENDDO
          ENDDO
          GOTO 999
         ENDIF

C         print*,'IGAS,IREC0 : ',IGAS,IREC0
C         print*,'delv = ',delv
         IF(delv.gt.0.0)then
          N1 = int((VWAVE - VMIN)/DELV + eps)	! find nearest point in table
C						  below current wavelength.
C						  Parameter eps is there 
C						  to prevent small numerical
C						  errors in VWAVE screwing
C					          things up between platforms
          IREC = IREC0+NP*NT*NG*N1
         ELSE
C         If irregularly gridded table, then it is assumed that the 
C         calculation wavelengths coincide with the central wavelengths
C         and thus IREC0 is assumed to hold the current record number, not
C         that of the start of the table
          irec = irec0 
          N1=-1   
C          print*,'get_k',iwave,vwave,irec
         ENDIF

C	 print*,'igas,lun0,irec0,irec',IGAS,' ',lun0,' ',irec0,' ',irec
         READ(LUN0,REC=IREC)KTEST

         IF(KTEST.EQ.0.0)THEN
          print*,'GET_K: Zero k-data for GAS: ',IGAS
          DO LAYER=1,NLAYER
           DO I=1,NG
            KOUT(LAYER,IGAS,I)=0.0
           ENDDO
          ENDDO
          GOTO 999
         ENDIF


         IF(delv.gt.0.0)then
          tmp = VMIN + N1*DELV	! Calculate wavelength in table below
C				  current wavelength

C         Does wavenumber coincide with a tabulated wavenumber?
          COINC=.FALSE.
          IF (abs(tmp - vwave).LT.eps) COINC=.TRUE.
         ELSE
          COINC=.TRUE.
         ENDIF

C         print*,abs(tmp-vwave)
C         print*,'Get_K. COINC,N1,VWAVE = ',COINC,N1,VWAVE
C         stop

         NTAB = NT*NP*NG
         if(NTAB.gt.MTAB)then
          print*,'Error in get_kg, NTAB>MTAB'
          print*,NTAB,MTAB
          stop
         endif

         DO I=1,NTAB
          READ(LUN0,REC=IREC)TABLE(I)
C          print*,IREC,TABLE(I)
          IREC=IREC+1
         ENDDO

C        If not in-line, read in k-coeffs for next bin also
         IF(.NOT.COINC)THEN
          DO I=1,NTAB
           READ(LUN0,REC=IREC)TABLE2(I)
           IREC=IREC+1
          ENDDO
         ENDIF

C	 Now interpolate k-coefficients for conditions in each layer
         DO 1050 LAYER=1,NLAYER
C          print*,'LAYER = ',LAYER,PRESS(LAYER),TEMP(LAYER)

C        IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C           Work out where P,T for each layer lies in the tables
C        ------------------------------------
	  IF(INTERP)THEN
C	    First check range is fine.
            T1 = TEMP(LAYER)
            P1 = LOG(PRESS(LAYER))

            IF(P1.LT.PMIN)P1=PMIN
            IF(P1.GT.PMAX)P1=PMAX
            IF(T1.LT.TMIN)T1=TMIN
            IF(T1.GT.TMAX)T1=TMAX

C            print*,'P1,T1',P1,T1
C	    Find position of temp and pressure values in k table, then
C	    interpolate
C-----------------------------------------------------------------------

 	    CALL locate(p, np, p1, cp)
	    CALL locate(t, nt, t1, ct)

 	    IF(cp.lt.1)cp = 1
	    IF(ct.lt.1)ct = 1
	    IF(cp.ge.np)cp=np-1
	    IF(ct.ge.nt)ct=nt-1

C            PRINT*,'CP,P(CP),P(CP+1),P1',CP,P(CP),P(CP+1),P1
C            PRINT*,'CT,T(CT),T(CT+1),T1',CT,T(CT),T(CT+1),T1
           
 	    VT(LAYER)=(P1-P(CP))/(P(CP+1)-P(CP))
	    UT(LAYER)=(T1-T(CT))/(T(CT+1)-T(CT))

 	   IOFF(LAYER,1) = (nt * ng * (cp-1)) + (ng * (ct-1))
	   IOFF(LAYER,2) = (nt * ng * cp) + (ng * (ct-1))
	   IOFF(LAYER,3) = (nt * ng * cp) + (ng * ct)
	   IOFF(LAYER,4) = (nt * ng * (cp-1)) + (ng * ct)

          ENDIF

C         IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	  U=UT(LAYER)
	  V=VT(LAYER)
	  I1 = IOFF(LAYER,1)
	  I2 = IOFF(LAYER,2)
	  I3 = IOFF(LAYER,3)
	  I4 = IOFF(LAYER,4)


C         If wavenumber lines up with table, just read in the k-data for
C                                                     that wavenumber:
          IF (COINC) THEN
C           print*,'Wavenumber coincides with tabulated value'

           DO I = 1, NG

                Y1=TABLE(I1+1)
                Y2=TABLE(I2+1)
                Y3=TABLE(I3+1)
                Y4=TABLE(I4+1)


                IF((Y1.LE.0.0).OR.(Y2.LE.0.0).OR.
     &		   (Y3.LE.0.0).OR.(Y4.LE.0.0))THEN

                        Y1=TABLE(I1+I)
                        Y2=TABLE(I2+I)
                        Y3=TABLE(I3+I)
                        Y4=TABLE(I4+I)
		        Z = -1.0
                ELSE
			Y1 = TABLE(I1+I)
			Y2 = TABLE(I2+I)
			Y3 = TABLE(I3+I)
			Y4 = TABLE(I4+I)
                        Y1=LOG(Y1)
                        Y2=LOG(Y2)
                        Y3=LOG(Y3)
                        Y4=LOG(Y4)
			Z = 1.0
                ENDIF

		X=(1.0-V)*(1.0-U)*Y1 + V*(1.0-U)*Y2 + 
     1			V*U*Y3 + (1.0-V)*U*Y4

		IF(Z.LE.0.0)THEN
			K_G(I)=X
		ELSE
			K_G(I)=EXP(X)
		ENDIF
           ENDDO

          ELSE
C           print*,'Wavenumber does not coincide with tabulated value'
           weight(2) = ((vwave - tmp)/delv)
           if(weight(2).gt.1.0)weight(2)=1.0
           if(weight(2).lt.0.0)weight(2)=0.0
           weight(1) = 1. - weight(2)

           NTAB = NP*NT*NG
           count = 0
           do loop = 1, 2
                DO I = 1, NG
                 if (loop.eq.1) then
                         Y1=TABLE(I1+I)
                         Y2=TABLE(I2+I)
                         Y3=TABLE(I3+I)
                         Y4=TABLE(I4+I)
                 else
                         Y1=TABLE2(I1+I)
                         Y2=TABLE2(I2+I)
                         Y3=TABLE2(I3+I)
                         Y4=TABLE2(I4+I)
		 endif

                 K_g(I)=(1.0-V)*(1.0-U)*Y1 + V*(1.0-U)
     1              *Y2 + V*U*Y3 + (1.0-V)*U*Y4
                 count = count + 1
                 cont(count) = k_g(I)
                 fac(count) = delg(I) * weight(loop)
                ENDDO
           enddo

           sum = 0.
           do I = 1, count
                sum = sum + fac(I)
           enddo
           call rank (delg, ng, cont, fac, count, k_g)

         ENDIF


         DO I=1,NG
           KOUT(LAYER,IGAS,I)=K_G(I)
C           print*,K_G(I)
         ENDDO

1050     CONTINUE

	 INTERP = .FALSE.

999      CONTINUE

1000    CONTINUE

C        print*,'kout',kout(1,1,1)
	RETURN

	END

************************************************************************
************************************************************************
