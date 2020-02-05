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

        INCLUDE '../includes/arrdef.f'

	INTEGER	NP, NT, NG, CP, CT, I, I1, I2, I3, I4, N1,IREC
	INTEGER IWAVE,NGAS,NLAYER,IGAS,LAYER,LUN0,IREC0,maxc
        INTEGER NTAB,loop,count,MTAB,J,IRECX,CT1,CT2, NX
        LOGICAL NTEST,ISNAN
	LOGICAL COINC
        parameter (maxc=2*maxg,MTAB=maxk*maxk*maxg)
        REAL TABLE(MTAB),TABLE2(MTAB),X1,X2,U2

	REAL	P1, T1,DELV,VMIN, tmp, eps, KTEST,
     1		Y1, Y2, Y3, Y4, U, V, pmax, pmin, tmax, tmin, X, Z
        real weight(2),fac(maxc),cont(maxc),sum
        REAL PRESS(NLAYER),TEMP(NLAYER),VWAVE,FRAC

C       Defines the maximum values for a series of variables (layers, 
C       bins, paths, etc.)
        integer lun(maxbin,maxgas), ireck(maxbin,maxgas)
        real xmin(maxbin,maxgas), delx(maxbin,maxgas)
	REAL KOUT(MAXLAY,MAXGAS,MAXG),K_G(MAXG),G_ORD(MAXG)
        REAL DELG(MAXG),fracx(maxbin,maxgas)
        REAL dkoutdt(maxlay,maxgas,maxg)
        REAL P(MAXK),T(MAXK),T2(MAXK,MAXK),TN(MAXK)
	REAL UT(MAXLAY),VT(MAXLAY),FWHMK,DELVK,UT2(MAXLAY)
        INTEGER IOFF(MAXLAY,4)
        LOGICAL INTERP

	COMMON /INTERPK/LUN,IRECK,XMIN,DELX,FRACX,P,NP,T,T2,NT,
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

        INTERP = .TRUE.

        DO 1000 IGAS=1,NGAS

         LUN0 = LUN(IWAVE,IGAS)
         VMIN = XMIN(IWAVE,IGAS)
         DELV = DELX(IWAVE,IGAS)
         FRAC = FRACX(IWAVE,IGAS)
         IREC0 = IRECK(IWAVE,IGAS)

C         print*,'IGAS,DELV,IWAVE,VMIN',IGAS,DELV,IWAVE,VMIN
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
C         print*,'delv,vmin,vwave = ',delv,vmin,vwave


         IF(delv.gt.0.0)then

C         Set minimum closeness to tabulated wavenumbers to be 1/50 of 
C         the separation
          eps = 0.02*delv

C         find nearest point in table below current wavelength/wavenumber
          N1 = int((VWAVE + EPS - VMIN)/DELV)
C	   print*,vmin,vwave,delv,(vwave-vmin)/delv,eps,N1
C						  below current wavelength.
C						  Parameter eps is there 
C						  to prevent small numerical
C						  errors in VWAVE screwing
C					          things up between platforms
          IREC = IREC0+NP*ABS(NT)*NG*N1
C          print*,'IREC = ',IREC
         ELSE       
C         For irregularly gridded tables IREC0 is assumed to hold the 
C         nearest record number in the table to the requested wavelength, not
C         that of the start of the table. This has already been allocated
C         by read_klist.f
          irec = irec0 
          N1=-1   
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

202      CONTINUE
         IREC=IRECX

         IF(KTEST.EQ.0.0)THEN
c          print*,'GET_K: Zero k-data for GAS: ',IGAS
c          print*,'Wavelength/waveumber = ',VWAVE
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
          frac = (vwave-tmp)/delv
          if(frac.lt.eps)then
C          delv > 0 and requested wavenumber close enough to tabulated to
C          be considered coincident (and thus need to interrogate k-table once
C          no twice.
           COINC=.TRUE.
          else
           COINC=.FALSE.
          endif

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


         NTAB = ABS(NT)*NP*NG
         if(NTAB.gt.MTAB)then
          print*,'Error in get_kg, NTAB>MTAB'
          print*,NTAB,MTAB
          stop
         endif

         DO I=1,NTAB
          READ(LUN0,REC=IREC)TABLE(I)
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
C          print*,'LAYER = ',LAYER,PRESS(LAYER),
C     1    TEMP(LAYER)

C        IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C           Work out where P,T for each layer lies in the tables
C        ------------------------------------
	  IF(INTERP)THEN
C	    First check range is fine.
            P1 = LOG(PRESS(LAYER))

C	    Find position of temp and pressure values in k table, then
C	    interpolate
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

              IOFF(LAYER,1) = (nt * ng * (cp-1)) + (ng * (ct-1))
              IOFF(LAYER,2) = (nt * ng * cp) + (ng * (ct-1))
              IOFF(LAYER,3) = (nt * ng * cp) + (ng * ct)
              IOFF(LAYER,4) = (nt * ng * (cp-1)) + (ng * ct)

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
              IOFF(LAYER,1) = (nx * ng * (cp-1)) + (ng * (ct1-1))
              IOFF(LAYER,2) = (nx * ng * cp) + (ng * (ct2-1))
              IOFF(LAYER,3) = (nx * ng * cp) + (ng * ct2)
              IOFF(LAYER,4) = (nx * ng * (cp-1)) + (ng * ct1)

            ENDIF            


          ENDIF

C         IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	  U=UT(LAYER)
	  U2=UT2(LAYER)
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
                        IF(Y1.LT.1E-37)Y1=1E-37
                        IF(Y2.LT.1E-37)Y2=1E-37
                        IF(Y3.LT.1E-37)Y3=1E-37
                        IF(Y4.LT.1E-37)Y4=1E-37
                        Y1=LOG(Y2)
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
			K_G(I)=X
		ELSE
			K_G(I)=EXP(X)
		ENDIF
           ENDDO

          ELSE
C           print*,'Wavenumber does not coincide with tabulated value'
           if(delv.gt.0)then
             weight(2) = ((vwave - tmp)/delv)
             if(weight(2).gt.1.0)weight(2)=1.0
             if(weight(2).lt.0.0)weight(2)=0.0
             weight(1) = 1. - weight(2)
C             print*,'weight',vwave,tmp,weight(1),weight(2)
           else
             weight(1)=frac
             weight(2)=1.0-frac

           endif
C           print*,'vwave,Weight : ',vwave,weight(1),weight(2)
           NTAB = NP*ABS(NT)*NG
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
                 IF(NT.GT.0)THEN
                  K_g(I)=(1.0-V)*(1.0-U)*Y1 + V*(1.0-U)
     1              *Y2 + V*U*Y3 + (1.0-V)*U*Y4
                 ELSE
                  X1=(1.0-U)*Y1 + U*Y4
                  X2=(1.0-U2)*Y2 + U2*Y3
                  K_G(I) = (1.0-V)*X1 + V*X2
                 ENDIF

                 count = count + 1
                 cont(count) = k_g(I)
                 fac(count) = delg(I) * weight(loop)
C                 print*,count,cont(count),fac(count)
                ENDDO
           enddo

           sum = 0.
           do I = 1, count
                sum = sum + fac(I)
           enddo
C           print*,'sum',sum

           if(weight(1).gt.0.and.weight(2).gt.0)then
             call rank (delg, ng, cont, fac, count, k_g)
           else
             do i=1,ng
              if(weight(2).eq.0)then
                k_g(i)=cont(i)
              else
                k_g(i)=cont(i+ng)
              endif
             enddo
           endif
         ENDIF

         DO I=1,NG
C           print*,igas,K_G(I)
           NTEST=ISNAN(K_G(I))
           IF(NTEST)THEN
            KOUT(LAYER,IGAS,I)=1e-37
            PRINT*,'Warning, NAN returned by get_k.f for gas',igas
            PRINT*,'Setting to zero here, but should check files'
            PRINT*,'VWAVE, IG = ',vwave,I
           ELSE
            KOUT(LAYER,IGAS,I)=K_G(I)
           ENDIF
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
