************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE LBLCONV1
C
C	Convolves input spectrum (vwave, y) with a bin of width
C	fwhm to produce an output spectrum (vconv, yout).
C
C-----------------------------------------------------------------------

	SUBROUTINE lblconv1(runname,fwhm, ishape,nwave, vwave, y, 
     1    nconv, vconv,yout)

	IMPLICIT NONE

C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)
	INCLUDE '../includes/arrdef.f'

	REAL		fwhm
	INTEGER		nstep,ishape
	PARAMETER	(nstep=20)

        INTEGER		nwave, nconv, nc, I, J,nconv1,nsub,k,nc1
        LOGICAL		FLAGNAN
	REAL		vwave(nwave), y(maxout), vconv(nconv),
     1			yout(maxout), xc(maxbin), yc(maxbin),
     2			y2(maxbin), x1, x2, delx, xi, dv, y1,
     3 			xc1(maxbin),yc1(maxbin)
	REAL		vfil(1000),fil(1000),yy,delv
	REAL		vcentral,ytmp(maxout)
        DOUBLE PRECISION sum,sumf,yi,yold
	CHARACTER*100	runname
C-----------------------------------------------------------------------


        IF(fwhm.gt.0.0)THEN

C        Set total width of Hamming/Hanning function window in terms of
C        numbers of FWHMs for ISHAPE=3 and ISHAPE=4
         NFW = 3.

         DO J=1,NCONV
C         Find limits of instrument width in wavenumbers
          IF(ISHAPE.EQ.0)THEN
           V1=VCONV(J)-0.5*FWHM
           V2=V1+FWHM
          ELSEIF(ISHAPE.EQ.1)THEN
           V1=VCONV(J)-FWHM
           V2=VCONV(J)+FWHM
          ELSEIF(ISHAPE.EQ.2)THEN
           SIG = 0.5*FWHM/SQRT(ALOG(2.0))
           V1=VCONV(J)-3.*SIG
           V2=VCONV(J)+3.*SIG
          ELSE
           V1=VCONV(J)-NFW*FWHM
           V2=VCONV(J)+NFW*FWHM
          ENDIF
          VCEN=VCONV(J)

C         Find relevant points in tabulated files.
          V1A=-1.
          DO I=1,NWAVE
           IF(VWAVE(I).GE.V1.AND.V1A.LT.0)THEN
            I1=I
            V1A=VWAVE(I)
           ENDIF
           IF(VWAVE(I).LE.V2)THEN
            I2=I
            V2A=VWAVE(I)
           ENDIF
          ENDIF
 
          DO I=I1,I2
           F1=0.
           IF(ISHAPE.EQ.0)THEN
C           Square Instrument Shape
            IF(VWAVE(I).GE.V1.AND.VWAVE(I).LE.V2)F1=1.0
           ELSEIF(ISHAPE.EQ.1)THEN
C           Triangular Instrument Shape
            IF(VWAVE(I).GE.V1.AND.VWAVE(I).LE.V2)THEN
             F1=1.0 - ABS(VWAVE-VCEN)/FWHM
            ENDIF
           ENDIF

           IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             F1=1.0 - ABS(VWAVE1-VCEN)/FWHM
            ELSE
             F1=1.0 - ABS(1E4/VWAVE1-VCEN)/FWHM
            ENDIF
           ENDIF
           IF(F2.LT.0.0)THEN
            PRINT*,'F2 gone slightly negative. Fixing at 0.',F2
            F2=0.0
           ENDIF
           IF(F1.LT.0.0)THEN
            PRINT*,'F1 gone slightly negative. Fixing at 0.',F1
            F1=0.0
           ENDIF

          ELSEIF(ISHAPE.EQ.2)THEN
C          Gaussian Instrument Shape
           IF(VWAVE.GE.V1.AND.VWAVE.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             F2=EXP(-((VWAVE-VCEN)/SIG)**2)
            ELSE
             F2=EXP(-((1E4/VWAVE-VCEN)/SIG)**2)
            ENDIF
           ENDIF
           IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             F1=EXP(-((VWAVE1-VCEN)/SIG)**2)
            ELSE
             F1=EXP(-((1E4/VWAVE1-VCEN)/SIG)**2)
            ENDIF
           ENDIF

          ELSEIF(ISHAPE.EQ.3)THEN
C          Hamming Instrument Shape
           IF(VWAVE.GE.V1.AND.VWAVE.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             XOFF = VWAVE-VCEN
            ELSE
             XOFF = 1E4/VWAVE - VCEN
            ENDIF
            F2 = HAMMING(FWHM,XOFF)
          ENDIF

           IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             XOFF = VWAVE1-VCEN
            ELSE
             XOFF = 1E4/VWAVE1-VCEN
            ENDIF
            F1 = HAMMING(FWHM,XOFF)
           ENDIF

          ELSE
C          Hanning Instrument Shape
           IF(VWAVE.GE.V1.AND.VWAVE.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             XOFF = VWAVE-VCEN
            ELSE
             XOFF = 1E4/VWAVE - VCEN
            ENDIF
            F2 = HANNING(FWHM,XOFF)
           ENDIF

           IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             XOFF = VWAVE1-VCEN
            ELSE
             XOFF = 1E4/VWAVE1-VCEN
            ENDIF
            F1 = HANNING(FWHM,XOFF)
           ENDIF

          ENDIF

          IF(ISHAPE.GT.0.AND.ISHAPE.LT.3)THEN
           IF(F1.LT.0.0.OR.F2.LT.0.0)THEN
            print*,'Warning: lblconv. F1 or F2 < 0.0'
            print*,F1,F2
            stop
           ENDIF
          ENDIF

          IF(F1.GT.0.0.OR.F2.GT.0.0)THEN
C           print*,J,V1,VWAVE1,VWAVE,V2,F1,F2
           DO I=1,NPATH
C            print*,I,Y(I,1),Y(I,2),0.5*(F1*Y(I,1)+F2*Y(I,2))*DELX
C            print*,0.5*(F1+F2)*DELX
            YOUT(I,J)=YOUT(I,J)+0.5*(F1*Y(I,1)+F2*Y(I,2))*DELX
            YNOR(I,J)=YNOR(I,J)+0.5*(F1+F2)*DELX
C            print*,J,YOUT(I,J),YNOR(I,J)
           ENDDO
          ENDIF
         ENDDO

        ELSE
C         Channel Integrator Mode: Slightly more advanced than previous

          CALL FILE(runname,runname,'fil')
          OPEN(12,FILE=runname,STATUS='old')
          READ(12,*)NCONV1
          DO 200 K=1,nconv1
           READ(12,*)vcentral
           READ(12,*)nsub
           do j=1,nsub
            read(12,*)vfil(j),fil(j)
           enddo

           do 205 j=1,nconv
C           Make sure you're using the right filter function for the
C           channel requested.
            dv = 100*abs(vcentral-vconv(j))/vconv(j)
            if(dv.lt.1.0)then


              if((vwave-delv).ge.vfil(1).and.(vwave-delv).
     1           le.vfil(nsub))then
              CALL interp(vfil,fil,nsub,f1,vwave)
              else
               f1=0.0
              endif

              if(vwave.ge.vfil(1).and.vwave.le.vfil(nsub))then
               CALL interp(vfil,fil,nsub,f2,vwave)
              else
               f2=0.0
              endif
              IF(F1.GT.0.0.OR.F2.GT.0.0)THEN
               DO I=1,NPATH
                YOUT(I,J)=YOUT(I,J)+0.5*(F1*Y(I,1)+F2*Y(I,2))*DELV
                YNOR(I,J)=YNOR(I,J)+0.5*(F1+F2)*DELV
               ENDDO
              ENDIF
            endif

205        continue

200       continue

          close(12)

        ENDIF




C       set the FWHM to be the same as the spacing of k-distribs
C       in look-up table

C	WRITE(*,*)'CIRSCONV: FWHM of boxcar = ',FWHM

	IF (nwave.eq.1.and.fwhm.ne.0.0) THEN
		WRITE(*,*)'CIRSCONV: Too few input points: nwave= ',nwave
		nconv = nwave
		vconv(1) = vwave(1)
		yout(1) = 1
		GOTO 10
	ENDIF

	nc = nwave
	DO I = 1, nc
		xc(I) = vwave(I)
		yc(I) = y(I)
C                print*,'raw',i,xc(I),yc(I)
	ENDDO

C-----------------------------------------------------------------------
C
C	Extrapolate if vwave range does not cover vconv range
C
C-----------------------------------------------------------------------
        IF(fwhm.gt.0.0)THEN

 	   IF(vconv(nconv).gt.vwave(nwave)-fwhm/2.)THEN
C		WRITE(*,*)'CIRSCONV: having to extrapolate vwave up'
C		WRITE(*,*)'vconv(nconv),vwave(nwave)-0.5*fwhm',
C     1			vconv(nconv),vwave(nwave)-0.5*fwhm

		nc = nc + 1
		xc(nc) = vconv(nconv) + fwhm
		yc(nc) = y(nwave) + ((y(nwave) - y(nwave-1)) /
     1			(vwave(nwave) - vwave(nwave-1))
     2			* (xc(nc) - vwave(nwave)))
	   ENDIF

	  IF (vconv(1).lt.vwave(1)+fwhm/2.) THEN
C		WRITE(*,*)'CIRSCONV: having to extrapolate vwave down' 
C		WRITE(*,*)'vconv(1),vwave(1)+0.5*fwhm',vconv(1),
C     1			vwave(1)+0.5*fwhm

		nc = nc + 1
		DO I = nc, 2, -1
			xc(I) = xc(I-1)
			yc(I) = yc(I-1)
		ENDDO
		xc(1) = vconv(1) - fwhm
		yc(1) = y(1) + ((y(2) - y(1))/(vwave(2) - vwave(1))
     1			* (xc(1) - vwave(1)))
	  ENDIF

C-----------------------------------------------------------------------
C
C	Robust integrator based on brute force.
C
C-----------------------------------------------------------------------

C         Check to make sure spectrum has no NaN's
          nc1=0
          do i=1,nc
           if(.not.isnan(yc(i)))then
            nc1=nc1+1
            xc1(nc1)=xc(i)
            yc1(nc1)=yc(i)
           endif
          enddo

          FLAGNAN=.FALSE.
          if(nc1.lt.nc)FLAGNAN=.TRUE.

          IF(FLAGNAN)THEN
           print*,'Warning from cirsconv.f: Input spectrum contains'
           print*,'a NaN'
           do i=1,nc
            print*,i,xc(i),yc(i)
           enddo
          ENDIF

C         Delete the NaNs and fit output spectrum to remaining points
          nc=nc1
          do i=1,nc
           xc(i)=xc1(i)
           yc(i)=yc1(i)
          enddo

	  CALL cspline (xc, yc, nc, 5.e30, 5.e30, y2)

	  DO I = 1, nconv
		x1 = vconv(I) - fwhm/2.
		x2 = vconv(I) + fwhm/2.

		delx = (x2-x1)/FLOAT(nstep-1)
		DO J = 1, nstep
			xi = x1 + (J-1) * delx

			CALL csplint(xc, yc, y2, nc, xi, y1)
                        yi=dble(y1)
			IF (J.eq.1) THEN
				sum = 0.
			ELSE
				sum = sum + (yi + yold)*dble(delx/2.)
			ENDIF
			yold = yi
		ENDDO
		yout(I) = sngl(sum/fwhm)
	  ENDDO

     	ELSEIF(FWHM.EQ.0.0)THEN
C         Channel Integrator mode where the k-tables have been previously
C         tabulated INCLUDING the filter profile. In which case all we 
C         need do is just transfer the outputs

          DO I=1,nconv
            	yout(I)=y(I)
          ENDDO

        ELSE
C         Channel Integrator Mode: Slightly more advanced than previous

          CALL FILE(runname,runname,'fil')
          OPEN(12,FILE=runname,STATUS='old')
          READ(12,*)NCONV1
          DO 200 K=1,nconv1
           READ(12,*)vcentral
           READ(12,*)nsub
           do j=1,nsub
            read(12,*)vfil(j),fil(j)   
           enddo

           do 205 i=1,nconv
C           Make sure you're using the right filter function for the
C           channel requested.
C            dv = 100*abs(vcentral-vconv(i))/vconv(i)
           dv = abs(vcentral-vconv(i))
C           print*,'averaged consistency',i,vconv(i),vcentral
            if(dv.lt.0.0001)then
             do j=1,nwave
              if(vwave(j).ge.vfil(1).and.vwave(j).le.vfil(nsub))then
               CALL interp(vfil,fil,nsub,yy,vwave(j))
               ytmp(j)=yy
              else
               ytmp(j)=0.0
              endif
             enddo

             sumf = 0.0
             sum = 0.0



             do j=1,nwave-1

              if(ytmp(j).ne.0)then
               delv = vwave(j+1)-vwave(j)
               sum=sum+ytmp(j)*y(j)*delv
               sumf=sumf+ytmp(j)*delv
              else
               sum=sum
               sumf=sumf
              endif
             enddo

             yout(I)=sngl(sum/sumf)
C             print*,'nconv,vconv,yout',nconv,vconv(i),yout(I)
            endif

205        continue

200       continue

          close(12)

        ENDIF
C-----------------------------------------------------------------------
C
C	Return and end.
C
C-----------------------------------------------------------------------

10	RETURN

	END

************************************************************************
************************************************************************
