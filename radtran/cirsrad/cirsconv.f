************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE CIRSCONV
C
C	Convolves input spectrum (vwave, y) with a bin of width
C	fwhm to produce an output spectrum (vconv, yout).
C
C-----------------------------------------------------------------------

	SUBROUTINE cirsconv(runname,fwhm, nwave, vwave, y, 
     1    nconv, vconv,yout,FWHMEXIST,NFWHM,VFWHM,XFWHM)

	IMPLICIT NONE

C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)
	INCLUDE '../includes/arrdef.f'

	REAL		fwhm
	INTEGER		nstep,NFWHM
	PARAMETER	(nstep=20)
        REAL		VFWHM(NFWHM),XFWHM(NFWHM),YFWHM
        INTEGER		nwave, nconv, nc, I, J,nconv1,nsub,k,nc1
        LOGICAL		FLAGNAN,FWHMEXIST
	REAL		vwave(nwave), y(maxout), vconv(nconv),
     1			yout(maxout), xc(maxbin), yc(maxbin),
     2			y2(maxbin), x1, x2, delx, xi, dv, y1,
     3 			xc1(maxbin),yc1(maxbin)
	REAL		vfil(1000),fil(1000),yy,delv
	REAL		vcentral,ytmp(maxout)
        DOUBLE PRECISION sum,sumf,yi,yold
	CHARACTER*100	runname
        LOGICAL		xnorm
        integer idiag,iquiet
        common/diagnostic/idiag,iquiet
C-----------------------------------------------------------------------
C
C	Do a quick check first. Extrapolate end points if needed.
C
C-----------------------------------------------------------------------

        xnorm = .true.

C       set the FWHM to be the same as the spacing of k-distribs
C       in look-up table


	IF (nwave.eq.1.and.fwhm.ne.0.0) THEN
	 if(idiag.gt.0)then
           WRITE(*,*)'CIRSCONV: Too few input points: nwave= ',nwave
	 endif
         	nconv = nwave
		vconv(1) = vwave(1)
		yout(1) = 1
		GOTO 10
	ENDIF

	nc = nwave
	DO I = 1, nc
		xc(I) = vwave(I)
		yc(I) = y(I)
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

          IF(FLAGNAN.AND.IDIAG.GT.0)THEN
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
                yfwhm=fwhm
                if(fwhmexist)then
                 call verint(vfwhm,xfwhm,nfwhm,yfwhm,vconv(i))
                endif

		x1 = vconv(I) - yfwhm/2.
		x2 = vconv(I) + yfwhm/2.

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
            if(dv.lt.0.00001)then
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


C	      Uncomment following line to use WASP43-b hybrid k-tables.
C             xnorm=.false.

             if(xnorm)then
              do j=1,nwave-1
 
               if(ytmp(j).ne.0)then
                delv = vwave(j+1)-vwave(j)
                sum=sum+ytmp(j)*y(j)*delv
                sumf=sumf+ytmp(j)*delv
               endif
              enddo
             else
              do j=1,nwave
 
               if(ytmp(j).ne.0)then
                sum=sum+ytmp(j)*y(j)
                sumf=sumf+ytmp(j)
               endif
              enddo

             endif

             yout(I)=sngl(sum/sumf)
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
