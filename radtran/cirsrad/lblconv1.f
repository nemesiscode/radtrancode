************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE LBLCONV1
C
C	Convolves input spectrum (vwave, y) with a bin of width
C	fwhm to produce an output spectrum (vconv, yout).
C       Operates on a single path.
C
C-----------------------------------------------------------------------

	SUBROUTINE lblconv1(runname,fwhm, ishape,nwave, vwave, y, 
     1    nconv, vconv,yout,FWHMEXIST,NFWHM,VFWHM,XFWHM)

	IMPLICIT NONE

C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)
	INCLUDE '../includes/arrdef.f'

	REAL		fwhm
	INTEGER		nstep,ishape,NFWHM
	PARAMETER	(nstep=20)
        REAL            VFWHM(NFWHM),XFWHM(NFWHM),YFWHM
        LOGICAL		FWHMEXIST
        INTEGER		nwave, nconv, nc, I, J,nconv1,nsub,k,nc1
        INTEGER		I1,I2
	REAL		vwave(nwave), y(maxout), vconv(nconv),
     1			yout(maxout), ynor(maxout), f1
	REAL		vfil(1000),fil(1000),yy,delv,hanning
        REAL		hamming,dv,xoff,NFW,V1,V2,SIG,VCEN
	REAL		vcentral
	CHARACTER*100	runname
C-----------------------------------------------------------------------

C        print*,'LBLCONV1 --> FWHM = ',FWHM

        IF(fwhm.gt.0.0)THEN

C        Set total width of Hamming/Hanning function window in terms of
C        numbers of FWHMs for ISHAPE=3 and ISHAPE=4
         NFW = 3.

    
 
         DO 101 J=1,NCONV
          YFWHM=FWHM
          IF(FWHMEXIST)THEN
           CALL VERINT(VFWHM,XFWHM,NFWHM,YFWHM,VCONV(J))
          ENDIF
C          print*,'J,VCONV(J),YFWHM',J,VCONV(J),YFWHM
C         Find limits of instrument width in wavenumbers
          IF(ISHAPE.EQ.0)THEN
           V1=VCONV(J)-0.5*YFWHM
           V2=V1+YFWHM
          ELSEIF(ISHAPE.EQ.1)THEN
           V1=VCONV(J)-YFWHM
           V2=VCONV(J)+YFWHM
          ELSEIF(ISHAPE.EQ.2)THEN
           SIG = 0.5*YFWHM/SQRT(ALOG(2.0))
           V1=VCONV(J)-3.*SIG
           V2=VCONV(J)+3.*SIG
          ELSE
           V1=VCONV(J)-NFW*YFWHM
           V2=VCONV(J)+NFW*YFWHM
          ENDIF
          VCEN=VCONV(J)

C         Find relevant points in tabulated files.
          I1=-1
          I2=-1
          DO I=1,NWAVE
           IF(VWAVE(I).GE.V1.AND.I1.LT.0)THEN
            I1=I
           ENDIF
           IF(VWAVE(I).LE.V2)THEN
            I2=I
           ENDIF
          ENDDO

          IF(I1.LT.0.OR.I2.LT.0)THEN
           PRINT*,'Error in lblconv1 - wavelength/wavenumber not'
           PRINT*,'covered by lbltables'
           STOP         
          ENDIF 

          YOUT(J)=0.
          YNOR(J)=0.

          DO 102 I=I1,I2
           F1=0.
           IF(ISHAPE.EQ.0)THEN
C           Square Instrument Shape
            F1=1.0
           ELSEIF(ISHAPE.EQ.1)THEN
C           Triangular Instrument Shape
            F1=1.0 - ABS(VWAVE(I)-VCEN)/YFWHM
           ELSEIF(ISHAPE.EQ.2)THEN
C           Gaussian Instrument Shape
            F1=EXP(-((VWAVE(I)-VCEN)/SIG)**2)
           ELSEIF(ISHAPE.EQ.3)THEN
C           Hamming Instrument Shape
            XOFF = VWAVE(I)-VCEN
            F1 = HAMMING(YFWHM,XOFF)
           ELSEIF(ISHAPE.EQ.4)THEN
C           Hanning Instrument Shape
            XOFF = VWAVE(I)-VCEN
            F1 = HANNING(YFWHM,XOFF)
           ELSE
            F1=0.
           ENDIF

           IF(ISHAPE.GT.0.AND.ISHAPE.LT.3)THEN
            IF(F1.LT.0.0)THEN
             print*,'Warning: lblconv1. F1 < 0.0'
             print*,F1
             stop
            ENDIF
           ENDIF

           IF(F1.GT.0.0)THEN
             YOUT(J)=YOUT(J)+F1*Y(I)
             YNOR(J)=YNOR(J)+F1
           ENDIF

102       CONTINUE

          YOUT(J)=YOUT(J)/YNOR(J)

101      CONTINUE

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

           do 205 J=1,NCONV
C           Make sure you're using the right filter function for the
C           channel requested.
            dv = 100*abs(vcentral-vconv(j))/vconv(j)
c            if(dv.lt.1.0)then
            if(dv.lt.0.0001)then
             V1 = VFIL(1)
             V2 = VFIL(NSUB)

C            Find relevant points in tabulated files.
             I1=-1
             I2=-1
             DO I=1,NWAVE
              IF(VWAVE(I).GE.V1.AND.I1.LT.0)THEN
               I1=I
              ENDIF
              IF(VWAVE(I).LE.V2)THEN
               I2=I
              ENDIF
             ENDDO

             IF(I1.LT.0.OR.I2.LT.0)THEN
              PRINT*,'Error in lblconv1 - wavelength/wavenumber not'
              PRINT*,'covered by lbltables'
              STOP         
             ENDIF 

             YOUT(J)=0.
             YNOR(J)=0.

             DO 202 I=I1,I2

              CALL interp(vfil,fil,nsub,f1,vwave(I))

              IF(F1.GT.0.0)THEN
                YOUT(J)=YOUT(J)+F1*Y(I)
                YNOR(J)=YNOR(J)+F1
              ENDIF

202          CONTINUE

             YOUT(J)=YOUT(J)/YNOR(J)

            endif

205        continue

200       CONTINUE

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
