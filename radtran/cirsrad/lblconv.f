************************************************************************
************************************************************************
C-----------------------------------------------------------------------
C
C			SUBROUTINE LBLCONV
C
C	Convolves input spectrum (vwave, y) with a bin of width
C	fwhm to produce an output spectrum (vconv, yout).
C
C-----------------------------------------------------------------------

	SUBROUTINE lblconv(runname,fwhm, ishape, npath, ispace, 
     1     vwave, delv, y, nconv, vconv,yout,ynor,FWHMEXIST,NFWHM,
     2     VFWHM,XFWHM,
     3     FINSTEXIST,minst,ninst,vinst,finst)

	IMPLICIT NONE

C       Defines the maximum values for a series of variables (layers,
C       bins, paths, etc.)
	INCLUDE '../includes/arrdef.f'

	REAL		fwhm,x,sig
	INTEGER		nstep,mconv,npath,ispace,ishape,NFWHM
	PARAMETER	(nstep=20, mconv=2000)
        REAL		VFWHM(NFWHM),XFWHM(NFWHM),YFWHM
        LOGICAL		FWHMEXIST

        INTEGER		nconv, nc, I, J,nconv1,nsub,k
	REAL		y(maxpat,2),
     1			yout(maxpat,mconv), ynor(maxpat,mconv),
     2			y2(maxbin), x1, x2, delx, xi, yi, yold, dv
        REAL*8		vwave,vconv(nconv),vwave1,vcentral,v1,v2,vcen,
     1			  vt
	REAL		vfil(1000),fil(1000),yy,sumf,delv
	REAL		ytmp(maxout)
	REAL		f1,f2,XOFF,HAMMING,NFW,HANNING
	CHARACTER*100	runname
c  ** instrument function (should be defined if ishape=5) **
        integer	minst
        logical     FINSTEXIST
        integer	ninst
        real		vinst(minst),finst(minst)
        integer idiag,iquiet
        common/diagnostic/idiag,iquiet

C        print*,'LBLCONV --> FWHM = ',FWHM


	IF((ishape.eq.5).and.FINSTEXIST)then
c  **     apply an instrument function read from .fin file in 
c  **	    lblrtf_wave.f. If a .fin file doesn't exist then ignore
c  **	    this block and move to next option

C        vwave, delv are in wavenumbers

         vwave1=vwave-delv

         delx=delv
         if(ispace.eq.1)then
          delx=1e4*sngl(delv/(vwave*vwave))
         endif

C         print*,'delv,delx = ',delv,delx

         DO J=1,NCONV

C         Find limits of instrument function in wavenumbers
          V1=VCONV(J)+vinst(1)
          V2=VCONV(J)+vinst(ninst)
          VCEN=VCONV(J)

c	    convert to wavelength if ispace=1
          IF(ISPACE.EQ.1)THEN
            VT=V1
            V1=1e4/V2
            V2=1e4/VT
          ENDIF

          F1=0.0
          F2=0.0
          IF(VWAVE.GE.V1.AND.VWAVE.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
              xoff = SNGL(VWAVE-VCEN)
            ELSE
              xoff = 1E4/SNGL(VWAVE-VCEN)
            ENDIF
            call verint(vinst,finst,ninst,F2,xoff)
          ENDIF
          IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
              xoff = SNGL(VWAVE1-VCEN)
            ELSE
              xoff = 1E4/SNGL(VWAVE1-VCEN)
            ENDIF
            call verint(vinst,finst,ninst,F1,xoff)
          ENDIF
          IF(F2.LT.0.0)THEN
C            if(idiag.gt.0)then
C             print*,'F2 gone slightly negative. Fixing at 0.',F2
C            endif
            F2=0.0
          ENDIF
          IF(F1.LT.0.0)THEN
C            if(idiag.gt.0)then
C             print*,'F1 gone slightly negative. Fixing at 0.',F1
C            endif
            F1=0.0
          ENDIF

          IF(F1.GT.0.0.OR.F2.GT.0.0)THEN 
C           print*,J,V1,VWAVE1,VWAVE,V2,F1,F2
           DO I=1,NPATH
            YOUT(I,J)=YOUT(I,J)+0.5*(F1*Y(I,1)+F2*Y(I,2))*DELX
            YNOR(I,J)=YNOR(I,J)+0.5*(F1+F2)*DELX
C            print*,'BB'
C            print*,I,Y(I,1),Y(I,2),0.5*(F1*Y(I,1)+F2*Y(I,2))*DELX
C            print*,0.5*(F1+F2)*DELX
C            print*,J,YOUT(I,J),YNOR(I,J)
           ENDDO
          ENDIF
         ENDDO

        
        ELSEIF(fwhm.gt.0.0)THEN

C         print*,'EE',vwave,delv,fwhm,npath,ishape,ispace
C         print*,runname
C         print*,nconv,(vconv(j),j=1,nconv)
C         do i=1,npath
C          print*,y(i,1),y(i,2)
C         enddo

C        vwave, delv are in wavenumbers

         vwave1=vwave-delv

         delx=delv
         if(ispace.eq.1)then
          delx=1e4*sngl(delv/(vwave*vwave))
         endif

C         print*,'CC',vwave,vwave1,delv,delx

C        Set total width of Hamming/Hanning function window in terms of
C        numbers of FWHMs for ISHAPE=3 and ISHAPE=4
         NFW = 3.

         DO J=1,NCONV

          YFWHM=FWHM
          if(fwhmexist)then
             call verint(vfwhm,xfwhm,nfwhm,yfwhm,sngl(vconv(j)))
          endif
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

          IF(ISPACE.EQ.1)THEN
            VT=V1
            V1=1e4/V2
            V2=1e4/VT
          ENDIF

          F1=0.0
          F2=0.0

          IF(ISHAPE.EQ.0)THEN
C          Square Instrument Shape
           IF(VWAVE.GE.V1.AND.VWAVE.LE.V2)F2=1.0
           IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)F1=1.0

          ELSEIF(ISHAPE.EQ.1)THEN
C          Triangular Instrument Shape
           IF(VWAVE.GE.V1.AND.VWAVE.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             F2=1.0 - ABS(SNGL(VWAVE-VCEN))/YFWHM
            ELSE
             F2=1.0 - ABS(SNGL(1E4/VWAVE-VCEN))/YFWHM          
            ENDIF
           ENDIF

           IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             F1=1.0 - ABS(SNGL(VWAVE1-VCEN))/YFWHM
            ELSE
             F1=1.0 - ABS(1E4/SNGL(VWAVE1-VCEN))/YFWHM
            ENDIF
           ENDIF           
           IF(F2.LT.0.0)THEN
C            if(idiag.gt.0)then
C             print*,'F2 gone slightly negative. Fixing at 0.',F2
C            endif
            F2=0.0
           ENDIF
           IF(F1.LT.0.0)THEN
C            if(idiag.gt.0)then
C             print*,'F1 gone slightly negative. Fixing at 0.',F1
C            endif
            F1=0.0
           ENDIF

          ELSEIF(ISHAPE.EQ.2)THEN
C          Gaussian Instrument Shape
           IF(VWAVE.GE.V1.AND.VWAVE.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             F2=EXP(-(SNGL(VWAVE-VCEN)/SIG)**2)
            ELSE
             F2=EXP(-(SNGL(1E4/VWAVE-VCEN)/SIG)**2)
            ENDIF
           ENDIF
           IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             F1=EXP(-(SNGL(VWAVE1-VCEN)/SIG)**2)
            ELSE
             F1=EXP(-(SNGL(1E4/VWAVE1-VCEN)/SIG)**2)
            ENDIF
           ENDIF           

          ELSEIF(ISHAPE.EQ.3)THEN
C          Hamming Instrument Shape
           IF(VWAVE.GE.V1.AND.VWAVE.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             XOFF = SNGL(VWAVE-VCEN)
            ELSE
             XOFF = SNGL(1E4/VWAVE - VCEN)
            ENDIF
            F2 = HAMMING(YFWHM,XOFF)
           ENDIF

           IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             XOFF = SNGL(VWAVE1-VCEN)
            ELSE
             XOFF = SNGL(1E4/VWAVE1-VCEN)
            ENDIF
            F1 = HAMMING(YFWHM,XOFF)
           ENDIF           

          ELSE
C          Hanning Instrument Shape
           IF(VWAVE.GE.V1.AND.VWAVE.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             XOFF = SNGL(VWAVE-VCEN)
            ELSE
             XOFF = SNGL(1E4/VWAVE - VCEN)
            ENDIF
            F2 = HANNING(YFWHM,XOFF)
           ENDIF

           IF(VWAVE1.GE.V1.AND.VWAVE1.LE.V2)THEN
            IF(ISPACE.EQ.0)THEN
             XOFF = SNGL(VWAVE1-VCEN)
            ELSE
             XOFF = SNGL(1E4/VWAVE1-VCEN)
            ENDIF
            F1 = HANNING(YFWHM,XOFF)
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
           DO I=1,NPATH
            YOUT(I,J)=YOUT(I,J)+0.5*(F1*Y(I,1)+F2*Y(I,2))*DELX
            YNOR(I,J)=YNOR(I,J)+0.5*(F1+F2)*DELX
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
            dv = 100*sngl(abs(vcentral-vconv(j))/vconv(j))
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
C-----------------------------------------------------------------------
C
C	Return and end.
C
C-----------------------------------------------------------------------

10	RETURN

	END

************************************************************************
************************************************************************
