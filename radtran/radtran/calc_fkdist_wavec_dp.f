      SUBROUTINE CALC_FKDIST_WAVEC(IWAVE,VMIN,DELV,NPOINT,
     1 NFIL,VFIL,TFIL,G_ORD,DEL_G,K_G,NGMAX,NG)
C     $Id:
C***********************************************************************
C_TITL:	CALC_KDIST_WAVEC.f
C
C_DESC:	Subroutine to calculate the cumulative K-distribution in a band
C	by binning a pre-calculated lbl absorption coefficient spectrum. 
C	LBL spectrum is analysed according to the equation given by Lacis
C	and Oinas (1991):
C		f(k) = (1/(V2-V1))*SUM(ABS(dV/DK))
C	f(k) is then summed to give the cumulative k distribution. 
C
C_ARGS:	Input variables:
C	OUTPUT(MPOINT)	REAL	Supplied lbl spectrum. Variable is passed
C				via the /SPECTRUM/ common block for speed.
C       VMIN		REAL	Lowest wavenumber.
C  	DELV		REAL	Wavenumber spacing of LBL spectrum.
C       NPOINT		INTEGER	Number of points in spectrum.
C       G_ORD(MAXG)	REAL	Gauss-Legendre ordinates for calculating the
C				k-distribution.
C       NG              INTEGER Number of ordinates in k-distribution.
C
C	Output variables:
C	K_G(MAXG)		REAL	Calculated k-distribution.
C
C_FILE: No files openned.
C
C_CALL: No calls made.
C
C_HIST:	
C	4/4/06	NT Del_g and NGMAX now passed as parameters from
C			lbl_fknew (avoids duplicating data statements)
C***************************** VARIABLES *******************************

      IMPLICIT NONE

      INCLUDE '../includes/arrdef.f'

C Input parameters ...
      INTEGER NPOINT,NG,NGMAX,IWAVE,ICH,NFIL
      DOUBLE PRECISION VMIN,VMAX,DELV,V1,V2
      REAL VV,G_ORD(NGMAX)
      REAL VFIL(NFIL),TFIL(NFIL),WFIL     

      REAL K_G(NGMAX)
      
      DOUBLE PRECISION OUTPUT(MPOINT),OUTPUT1(MPOINT)
      
      DOUBLE PRECISION X(MPOINT),VMIN1,VMAX1,DELV1,Y1,X1(MPOINT)


      COMMON /SPECTRUM/ OUTPUT


C General parameters ...
      INTEGER I,J,K,NKINT

      INTEGER MAXKK
      PARAMETER (MAXKK=5000)

      DOUBLE PRECISION YY(MPOINT),YMIN,YMAX
      REAL KMAX,KMIN,DELK,XFR,YK,SUMK
      REAL F(MAXKK),ORD(MAXKK),G(MAXKK),DK(MAXKK)
      REAL FRAC1,DG,XMINK
      REAL GG,SUM
      REAL DEL_G(NGMAX)
      PARAMETER (XMINK=-35.0)

C******************************** CODE *********************************

C     Spectrum is calculated in wavenumber space. If we want
C     a k-table for wavelength space we need to interpolate onto a grid
C     of equally spaced wavelengths

      PRINT*,'CALC_FKDIST_WAVEC: IWAVE = ',IWAVE

      DO I=1,NPOINT
        X(I)=VMIN+(I-1)*DELV
      ENDDO
      
      
      open(12,file='kdist_wave_spec1.dat',status='unknown')
       write(12,*)npoint
       do i=1,npoint
        write(12,*)x(i),output(i)
       enddo
      close(12)

      

      IF(IWAVE.EQ.0)THEN
       VMAX=VMIN+(NPOINT-1)*DELV
       DO I=1,NPOINT
        X(I)=VMIN+(I-1)*DELV
       ENDDO

       print*,'NPOINT,DELV=',NPOINT,DELV

       VMIN1 = 1E4/VMAX
       VMAX1 = 1E4/VMIN
       DELV1 = (VMAX1-VMIN1)/FLOAT(NPOINT-1)

       print*,'VMIN1,VMAX1=',VMIN1,VMAX1,VMIN,VMAX

       DO I=1,NPOINT
        V1 = VMIN1+(I-1)*DELV1
        V2 = 1E4/V1
c        CALL VERINT(X,OUTPUT,NPOINT,Y1,V2)
	CALL VERINT_DP(X,OUTPUT,NPOINT,Y1,V2)
        OUTPUT1(I)=Y1
c Fletcher addition:
        X1(I)=V2
	
       ENDDO
       print*,'End Loop:',v1,v2,y1,vmin1,delv1
       
       DO I=1,NPOINT
        OUTPUT(I)=OUTPUT1(I)
c	X(I)=X1(I)
c	print*,'Loop 2:',VMIN1+(I-1)*DELV1,output(i),x(i),1e4/x(i)
       ENDDO
       
       
       
      ENDIF
      
      
      open(12,file='kdist_wave_spec2.dat',status='unknown')
       write(12,*)npoint
       do i=1,npoint
        write(12,*)x1(i),output1(i)
       enddo
      close(12)

C=======================================================================
C
C	Calculate the cumulative k-distribution.
C
C=======================================================================
      YMIN = 10000.0
      ICH=0
      YMAX = XMINK-10.
      DO 301 I=1,NPOINT
       IF(OUTPUT(I).GT.0)THEN
        YMIN = MIN(YMIN,OUTPUT(I))
        ICH=1
       ENDIF
301   CONTINUE
      IF(ICH.EQ.1)THEN
       YMIN=DLOG10(YMIN)
      ELSE
       YMIN=XMINK
      ENDIF
      DO 304 I=1,NPOINT
        IF(OUTPUT(I).LE.0.)THEN
          YY(I) = YMIN
        ELSE
          YY(I) = DLOG10(OUTPUT(I))
        ENDIF
        YMAX = MAX(YMAX,YY(I))
304   CONTINUE

      SUMK = 0.0
      DO I=1,NPOINT-1
        SUMK = SUMK + 0.5*(OUTPUT(I) + OUTPUT(I+1))
      ENDDO
      KMIN = FLOOR(YMIN)
      KMAX = CEILING(YMAX)

      IF(YMIN.EQ.YMAX)THEN
        DO I=1,NG
          IF(KMIN.EQ.XMINK)THEN
            K_G(I) = 0.0
          ELSE
            K_G(I) = 10**KMIN
          ENDIF
        ENDDO
C        print*,'calc_fkdist_wavec - zero spec - returning'
        RETURN
      ENDIF

      NKINT = MAXKK
      DELK = (KMAX - KMIN)/(1.0*NKINT-1.0)
      DO 315 K=1,NKINT
        F(K) = 0.0
        G(K) = 0.0
        ORD(K) = 10**(KMIN + (K - 1)*DELK)
315   CONTINUE

      DO 306 I=1,NPOINT
        VV = VMIN + (I-1)*DELV
        WFIL=1.0
        IF(NFIL.GT.0)THEN
         IF(VV.GE.VFIL(1).AND.VV.LE.VFIL(NFIL))THEN
          CALL INTERP(VFIL,TFIL,NFIL,WFIL,VV)
         ELSE
           WFIL=0.0
         ENDIF
        ENDIF

        K = 1 + INT((YY(I) - KMIN)/DELK)
        YK = KMIN + (K - 1)*DELK
        XFR = (YY(I) - YK)/DELK
        IF(XFR.GE.0.5.AND.K.LT.NKINT)K = K + 1
        F(K) = F(K) + WFIL
306   CONTINUE

      SUM = 0.0
      DO 307 K=1,NKINT
        SUM = SUM + F(K)
        G(K) =  SUM
307   CONTINUE

      DO 311 K=1,NKINT
       G(K)=G(K)/SUM
311   CONTINUE

      SUMK = 0.0       
      DO K=1,NKINT-1
        DG = G(K+1) - G(K)
        SUMK = SUMK + 0.5*(ORD(K) + ORD(K+1))*DG
      ENDDO

      SUMK = 0.0
      DO 308 J=1,NG
        GG = G_ORD(J)

        K = 1
243     IF(K.LT.NKINT)THEN
          IF(G(K).LE.GG.AND.G(K+1).GT.GG)THEN
            FRAC1 = (GG - G(K))/(G(K+1) - G(K))
            K_G(J) = (1 - FRAC1)*ORD(K) + FRAC1*ORD(K+1)
          ELSE
            K = K + 1
            GOTO 243
          ENDIF
        ENDIF
        SUMK = SUMK + DEL_G(J)*K_G(J)
308   CONTINUE

      open(12,file='calc_kdist_wave.dat',status='unknown')
       write(12,*)npoint
       do i=1,npoint
        write(12,*)x(i),output(i)
       enddo
       write(12,*)nkint
       do i=1,nkint
        write(12,*)f(i),g(i),ord(i)
       enddo
      close(12)



      RETURN      

      END
