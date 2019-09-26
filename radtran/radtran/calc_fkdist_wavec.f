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
C       IWAVE		INTEGER	0=Wavelength, 1=Wavenumber
C       VMIN		REAL	Lowest wavenumber.
C  	DELV		REAL	Wavenumber spacing of LBL spectrum.
C       NPOINT		INTEGER	Number of points in spectrum.
C       NFIL		INTEGER	Number of points in filter function
C	VFIL(NFIL)	REAL	Filter function wavelengths/wavenumbers
C	TFIL(NFILE)	REAL	Filter function transmissions
C       G_ORD(MAXG)	REAL	Gauss-Legendre ordinates for calculating the
C				k-distribution.
C       NG              INTEGER Number of ordinates in k-distribution.
C	NGMAX		INTEGER Maximum allowable size of NG
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
      INTEGER NPOINT,NG,NGMAX,IWAVE,ICH,NFIL,ITYPE
      REAL VMIN,VMAX,DELV,VV,G_ORD(NGMAX),V1,V2
      REAL VFIL(NFIL),TFIL(NFIL),WFIL     
      CHARACTER*1 ANS

      REAL OUTPUT(MPOINT),K_G(NGMAX),OUTPUT1(MPOINT)
      REAL X(MPOINT),VMIN1,VMAX1,DELV1,Y1


      COMMON /SPECTRUM/ OUTPUT


C General parameters ...
      INTEGER I,J,K,NKINT

      INTEGER MAXKK
      PARAMETER (MAXKK=5000)

      REAL YY(MPOINT),YMIN,YMAX
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

1     FORMAT(A)
      PRINT*,'CALC_FKDIST_WAVEC: IWAVE = ',IWAVE

      DO I=1,NPOINT
        X(I)=VMIN+(I-1)*DELV
      ENDDO

      IF(IWAVE.EQ.0)THEN
       VMAX=VMIN+(NPOINT-1)*DELV
       DO I=1,NPOINT
        X(I)=VMIN+(I-1)*DELV
       ENDDO
       
       print*,'NPOINT,VMIN,DELV=',NPOINT,DELV

       VMIN1 = 1E4/VMAX
       VMAX1 = 1E4/VMIN
       DELV1 = (VMAX1-VMIN1)/FLOAT(NPOINT-1)

       print*,'NPOINT,VMIN1,VMAX1,DELV1=',NPOINT,VMIN1,
     1  VMAX1,DELV1

       DO I=1,NPOINT
        V1 = VMIN1+(I-1)*DELV1
        V2 = 1E4/V1
        CALL VERINT(X,OUTPUT,NPOINT,Y1,V2)
        OUTPUT1(I)=Y1
       ENDDO

       DO I=1,NPOINT
        OUTPUT(I)=OUTPUT1(I)
       ENDDO

      ENDIF

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
       YMIN=ALOG10(YMIN)
      ELSE
       YMIN=XMINK
      ENDIF
      DO 304 I=1,NPOINT
        IF(OUTPUT(I).LE.0.)THEN
          YY(I) = YMIN
        ELSE
          YY(I) = ALOG10(OUTPUT(I))
          IF(ISNAN(YY(I)))THEN
           print*,'NAN',I,OUTPUT(I)
           YY(I)=YMIN
          ENDIF
        ENDIF
        YMAX = MAX(YMAX,YY(I))
304   CONTINUE

      SUMK = 0.0
      DO I=1,NPOINT-1
        SUMK = SUMK + 0.5*(OUTPUT(I) + OUTPUT(I+1))
      ENDDO
      KMIN = FLOOR(YMIN)
      KMAX = CEILING(YMAX)
      print*,KMIN,KMAX,YMIN,YMAX
      IF(YMIN.EQ.YMAX)THEN
        DO I=1,NG
          IF(KMIN.EQ.XMINK)THEN
            K_G(I) = 0.0
          ELSE
            K_G(I) = 10**KMIN
          ENDIF
        ENDDO
        print*,'calc_fkdist_wavec - zero spec - returning'
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
C        print*,VV,NFIL,WFIL

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


      ITYPE=0

      IF(ITYPE.EQ.0)THEN
C      Technically correct method of sampling ORD at the 
C      Gaussian ordinates G_ORD
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
308    CONTINUE

      ELSE
C      Empirical method of sampling ORD by ensuring 
C      Integral(ORD.DG) = Sum(K_G.DEL_G)

       CALL KINTEGRATE(NKINT,ORD,G,NG,K_G,DEL_G)

      ENDIF

C      open(12,file='calc_kdist_wavex.dat',status='unknown')
C       write(12,*)npoint
C       do i=1,npoint
C        write(12,*)x(i),output(i)
C       enddo
C       write(12,*)nkint
C       do i=1,nkint
C        write(12,*)f(i),g(i),ord(i)
C       enddo
C       write(12,*)ng
C       do i=1,ng
C        write(12,*)g_ord(i),k_g(i)
C       enddo
C      close(12)

C      print*,'press a key to continue'
C      read(5,1)ans

      print*,'CALC_FKDIST_WAVEC called OK'

      RETURN      

      END


      SUBROUTINE KINTEGRATE(NKINT,ORD,G,NG,K_G,DEL_G)
C     *************************************************************
C     Subroutine to integrate ord.dg and set K_G to be equal to average of
C     integrals over subranges of width DEL_G
C
C     Pat Irwin 17/11/18
C
C     *************************************************************
      IMPLICIT NONE
      INTEGER NKINT,NG,I,IG
      REAL ORD(NKINT),G(NKINT),DG,F,ORDX
      REAL K_G(NG),DEL_G(NG),SUM,GG,DG1
      SUM=0.
      DO I=1,NG
       K_G(I)=0.
       SUM=SUM+DEL_G(I)
      ENDDO
C      print*,SUM
      DO I=1,NG
       DEL_G(I)=DEL_G(I)/SUM
      ENDDO

      IG=1
      GG=DEL_G(IG)
      SUM=0.
      DO 101 I=2,NKINT
       IF(G(I).LT.GG)THEN
        K_G(IG)=K_G(IG)+0.5*(ORD(I)+ORD(I-1))*(G(I)-G(I-1))
       ELSE
         DG = (GG-G(I-1))
         F = DG/(G(I)-G(I-1))
         ORDX = ORD(I-1)*(1.0-F) + ORD(I)*F
         K_G(IG)=K_G(IG)+0.5*(ORDX+ORD(I-1))*DG
C         print*,IG,GG,K_G(IG)
         IF(IG.LT.NG)THEN
          DG1 = G(I)-GG
          K_G(IG+1)=K_G(IG+1)+0.5*(ORD(I)+ORDX)*DG1
          IG=IG+1
          GG=GG+DEL_G(IG) 
C          print*,'A',IG,GG,K_G(IG)
         ENDIF
       ENDIF 
101   CONTINUE

      SUM=0
      DO 102 I=1,NG
       K_G(I)=K_G(I)/DEL_G(I)
       SUM=SUM+K_G(I)*DEL_G(I)
102   CONTINUE
C      print*,'Rough Sum = ',SUM

      SUM=0.
      DO 103 I=2,NKINT
       SUM=SUM+0.5*(ORD(I-1)+ORD(I))*(G(I)-G(I-1))
103   CONTINUE
C      print*,'Fine Sum = ',SUM
C      print*,'Last elements = ',G(NKINT),ORD(NKINT)

      RETURN

      END
