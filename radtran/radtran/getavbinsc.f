      SUBROUTINE GETAVBINSC(NPOINT,KWAVE,WAVEIN,FWHMIN,DELV,WCENTRAL,
     1 WFWHM,ISHAPE,NFILBIN,WFIL,TFIL,IWAVE,NAV,IAV,FBIN)
C     ***************************************************************
C     Subroutine to find bins needed to average over a prescribed 
C     wavelength/wavenumber interval
C
C     Input variables
C	NPOINT	      INTEGER	Number of wavelengths in table to be averaged
C	KWAVE	      INTEGER	wavelength space of table(0=wavenumber,1=wavelength)
C	WAVEIN(MBIN)  REAL      Wavelengths/wavenumbers
C	FWHMIN(MBIN)  REAL      FWHM of table
C	DELV	      REAL	Step of table
C       WCENTRAL      REAL	centre of averaging bin
C	WFWHM	      REAL	FWHM of averaging bin
C	ISHAPE	      INTEGER	Shape of averaging bin (0=square, 1=triangle)
C				 2= gaussian, 3 = filter file
C       NFILBIN       INTEGER   Number of bins for the filter
C       WFIL(NFILBIN) REAL      Wavelengths of the filter bins
C       TFIL(NFILBIN) REAL      Transmission of each filter bin
C	IWAVE	      INTEGER   wavelength space of bin
C
C     Output variables
C	NAV	      INTEGER   Number of bins needed
C	IAV(MBIN)     INTEGER   Location of bins
C	FBIN(MBIN)    REAL	Required weights of bins
C
C     Pat Irwin		23/01/04
C     SPH               15/08/13 Created from getavbinsc
C     ***************************************************************
      IMPLICIT NONE
      INTEGER NPOINT,IS1,IS2,MBIN,KWAVE,I,IWAVE
      INTEGER NAV,NODD,NEVEN
      REAL SUMODD,SUMEVEN,SUM
      LOGICAL IODD
      PARAMETER(MBIN=8000)
      REAL WCENTRAL,W1,W2,X1,X2,WAVEIN(MBIN),FWHMIN(MBIN),FRAC(MBIN)
      INTEGER IAV(MBIN),ISHAPE
      REAL DELV,V1,V2,FBIN(MBIN),WFWHM,MAXDV,DV,WCEN1
      INTEGER NFILBIN
      REAL WFIL(MBIN), TFIL(MBIN), WAVE, VFIL(MBIN), XFRAC

      IF(ISHAPE.EQ.0)THEN
       W1 = WCENTRAL-0.5*WFWHM
       W2 = W1+WFWHM
      ELSEIF(ISHAPE.EQ.1)THEN
       W1 = WCENTRAL-WFWHM
       W2 = W1+2*WFWHM
      ELSEIF(ISHAPE.EQ.2)THEN
       W1 = WCENTRAL-3*WFWHM
       W2 = W1 + 6*WFWHM
      ELSE
       W1 = WFIL(1)
       W2 = WFIL(NFILBIN)
       DO I=1,NFILBIN
        VFIL(I)=WFIL(I)-WCENTRAL
       ENDDO
      ENDIF
      
      print*,'AA',ISHAPE,W1,W2,IWAVE,KWAVE
      CALL TRANSLATE(IWAVE,W1,W2,KWAVE,X1,X2)
      print*,'W1,W2',W1,W2
      print*,'X1,X2,DELX',X1,X2,X2-X1  
      CALL FINDLOC(WAVEIN,NPOINT,X1,X2,IS1,IS2)
      IS1 = MAX(1,IS1-2)
      IS2 = MIN(NPOINT,IS2+2)

      print*,'IS1,IS2',IS1,IS2       
      print*,'table wavenumber range : ',WAVEIN(IS1),WAVEIN(IS2)
      DO 15 I=IS1,IS2 

C      Find offset from centre of averaging function in wavenumber space
C      of averaging function
       IF(IWAVE.EQ.0)THEN
         IF(KWAVE.EQ.0)THEN
          DV = ABS(WAVEIN(I)-WCENTRAL)
         ELSE
          DV = ABS(1E4/WAVEIN(I) - WCENTRAL)
         ENDIF
       ELSE
         IF(KWAVE.EQ.0)THEN
          DV = ABS(1E4/WAVEIN(I)-WCENTRAL)
         ELSE
          DV = ABS(WAVEIN(I)-WCENTRAL)
         ENDIF
       ENDIF

        PRINT*,'ISHAPE,DV = ',ISHAPE,DV
        IF(ISHAPE.EQ.0)THEN
         IF(DV.LE.0.5*WFWHM)THEN
          FRAC(I)=1.0
         ELSE
          FRAC(I)=0.0
         ENDIF
        ELSEIF(ISHAPE.EQ.1)THEN
         FRAC(I)=1.0-DV/WFWHM
        ELSEIF(ISHAPE.EQ.2)THEN
         FRAC(I)=EXP(-0.693*(2*DV/WFWHM)**2)
        ELSE
         CALL INTERP(VFIL,TFIL,NFILBIN,XFRAC,DV)
         FRAC(I)=XFRAC
        ENDIF

        IF(FRAC(I).GT.1)FRAC(I)=1
        IF(FRAC(I).LT.0)FRAC(I)=0

        print*,'I,DV,FRAC',I,DV,FRAC(I)


C       Now tidy up the ends if the bins of the band table are of a similar
C       size to the size of the averaging window.
        V1 = WAVEIN(I)-0.5*FWHMIN(I)  
        V2 = WAVEIN(I)+0.5*FWHMIN(I)
        IF(V1.LT.X1)THEN
         FRAC(I)=FRAC(I)*(V2-X1)/FWHMIN(I)
         IF(FRAC(I).LT.0.0)FRAC(I)=0.0
        ENDIF
        IF(V2.GT.X2)THEN 
         FRAC(I)=FRAC(I)*(X2-V1)/FWHMIN(I)
         IF(FRAC(I).LT.0.0)FRAC(I)=0.0
        ENDIF
        print*,'Mod',I,FRAC(I)


15    CONTINUE



C     Now average band table over the bins.

      IF(FWHMIN(1).EQ.2*DELV)THEN
C       NYQUIST SAMPLED BAND DATA
        NAV = 0
        IODD = .FALSE.
        NODD = 0
        NEVEN = 0
        SUMODD = 0.0
        SUMEVEN = 0.0
        DO I=IS1,IS2
         IF(FRAC(I).GT.0)THEN
          NAV=NAV+1
          IAV(NAV)=I
          FBIN(NAV)=FRAC(I)
          IF(KWAVE.EQ.0.AND.IWAVE.EQ.1)THEN
           FBIN(NAV)=FBIN(NAV)*(WAVEIN(IS1)/WAVEIN(I))**2
          ENDIF
          IODD = .NOT.IODD
          IF(IODD) THEN
           SUMODD = SUMODD+FBIN(NAV)
           NODD = NODD + 1
          ELSE
           SUMEVEN = SUMEVEN+FBIN(NAV)
           NEVEN = NEVEN+1
          ENDIF
         ENDIF
        ENDDO
        IF(SUMODD.GT.SUMEVEN)THEN
         DO I=1,NODD
          IAV(I)=IAV(2*I-1)
          FBIN(I)=FBIN(2*I-1)/SUMODD
         ENDDO
         NAV = NODD
        ENDIF

        IF(SUMODD.LT.SUMEVEN)THEN
         DO I=1,NEVEN
          IAV(I)=IAV(2*I)
          FBIN(I)=FBIN(2*I)/SUMEVEN
         ENDDO
         NAV = NEVEN
        ENDIF
        IF(SUMODD.EQ.SUMEVEN)THEN

         IF(NEVEN.LT.NODD)THEN
          DO I=1,NEVEN
           IAV(I)=IAV(2*I)
           FBIN(I)=FBIN(2*I)/SUMEVEN
          ENDDO
          NAV = NEVEN
         ELSE
          DO I=1,NODD
           IAV(I)=IAV(2*I-1)
           FBIN(I)=FBIN(2*I-1)/SUMODD
          ENDDO
          NAV = NODD
         ENDIF
        ENDIF
          
      ELSEIF(FWHMIN(1).EQ.DELV)THEN
        NAV = 0
        SUM = 0.0
        DO I=IS1,IS2
         IF(FRAC(I).GT.0)THEN
          NAV=NAV+1
          IAV(NAV)=I
          FBIN(NAV)=FRAC(I)
          IF(KWAVE.EQ.0.AND.IWAVE.EQ.1)THEN
           FBIN(NAV)=FBIN(NAV)*(WAVEIN(IS1)/WAVEIN(I))**2
          ENDIF
         ENDIF
        ENDDO
      ELSE
        PRINT*,'GETAVBINSC - STEP OF BAND MODEL IS NOT RECOGNISED'
        print*,'FWHMIN(1),DELV',FWHMIN(1),DELV
        STOP
      ENDIF

      IF(NAV.EQ.0)THEN
       IF(KWAVE.EQ.0.AND.IWAVE.EQ.1)THEN
        WCEN1 = 1E4/WCENTRAL
       ELSE
        WCEN1 = WCENTRAL
       ENDIF

       IF(WCEN1.GE.WAVEIN(1).AND.WCEN1.LE.WAVEIN(NPOINT))THEN
        PRINT*,'Find closest point if wavelength in range'
        MAXDV = WAVEIN(NPOINT)-WAVEIN(1)
        DO 16 I=IS1,IS2 

         IF(IWAVE.EQ.0)THEN
          IF(KWAVE.EQ.0)THEN
           DV = ABS(WAVEIN(I)-WCENTRAL)
          ELSE
           DV = ABS(1E4/WAVEIN(I) - WCENTRAL)
          ENDIF
         ELSE
          IF(KWAVE.EQ.0)THEN
           DV = ABS(1E4/WAVEIN(I)-WCENTRAL)
          ELSE
           DV = ABS(WAVEIN(I)-WCENTRAL)
          ENDIF
         ENDIF
         IF(DV.LT.MAXDV)THEN
          MAXDV=DV
          IAV(1)=I
          FBIN(1)=1
         ENDIF
16      CONTINUE
        NAV=1       
       ENDIF
      ENDIF

      PRINT*,'NAV = ',NAV
      
      SUM=0.0
      DO I=1,NAV
        PRINT*,I,WAVEIN(IAV(I)),FBIN(I)
        SUM=SUM+FBIN(I)
      ENDDO

      PRINT*,'SUM = ',SUM


      RETURN

      END
