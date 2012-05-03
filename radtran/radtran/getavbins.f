      SUBROUTINE GETAVBINS(NPOINT,KWAVE,WAVEIN,FWHMIN,DELV,W1,W2,
     1 IWAVE,NAV,IAV,FBIN)
C     ***************************************************************
C     Subroutine to find bins needed to average over a prescribed 
C     wavelength/wavenumber interval
C
C     Input variables
C	NPOINT	INTEGER	Number of wavelengths in table to be averaged
C	KWAVE	INTEGER	wavelength space of table(0=wavenumber,1=wavelength)
C	WAVEIN(MBIN) REAL Wavelengths/wavenumbers
C	FWHMIN(MBIN) REAL FWHM of table
C	DELV	REAL	Step of table
C       W1	REAL	Lower limit of averaging bin
C	W2	REAL	Upper limit of averaging bin
C	IWAVE	INTEGER	wavelength space of bin
C
C     Output variables
C	NAV	INTEGER	Number of bins needed
C	IAV(MBIN) INTEGER Location of bins
C	FBIN(MBIN) REAL	Required weights of bins
C
C     Pat Irwin		23/1/04
C
C     ***************************************************************
      IMPLICIT NONE
      INTEGER NPOINT,IS1,IS2,MBIN,KWAVE,I,IWAVE
      INTEGER NAV,NODD,NEVEN
      REAL SUMODD,SUMEVEN,SUM
      LOGICAL IODD
      PARAMETER(MBIN=8000)
      REAL W1,W2,X1,X2,WAVEIN(MBIN),FWHMIN(MBIN),FRAC(MBIN)
      INTEGER IAV(MBIN)
      REAL DELV,V1,V2,FBIN(MBIN)


      CALL TRANSLATE(IWAVE,W1,W2,KWAVE,X1,X2)

      print*,'X1,X2',X1,X2  
      CALL FINDLOC(WAVEIN,NPOINT,X1,X2,IS1,IS2)
      IS1 = MAX(1,IS1-2)
      IS2 = MIN(NPOINT,IS2+2)

      print*,'IS1,IS2',IS1,IS2       
      DO 15 I=IS1,IS2 

        V1 = WAVEIN(I)-0.5*FWHMIN(I)  
        V2 = WAVEIN(I)+0.5*FWHMIN(I)
        FRAC(I)=1.0
        IF(V1.LT.X1)THEN
         FRAC(I)=(V2-X1)/FWHMIN(I)
         IF(FRAC(I).LT.0.0)FRAC(I)=0.0
        ENDIF
        IF(V2.GT.X2)THEN 
         FRAC(I)=(X2-V1)/FWHMIN(I)
         IF(FRAC(I).LT.0.0)FRAC(I)=0.0
        ENDIF
C        print*,I,FRAC(I)
15    CONTINUE

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
C        print*,'sub-Nyquist'
        NAV = 0
        SUM = 0.0
        DO I=IS1,IS2
         IF(FRAC(I).GT.0)THEN
          NAV=NAV+1
          IAV(NAV)=I
          FBIN(NAV)=FRAC(I)
C          print*,NAV,FBIN(NAV)
          SUM=SUM+FBIN(NAV)
          IF(KWAVE.EQ.0.AND.IWAVE.EQ.1)THEN
           FBIN(NAV)=FBIN(NAV)*(WAVEIN(IS1)/WAVEIN(I))**2
          ENDIF
         ENDIF
        ENDDO

        
        DO I=1,NAV
         FBIN(I)=FBIN(I)/SUM
        ENDDO

      ELSE
        PRINT*,'GETAVBINS - STEP OF BAND MODEL IS NOT RECOGNISED'
        print*,'FWHMIN(1),DELV',FWHMIN(1),DELV
        STOP
      ENDIF

      PRINT*,'NAV = ',NAV
      SUM=0.0
      DO I=1,NAV
C        PRINT*,I,IAV(I),FBIN(I)
        SUM=SUM+FBIN(I)
      ENDDO
      PRINT*,'SUM = ',SUM

      IF(SUM.EQ.0.0)THEN
       PRINT*,'Wavelength not covered'
      ELSE
       IF(ABS(SUM-1.0).GT.0.01)THEN
        PRINT*,'Error in Getavbins - sum of weights does not add to 1.0'
        STOP
       ENDIF
      ENDIF
      RETURN

      END
