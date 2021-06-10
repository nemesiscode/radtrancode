      SUBROUTINE RAMANJSOURCE(VV,NLAY,FMEAN)
C     ***********************************************************************
C     Subroutine to calculate Raman source functions using outline
C     method suggested by Sromovsky (2005)
C
C     Input variables
C	VV	REAL	Calculation wavenumber
C	NLAY	INTEGER	Number of atmospheric layers
C	FMEAN(NLAY) REAL	density-weight sphere-integrated flux
C				for each layer (W cm-2 micron-1)
C    
C     Output variable
C	None returned, but Ramin-scattered light added to source function
C       matrix JRAMAN(1000,MAXSCATLAY) at longer wavelengths, held in a
C       COMMON block, which is later added to the thermal emission 
C       source function by double1.f.
C
C     Pat Irwin	7/6/21
C
C     ***********************************************************************

      IMPLICIT NONE
      INCLUDE '../includes/arrdef.f'
     
      REAL FMEAN(MAXSCATLAY),VV,LAMBDA,H2ABUND(MAXSCATLAY)
      integer idiag,iquiet,NLAY,I
      common/diagnostic/idiag,iquiet
      INTEGER ILAMBDA,IRAMAN,NRAMAN,ILAMBDA0
      REAL DENS(MAXSCATLAY),FPRAMAN(MAXSCATLAY)
      REAL JRAMAN(1000,MAXSCATLAY),VRAM0,VRAMST
      REAL PI,KWT
      REAL LAMBDA0,DNU,FPARA
      DOUBLE PRECISION SPECS1,SPECQ,SPECS0,SPEC,SIG
      LOGICAL RAMAN
      COMMON/RAMAN/RAMAN,DENS,JRAMAN,VRAM0,VRAMST,NRAMAN,FPRAMAN,
     1 H2ABUND,IRAMAN,KWT

      REAL DNUS0,DNUS1,DNUQ

      PARAMETER(PI=3.1415927)
C     Wavenumber shifts of S0, S1 and combined 'Q' transitions 
C     See Sromovsky (2005)
      PARAMETER(DNUS0=354.69,DNUS1=587.07,DNUQ=4161.0)

      LAMBDA0=1e4/VV
      ILAMBDA0 = 1+INT((LAMBDA0-VRAM0)/VRAMST)

      CALL RAMANXSEC(LAMBDA0,SPECS0,SPECS1,SPECQ)

C     Go three the three possible scattered terms to add

C     1 - S0 transition
      DNU = DNUS0
      LAMBDA = 1E4/(VV-DNU)
      ILAMBDA = 1+INT((LAMBDA-VRAM0)/VRAMST)
      IF(ILAMBDA.EQ.ILAMBDA0)GOTO 1000
      IF(ILAMBDA.GE.1.AND.ILAMBDA.LE.NRAMAN)THEN
       DO 100 I=1,NLAY
        FPARA=FPRAMAN(I)
        spec = fpara*SPECS0
        SIG = FMEAN(I)*H2ABUND(I)*spec/4*PI
C       photons re-radiated at lower frequency, so energy is lower (photons conserved)
        SIG=(VV-DNU)*SIG/VV
C       Assume triangular ILS so average J over target bin and one either side.
        IF((ILAMBDA-ILAMBDA0).GT.1.AND.ILAMBDA.GT.1)THEN
         JRAMAN(ILAMBDA-1,I)=JRAMAN(ILAMBDA-1,I)+0.25*SNGL(SIG*KWT)
         JRAMAN(ILAMBDA,I)=JRAMAN(ILAMBDA,I)+0.5*SNGL(SIG*KWT)
         JRAMAN(ILAMBDA+1,I)=JRAMAN(ILAMBDA+1,I)+0.25*SNGL(SIG*KWT)
        ELSE
         JRAMAN(ILAMBDA,I)=JRAMAN(ILAMBDA,I)+SNGL(SIG*KWT)
        ENDIF
100    CONTINUE
      ENDIF

C     1 - S1 transition
      DNU = DNUS1
      LAMBDA = 1E4/(VV-DNU)
      ILAMBDA = 1+INT((LAMBDA-VRAM0)/VRAMST)
      IF(ILAMBDA.EQ.ILAMBDA0)GOTO 1000
      IF(ILAMBDA.GE.1.AND.ILAMBDA.LE.NRAMAN)THEN
       DO 101 I=1,NLAY       
        FPARA=FPRAMAN(I)
        spec = (1-fpara)*SPECS1
        SIG = FMEAN(I)*H2ABUND(I)*spec/4*PI
C       photons re-radiated at lower frequency, so energy is lower (photons conserved)
        SIG=(VV-DNU)*SIG/VV
C       Assume triangular ILS so average J over target bin and one either side.
        IF((ILAMBDA-ILAMBDA0).GT.1.AND.ILAMBDA.GT.1)THEN
         JRAMAN(ILAMBDA-1,I)=JRAMAN(ILAMBDA-1,I)+0.25*SNGL(SIG*KWT)
         JRAMAN(ILAMBDA,I)=JRAMAN(ILAMBDA,I)+0.5*SNGL(SIG*KWT)
         JRAMAN(ILAMBDA+1,I)=JRAMAN(ILAMBDA+1,I)+0.25*SNGL(SIG*KWT)
        ELSE
         JRAMAN(ILAMBDA,I)=JRAMAN(ILAMBDA,I)+SNGL(SIG*KWT)
        ENDIF
101    CONTINUE
      ENDIF       

C     1 - 'Q' transitions
      DNU = DNUQ
      LAMBDA = 1E4/(VV-DNU)
      ILAMBDA = 1+INT((LAMBDA-VRAM0)/VRAMST)
      IF(ILAMBDA.EQ.ILAMBDA0)GOTO 1000
      IF(ILAMBDA.GE.1.AND.ILAMBDA.LE.NRAMAN)THEN
       DO 102 I=1,NLAY
        spec = SPECQ
        SIG = FMEAN(I)*H2ABUND(I)*spec/4*PI
C       photons re-radiated at lower frequency, so energy is lower (photons conserved)
        SIG=(VV-DNU)*SIG/VV
C       Assume triangular ILS so average J over target bin and one either side.
        IF((ILAMBDA-ILAMBDA0).GT.1.AND.ILAMBDA.GT.1)THEN
         JRAMAN(ILAMBDA-1,I)=JRAMAN(ILAMBDA-1,I)+0.25*SNGL(SIG*KWT)
         JRAMAN(ILAMBDA,I)=JRAMAN(ILAMBDA,I)+0.5*SNGL(SIG*KWT)
         JRAMAN(ILAMBDA+1,I)=JRAMAN(ILAMBDA+1,I)+0.25*SNGL(SIG*KWT)
        ELSE
         JRAMAN(ILAMBDA,I)=JRAMAN(ILAMBDA,I)+SNGL(SIG*KWT)
        ENDIF
102    CONTINUE
      ENDIF

      RETURN

1000  CONTINUE
      print*,'Error in ramanjsource.f'
      print*,'Raman source function computed to lie within'
      print*,'current bin. Need to use JRAMAN array with smaller'
      print*,'step-size'
      STOP

      END
