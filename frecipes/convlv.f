      SUBROUTINE CONVLV(DATA,N,RESPNS,M,ISIGN,ANS)
      PARAMETER(NMAX=8192)
      DIMENSION DATA(N),RESPNS(N)
      COMPLEX FFT(NMAX),ANS(N)
      DO 11 I=1,(M-1)/2
        RESPNS(N+1-I)=RESPNS(M+1-I)
11    CONTINUE
      DO 12 I=(M+3)/2,N-(M-1)/2
        RESPNS(I)=0.0
12    CONTINUE
      CALL TWOFFT(DATA,RESPNS,FFT,ANS,N)
      NO2=N/2
      DO 13 I=1,NO2+1
        IF (ISIGN.EQ.1) THEN
          ANS(I)=FFT(I)*ANS(I)/NO2
        ELSE IF (ISIGN.EQ.-1) THEN
          IF (CABS(ANS(I)).EQ.0.0) PAUSE 'deconvolving at a response zer
     *o'
          ANS(I)=FFT(I)/ANS(I)/NO2
        ELSE
          PAUSE 'no meaning for ISIGN'
        ENDIF
13    CONTINUE
      ANS(1)=CMPLX(REAL(ANS(1)),REAL(ANS(NO2+1)))
      CALL REALFT(ANS,NO2,-1)
      RETURN
      END
