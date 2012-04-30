      SUBROUTINE SMOOFT(Y,N,PTS)
      PARAMETER(MMAX=1024)
      DIMENSION Y(MMAX)
      M=2
      NMIN=N+2.*PTS
1     IF(M.LT.NMIN)THEN
        M=2*M
      GO TO 1
      ENDIF
      IF(M.GT.MMAX) PAUSE 'MMAX too small'
      CONST=(PTS/M)**2
      Y1=Y(1)
      YN=Y(N)
      RN1=1./(N-1.)
      DO 11 J=1,N
        Y(J)=Y(J)-RN1*(Y1*(N-J)+YN*(J-1))
11    CONTINUE
      IF(N+1.LE.M)THEN
        DO 12 J=N+1,M
          Y(J)=0.
12      CONTINUE
      ENDIF
      MO2=M/2
      CALL REALFT(Y,MO2,1)
      Y(1)=Y(1)/MO2
      FAC=1.
      DO 13 J=1,MO2-1
        K=2*J+1
        IF(FAC.NE.0.)THEN
          FAC=AMAX1(0.,(1.-CONST*J**2)/MO2)
          Y(K)=FAC*Y(K)
          Y(K+1)=FAC*Y(K+1)
        ELSE
          Y(K)=0.
          Y(K+1)=0.
        ENDIF
13    CONTINUE
      FAC=AMAX1(0.,(1.-0.25*PTS**2)/MO2)
      Y(2)=FAC*Y(2)
      CALL REALFT(Y,MO2,-1)
      DO 14 J=1,N
        Y(J)=RN1*(Y1*(N-J)+YN*(J-1))+Y(J)
14    CONTINUE
      RETURN
      END
