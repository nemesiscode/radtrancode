      SUBROUTINE COSFT(Y,N,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION Y(N)
      THETA=3.14159265358979D0/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      SUM=Y(1)
      M=N/2
      DO 11 J=1,M-1
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
        Y1=0.5*(Y(J+1)+Y(N-J+1))
        Y2=(Y(J+1)-Y(N-J+1))
        Y(J+1)=Y1-WI*Y2
        Y(N-J+1)=Y1+WI*Y2
        SUM=SUM+WR*Y2
11    CONTINUE
      CALL REALFT(Y,M,+1)
      Y(2)=SUM
      DO 12 J=4,N,2
        SUM=SUM+Y(J)
        Y(J)=SUM
12    CONTINUE
      IF (ISIGN.EQ.-1) THEN
        EVEN=Y(1)
        ODD=Y(2)
        DO 13 I=3,N-1,2
          EVEN=EVEN+Y(I)
          ODD=ODD+Y(I+1)
13      CONTINUE
        ENF0=2.0*(EVEN-ODD)
        SUMO=Y(1)-ENF0
        SUME=(2.0*ODD/FLOAT(N))-SUMO
        Y(1)=0.5*ENF0
        Y(2)=Y(2)-SUME
        DO 14 I=3,N-1,2
          Y(I)=Y(I)-SUMO
          Y(I+1)=Y(I+1)-SUME
14      CONTINUE
      ENDIF
      RETURN
      END
