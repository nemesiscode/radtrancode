      SUBROUTINE SHELL(N,ARR)
      PARAMETER (ALN2I=1.4426950, TINY=1.E-5)
      DIMENSION ARR(N)
      LOGNB2=INT(ALOG(FLOAT(N))*ALN2I+TINY)
      M=N
      DO 12 NN=1,LOGNB2
        M=M/2
        K=N-M
        DO 11 J=1,K
          I=J
3         CONTINUE
          L=I+M
          IF(ARR(L).LT.ARR(I)) THEN
            T=ARR(I)
            ARR(I)=ARR(L)
            ARR(L)=T
            I=I-M
            IF(I.GE.1)GO TO 3
          ENDIF
11      CONTINUE
12    CONTINUE
      RETURN
      END
