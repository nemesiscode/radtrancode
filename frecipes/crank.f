      SUBROUTINE CRANK(N,W,S)
      DIMENSION W(N)
      S=0.
      J=1
1     IF(J.LT.N)THEN
        IF(W(J+1).NE.W(J))THEN
          W(J)=J
          J=J+1
        ELSE
          DO 11 JT=J+1,N
            IF(W(JT).NE.W(J))GO TO 2
11        CONTINUE
          JT=N+1
2         RANK=0.5*(J+JT-1)
          DO 12 JI=J,JT-1
            W(JI)=RANK
12        CONTINUE
          T=JT-J
          S=S+T**3-T
          J=JT
        ENDIF
      GO TO 1
      ENDIF
      IF(J.EQ.N)W(N)=N
      RETURN
      END
