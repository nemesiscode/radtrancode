      SUBROUTINE ELMHES(A,N,NP)
      DIMENSION A(NP,NP)
      IF(N.GT.2)THEN
        DO 17 M=2,N-1
          X=0.
          I=M
          DO 11 J=M,N
            IF(ABS(A(J,M-1)).GT.ABS(X))THEN
              X=A(J,M-1)
              I=J
            ENDIF
11        CONTINUE
          IF(I.NE.M)THEN
            DO 12 J=M-1,N
              Y=A(I,J)
              A(I,J)=A(M,J)
              A(M,J)=Y
12          CONTINUE
            DO 13 J=1,N
              Y=A(J,I)
              A(J,I)=A(J,M)
              A(J,M)=Y
13          CONTINUE
          ENDIF
          IF(X.NE.0.)THEN
            DO 16 I=M+1,N
              Y=A(I,M-1)
              IF(Y.NE.0.)THEN
                Y=Y/X
                A(I,M-1)=Y
                DO 14 J=M,N
                  A(I,J)=A(I,J)-Y*A(M,J)
14              CONTINUE
                DO 15 J=1,N
                  A(J,M)=A(J,M)+Y*A(J,I)
15              CONTINUE
              ENDIF
16          CONTINUE
          ENDIF
17      CONTINUE
      ENDIF
      RETURN
      END
