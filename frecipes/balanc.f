      SUBROUTINE BALANC(A,N,NP)
      PARAMETER (RADIX=2.,SQRDX=4.)
      DIMENSION A(NP,NP)
1     CONTINUE
        LAST=1
        DO 14 I=1,N
          C=0.
          R=0.
          DO 11 J=1,N
            IF(J.NE.I)THEN
              C=C+ABS(A(J,I))
              R=R+ABS(A(I,J))
            ENDIF
11        CONTINUE
          IF(C.NE.0..AND.R.NE.0.)THEN
            G=R/RADIX
            F=1.
            S=C+R
2           IF(C.LT.G)THEN
              F=F*RADIX
              C=C*SQRDX
            GO TO 2
            ENDIF
            G=R*RADIX
3           IF(C.GT.G)THEN
              F=F/RADIX
              C=C/SQRDX
            GO TO 3
            ENDIF
            IF((C+R)/F.LT.0.95*S)THEN
              LAST=0
              G=1./F
              DO 12 J=1,N
                A(I,J)=A(I,J)*G
12            CONTINUE
              DO 13 J=1,N
                A(J,I)=A(J,I)*F
13            CONTINUE
            ENDIF
          ENDIF
14      CONTINUE
      IF(LAST.EQ.0)GO TO 1
      RETURN
      END
