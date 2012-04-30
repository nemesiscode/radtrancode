      SUBROUTINE SPCTRM(P,M,K,OVRLAP,W1,W2)
      LOGICAL OVRLAP
      DIMENSION P(M),W1(*),W2(M)
      WINDOW(J)=(1.-ABS(((J-1)-FACM)*FACP))
C     WINDOW(J)=1.
C     WINDOW(J)=(1.-(((J-1)-FACM)*FACP)**2)
      MM=M+M
      M4=MM+MM
      M44=M4+4
      M43=M4+3
      DEN=0.
      FACM=M-0.5
      FACP=1./(M+0.5)
      SUMW=0.
      DO 11 J=1,MM
        SUMW=SUMW+WINDOW(J)**2
11    CONTINUE
      DO 12 J=1,M
        P(J)=0.
12    CONTINUE
      IF(OVRLAP)THEN
        READ (9,*) (W2(J),J=1,M)
      ENDIF
      DO 18 KK=1,K
        DO 15 JOFF=-1,0,1
          IF (OVRLAP) THEN
            DO 13 J=1,M
              W1(JOFF+J+J)=W2(J)
13          CONTINUE
            READ (9,*) (W2(J),J=1,M)
            JOFFN=JOFF+MM
            DO 14 J=1,M
              W1(JOFFN+J+J)=W2(J)
14          CONTINUE
          ELSE
            READ (9,*) (W1(J),J=JOFF+2,M4,2)
          ENDIF
15      CONTINUE
        DO 16 J=1,MM
          J2=J+J
          W=WINDOW(J)
          W1(J2)=W1(J2)*W
          W1(J2-1)=W1(J2-1)*W
16      CONTINUE
        CALL FOUR1(W1,MM,1)
        P(1)=P(1)+W1(1)**2+W1(2)**2
        DO 17 J=2,M
          J2=J+J
          P(J)=P(J)+W1(J2)**2+W1(J2-1)**2
     *        +W1(M44-J2)**2+W1(M43-J2)**2
17      CONTINUE
        DEN=DEN+SUMW
18    CONTINUE
      DEN=M4*DEN
      DO 19 J=1,M
        P(J)=P(J)/DEN
19    CONTINUE
      RETURN
      END
