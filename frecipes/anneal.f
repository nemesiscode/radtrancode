      SUBROUTINE ANNEAL(X,Y,IORDER,NCITY)
C Rev. 12/13/85
      DIMENSION X(NCITY),Y(NCITY),IORDER(NCITY),N(6)
      LOGICAL ANS
      ALEN(X1,X2,Y1,Y2)=SQRT((X2-X1)**2+(Y2-Y1)**2)
      NOVER=100*NCITY
      NLIMIT=10*NCITY
      TFACTR=0.9
      PATH=0.0
      T=0.5
      DO 11 I=1,NCITY-1
          I1=IORDER(I)
          I2=IORDER(I+1)
          PATH=PATH+ALEN(X(I1),X(I2),Y(I1),Y(I2))
11    CONTINUE
      I1=IORDER(NCITY)
      I2=IORDER(1)
      PATH=PATH+ALEN(X(I1),X(I2),Y(I1),Y(I2))
      IDUM=-1
      ISEED=111
      DO 15 J=1,100
          NSUCC=0
          DO 13 K=1,NOVER
12            N(1)=1+INT(NCITY*RAN3(IDUM))
              N(2)=1+INT((NCITY-1)*RAN3(IDUM))
              IF (N(2).GE.N(1)) N(2)=N(2)+1
              IDEC=IRBIT1(ISEED)
              NN=1+MOD((N(1)-N(2)+NCITY-1),NCITY)
              IF (NN.LT.3) GOTO 12
              IF (IDEC.EQ.0) THEN
                  N(3)=N(2)+INT(ABS(NN-2)*RAN3(IDUM))+1
                  N(3)=1+MOD(N(3)-1,NCITY)
                  CALL TRNCST(X,Y,IORDER,NCITY,N,DE)
                  CALL METROP(DE,T,ANS)
                  IF (ANS) THEN
                      NSUCC=NSUCC+1
                      PATH=PATH+DE
                      CALL TRNSPT(IORDER,NCITY,N)
                  ENDIF
              ELSE
                  CALL REVCST(X,Y,IORDER,NCITY,N,DE)
                  CALL METROP(DE,T,ANS)
                  IF (ANS) THEN
                      NSUCC=NSUCC+1
                      PATH=PATH+DE
                      CALL REVERS(IORDER,NCITY,N)
                  ENDIF
              ENDIF
              IF (NSUCC.GE.NLIMIT) GOTO 14
13        CONTINUE
14        WRITE(*,*)
          WRITE(*,'(1X,A,F10.6,A,F10.6)') 'T =',T,
     *            '   Path Length =',PATH
          WRITE(*,'(1X,A,I6)') 'Successful Moves: ',NSUCC
          T=T*TFACTR
          IF (NSUCC.EQ.0) RETURN
15    CONTINUE
      RETURN
      END
