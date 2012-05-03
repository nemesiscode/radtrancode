      SUBROUTINE FGAUSS(X,A,Y,DYDA,NA)
      DIMENSION A(NA),DYDA(NA)
      Y=0.
      DO 11 I=1,NA-1,3
        ARG=(X-A(I+1))/A(I+2)
        EX=EXP(-ARG**2)
        FAC=A(I)*EX*2.*ARG
        Y=Y+A(I)*EX
        DYDA(I)=EX
        DYDA(I+1)=FAC/A(I+2)
        DYDA(I+2)=FAC*ARG/A(I+2)
11    CONTINUE
      RETURN
      END
