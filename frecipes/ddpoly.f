      SUBROUTINE DDPOLY(C,NC,X,PD,ND)
      DIMENSION C(NC),PD(ND)
      PD(1)=C(NC)
      DO 11 J=2,ND
        PD(J)=0.
11    CONTINUE
      DO 13 I=NC-1,1,-1
        NND=MIN(ND,NC+1-I)
        DO 12 J=NND,2,-1
          PD(J)=PD(J)*X+PD(J-1)
12      CONTINUE
        PD(1)=PD(1)*X+C(I)
13    CONTINUE
      CONST=2.
      DO 14 I=3,ND
        PD(I)=CONST*PD(I)
        CONST=CONST*I
14    CONTINUE
      RETURN
      END
