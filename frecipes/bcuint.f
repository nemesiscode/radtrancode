      SUBROUTINE BCUINT(Y,Y1,Y2,Y12,X1L,X1U,X2L,X2U,X1,X2,ANSY,ANSY1,ANS
     *Y2)
      DIMENSION Y(4),Y1(4),Y2(4),Y12(4),C(4,4)
      CALL BCUCOF(Y,Y1,Y2,Y12,X1U-X1L,X2U-X2L,C)
      IF(X1U.EQ.X1L.OR.X2U.EQ.X2L)PAUSE 'bad input'
      T=(X1-X1L)/(X1U-X1L)
      U=(X2-X2L)/(X2U-X2L)
      ANSY=0.
      ANSY2=0.
      ANSY1=0.
      DO 11 I=4,1,-1
        ANSY=T*ANSY+((C(I,4)*U+C(I,3))*U+C(I,2))*U+C(I,1)
        ANSY2=T*ANSY2+(3.*C(I,4)*U+2.*C(I,3))*U+C(I,2)
        ANSY1=U*ANSY1+(3.*C(4,I)*T+2.*C(3,I))*T+C(2,I)
11    CONTINUE
      ANSY1=ANSY1/(X1U-X1L)
      ANSY2=ANSY2/(X2U-X2L)
      RETURN
      END
