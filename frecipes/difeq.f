      SUBROUTINE DIFEQ(K,K1,K2,JSF,IS1,ISF,INDEXV,NE,S,NSI,NSJ,Y,NYJ,NYK
     *)
      PARAMETER(M=41)
      COMMON X(M),H,MM,N,C2,ANORM
      DIMENSION Y(NYJ,NYK),S(NSI,NSJ),INDEXV(NYJ)
      IF(K.EQ.K1) THEN
        IF(MOD(N+MM,2).EQ.1)THEN
          S(3,3+INDEXV(1))=1.
          S(3,3+INDEXV(2))=0.
          S(3,3+INDEXV(3))=0.
          S(3,JSF)=Y(1,1)
        ELSE
          S(3,3+INDEXV(1))=0.
          S(3,3+INDEXV(2))=1.
          S(3,3+INDEXV(3))=0.
          S(3,JSF)=Y(2,1)
        ENDIF
      ELSE IF(K.GT.K2) THEN
        S(1,3+INDEXV(1))=-(Y(3,M)-C2)/(2.*(MM+1.))
        S(1,3+INDEXV(2))=1.
        S(1,3+INDEXV(3))=-Y(1,M)/(2.*(MM+1.))
        S(1,JSF)=Y(2,M)-(Y(3,M)-C2)*Y(1,M)/(2.*(MM+1.))
        S(2,3+INDEXV(1))=1.
        S(2,3+INDEXV(2))=0.
        S(2,3+INDEXV(3))=0.
        S(2,JSF)=Y(1,M)-ANORM
      ELSE
        S(1,INDEXV(1))=-1.
        S(1,INDEXV(2))=-.5*H
        S(1,INDEXV(3))=0.
        S(1,3+INDEXV(1))=1.
        S(1,3+INDEXV(2))=-.5*H
        S(1,3+INDEXV(3))=0.
        TEMP=H/(1.-(X(K)+X(K-1))**2*.25)
        TEMP2=.5*(Y(3,K)+Y(3,K-1))-C2*.25*(X(K)+X(K-1))**2
        S(2,INDEXV(1))=TEMP*TEMP2*.5
        S(2,INDEXV(2))=-1.-.5*TEMP*(MM+1.)*(X(K)+X(K-1))
        S(2,INDEXV(3))=.25*TEMP*(Y(1,K)+Y(1,K-1))
        S(2,3+INDEXV(1))=S(2,INDEXV(1))
        S(2,3+INDEXV(2))=2.+S(2,INDEXV(2))
        S(2,3+INDEXV(3))=S(2,INDEXV(3))
        S(3,INDEXV(1))=0.
        S(3,INDEXV(2))=0.
        S(3,INDEXV(3))=-1.
        S(3,3+INDEXV(1))=0.
        S(3,3+INDEXV(2))=0.
        S(3,3+INDEXV(3))=1.
        S(1,JSF)=Y(1,K)-Y(1,K-1)-.5*H*(Y(2,K)+Y(2,K-1))
        S(2,JSF)=Y(2,K)-Y(2,K-1)-TEMP*((X(K)+X(K-1))
     *      *.5*(MM+1.)*(Y(2,K)+Y(2,K-1))-TEMP2*
     *      .5*(Y(1,K)+Y(1,K-1)))
        S(3,JSF)=Y(3,K)-Y(3,K-1)
      ENDIF
      RETURN
      END
