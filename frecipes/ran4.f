      FUNCTION RAN4(IDUM)
      PARAMETER (IM=11979,IA=430,IC=2531,NACC=24)
      DIMENSION INP(64),JOT(64),KEY(64),POW(65)
      DATA IFF/0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IDUM,IM)
        POW(1)=0.5
        DO 11 J=1,64
          IDUM=MOD(IDUM*IA+IC,IM)
          KEY(J)=(2*IDUM)/IM
          INP(J)=MOD((4*IDUM)/IM,2)
          POW(J+1)=0.5*POW(J)
11      CONTINUE
        NEWKEY=1
      ENDIF
      ISAV=INP(64)
      IF(ISAV.NE.0)THEN
        INP(4)=1-INP(4)
        INP(3)=1-INP(3)
        INP(1)=1-INP(1)
      ENDIF
      DO 12 J=64,2,-1
        INP(J)=INP(J-1)
12    CONTINUE
      INP(1)=ISAV
      CALL DES(INP,KEY,NEWKEY,0,JOT)
      RAN4=0.0
      DO 13 J=1,NACC
        IF(JOT(J).NE.0)RAN4=RAN4+POW(J)
13    CONTINUE
      RETURN
      END
