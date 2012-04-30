      FUNCTION POIDEV(XM,IDUM)
      PARAMETER (PI=3.141592654)
      DATA OLDM /-1./
      IF (XM.LT.12.)THEN
        IF (XM.NE.OLDM) THEN
          OLDM=XM
          G=EXP(-XM)
        ENDIF
        EM=-1
        T=1.
2       EM=EM+1.
        T=T*RAN1(IDUM)
        IF (T.GT.G) GO TO 2
      ELSE
        IF (XM.NE.OLDM) THEN
          OLDM=XM
          SQ=SQRT(2.*XM)
          ALXM=ALOG(XM)
          G=XM*ALXM-GAMMLN(XM+1.)
        ENDIF
1       Y=TAN(PI*RAN1(IDUM))
        EM=SQ*Y+XM
        IF (EM.LT.0.) GO TO 1
        EM=INT(EM)
        T=0.9*(1.+Y**2)*EXP(EM*ALXM-GAMMLN(EM+1.)-G)
        IF (RAN1(IDUM).GT.T) GO TO 1
      ENDIF
      POIDEV=EM
      RETURN
      END
