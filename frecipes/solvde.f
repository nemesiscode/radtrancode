      SUBROUTINE SOLVDE(ITMAX,CONV,SLOWC,SCALV,INDEXV,NE,NB,M,
     *    Y,NYJ,NYK,C,NCI,NCJ,NCK,S,NSI,NSJ)
      PARAMETER (NMAX=10)
      DIMENSION Y(NYJ,NYK),C(NCI,NCJ,NCK),S(NSI,NSJ),SCALV(NYJ),INDEXV(N
     *YJ)
      DIMENSION ERMAX(NMAX),KMAX(NMAX)
      K1=1
      K2=M
      NVARS=NE*M
      J1=1
      J2=NB
      J3=NB+1
      J4=NE
      J5=J4+J1
      J6=J4+J2
      J7=J4+J3
      J8=J4+J4
      J9=J8+J1
      IC1=1
      IC2=NE-NB
      IC3=IC2+1
      IC4=NE
      JC1=1
      JCF=IC3
      DO 16 IT=1,ITMAX
        K=K1
        CALL DIFEQ(K,K1,K2,J9,IC3,IC4,INDEXV,NE,S,NSI,NSJ,Y,NYJ,NYK)
        CALL PINVS(IC3,IC4,J5,J9,JC1,K1,C,NCI,NCJ,NCK,S,NSI,NSJ)
        DO 11 K=K1+1,K2
          KP=K-1
          CALL DIFEQ(K,K1,K2,J9,IC1,IC4,INDEXV,NE,S,NSI,NSJ,Y,NYJ,NYK)
          CALL RED(IC1,IC4,J1,J2,J3,J4,J9,IC3,JC1,JCF,KP,
     *        C,NCI,NCJ,NCK,S,NSI,NSJ)
          CALL PINVS(IC1,IC4,J3,J9,JC1,K,C,NCI,NCJ,NCK,S,NSI,NSJ)
11      CONTINUE
        K=K2+1
        CALL DIFEQ(K,K1,K2,J9,IC1,IC2,INDEXV,NE,S,NSI,NSJ,Y,NYJ,NYK)
        CALL RED(IC1,IC2,J5,J6,J7,J8,J9,IC3,JC1,JCF,K2,
     *      C,NCI,NCJ,NCK,S,NSI,NSJ)
        CALL PINVS(IC1,IC2,J7,J9,JCF,K2+1,C,NCI,NCJ,NCK,S,NSI,NSJ)
        CALL BKSUB(NE,NB,JCF,K1,K2,C,NCI,NCJ,NCK)
        ERR=0.
        DO 13 J=1,NE
          JV=INDEXV(J)
          ERMAX(J)=0.
          ERRJ=0.
          KMAX(J)=0
          VMAX=0.
          DO 12 K=K1,K2
            VZ=ABS(C(J,1,K))
            IF(VZ.GT.VMAX) THEN
               VMAX=VZ
               KM=K
            ENDIF
            ERRJ=ERRJ+VZ
12        CONTINUE
          ERR=ERR+ERRJ/SCALV(JV)
          ERMAX(J)=C(J,1,KM)/SCALV(JV)
          KMAX(J)=KM
13      CONTINUE
        ERR=ERR/NVARS
        FAC=SLOWC/MAX(SLOWC,ERR)
        DO 15 JV=1,NE
          J=INDEXV(JV)
          DO 14 K=K1,K2
            Y(J,K)=Y(J,K)-FAC*C(JV,1,K)
14        CONTINUE
15      CONTINUE
        WRITE(*,100) IT,ERR,FAC,(KMAX(J),ERMAX(J),J=1,NE)
        IF(ERR.LT.CONV) RETURN
16    CONTINUE
      PAUSE 'ITMAX exceeded'
100   FORMAT(1X,I4,2F12.6,(/5X,I5,F12.6))
      RETURN
      END
