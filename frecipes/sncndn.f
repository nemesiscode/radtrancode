      SUBROUTINE SNCNDN(UU,EMMC,SN,CN,DN)
      PARAMETER (CA=.0003)
      LOGICAL BO
      DIMENSION EM(13),EN(13)
      EMC=EMMC
      U=UU
      IF(EMC.NE.0.)THEN
        BO=(EMC.LT.0.)
        IF(BO)THEN
          D=1.-EMC
          EMC=-EMC/D
          D=SQRT(D)
          U=D*U
        ENDIF
        A=1.
        DN=1.
        DO 11 I=1,13
          L=I
          EM(I)=A
          EMC=SQRT(EMC)
          EN(I)=EMC
          C=0.5*(A+EMC)
          IF(ABS(A-EMC).LE.CA*A)GO TO 1
          EMC=A*EMC
          A=C
11      CONTINUE
1       U=C*U
        SN=SIN(U)
        CN=COS(U)
        IF(SN.EQ.0.)GO TO 2
        A=CN/SN
        C=A*C
        DO 12 II=L,1,-1
          B=EM(II)
          A=C*A
          C=DN*C
          DN=(EN(II)+A)/(B+A)
          A=C/B
12      CONTINUE
        A=1./SQRT(C**2+1.)
        IF(SN.LT.0.)THEN
          SN=-A
        ELSE
          SN=A
        ENDIF
        CN=C*SN
2       IF(BO)THEN
          A=DN
          DN=CN
          CN=A
          SN=SN/D
        ENDIF
      ELSE
        CN=1./COSH(U)
        DN=CN
        SN=TANH(U)
      ENDIF
      RETURN
      END
