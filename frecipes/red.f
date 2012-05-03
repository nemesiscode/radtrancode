      SUBROUTINE RED(IZ1,IZ2,JZ1,JZ2,JM1,JM2,JMF,IC1,JC1,JCF,KC,
     *    C,NCI,NCJ,NCK,S,NSI,NSJ)
      DIMENSION C(NCI,NCJ,NCK),S(NSI,NSJ)
      LOFF=JC1-JM1
      IC=IC1
      DO 14 J=JZ1,JZ2
        DO 12 L=JM1,JM2
          VX=C(IC,L+LOFF,KC)
          DO 11 I=IZ1,IZ2
            S(I,L)=S(I,L)-S(I,J)*VX
11        CONTINUE
12      CONTINUE
        VX=C(IC,JCF,KC)
        DO 13 I=IZ1,IZ2
          S(I,JMF)=S(I,JMF)-S(I,J)*VX
13      CONTINUE
        IC=IC+1
14    CONTINUE
      RETURN
      END
