      program testfracpara
      INTEGER NORMAL,JROT
      REAL CALCFRACPARA,F,X,TMOD(11)
      REAL OUTPUT1(6,11),OUTPUT2(6,11)
      DOUBLE PRECISION TEMP      

      F=0.25
      DO I=1,11
       TMOD(I)=50.+25.*(I-1)
      ENDDO

      DO 10 I=1,11
       TEMP=TMOD(I)
       OUTPUT1(1,I)=TMOD(I)
       OUTPUT2(1,I)=TMOD(I)
       DO 20 JROT=0,4
        NORMAL=0
        OUTPUT1(JROT+2,I)=CALCFRACPARA(NORMAL,F,JROT,TEMP)
        NORMAL=1
        OUTPUT2(JROT+2,I)=CALCFRACPARA(NORMAL,F,JROT,TEMP)
20     CONTINUE
10    CONTINUE

      print*,'Normal = 0'
      DO I=1,11
       print*,(output1(J,I),J=1,6)
      ENDDO
      print*,'Normal = 1'
      DO I=1,11
       print*,(output2(J,I),J=1,6)
      ENDDO
   
      end
