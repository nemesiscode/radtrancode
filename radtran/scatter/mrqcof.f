      subroutine mrqcof(nphase,theta,phase,x,alpha,beta,chisq)

      integer max_thet,nphase,mx,i,j,MY
      parameter (max_thet=100,mx=3,MY=100)
      real theta(max_thet),phase(max_thet),x(3)
      REAL VCONV(MY),Y(MY),SEI(MY,MY),ALPHA(MX,MX),BETA(MX),
     1 CHISQ,YMOD,SIG2I,DY,WT,KK(MY,MX),CPHASE(MY)

      DO 12 J=1,MX
        DO 11 K=1,J
          ALPHA(J,K)=0.
11      CONTINUE
        BETA(J)=0.
12    CONTINUE
      CHISQ=0.

      call subhgphas(nphase,theta,x,cphase,kk)

C      print*,nphase
C      print*,theta
C      print*,x
      DO 15 I=1,NPHASE
        DY=PHASE(I)-CPHASE(I)
C        print*,phase(i),cphase(i),kk(i,1),kk(i,2),kk(i,3)
        DO 14 J=1,MX
          WT=KK(I,J)
          DO 13 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+WT*KK(I,K)
13        CONTINUE
          BETA(J)=BETA(J)+DY*WT
14      CONTINUE
        CHISQ=CHISQ+DY*DY
15    CONTINUE
      DO 17 J=2,MX
        DO 16 K=1,J-1
          ALPHA(K,J)=ALPHA(J,K)
16      CONTINUE
17    CONTINUE

C      print*,'alpha'
C      print*,(alpha(1,i),i=1,3)
C      print*,(alpha(2,i),i=1,3)
C      print*,(alpha(3,i),i=1,3)
C      print*,'beta'
C      print*,beta
      
      RETURN
      END
