      subroutine mrqminl(nphase,theta,phase,x,covar,alpha,beta,chisq,
     1ochisq,alamda)
C     Fit phase function with henyey-greenstein parameters in log space

      implicit none
      integer max_thet,nphase,mx,i,j,k,MY
      parameter (max_thet=100,mx=3,MY=100)
      real theta(max_thet),phase(max_thet),x(3),xt(3)
      REAL VCONV(MY),Y(MY),SEI(MY,MY),ALPHA(MX,MX),BETA(MX),
     1 CHISQ,YMOD,SIG2I,DY,WT,KK(MY,MX),CPHASE(MY),COVAR(MX,MX),
     2 ALAMDA,OCHISQ,DA(MX)

      
      IF(ALAMDA.LT.0.)THEN

        call mrqcofl(nphase,theta,phase,x,alpha,beta,chisq)

        OCHISQ=CHISQ
        ALAMDA=1000.0

      ENDIF


      DO 15 J=1,MX
        DO 14 K=1,MX
          COVAR(J,K)=ALPHA(J,K)
14      CONTINUE
        COVAR(J,J)=ALPHA(J,J)*(1.0 + ALAMDA)
        DA(J)=BETA(J)
15    CONTINUE
      CALL GAUSSJ(COVAR,MX,MX,DA,1,1)
      IF(ALAMDA.EQ.0.0)THEN
        RETURN
      ENDIF


      do i=1,3
       xt(i)=x(i)+da(i)
       if(i.eq.1)then
C        if(xt(i).gt.0.99)xt(i)=0.99
C        if(xt(i).lt.0.01)xt(i)=0.01
        if(xt(i).gt.0.999999)xt(i)=0.999999
        if(xt(i).lt.0.000001)xt(i)=0.000001
       else if(i.eq.2)then
        if(xt(i).gt.0.98)xt(i)=0.98
        if(xt(i).lt.0.0)xt(i)=0.0
       else if(i.eq.3)then
        if(xt(i).lt.-0.98)xt(i)=-0.98
C        if(xt(i).gt.0.0)xt(i)=0.0
        if(xt(i).gt.-0.1)xt(i)=-0.1
       endif
      end do

      call mrqcofl(nphase,theta,phase,xt,covar,da,chisq)

      IF(CHISQ.LE.OCHISQ)THEN
        ALAMDA=0.9*ALAMDA
        OCHISQ=CHISQ
        DO 18 J=1,MX
          DO 17 K=1,MX
            ALPHA(J,K)=COVAR(J,K)
17        CONTINUE
          BETA(J)=DA(J)
          X(J)=XT(J)
18      CONTINUE
      ELSE
        ALAMDA=1.5*ALAMDA
        CHISQ=OCHISQ
        IF(ALAMDA.GT.1.0E36)ALAMDA=1.0E36
      ENDIF
      RETURN
      END
