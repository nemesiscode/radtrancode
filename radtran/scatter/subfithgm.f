      subroutine subfithgm(nphase,theta,phase,f,g1,g2,rms)
C     $Id: subfithgm.f,v 1.3 2011-06-17 15:57:55 irwin Exp $
C     ******************************************************************
C     Fits a combined Henyey-Grrenstein function to an ASCII phase function
C     file using Levenburg-Marquardt technique
C
C     Pat Irwin		13/5/97
C
C     ******************************************************************
      implicit none
      real f,g1,g2

      integer max_thet,nphase,mx,i,j,MY,itemp,nover,nc
      parameter (max_thet=100,mx=3,MY=100)
      real theta(max_thet),phase(max_thet),x(3),rms
      REAL VCONV(MY),Y(MY),SEI(MY,MY),ALPHA(MX,MX),BETA(MX),
     1 CHISQ,YMOD,SIG2I,DY,WT,KK(MY,MX),CPHASE(MY),COVAR(MX,MX),
     2 ALAMDA,OCHISQ



      x(1)=0.5
      x(2)=0.5
      x(3)=-0.5

      alamda = -1
      


      nover=1000
      nc = 0
      ochisq = 0.0
      do 1000 itemp=1,nover

         call mrqmin(nphase,theta,phase,x,covar,alpha,beta,chisq,
     1ochisq,alamda)
         
         
          if(chisq.eq.ochisq)then
            nc=nc+1
          else
            ochisq=chisq
            nc=0
          endif

1000  continue

      f=x(1)
      g1=x(2)
      g2=x(3)

      rms = sqrt(chisq)

      return
      end

