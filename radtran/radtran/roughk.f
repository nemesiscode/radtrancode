      subroutine roughk(ndata,x,y,ng,del_g,k_g)
C     ****************************************************************
C     Subroutine to fit a rough k-distribution to a transmission curve
C     prior to more elaborate fitting algorithms
C     ****************************************************************
      implicit none
      integer ndata,ng,i,j,k
      real x(20),y(20),del_g(10),k_g(10),sum1,sum2
      real sum,tmp,dpexp,y1,y2

C     Simple derivation of rough k-distribution
      do 10 i=1,ng      
       j = 20-2*(i-1)
       k_g(i)=-alog(y(j))/x(j) 
10    continue

      return

      end
