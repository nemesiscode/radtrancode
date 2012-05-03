      subroutine roughk(ndata,x,y,ng,del_g,k_g)
C     ****************************************************************
C     Subroutine to fit a rough k-distribution to a transmission curve
C     prior to more elaborate fitting algorithms
C     ****************************************************************
      implicit none
      integer ndata,ng,i,j,k
      real x(20),y(20),del_g(10),k_g(10),sum1,sum2
      real k1(10),k2(10),sum,tmp,dpexp,y1,y2

C     Method 1
      do 10 i=1,ng      
       j = 20-2*(i-1)
       k1(i)=-alog(y(j))/x(j) 
10    continue

C     Method 2
C      do 20 i=1,ng
C       j = 20-2*(i-1)
C       if(i.eq.1)then
C        k2(i)=-alog(y(j)/del_g(i))/x(j)
C       else
C        sum=0.0
C        do 15 k=1,i-1
C         sum = sum + del_g(k)*dpexp(-x(j)*k2(k))
C15      continue
C        tmp = (y(j)-sum)/del_g(i)
C        if(tmp.gt.0.0)then
C          k2(i) = -alog(tmp)/x(j) 
C        else
C          k2(i) = k1(i)
C        endif
C       endif
C20    continue
C     Sort just in case
C      call sort(ng,k2)

C      do i=1,ng
C       print*,i,k1(i),k2(i)
C      enddo

C      sum1 = 0.0
C      sum2 = 0.0
C      do i=1,ndata
C       y1 = 0.0
C       y2 = 0.0
C       do k=1,ng
C        y1 = y1+del_g(k)*dpexp(-k1(k)*x(i))
C        y2 = y2+del_g(k)*dpexp(-k2(k)*x(i))
C       enddo
C       sum1 = sum1 + (y1-y(i))**2
C       sum2 = sum2 + (y2-y(i))**2
C       print*,i,x(i),y(i),y1,y2
C      enddo

C      sum1 = sqrt(sum1/float(ndata))
C      sum2 = sqrt(sum2/float(ndata))
C      print*,'roughk',sum1,sum2
C      if(sum1.le.sum2)then
C       print*,'k1'
C      else
C       print*,'k2'
C      endif
C      do i=1,ng
C       if(sum1.le.sum2)then
C        k_g(i)=k1(i)
C       else
C        k_g(i)=k2(i)
C       endif
C       print*,k_g(i)
C      enddo

      do i=1,ng
        k_g(i)=k1(i)
      enddo      

      return

      end
