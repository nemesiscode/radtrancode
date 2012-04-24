      double precision function dslvf(k,aa,bb,B,gg)
C     $Id: dslvf.f,v 1.2 2011-06-17 15:40:26 irwin Exp $

      implicit none
      double precision aa,bb,B,k,gg,pi,xx,yy,tmp,derf
      external derf

	pi = 2.0d0 * dasin(1.0d0)


      xx=aa/dsqrt(k)
      yy=bb*dsqrt(k)


      if(derf(xx+yy).lt.1.0D0)then
       tmp=0.5*(1.0D0-derf(xx-yy)) + 0.5*(1.0D0-derf(xx+yy))*dexp(pi*B)
     1 - dble(gg)
      else
       tmp=0.5*(1.0D0-derf(xx-yy)) - dble(gg)
      end if

      dslvf=tmp
C      print*,'dslvf = ',dslvf
      return

      end

