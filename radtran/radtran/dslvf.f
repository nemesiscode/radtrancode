      double precision function dslvf(k,aa,bb,B,gg)
C     $Id:
C     **************************************************************
C     Double precision k-distribution utility routine used by kml, intrpk,
C     calc_mlk1_k, and calc_mlk2_k.
C     The function dslvf calculates the function (g(ke1) - gg)
C     given the LO parameters 'a', 'b' and 'B'
C
C     Pat Irwin	?/?/??	Original
C     Pat Irwin 26/4/12	Commented
C
C     **************************************************************

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

