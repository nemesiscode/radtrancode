      real function slvf(k,aa,bb,B,gg)
C     $Id:
C     **************************************************************
C     Single precision k-distribution utility routine used by kml, intrpk,
C     calc_mlk1_k, and calc_mlk2_k.
C     The function dslvf calculates the function (g(ke1) - gg)
C     given the LO parameters 'a', 'b' and 'B'
C
C     Pat Irwin ?/?/??  Original
C     Pat Irwin 26/4/12 Commented
C
C     **************************************************************

      implicit none
      real k,gg
      double precision aa,bb,B,pi
      double precision xx,yy,tmp,derf
      parameter (pi=3.1415927)
      external derf

C      print*,'slvf,k,aa,bb,B,gg',k,aa,bb,B,gg

      xx=aa/dble(sqrt(k))
      yy=bb*dble(sqrt(k))

C      print*,'xx,yy,pi*B,gg',xx,yy,pi*B,gg
C      print*,'xx+yy,derf(xx+yy)',xx+yy,derf(xx+yy)
C      print*,'xx-yy,derf(xx-yy)',xx-yy,derf(xx-yy)
C      print*,'dexp(pi*B)',dexp(pi*B)

      if(derf(xx+yy).lt.1.0D0)then
       tmp=0.5*(1.0D0-derf(xx-yy)) + 0.5*(1.0D0-derf(xx+yy))*dexp(pi*B)
     1 - dble(gg)
      else
       tmp=0.5*(1.0D0-derf(xx-yy)) - dble(gg)
      end if

      slvf=sngl(tmp)
C      print*,'slvf = ',slvf
      return

      end

