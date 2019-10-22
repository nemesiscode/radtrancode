      real function bilinear(y,mx,my,j,k,t,u)
C     *********************************************************
C     Simple bilinear interpolation code.
C
C     Pat Irwin   11/10/18
C
C     *********************************************************

      implicit none
      integer mx,my,j,k
      real y(mx,my),yout,t,u

      yout = (1.0-t)*(1.0-u)*y(j,k) + t*(1.0-u)*y(j+1,k)+
     1 t*u*y(j+1,k+1) + (1.0-t)*u*y(j,k+1)

      bilinear = yout

      return

      end
      
