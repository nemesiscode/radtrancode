      subroutine polar2cartesian(R,theta,phi,v)
      implicit none
      real R,theta,phi,dtr,theta1,phi1,v(3),pi
      parameter (pi=3.1415927)

      dtr = pi/180.0
      theta1 = theta*dtr
      phi1 = phi*dtr

      v(1) = R*sin(theta1)*cos(phi1)
      v(2) = R*sin(theta1)*sin(phi1)
      v(3) = R*cos(theta1)

      return

      end
