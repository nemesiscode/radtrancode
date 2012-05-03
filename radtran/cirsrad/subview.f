      subroutine subview(theta_sun,phi_sun,R_sat,theta_sat,phi_sat,
     1 R_loc,theta_loc,phi_loc,zen_sat,zen_sun,phi)
C     *****************************************************************
C     Subroutine to compute angles of beam along a Line of sight in
C     a spherical atmosphere.
C  
C     Input variables
C	theta_sun	real	Solar zenith angle (as seen from centre
C				of planet relative to either the tangent point
C				of the beam or the observation point) 
C	phi_sun		real	Azimuth angle position of sun
C	R_sat		real	Distance from satellite to centre of planet
C	theta_sat	real	Zenith angle position of satellite
C	phi_sat		real	Azimuth angle position of satellite
C	R_loc		real	Radial position of point being computed
C	theta_loc	real	Zenith angle position of local point
C	phi_loc		real	Azimuth angle position of local point
C
C     Output variables
C	zen_sat		real	Local zenith angle of satellite
C	zen_sun		real	Local zenith angle of Sun
C	phi		real	Local azimuth between Sun and satellite.
C
C     Pat Irwin	17/6/11	Original version
C     Pat Irwin	29/2/12	Updated for Radtrans2.0
C
C     *****************************************************************
      implicit none
      real theta_sun,phi_sun,R_sat,theta_sat,phi_sat
      real R_loc,theta_loc,phi_loc,zen_sat,zen_sun,phi,dtr,pi
      real solvec(3),viewloc(3),satloc(3),satvec(3)
      real d,zsat,zsun,x,tmp,dotprod
      integer i
      parameter(pi=3.1415927)

      dtr = pi/180.0

C     direction vector of sunlight
      call polar2cartesian(1.0,theta_sun,phi_sun,solvec)

C     position vector of point in planet's atmosphere 
      call polar2cartesian(R_loc,theta_loc,phi_loc,viewloc)

C     position vector of satellite
      call polar2cartesian(R_sat,theta_sat,phi_sat,satloc)

C     relative vector of satellite w.r.t. viewing position
      do i=1,3 
       satvec(i) = satloc(i) - viewloc(i)
      enddo
      d =sqrt(dotprod(satvec,satvec))

      do i=1,3
       satvec(i)=satvec(i)/d
      enddo

      zen_sat = acos(dotprod(satvec,viewloc)/R_loc)/dtr
      zen_sun = acos(dotprod(solvec,viewloc)/R_loc)/dtr

      zsat = zen_sat*dtr 
      zsun = zen_sun*dtr

C     Now need to find azimuth angle
      x = dotprod(satvec,solvec)
      tmp = (x - cos(zsat)*cos(zsun))/(sin(zsat)*sin(zsun))
      if(tmp.gt.1.0) tmp=1.0
      if(tmp.lt.-1.0) tmp=-1.0
      phi = 180-acos(tmp)/dtr

      return

      end

      real function dotprod(v1,v2)
      real v1(3),v2(3),s
      integer i

      s=0.0
      do i=1,3
       s = s +v1(i)*v2(i)
      enddo

      dotprod=s

      return

      end
