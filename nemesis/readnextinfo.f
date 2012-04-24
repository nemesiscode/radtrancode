      subroutine readnextinfo(linfo,altbore,marsradius,satrad,thetrot)
C     $Id:
C     ****************************************************************
C       
C     Subroutine to read in MCS pointing infor file.
C
C     Input variables
C	linfo		integer		file unit number (already open)
C
C     Output variables
C	altbore		real		altitude of bore sight above surface
C	marsradius	real		radius of mars at tangent point
C	satrad		real		radius of satellite orbit at time
C					of observation
C	thetrot		real		rotation angle of pixel array
C
C     Pat Irwin		1/11/07
C
C     ****************************************************************
      implicit none
      integer linfo
      real altbore,marsradius,satrad,thetrot

      integer iblock,ilimb,inad
      real sclk,p_lat,p_lon,p_rad,p_alt,sc_rad,p_sol_ang,p_time
      real p_rot,p_latnad,p_lonnad,p_radnad,p_altnad,sc_radnad
      real p_gssep 
   
      read(linfo,*) iblock,sclk,ilimb,p_lat,p_lon,p_rad,p_alt,
     1  sc_rad,p_sol_ang,p_time,p_rot,inad,p_latnad,p_lonnad,
     2  p_radnad,p_altnad,sc_radnad,p_gssep


      marsradius = p_rad
      altbore = p_alt
      satrad = sc_rad
      thetrot = p_rot

      print*,'readnextinfo : lat,long = ',p_lat,p_lon
      print*,'readnextinfo : radius,satrad,altbore,thetrot : ',
     1 marsradius,satrad,altbore,thetrot

      return

      end
