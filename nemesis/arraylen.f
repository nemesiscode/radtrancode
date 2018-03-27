C     Lengths of measurement and vector arrays included by all other routines
C     in Nemesis
C     $Id:
C     ****************************************************************
      integer mx,mvar,my,mwave,mparam,mgeom,mav,mconv,jmod
      real lat_tolerance
C     mx is the maximum size of the state vector
C     my is the maximum size of the measurement vector
C     mparam is number of extra elements in varparam
C     mvar is the maximum number of variable profiles
C     mgeom is the maximum number of observation geometries per retrieval
C     mav is the maximum number of ordinates in FOV calculation
C     mwave is the maximum number of calculation wavelengths
C     mconv is the maximum number of convolution wavelengths
C	lat_tolerance is the maximum allowable lat diff between subsequent
C	  retrieval stages
 
C     Set mx to equal maxv, defined in arrdef.f
C     Set my to idim, defined in arrdef.f
C     Set mwave to equal maxbin, defined in arrdef.f
C     Set mconv to idim, defined in arrdef.f
      parameter (my=idim,mconv=idim,mx=maxv,mwave=maxbin)

C     Set other parameters
      parameter (mparam=20,mvar=12,mgeom=100,mav=61)
      parameter (lat_tolerance=3.0)

C     ****************************************************************

