      real function calcsoltau(Radius,R_loc,ilayer,nlayer,baseH,
     1  delH,tau,zen_sun)
C     *****************************************************************
C     Function to calculate the optical depth to the Sun.
C
C     Input variables
C    	Radius		real	Radius of planet at 0km altitude
C	R_loc		real	Local altitude
C	ilayer		integer	Local layer number
C	nlayer		integer	Total number of vertical layers in model
C	baseH(maxlay)	real	Base altitudes of layers
C	delH(maxlay)	real	Layer heights
C	tau(maxlay)	real	nadir optical depth of each layer
C	zen_sun		real	Local solar zenith angle
C
C     Output variables
C	 calcsoltau	real	Optical path to Sun.
C
C     Pat Irwin	17/6/11	Original Version
C     Pat Irwin	29/2/12	Updated for Radtrans2.0
C
C     *****************************************************************

      include '../includes/arrdef.f'
      integer i,ilayer,nlayer
      real pi
      parameter (pi = 3.1415927)
      real Radius, R_loc, baseH(maxlay),tau(maxlay),zen_sun
      real dist(maxlay), delH(maxlay), sf(maxlay),tausum,R1,R0
      dtr = pi/180.0

      R0 = R_loc
      tausum=0.0
      do 10 ilev=ilayer,nlayer 
       R1 = Radius + baseH(ilev)+delH(ilev)
       dist(ilev) = sqrt(R1**2-(R0*sin(zen_sun*dtr))**2) -
     &     R0*cos(zen_sun*dtr)
       if(ilev.eq.ilayer)then
           sf(ilev)=dist(ilev)/delH(ilev)
       else
           sf(ilev)=(dist(ilev)-dist(ilev-1))/delH(ilev) 
       endif
       tausum = tausum+sf(ilev)*tau(ilev)
10    continue 
     
      calcsoltau = tausum

      return 

      end



