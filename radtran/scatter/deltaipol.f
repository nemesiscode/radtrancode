
      real function deltaipol(theta,omega,tau)
C      ****************************************
C      Routine to calculation delta_I after Eq A.1 of 
C      Sromovsky (2005) - polarisation
C
C      Input variables
C	theta	real	zenith angle (radians)
C	omega	real	single-scattering albedo
C	tau	real	opacity
C
C      Output variables
C	deltaipol	Required correction
C
C      Pat Irwin 9/9/21
C
C      ****************************************
       implicit none
       real theta,omega,tau
       real thetmax,thetwid,polimax,poli0
       real thmax,thwid,imax,i0,x1,x2,x,pi
       parameter (pi=3.1415927)

       thmax=thetmax(omega,tau)
       thwid=thetwid(omega,tau)
       imax=polimax(omega,tau)
       i0=poli0(omega,tau)

       if(theta.ge.0.0.and.theta.lt.thmax) then
         x1 = (1-exp(-((theta-thmax)/thwid)**2))
         x2 = 1.0/(1-exp(-(thmax/thwid)**2))
         deltaipol = imax - (imax-i0)*x1*x2
       else
         x = abs((theta-thmax)/(0.5*pi-thmax))
         x1 = (cos(0.5*pi*(x**1.15)))**0.6
         x2 = cos(0.5*pi*(x**1.4))
         deltaipol = imax*(x1*(1-omega**4) + x2*omega**4)
       endif

       return

       end

       real function poli0(omega,tau)
C      ****************************************
C      Routine to calculation I0 after Eq A.4 of 
C      Sromovsky (2005) - polarisation
C
C      Input variables
C	omega	real	single-scattering albedo
C	tau	real	opacity
C
C      Output variables
C	poli0	real	Reflectivity correction at theta=0
C
C      Pat Irwin 9/9/21
C
C      ****************************************
       real omega,tau,y,y1,y2,arg,wf

       arg=1.0-omega
       if(arg.ge.0.0)then
        y = abs((alog10(arg)))**1.789
        y1 = 1.0-exp(-y/0.458)
        if(isnan(y1))then
         y1=1.0
        endif
       else
        y1=1.0
       endif
C       wf=1.01
       wf=1.01*0.85
       y2 = (1-exp(-(tau**wf)/0.493))**(2.0/wf)

       poli0 = 0.0403*y1*y2

       return

       end


       real function polimax(omega,tau)
C      ****************************************
C      Routine to calculation f(omega) after Eq A.5 of 
C      Sromovsky (2005) - polarisation
C
C      Input variables
C	omega	real	single-scattering albedo
C	tau	real	opacity
C
C      Output variables
C	polimax	real	Reflectivity correction at theta=thetmax
C
C      Pat Irwin 9/9/21
C
C      ****************************************
       real tau,omega,y1,y2,fomega

       y1 = (1-exp(-(tau**1.617)/1.015))**0.6
       y2 = fomega(omega)

       polimax = 0.0457*y1*y2

       return

       end

       real function fomega(omega)
C      ****************************************
C      Routine to calculation f(omega) after Eq A.6 of 
C      Sromovsky (2005) - polarisation.
C
C      Input variables
C	omega	real	single-scattering albedo
C
C      Output variables
C	fomega	real	Sromovsky polarisation function
C
C      Pat Irwin 9/9/21
C
C      ****************************************
       real omega,y

       if(omega.gt.0.95) then
        fomega = omega**3 
       else
        y = abs((alog10(1.0-omega)))**1.789
        fomega = 0.882*(1-exp(-y/0.458))
       endif 

       return

       end

       real function thetwid(omega,tau)
C      ****************************************
C      Routine to calculation theta_wid after Eq A.3 of 
C      Sromovsky (2005) - polarisation
C
C      Input variables
C	omega	real	single-scattering albedo
C	tau	real	opacity
C
C      Output variables
C	thetwid	real	Angular width of polarisation correction maximum
C			at small angles (radians)
C      Pat Irwin 9/9/21
C
C      ****************************************
       real omega,tau,thet0,thet1,pi,x1,x2,thetmax
       parameter (pi=3.1415927)

       thet0 = 36.8*pi/180.0
       thet1 = 44.0*pi/180.0

       x1 = thet0*(1-thetmax(omega,tau)/(0.5*pi))**0.75
       x2 = (thetmax(omega,tau)/thet1)**1.3

       thetwid = x1*x2

       return

       end

       real function thetmax(omega,tau)
C      ****************************************
C      Routine to calculation theta_max after Eq A.2
C      of Sromovsky (2005) - polarisation
C
C      Input variables
C	omega	real	single-scattering albedo
C	tau	real	opacity
C
C      Output variables
C	thetmax	real	Angle at which polarisation correction is
C			maximum (radians)
C      Pat Irwin 9/9/21
C
C      ****************************************
       real omega,tau,pi,thetinf
       parameter (pi=3.1415927)

       thetinf = (44.0*pi/180.0)*omega**4
       thetmax = thetinf + (0.5*pi-thetinf)*exp(-tau/0.935)

       return

       end



