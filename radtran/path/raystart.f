      subroutine raystart(inp_angle,inp_solangle,inp_phi,ipzen,
     1  nlay,rad,r0,vec,svec)
C     ****************************************************************
C     Subroutine to compute starting position and direction vector of 
C     line of sight and also the direction to the Sun given a requested 
C     viewing geometry.
C
C     Input variables
C	inp_angle	real	input viewing zenith angle
C	inp_solangle	real	input solar zenith angle
C	inp_phi		real	input azimuth angle
C	ipzen		integer	ipzen=0 for angles defined at surface
C				ipzen=1 for angles refined at TOA
C	nlay		integer	Number of layers in atmosphere
C	rad(maxlay+1)	real	Radii of layer boundaries
C
C     Output variables
C	r0(3)		real	Starting position vector of LOS
C	vec(3)		real	Direction vector of LOS
C	svec(3)		real	Direction vector of vector towards Sun
C
C     Pat Irwin	Original	20/3/15
C
C     ****************************************************************
      implicit none
      include '../includes/arrdef.f'
      include '../includes/constdef.f'
      real inp_angle,inp_solangle,inp_phi,rad(maxlay+1),dtr
      real angle,rmin,rbot,rtop,solangle,phi,xs1,ys1,zs1
      real xs2,ys2,zs2,vec(3),svec(3),r0(3),sur_angle,a,b
      real x1,throt,xnow,znow,aa,bb,cc,xx(2)
      parameter (dtr=pi/180.)
      integer ipzen,ilimb,nlay,ierr,igrad
 
      rbot=rad(1)
      rtop=rad(nlay+1)

C     Need solar direction to compute solar arrays.
      solangle=inp_solangle
      phi=inp_phi
C     compute vector towards Sun.
      xs1=sin(solangle*dtr)*cos((180.0-phi)*dtr)
      ys1=sin(solangle*dtr)*sin((180.0-phi)*dtr)
      zs1=cos(solangle*dtr)
      svec(1)=xs1
      svec(2)=ys1
      svec(3)=zs1      

C     See if ray hits surface
      ilimb=0
      if(ipzen.eq.2) then
        angle=inp_angle
        rmin=rtop*sin(angle*dtr)
        if(rmin.gt.rbot)then
         ilimb=1
        else
         sur_angle=asin(rtop*sin(angle*dtr)/rbot)/dtr
        endif
      else
        if(ipzen.eq.0)then
         sur_angle=inp_angle
         angle = asin(rbot*sin(pi-inp_angle*dtr)/rtop)/dtr
        else
         print*,'raystart: ipzen must be 0 or 2'
         stop
        endif
      endif

C     turn ray into form z=a+bx
      a=rtop
C     If angle=0, then ray is straight down with a gradient of infinity.
      igrad=0
      if(angle.gt.0.)then
        b=1.0/tan(angle*dtr)
      else
        igrad=1
      endif

      if(ilimb.eq.0)then
C       get ray to hit surface at x=0.
        if(angle.gt.0)then
         x1=(rbot-a)/b
        else
         x1=0.
        endif
        if(sur_angle.gt.0.)then
         a=rbot
         b=1.0/tan(sur_angle*dtr)
        else
         igrad=1
        endif
C       if IPZEN=2 then rotate frame to correct solar vector
        if(ipzen.eq.2)then
           throt=(sur_angle-angle)*dtr
           zs2=zs1*cos(throt)-xs1*sin(throt)
           xs2=xs1*cos(throt)+zs1*sin(throt)
           zs1=zs2
           xs1=xs2
           svec(1)=xs2
           svec(3)=zs2
        endif
      endif

      xnow=0.
      znow=rtop

      if(igrad.eq.0)then
       aa=1.0+b**2
       bb=2*a*b

       if(ilimb.eq.0)then
        cc=a**2-rtop**2
        call quadsolve(aa,bb,cc,xx,ierr)
        if(xx(1).gt.xx(2))then
          xnow=xx(1) 
        else 
          xnow=xx(2)
        endif
        znow=a+b*xnow
       endif
      endif

      r0(1)=xnow
      r0(2)=0.
      r0(3)=znow

      vec(1)=-sin(angle*dtr)
      vec(2)=0.
      vec(3)=-cos(angle*dtr)

      return

      end

