      subroutine traceray(nlay,rad,refrac,startlay,r0,vec,ethick,
     1 scaleH,npath,trace,jtrace,disttrace,scaletrace,vectrace)
C     ***************************************************************
C     Subroutine to calculate ray path through a refracting, spherical 
C     atmosphere.
C 
C     Input variables
C	nlay		integer	Number of layers in atmosphere
C	rad(maxlay+1)	real	Radii of layer boundaries (km)
C	refrac(maxlay)	real	Refractive index of each layer
C	startlay	integer	Ordinate of starting radius
C	r0(3)		real	Initial position vector
C	vec(3)		real	Initial direction vector 
C	ethick		real	Minimum distance to travel before registering
C				change of radius
C	scaleH		real	Mean scale height of atmosphere (km)
C
C     Output variables
C	vec(3)		real	Final direction vector
C	npath		integer	Number of points along ray path
C	trace(maxpat,3)	real	Position vectors of intersection with radii
C				 boundaries
C	jtrace		integer	Ordinates of radii intersections
C	disttrace(maxpat) real	Distance from previous point in path
C	scaletrace(maxpat) real 	Airmass of layers along path.
C	vectrace(maxpat,3) real		Path vector at start of each
C					intersection 
C     Pat Irwin	Original	20/3/15
C	
C     ***************************************************************
      implicit none
      include '../includes/arrdef.f'
      integer ioff,i
      integer nlay,startlay,npath,jtrace(maxpat),jlay,ilay
      real rad(maxlay+1),refrac(maxlay),r0(3),vec(3),ethick,scaleH      
      real trace(maxpat,3),disttrace(maxpat),scaletrace(maxpat)
      real dtr,rnow,sproduct,theta,vec1(3),p0,thick,p1,thick1,f
      real r1(3),r2(3),dist,scale,n1,n2,theta1,theta2,xx,dtheta
      real rmin,p,vectrace(maxpat,3)
      parameter (dtr=3.1415927/180.)

C      ethick=1.
      rnow=sqrt(sproduct(r0,r0))
      theta=acos(sproduct(r0,vec)/rnow)/dtr

      jlay=startlay
      ioff=0
      do i=1,3
       r1(i)=r0(i)
      enddo
      dist=0.
      scale=0.

10    continue

      ioff=ioff+1
      jtrace(ioff)=jlay
      disttrace(ioff)=dist
      scaletrace(ioff)=scale
      do i=1,3 
       trace(ioff,i)=r0(i)
       vectrace(ioff,i)=vec(i)
      enddo

C     See if ray is going up or down and determine if it needs refracting or not
      if(theta.gt.90.0)then
C      ray going down
       if(jlay.eq.nlay+1) then
C       starting at top.
        n1=1.0
        n2=refrac(nlay)
       else
        if(jlay.eq.1) then
C        we have reached the ground
         goto 100
        else
         n1=refrac(jlay)
         n2=refrac(jlay-1)
        endif
       endif
      else
C      ray going up
       if(jlay.eq.1) then
        n1=refrac(jlay)
        n2=refrac(jlay)
       else
        if(jlay.eq.nlay+1) then
C        we have reached space
         goto 100
        else
         n1=refrac(jlay-1)
         n2=refrac(jlay)
        endif
       endif
      endif

C     deflect ray
      theta1=theta
      xx = n1*sin(theta1*dtr)/n2
      if(xx.le.1.0) then
        theta2=asin(n1*sin(theta1*dtr)/n2)/dtr
      else
C        print*,'Critical angle exceeded.'
C        print*,'Assume ray goes straight on'
        theta2=theta1
      endif
      if(theta.gt.90.0)theta2=180.0-theta2
      dtheta=theta2-theta1
      call deflectray(vec,r0,dtheta,vec1)
      do i=1,3
       vec(i)=vec1(i)
      enddo

C     find next intersection with concentric spheres
      call findnear(jlay,nlay,rad,r0,vec,ethick,ilay,r1,dist,
     &  theta)

      if(ilay.ne.jlay)then
       scale=dist/(abs(rad(ilay)-rad(jlay)))
      else
C      find midpoint
       do i=1,3
        r2(i)=0.5*(r0(i)+r1(i))
       enddo
       rmin=sqrt(sproduct(r2,r2))
       p0=1.
       thick=rad(ilay)-rad(ilay-1)
       p1=exp(-thick/scaleH)
       thick1 = rmin-rad(ilay-1)
       p=exp(-thick1/scaleH)
       f=(p-p1)/(p0-p1)
       scale=f*dist/thick
      endif

      jlay=ilay
      do i=1,3
       r0(i)=r1(i)
      enddo

C     Find next intersection point of ray
      goto 10

C     Finish up and exit
100   continue

      npath=ioff

      return

      end
