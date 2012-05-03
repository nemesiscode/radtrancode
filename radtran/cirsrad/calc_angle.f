      subroutine calc_angle(sol_ang,emiss_ang,aphi,Radius,nlev,baseHS,
     1 ibaseH,nlayer,baseH,delH,tau,tzsat,tzsun,tphi,tausun)
C     ***************************************************************
C     Subroutine to calculate the local angles of a beam in a 
C     spherical atmosphere
C
C     Input variables:
C	sol_ang		real	solar zenith angle 
C	emiss_ang	real	viewing zenith angle
C	aphi		real	azimuth angle
C	Radius		real	Planet radius (at 0km altitude)
C	nlev		integer	Number of layers along LOS	
C	baseHS(maxlay)	real	base height of layers LOS
C	ibaseH(maxlay)	integer	level number of LOS layer
C	nlayer		integer number of vertical layers in atmosphere
C	baseH(maxlay)	real	Base height of layers (km)
C	delH(maxlay)	real	Vertical height of layers (km)
C       tau(maxlay)	real 	Opacity of layers
C
C     Output variables
C	tzsat(maxlay)	real	Local satellite zenith angle at each layer
C				along LOS.
C	tzsun(maxlay)	real	Local solar zenith angle at each layer
C				along LOS.
C	tphi(maxlay)	real	Local azimuth angle between satellite and
C				sun along LOS
C	tausun(maxlay)	real	solar opacity at each point along the LOS.
C
C 
C     Pat Irwin	17/6/11	Original version
C     Pat Irwin	29/2/12	Updated for Radtrans2.0
C
C     ***************************************************************
      implicit none
      include '../includes/arrdef.f'
 
      integer i,nlev,nlayer,ilayer
      real sol_ang,emiss_ang,aphi,Radius,baseHS(maxlay)
      real tau(maxlay),tausun(maxlay),delH(maxlay),baseH(maxlay)
      real tzsat(maxlay),tzsun(maxlay),tphi(maxlay),pi,dotprod
      real delH1(maxlay),arg
      integer ibaseH(maxlay)
      parameter (pi=3.1415927)
C     New Nemesis coding is that if emiss_ang is -ve, then a limb path
C     is defined where sol_ang holds the tangent altitude and sol_ang
C     is equal to -emiss_ang
      real th1,th2,th3,phi,satvec(3),zen_sun,zen_sat,theta_sat
      real theta_loc,phi_loc,R_loc,R_sat,phi_sat
      real theta_sun,phi_sun,htan,dtr,calcsoltau

      dtr = pi/180.0 
 
      do i=1,nlayer-1
       delH1(i)=baseH(i+1)-baseH(i)
      enddo
      i=nlayer
      delH1(i)=delH(i)

     
      if(emiss_ang.lt.0.0) then 
C      Limb geometry
       theta_sun = -emiss_ang
       phi_sun = aphi
       htan = sol_ang
       R_sat = sqrt((Radius+htan)**2 + (3*Radius)**2)
       phi_sat = 180.0
       theta_sat = atan((3*Radius)/(Radius+htan))/dtr 

       R_loc = Radius
       theta_loc = 0.0
       phi_loc = 0.0

       do 10 i=1,nlev
        R_loc = Radius+baseHS(i)
        arg = (Radius+htan)/(Radius+baseHS(i))
        if(arg.ge.1.0)then
         theta_loc = 0.0
        else
         theta_loc = acos(arg)/dtr
        endif
        phi_loc = 180.0
        if(i.gt.nlev/2) phi_loc = 0.0
 
        call subview(theta_sun,phi_sun,R_sat,theta_sat,phi_sat,R_loc,
     1   theta_loc,phi_loc,zen_sat,zen_sun,phi)

        tzsat(i) = zen_sat
        tzsun(i) = zen_sun
        tphi(i) = phi
        ilayer = ibaseH(i)
        tausun(i)= calcsoltau(Radius,R_loc,ilayer,nlayer,baseH,
     1    delH1,tau,zen_sun)

10     continue

      else
C      Nadir and off-nadir geometry

       th1 = emiss_ang*dtr
       satvec(1)=0.0-10*Radius*sin(th1)
       satvec(2)=0.0
       satvec(3)=Radius + 10*Radius*cos(th1)

       R_sat=sqrt(dotprod(satvec,satvec))
       theta_sat = asin(10*Radius*sin(th1)/R_sat)/dtr
       phi_sat = 180.0

       theta_sun = sol_ang
       phi_sun = aphi

       do 20 i=1,nlev

        th2 = asin(Radius*sin(th1)/(Radius+baseHS(i)))
        th3 = th1-th2
        R_loc = Radius+baseHS(i)
        ilayer = ibaseH(i)
        theta_loc=th3/dtr
        phi_loc = 180.0

        call subview(theta_sun,phi_sun,R_sat,theta_sat,phi_sat,
     1   R_loc,theta_loc,phi_loc,zen_sat,zen_sun,phi)

        tzsat(i) = zen_sat
        tzsun(i) = zen_sun
        tphi(i) = phi

        tausun(i)= calcsoltau(Radius,R_loc,ilayer,nlayer,baseH,
     1    delH1,tau,zen_sun)

20     continue

      endif

      return

      end
