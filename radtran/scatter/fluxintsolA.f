      subroutine fluxintsolA(nmu,mu,wt,nlayerf,basehf,radmi,radpl,
     2  baseh,f,g1,g2,theta0,Stheta0,phi0,radsol0,rad)
C     **************************************************************
C     Subroutine to perform the integral
C
C      /+1  /2pi
C    1 |   |
C    - |   | P(theta0,phi0,mu,phi)*I(mu,phi) d(mu)d(phi)  + Fsol(loc)*p
C   4pi|   |
C     /-1 /0
C
C
C     Phase function P(alpha), calculated from Henyey-Greenstein 
C     coefficients.
C    
C     Input variables
C      	nmu	integer	Number of steps in zenith integration
C	mu(maxmu)  double precision  cos(zenith) values
C	wt(maxmu)  double precision  Zenith quadrature weights
C	nlayerf	integer	Number of layers in internal rad field
C	basehf(nlayerf)	real  base heights for internal field
C	radmi(maxmu,maxscatlay,181) real Upwards internal radiation
C	radpl(maxmu,maxscatlay,181) real Downwards internal radiation
C       baseh	real	Height of layer for calculation
C	f	real	Henyey-Greenstein parameter
C	g1	real	Henyey-Greenstein parameter
C	g2	real	Henyey-Greenstein parameter
C       theta0  real	zenith angle of observer
C	Stheta0 real	zenith angle of Sun
C	phi0    real    azimuth angle of observer wrt sun
C	radsol0  real	Solar irradiance at local position in atmosphere
C
C     Output variable
C 	rad	real	Integrated scattered radiance along (theta0,phi0)
C
C     Pat Irwin		Original	3/8/05
C               	Revised		4/4/09
C     **************************************************************
      implicit none

      include '../includes/arrdef.f'
      real radmi(maxmu,maxscatlay,181),radpl(maxmu,maxscatlay,181)
      integer nmu,ilay,ifou,imu,iphi,nlayerf,i,jf,mphi,nphi,j0
      real f,g1,g2,theta0,phi0,pi,baseh,basehf(maxlay),dtr
      double precision mu(maxmu),wt(maxmu)
      real cangscat,henyey,aphi,up,down,rad,trad,xaphi
      real uptheta,downtheta,cupscat,cdownscat,dphi
      real sum,test,ff,xmif,xplf,dphi1,Stheta0,radsol0,cphase
      parameter (pi=3.1415927,mphi=10)

      dtr = pi/180.0
      dphi = 2.0*pi/float(mphi)
      rad=0.0
      sum=0.0

      do i=1,nlayerf-1
       if(baseh.ge.basehf(i).and.baseh.lt.basehf(i+1))then
        jf=i
        ff = (baseh-basehf(i))/(basehf(i+1)-basehf(i))
       endif
      enddo
      if(baseh.ge.basehf(nlayerf))then
       jf=nlayerf-1
       ff=1.0
      endif
      if(baseh.lt.basehf(1))then
       jf=1
       ff=0.0
      endif


      do 10 imu = 1,nmu
       uptheta = acos(sngl(mu(imu)))
       downtheta = pi-uptheta


       trad = 0.0
       test = 0.0

       nphi=360
       dphi = 2.0*pi/float(nphi)


       do 35 iphi=1,nphi
         aphi = dphi*(iphi-1)

         xaphi = aphi/dtr
         if(xaphi.gt.180.0)xaphi=360.0-xaphi
         if(abs(xaphi).gt.181) then
          print*,'fluxintsol - serious problem'
          print*,aphi,xaphi
          stop
         endif

         j0 = int(xaphi)
         up = (1.0-ff)*radmi(imu,jf,j0)+ff*radmi(imu,jf+1,j0)
         down = (1.0-ff)*radpl(imu,jf,j0)+ff*radpl(imu,jf+1,j0)

         cupscat = cangscat(theta0,phi0,uptheta,aphi)
         trad = trad + up*henyey(cupscat,f,g1,g2)*dphi
         test = test + henyey(cupscat,f,g1,g2)*dphi
         cdownscat = cangscat(theta0,phi0,downtheta,aphi)
         trad = trad + down*henyey(cdownscat,f,g1,g2)*dphi      
         test = test + henyey(cdownscat,f,g1,g2)*dphi
35     continue
       rad = rad+trad*sngl(wt(imu))/(4.0*pi)
       sum = sum+test*sngl(wt(imu))/(4.0*pi)
10    continue

      if(sum.ne.1.0)then
       rad = rad/sum
      endif

C     Now add single-scattered component of Sunlight

      aphi = 0.0
      cphase = -cangscat(theta0,phi0,Stheta0,aphi)

      rad = rad + radsol0*henyey(cphase,f,g1,g2)/(4.0*pi)

      return

      end    

