      subroutine fluxintf(nmu,mu,wt,nlayerf,basehf,umif,uplf,
     2  nf,baseh,f,g1,g2,theta0,phi0,rad)
C     **************************************************************
C     Subroutine to perform the integral
C
C      /+1  /2pi
C      |   |
C      |   | P(theta0,phi0,mu,phi)*I(mu,phi) d(mu)d(phi)
C      |   |
C     /-1 /0
C
C     Phase function P(alpha), calculated from Henyey-Greenstein 
C     coefficients.
C    
C     Input variables
C      	nmu	integer	Number of steps in zenith integration
C	mu(maxmu)  double precision  cos(zenith) values
C	wt(maxmu)  double precision  Zenith quadrature weights
C       umif(maxmu,maxscatlay,maxf) real Calculated upward internal 
C                                     radiances from each layer, where 
C				      top of atmosphere is layer 1
C	uplf(maxmu,maxscatlay,maxf) real Calculated downwards internal 
C				      radiances
C	nf	integer	Number of Fourier components in azimuth 
C			integration
C	ilay	integer	Layer ABOVE where integration to be performed
C	f	real	Henyey-Greenstein parameter
C	g1	real	Henyey-Greenstein parameter
C	g2	real	Henyey-Greenstein parameter
C       theta0  real	zenith angle of observer
C	phi0    real    azimuth angle of observer
C
C     Output variable
C 	rad	real	Integrated scattered radiance along (theta0,phi0)
C
C     Pat Irwin		Original	3/8/05
C               	
C     **************************************************************
      implicit none

      include '../includes/arrdef.f'
      real umif(maxmu,maxscatlay,maxf),uplf(maxmu,maxscatlay,maxf)
      integer nf,nmu,ilay,ifou,imu,iphi,nlayerf,i,jf,mphi,nphi
      real f,g1,g2,theta0,phi0,pi,baseh,basehf(maxlay),dtr
      double precision mu(maxmu),wt(maxmu)
      real cangscat,henyey,aphi,up,down,rad,trad
      real uptheta,downtheta,cupscat,cdownscat,dphi
      real sum,test,ff,xmif,xplf,dphi1
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

       nphi = 1+int((uptheta/dtr)/10.0)
       dphi = 2.0*pi/float(nphi)

       do 35 iphi=1,nphi
         aphi = dphi*(iphi-1)
         up=0.0
         down=0.0

         do 25 ifou=0,nf
          xmif = (1.0-ff)*umif(imu,jf,ifou+1)+ff*umif(imu,jf+1,ifou+1)
          xplf = (1.0-ff)*uplf(imu,jf,ifou+1)+ff*uplf(imu,jf+1,ifou+1)
          if(ifou.eq.0) then
           up = up + xmif
           down = down + xplf
          else
           up = up + xmif*2*cos(ifou*aphi)
           down = down + xplf*2*cos(ifou*aphi)
          endif
25       continue
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

      return

      end    

