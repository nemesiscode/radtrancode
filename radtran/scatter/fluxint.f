      subroutine fluxint(nmu,mu,wt,umif,uplf,nf,ilay,f,g1,g2,
     1  theta0,phi0,rad)
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
      integer nf,nmu,ilay,ifou,imu,iphi
      real f,g1,g2,theta0,phi0,pi
      double precision mu(maxmu),wt(maxmu)
      real cangscat,henyey,aphi,up,down,rad,trad
      real uptheta,downtheta,cupscat,cdownscat,dphi
      real sum,test
      parameter (pi=3.1415927)

      dphi = 2*pi/100.0
      rad=0.0
      sum=0.0
      do 10 imu = 1,nmu
       uptheta = acos(sngl(mu(imu)))
       downtheta = pi-uptheta


       trad = 0.0
       test = 0.0
       do 35 iphi=1,100
         aphi = 0.01*2*pi*(iphi-1)
         up=0.0
         down=0.0
         do 25 ifou=0,nf
          up = up + umif(imu,ilay+1,ifou+1)*cos(ifou*aphi)
          down = down + uplf(imu,ilay,ifou+1)*cos(ifou*aphi)
25       continue
         cupscat = cangscat(theta0,phi0,uptheta,aphi)
         trad = trad + up*henyey(cupscat,f,g2,g2)*dphi
         test = test + henyey(cupscat,f,g2,g2)*dphi
         cdownscat = cangscat(theta0,phi0,downtheta,aphi)
         trad = trad + down*henyey(cdownscat,f,g1,g2)*dphi      
         test = test + henyey(cdownscat,f,g1,g2)*dphi
35     continue
       rad = rad+trad*sngl(wt(imu))/(4.0*pi)
10    continue

      return

      end    

      real function cangscat(theta0,phi0,theta,phi)
      real theta0,phi0,theta,phi
      real x(3),y(3)

      x(1) = sin(theta0)*cos(phi0)
      x(2) = sin(theta0)*sin(phi0)
      x(3) = cos(theta0)

      y(1) = sin(theta)*cos(phi)
      y(2) = sin(theta)*sin(phi)
      y(3) = cos(theta)

      sum = 0
      do i=1,3
       sum = sum+x(i)*y(i)
      enddo

      cangscat = sum

      return

      end
