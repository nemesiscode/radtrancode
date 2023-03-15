      subroutine streamflux(nlays,nmu,mu1,wt1,radg,lowbc,galb,ig,
     1 umif,uplf,fup,fdown,ftop)
C     ***************************************************************
C     Subroutine to format output of scloud11flux into a more usable
C     internal radiation field and output to a meta file.
C
C     Input variables
C	nlays	integer	Number of scattering layers
C	nmu	integer	Number of zenith quadrature angles
C	mu1(maxmu) double precision	zenith angle values
C	wt1(maxmu) double precision	zenith angle weights
C	radg(maxmu)	real	Ground radiance field
C	lowbc	integer	Set to 1 if Extra lambertian layer added
C       galb	double precision	Ground albedo	
C	ig	integer	G-ordinate
C	umif(maxmu,maxscatlay,maxf) real Radiance field going
C                                        upwards out of the top of each
C					 layer
C	uplf(maxmu,maxscatlay,maxf) real Radiance field going
C                                        downwards out of the bottom of each
C					 layer
C
C     Output variables
C	fup(maxlay,maxg) real	Upward flux at the bottom of each layer
C	fdown(maxlay,maxg) real Downward flux at the bottom of each layer
C       ftop(maxg)	real	Upward flux at top of top layer (i.e. to space)
C
C     Documented by Pat Irwin	13/3/14
C     Converted from dumpflux by Pat Irwin 13/3/14
C     Major overhaul and simplification by Pat Irwin 15/3/23
C
C     ***************************************************************
      implicit none
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'
      integer nlays,nmu,i,ilay,jlay,ig,lowbc
      real x,umif(maxmu,maxscatlay,maxf),uplf(maxmu,maxscatlay,maxf)
      real radg(maxmu),fup(maxlay,maxg),fdown(maxlay,maxg)
      real tup(maxlay,maxg),tdown(maxlay,maxg)
      real ftop(maxg),xfac,xnorm
      double precision mu1(maxmu),wt1(maxmu),galb


C     Add in bottom radiance to complete radiance field at base of
C     atmosphere if no Lambertian ground layer included in calculation
C     NB. Scattering flux defines uplf to be radiance going 
C     downwards out of bottom of layer and umif to be the radiance going
C     upwards out of the top. We change the meaning here to be the 
C     total radiance field (up and down) at the bottom of each layer.  
C
C     In addition, the order of layers of the umif,uplf arrays is from 
C     top to bottom. We here reverse the order to be from bottom to top
C     so that the alignment with the baseh array is correct.


C     First integrate fluxes over zenith solid angle


C     Find correction for any quadrature errors
C      xfac=0.
C      do i=1,nmu
C       xfac=xfac+sngl(mu1(i)*wt1(i))
C      enddo
C      print*,'xfac (should be 0.5) = ',xfac
C      xnorm=PI/xfac
      xnorm=2*PI

C     First use scattering convention of height order where we start at top and go
C     down.
      do 10 ilay=1,nlays
        tdown(ilay,ig)=0.
        tup(ilay,ig)=0.
        do i=1,nmu
         tdown(ilay,ig)=tdown(ilay,ig)+
     &		xnorm*sngl(mu1(i)*wt1(i))*uplf(i,ilay,1)

         tup(ilay,ig)=tup(ilay,ig)+sngl(mu1(i)*wt1(i))*
     &		xnorm*umif(i,ilay+1,1)
        enddo
10    continue 

C     If no extra layer added, assume upwards radiation at bottom of bottom layer is just the 
C     ground radiation.
      if(lowbc.eq.0)then
        tup(nlays,ig)=0.
        do i=1,nmu
         tup(nlays,ig)=tup(nlays,ig)+xnorm*sngl(mu1(i)*wt1(i))*radg(i)
        enddo
      endif

C     Set top upwards radiation to upwards flux from top layer
       
      ftop(ig)=0.
      do i=1,nmu
        ftop(ig)=ftop(ig)+xnorm*sngl(mu1(i)*wt1(i))*umif(i,1,1)
      enddo

C     Reverse order of final output arrays
      do 15 ilay=1,nlays
       jlay=nlays+1-ilay
       fdown(jlay,ig)=tdown(ilay,ig)
       fup(jlay,ig)=tup(ilay,ig)
15    continue

      return

      end
