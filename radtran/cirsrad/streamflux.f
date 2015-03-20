      subroutine streamflux(nlays,nmu,mu1,wt1,radg,galb,ig,umif,uplf,
     1 fup,fdown,ftop)
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
C
C     ***************************************************************
      implicit none
      INCLUDE '../includes/arrdef.f'
      INCLUDE '../includes/constdef.f'
      integer nlays,nmu,i,ilay,jlay,ig
      real x,umif(maxmu,maxscatlay,maxf),uplf(maxmu,maxscatlay,maxf)
      real radg(maxmu),fup(maxlay,maxg),fdown(maxlay,maxg)
      real ftop(maxg),xfac,xnorm
      double precision mu1(maxmu),wt1(maxmu),galb


C     Add in bottom radiance to complete radiance field at base of
C     atmosphere. NB. Scattering flux defines uplf to be radiance going 
C     downwards out of bottom of layer and umif to be the radiance going
C     upwards out of the top. We change the meaning here to be the 
C     total radiance field (up and down) at the bottom of each layer.  
C
C     In addition, the order of layers of the umif,uplf arrays is from 
C     top to bottom. We here reverse the order to be from bottom to top
C     so that the alignment with the baseh array is correct.


C     First integrate fluxes over zenith solid angle


C     Find correction for any quadrature errors
      xfac=0.
      do i=1,nmu
       xfac=xfac+sngl(mu1(i)*wt1(i))
      enddo
C      print*,'xfac = ',xfac
      xnorm=PI/xfac
      xnorm=2*PI

C      print*,'Flux down out of bottom, flux up out of top'
C      print*,'Layer, Fdown, Fup'
C      do ilay=1,nlays+1
C       fup(ilay,ig)=0.
C       fdown(ilay,ig)=0.
C       do i=1,nmu
C        fdown(ilay,ig)=fdown(ilay,ig)+
C     &          sngl(mu1(i)*wt1(i))*uplf(i,ilay,1)
C        fup(ilay,ig)=fup(ilay,ig)+
C     &          sngl(mu1(i)*wt1(i))*umif(i,ilay,1)
C       enddo
C       print*,ilay,xnorm*fdown(ilay,ig),xnorm*fup(ilay,ig)
C      enddo

      do 10 ilay=1,nlays

C      reverse order of output so that 1 is for the bottom layer
       jlay=nlays+1-ilay

       if(ilay.eq.1)then
C       Bottom layer
        fdown(ilay,ig)=0.
        do i=1,nmu
         fdown(ilay,ig)=fdown(ilay,ig)+
     &		sngl(mu1(i)*wt1(i))*uplf(i,jlay,1)
        enddo

        fup(ilay,ig)=0.
        do i=1,nmu
C        Thermal emission from ground
        fup(ilay,ig)=fup(ilay,ig)+sngl(mu1(i)*wt1(i))*radg(i)*
     &		(1.0-sngl(galb))
C	 Reflected radiance from ground
        fup(ilay,ig)=fup(ilay,ig)+sngl(mu1(i)*wt1(i))*
     &        uplf(i,jlay,1)*sngl(galb)
        enddo

       else

        fdown(ilay,ig)=0
        fup(ilay,ig)=0

        do i=1,nmu
         fdown(ilay,ig)=fdown(ilay,ig)+sngl(mu1(i)*wt1(i))*
     &		uplf(i,jlay,1)
C         fup(ilay,ig)=fup(ilay,ig)+sngl(mu1(i)*wt1(i))*
C     &		umif(i,jlay,1)
         fup(ilay,ig)=fup(ilay,ig)+sngl(mu1(i)*wt1(i))*
     &		umif(i,jlay+1,1)
        enddo
C        print*,ilay,xnorm*fdown(ilay,ig),xnorm*fup(ilay,ig)
       endif
       
       if(ilay.eq.nlays)then

        ftop(ig)=0.
        do i=1,nmu
         ftop(ig)=ftop(ig)+xnorm*sngl(mu1(i)*wt1(i))*umif(i,jlay,1)
        enddo

       endif

       fdown(ilay,ig)=fdown(ilay,ig)*xnorm
       fup(ilay,ig)=fup(ilay,ig)*xnorm

10    continue

    

      return

      end
