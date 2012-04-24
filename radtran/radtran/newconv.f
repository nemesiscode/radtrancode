      subroutine newconv(npoint,y,iconv,ispace,wid,vmin,delv,fwhm,
     1ispec)
C     $Id: newconv.f,v 1.4 2011-07-12 11:54:10 irwin Exp $
C     ************************************************************************
C
C     Subroutine to smooth a spectrum by a given function
C
C     Input variables:
C	npoint		integer	Number of points in spectrum
C	y(npoint)	real	Spectrum
C	iconv		integer	Shape of convolution function:
C				1 = square
C				2 = triangle
C				3 = gaussian
C				4 = sinc
C	ispace		integer	Convolution wavespace:
C				0 = wavenumber
C				1 = wavelength
C	wid		real	Required smoothing width
C	vmin		real	minimum wavenumber
C	delv		real 	wavenumber spacing of spectrum
C	fwhm		real	wavenumber resolution of spectrum
C	ispec		integer	Spectrum type:
C				0 = transmission
C				1 = radiance
C
C     Output variables
C	y(npoint)	real	Smoothed spectrum
C
C     Pat Irwin		7/3/97
C
C     ************************************************************************
      implicit none
      integer npoint,i,j,iconv,ispace,npoint1
      real delv,vmin,vmin1,fwhm,x1,x2,v1,v2,dx,sum
      integer mpoint,ispec,iend,istart
      parameter (mpoint = 90000)
      real wid,xm,y(npoint),y1(mpoint),smooth,x(mpoint),xfac


      if(npoint.gt.mpoint)then
        print*,'newconv. npoint > mpoint = ',npoint,mpoint
        stop
      end if

      if(iconv.lt.3)then
        xfac = 2.0
      else
        xfac = 5.0
      endif

      if(ispace.eq.0)then
C     Easy bit! Do calculation in wavenumber space

      do 10 i=1,npoint
       x(i) = vmin + (i-1)*delv
10    continue


      do 20 i=1,npoint
       xm = x(i)
       x1 = xm-xfac*wid
       x2 = xm+xfac*wid
       istart = int(( x1 - vmin)/delv)
       iend = 1+ int(( x2 - vmin)/delv)
       if(istart.lt.1)istart = 1
       if(iend.gt.npoint)iend = npoint
       y1(i)=0.0
       sum=0.0
       do 30 j=istart,iend
        dx = x(j)-xm
        y1(i)=y1(i)+y(j)*smooth(iconv,dx,wid)
        sum=sum+smooth(iconv,dx,wid)
30     continue
       y1(i)=y1(i)/sum
20    continue
      
      do 40 i=1,npoint
       y(i)=y1(i)
40    continue
      
      else
C     Hard bit! Do calculation in wavelength space

      do 15 i=1,npoint
       x1 = vmin + (i-1)*delv
       x(i) = 1e4/x1
       if(ispec.eq.1)then
        y(i)=y(i)*1e4/(x(i)*x(i))
       endif
15    continue



      do 25 i=1,npoint
       xm = x(i)

       x1 = xm-xfac*wid
       x2 = xm+xfac*wid

       v1 = 1e4/x1
       v2 = 1e4/x2

       istart = int(( v2 - vmin)/delv)
       iend = 1+ int(( v1 - vmin)/delv)
       if(istart.lt.1)istart = 1
       if(iend.gt.npoint)iend = npoint

       y1(i)=0.0
       sum=0.0
       do 35 j=istart,iend
        dx = x(j)-xm
        y1(i)=y1(i)+y(j)*smooth(iconv,dx,wid)
        sum=sum+smooth(iconv,dx,wid)
35     continue
       y1(i)=y1(i)/sum
25    continue
      
      do 45 i=1,npoint
       if(ispec.eq.1)then
        y(i)=y1(i)*x(i)*x(i)/1e4
       else
        y(i)=y1(i)
       endif
45    continue

      endif

      return
      end

      real function smooth(iconv,dx,fwhm)
C     ****************************************************************
C     Subroutine which returns a normalised convolution function
C
C     Input variables
C	iconv	integer	convolution type
C	dx	real	distance from centre of convolution function
C	fwhm	real	Full width, half maximum of convolution function
C
C     ****************************************************************


      integer iconv
      real dx,fwhm,x,pi,xD,sig
      parameter (pi=3.1415927)

      if(iconv.eq.1)then
       if(abs(dx).le.0.5*fwhm)then
        smooth = 1.0/fwhm
       else
        smooth = 0.0
       endif

      elseif(iconv.eq.2)then

       x1 = 1.0 - abs(dx/fwhm)
     
       if(x1.ge.0)then
        smooth = x1/fwhm
       else
        smooth = 0.0
       endif

      elseif(iconv.eq.3)then

       xD = 2.*sqrt(log(2.0))
       sig = fwhm/xD
       smooth = exp(-(dx/sig)**2)


C       smooth = exp(-(dx/fwhm)**2)/(1.772*fwhm)

      elseif(iconv.eq.4)then
       if(dx.eq.0.0)then
        smooth = 1.0/fwhm
       else
        x = pi*dx/fwhm
        smooth = (sin(x)/x)**2
        smooth = smooth/fwhm
       endif
      end if


      return
      end
 
