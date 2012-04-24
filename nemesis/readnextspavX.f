      subroutine readnextspavX(lspec,woff,xlat,xlon,ngeom,nav,ny,y,
     1 se,fwhm,nconv,vconv,angles,wgeom,flat,flon)
C     $Id:
C     ****************************************************************
C
C     Subroutine to read in measured CIRS spectrum vector and covariance 
C     matrix from the associated .spe file (assumed to be already open)
C
C     Input variables
C	lspec		integer		file unit number (already open)
C       woff            real            Additional wavlength offset (if
C                                               required)
C     Output variables
C	xlat		real		Mean latitude
C	xlon		real		Mean longitude
C	ngeom		integer		Number of viewing geometries
C	nav(mgeom)	integer		Number of synthetic spectra required
C                                       to simulate each FOV-averaged
C                                       measurement spectrum.
C	ny		integer		Number of measured points in spectrum
C	y(my)		real		Measured spectra
C	se(my)		real		Covariance matrix (assumed diagonal)
C       fwhm            real            FWHM of measured spectrum
C       nconv(mgeom)    integer         Number of points in each
C                                       convoluted spectrum (there may be
C                                       several spectra in y(my) at
C                                       different observation angles
C       vconv(mgeom,mconv) real         Convolution wavenumbers
C       angles(mgeom,mav,3) real        Observation geometry
C       wgeom(mgeom,mav)real            Averaging weights
C       flat(mgeom,mav) real            FOV point latitudes (or nearest limb)
C       flon(mgeom,mav) real            FOV point longitudes ( "          )
C
C     Original	Pat Irwin		19/5/97
C     Converted from NIMS PGJI		21/3/00
C     Tidied for Nemesis PGJI           22/10/03
C     Modified to read integration weights PGJI 8/2/04
C     Modified to new file name PGJI 	9/2/04
C     Modified for new treatment of FOV. PGJI  10/12/04

C     ****************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'

      integer i,j,ny

      real y(my),se(my),err,vconv(mgeom,mconv),yx,angles(mgeom,mav,3)
      real xlat,xlon,woff,x,wgeom(mgeom,mav),flat(mgeom,mav),fwhm
      real flon(mgeom,mav)
      integer nconv(mgeom),lspec,ngeom,igeom,nav(mgeom)

1     format(a)

      do i=1,my
        se(i)=0.0
      end do

      read(lspec,*)fwhm,xlat,xlon,ngeom
      if(ngeom.gt.mgeom)then
       print*,'Error in readnextspavX.f - ngeom > mgeom'
       print*,ngeom,mgeom
       stop
      endif

      ny = 0

      do 20 igeom=1,ngeom
       read(lspec,*)nconv(igeom)
       if(nconv(igeom).gt.mconv)then
        print*,'Error in readnextspavX.f - nconv(igeom) > mconv'
        print*,igeom,nconv(igeom),mconv
        stop
       endif
       read(lspec,*)nav(igeom)
       if(nav(igeom).gt.mav)then
        print*,'Error in readnextspavX.f - nav(igeom)>mav'
        print*,igeom,nav(igeom),mav
        stop
       endif

       do i=1,nav(igeom)
        read(lspec,*)flat(igeom,i),flon(igeom,i),
     1    (angles(igeom,i,j),j=1,3),wgeom(igeom,i)
       enddo
 
       do 10 i=1,nconv(igeom)
        read(lspec,*)x,yx,err
        vconv(igeom,i)=x+woff               ! Shift wavelengths if required

        ny=ny+1

        if(ny.gt.my)then
         print*,'Error in readnextspavX.f - ny > my'
         print*,igeom,i,ny,my
         stop
        endif
        y(ny)=yx
        se(ny)=err**2    

10     continue
20    continue
 
      return

      end
