      real function interpem(nem,vem,emissivity,vconv)
C     $Id:
C     *************************************************
C     Routine to interpolate surface emissivity spectrum
C 
C     Input variables
C	nem	integer	Number of points in emissivity spectrum
C	vem(nem) real	Wavenumbers of table
C	emissivity(nem) real Emissivities
C 	vconv	real	Required output wavenumber
C
C     Pat Irwin	25/11/03
C     *************************************************
      integer nem
      real vem(nem),emissivity(nem),vconv

      if(vconv.lt.vem(1).or.vconv.gt.vem(nem))then
        print*,'Surface emissivity file not consistent with'
        print*,'convolution wavelength : ',vconv
        print*,'vem : ',vem(1),vem(nem)
        stop
      endif

      do i=1,nem-1
       if(vconv.ge.vem(i).and.vconv.lt.vem(i+1))then
        f = (vconv-vem(i))/(vem(i+1)-vem(i))
        x = (1-f)*emissivity(i)+f*emissivity(i+1)
       endif
      enddo

      if(vconv.eq.vem(nem))x=emissivity(nem)

      interpem = x

      return

      end
