      real function interpO3(nu,temp)
C     *************************************************************
C     Function to read in and interpolate ozone absorption data.
C
C     Pat Irwin		Original	1999?
C     Pat Irwin		Checked		6/6/12
C     Dan Dawson    Updated to Serdyuchenko 2012 o3 data   9/4/13
C
C     *************************************************************
      implicit none
      integer i,k,iread,icheck
      real f,xabs1,xabs2,T
      real o3k(5,88668),xl,temp,nu

      common /o3table/o3k,iread

C     Read in O3 table if not already read in before.
      if(iread.ne.-1) then
       print*,'Reading in o3  UV absorption data'
       call readO3(icheck)
      endif

C     Convert wavenumbers to wavelengths (nm)
      xl = 1e3*1e4/nu

      do i=1,88668
       if(xl.ge.o3k(1,i).and.xl.lt.o3k(1,i+1))then
        k=i
        f = (xl - o3k(1,i))/(o3k(1,i+1)-o3k(1,i))
        goto 10
       endif
      enddo

      if(xl.lt.o3k(1,1))then
C       print*,'interpO3 : xl too small',xl,o3k(1,1)
       interpO3=0.0
       return
      else
C       print*,'interpO3 : xl too big',xl,o3k(1,88668)
       interpO3=0.0
       return
      endif

10    continue

      
      T = temp-273.15

      xabs1 = o3k(3,k)+o3k(4,k)*T+o3k(5,k)*T**2
      xabs2 = o3k(3,k+1)+o3k(4,k+1)*T+o3k(5,k+1)*T**2

      interpO3 = ((1.0-f)*xabs1 + f*xabs2)*1e-20

      return
      
      end
