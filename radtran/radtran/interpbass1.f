      real function interpbass1(nu,temp)
C     *************************************************************
C     Function to read in and interpolate Bass UV ozone absorption data.
C
C     Pat Irwin		Original	1999?
C     Pat Irwin		Checked		6/6/12
C
C     *************************************************************
      implicit none
      integer i,k,iread,icheck
      real f,xabs1,xabs2,T
      real bassk(5,1900),xl,temp,nu

      common /basstable/bassk,iread

C     Read in Bass table if not already read in before.
      if(iread.ne.-1) then
       print*,'Reading in Bass O3 UV absorption data'
       call readbass1(icheck)
      endif

C     Convert wavenumbers to wavelengths (nm)
      xl = 1e3*1e4/nu

      do i=1,1899
       if(xl.ge.bassk(1,i).and.xl.lt.bassk(1,i+1))then
        k=i
        f = (xl - bassk(1,i))/(bassk(1,i+1)-bassk(1,i))
        goto 10
       endif
      enddo

      if(xl.lt.bassk(1,1))then
C       print*,'interpbass1 : xl too small',xl,bassk(1,1)
       interpbass1=0.0
       return
      else
C       print*,'interpbass1 : xl too big',xl,bassk(1,1900)
       interpbass1=0.0
       return
      endif

10    continue

      
      T = temp-273.15

      xabs1 = bassk(3,k)+bassk(4,k)*T+bassk(5,k)*T**2
      xabs2 = bassk(3,k+1)+bassk(4,k+1)*T+bassk(5,k+1)*T**2

      interpbass1 = ((1.0-f)*xabs1 + f*xabs2)*1e-20

      return
      
      end
