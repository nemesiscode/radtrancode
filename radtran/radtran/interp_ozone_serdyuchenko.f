      real function interp_ozone_serdyuchenko(nu,temp)
C     *************************************************************
C     Function to read in and interpolate ozone absorption data.
C
C     Input variables
C	nu	real	Required wavenumber
C	temp	real	Required temperature
C
C     Output variable
C	interp_ozone_serdyuchenko	real x-section (cm2/mol)
C 
C     Pat Irwin		Original	14/4/16
C
C     *************************************************************
      implicit none
      integer i,k,iread,icheck,nozone,ioff,joff
      real xioff,xjoff,temp1,temp,nu,xl
      parameter (nozone=88668)
      real wozone(nozone),kozone(nozone,11),tempozone(11)
      real x,y1,y2,y3,y4
      common /serdozonetable/wozone,tempozone,kozone      
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet


C     Read in O3 table if not already read in before.
      if(int(wozone(1)).ne.213) then
       if(idiag.gt.0)then
        print*,'Reading in Serdyuchenko ozone absorption data'
       endif
       call read_ozone_serdyuchenko(icheck)
      endif

C     Convert wavenumbers to wavelengths (nm)
      xl = 1e3*1e4/nu
C      print*,'wavelength = ',xl
      if(xl.lt.wozone(1).or.xl.gt.wozone(nozone))then
C       print*,'Wavelength not covered by Serdyuchenko table'
C       print*,xl,wozone(1),wozone(nozone)
       interp_ozone_serdyuchenko=0.0
       return
      endif

      ioff = 1+int((xl-wozone(1))/0.01)
      if(ioff.eq.nozone) ioff=nozone-1
      xioff=(xl-wozone(ioff))/0.01
C      print*,'ioff,xioff',ioff,xioff

      temp1=temp
      if(temp.lt.193.0)then
C        print*,'Temperature < 193K. Snapping to 193K'
        temp1=193.0
      endif
      if(temp.gt.293.0)then
C        print*,'Temperature > 293K. Snapping to 293K'
        temp1=293.0
      endif
      joff=1+int((temp1-193.)/10.)
      if(joff.eq.11) joff=10
      xjoff = (temp1-tempozone(joff))/10.0

C      print*,'joff,xjoff',joff,xjoff

      y1 = kozone(ioff,joff)
      y2 = kozone(ioff+1,joff)
      y3 = kozone(ioff+1,joff+1)
      y4 = kozone(ioff,joff+1)

      x = (1.0-xioff)*(1-xjoff)*y1 + xioff*(1.0-xjoff)*y2
      x = x + xioff*xjoff*y3 + (1.0-xioff)*xjoff*y4

      interp_ozone_serdyuchenko = x

      return
      
      end
