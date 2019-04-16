    
      real function interp_ph3_xsections(nu,temp)
C*************************************************************
C     Function to read in and interpolate ph3 absorption data.
C
C     Input variables
C	  nu	real	Required wavenumber
C	  temp	real	Required temperature
C
C     Output variable
C	  interp_ph3_xsections	real x-section (cm2/molecule)
C 
C     Henrik Melin		Original	10/6/18 
C     Based on interp_ozone_serdyuchenko.f
C
C*************************************************************
      implicit none
      integer i,k,iread,icheck,nwave, ntemps, ioff, joff, offset
      character*100 ipfile,buffer
      real xioff,xjoff,temp1,temp,nu,xl
      real x,y1,y2,y3,y4, wres, tstart, tend, tres, tmin, tmax

C     These papameters define the contents of the x section file
      parameter (nwave=101)
      parameter (ntemps=2)
      parameter (offset=10)
      real wph3x(nwave),kph3x(nwave,ntemps),tempph3x(ntemps)
      common /ph3xtable/wph3x,tempph3x,kph3x      

      ipfile = 'uvXsect_PH3_Chen_1991.dat'

C     Define the temperatures of the x section columns
      tempph3x(1:2) = (/ 155.0, 295.0 /)

C     Check if the array has been populated
      if(int(wph3x(1)).ne.130) then
       print*,'Reading in PH3 absorption data'
       
C*************************************************************
C No need to edit below this line

       call datarchive(ipfile)

       open(12,file=ipfile,status='old')

       do 10 i=1,offset
        read(12,1)buffer
10     continue
1      format(a)

       do 20 i=1,nwave
        read(12,*)wph3x(i),(kph3x(i,k),k=ntemps,1,-1)
20     continue

       close(12)

      endif
      
      wres = wph3x(2) - wph3x(1)
      tres = tempph3x(2) - tempph3x(1)
      tmin = tempph3x(1) 
      tmax = tempph3x(1) 

C     Check that we're inside the wavelength range
      xl = 1e3*1e4/nu
      if(xl.lt.wph3x(1).or.xl.gt.wph3x(nwave))then
       interp_ph3_xsections=0.0
       return
      endif

      ioff = 1+int((xl-wph3x(1))/wres)
      if(ioff.eq.nwave) ioff=nwave-1
      xioff=(xl-wph3x(ioff))/wres

C     Check that we're inside the temperature range
      temp1=temp
      if(temp.lt.tmin)then
        temp1=tmin
      endif
      if(temp.gt.tmax)then
        temp1=tmax
      endif
      joff=1+int((temp1-tmin)/tres)
      if(joff.eq.size(tempph3x)) joff=1
      xjoff = (temp1-tempph3x(joff))/tres

C      print*,'joff,xjoff',joff,xjoff

      y1 = kph3x(ioff,joff)
      y2 = kph3x(ioff+1,joff)
      y3 = kph3x(ioff+1,joff+1)
      y4 = kph3x(ioff,joff+1)

      x = (1.0-xioff)*(1-xjoff)*y1 + xioff*(1.0-xjoff)*y2
      x = x + xioff*xjoff*y3 + (1.0-xioff)*xjoff*y4

      interp_ph3_xsections = x

      return
      
      end
