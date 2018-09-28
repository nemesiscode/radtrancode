    
      real function interp_co_xsections(nu,temp)
C*************************************************************
C     Function to read in and interpolate co absorption data.
C
C     Input variables
C	  nu	real	Required wavenumber
C	  temp	real	Required temperature
C
C     Output variable
C	  interp_co_xsections	real x-section (cm2/molecule)
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
      parameter (nwave=35)
      parameter (ntemps=1)
      parameter (offset=9)
      real wcox(nwave),kcox(nwave,ntemps),tempcox(ntemps)
      common /coxtable/wcox,tempcox,kcox      

      ipfile = 'uvXsect_CO_Chan_1993.dat'

C     Define the temperatures of the x section columns
      tempcox(1:1) = (/ 298 /)

C     Check if the array has been populated
      if(int(wcox(1)).ne.10) then
       print*,'Reading in CO absorption data'
       
C*************************************************************
C No need to edit below this line

       call datarchive(ipfile)

       open(12,file=ipfile,status='old')

       do 10 i=1,offset
        read(12,1)buffer
10     continue
1      format(a)

       do 20 i=1,nwave
        read(12,*)wcox(i),(kcox(i,k),k=ntemps,1,-1)
20     continue

       close(12)

      endif
      
      wres = wcox(2) - wcox(1)
      tres = tempcox(2) - tempcox(1)
      tmin = tempcox(1) 
      tmax = tempcox(1) 

C     Check that we're inside the wavelength range
      xl = 1e3*1e4/nu
      if(xl.lt.wcox(1).or.xl.gt.wcox(nwave))then
       interp_co_xsections=0.0
       return
      endif

      ioff = 1+int((xl-wcox(1))/wres)
      if(ioff.eq.nwave) ioff=nwave-1
      xioff=(xl-wcox(ioff))/wres

C     Check that we're inside the temperature range
      temp1=temp
      if(temp.lt.tmin)then
        temp1=tmin
      endif
      if(temp.gt.tmax)then
        temp1=tmax
      endif
      joff=1+int((temp1-tmin)/tres)
      if(joff.eq.size(tempcox)) joff=1
      xjoff = (temp1-tempcox(joff))/tres

C      print*,'joff,xjoff',joff,xjoff

      y1 = kcox(ioff,joff)
      y2 = kcox(ioff+1,joff)
      y3 = kcox(ioff+1,joff+1)
      y4 = kcox(ioff,joff+1)

      x = (1.0-xioff)*(1-xjoff)*y1 + xioff*(1.0-xjoff)*y2
      x = x + xioff*xjoff*y3 + (1.0-xioff)*xjoff*y4

      interp_co_xsections = x

      return
      
      end
