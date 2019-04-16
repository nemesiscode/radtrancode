    
      real function interp_c2h6_xsections(nu,temp)
C*************************************************************
C     Function to read in and interpolate c2h6 absorption data.
C
C     Input variables
C	  nu	real	Required wavenumber
C	  temp	real	Required temperature
C
C     Output variable
C	  interp_c2h6_xsections	real x-section (cm2/molecule)
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
      parameter (nwave=41)
      parameter (ntemps=1)
      parameter (offset=9)
      real wc2h6x(nwave),kc2h6x(nwave,ntemps),tempc2h6x(ntemps)
      common /c2h6xtable/wc2h6x,tempc2h6x,kc2h6x      

      ipfile = 'uvXsect_C2H6_Lee_2001.dat'

C     Define the temperatures of the x section columns
      tempc2h6x(1:1) = (/ 295 /)

C     Check if the array has been populated
      if(int(wc2h6x(1)).ne.120) then
       print*,'Reading in C2H6 absorption data'
       
C*************************************************************
C No need to edit below this line

       call datarchive(ipfile)

       open(12,file=ipfile,status='old')

       do 10 i=1,offset
        read(12,1)buffer
10     continue
1      format(a)

       do 20 i=1,nwave
        read(12,*)wc2h6x(i),(kc2h6x(i,k),k=ntemps,1,-1)
20     continue

       close(12)

      endif
      
      wres = wc2h6x(2) - wc2h6x(1)
      tres = tempc2h6x(2) - tempc2h6x(1)
      tmin = tempc2h6x(1) 
      tmax = tempc2h6x(1) 

C     Check that we're inside the wavelength range
      xl = 1e3*1e4/nu
      if(xl.lt.wc2h6x(1).or.xl.gt.wc2h6x(nwave))then
       interp_c2h6_xsections=0.0
       return
      endif

      ioff = 1+int((xl-wc2h6x(1))/wres)
      if(ioff.eq.nwave) ioff=nwave-1
      xioff=(xl-wc2h6x(ioff))/wres

C     Check that we're inside the temperature range
      temp1=temp
      if(temp.lt.tmin)then
        temp1=tmin
      endif
      if(temp.gt.tmax)then
        temp1=tmax
      endif
      joff=1+int((temp1-tmin)/tres)
      if(joff.eq.size(tempc2h6x)) joff=1
      xjoff = (temp1-tempc2h6x(joff))/tres

C      print*,'joff,xjoff',joff,xjoff

      y1 = kc2h6x(ioff,joff)
      y2 = kc2h6x(ioff+1,joff)
      y3 = kc2h6x(ioff+1,joff+1)
      y4 = kc2h6x(ioff,joff+1)

      x = (1.0-xioff)*(1-xjoff)*y1 + xioff*(1.0-xjoff)*y2
      x = x + xioff*xjoff*y3 + (1.0-xioff)*xjoff*y4

      interp_c2h6_xsections = x

      return
      
      end
