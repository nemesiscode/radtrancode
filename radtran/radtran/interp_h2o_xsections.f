    
      real function interp_h2o_xsections(nu,temp)
C*************************************************************
C     Function to read in and interpolate h2o absorption data.
C
C     Input variables
C	  nu	real	Required wavenumber
C	  temp	real	Required temperature
C
C     Output variable
C	  interp_h2o_xsections	real x-section (cm2/molecule)
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
      parameter (nwave=81)
      parameter (ntemps=1)
      parameter (offset=9)
      real wh2ox(nwave),kh2ox(nwave,ntemps),temph2ox(ntemps)
      common /h2oxtable/wh2ox,temph2ox,kh2ox      

      ipfile = 'uvXsect_H2O_JPL_2011.dat'

C     Define the temperatures of the x section columns
      temph2ox(1:1) = (/ 298.0 /)

C     Check if the array has been populated
      if(int(wh2ox(1)).ne.120) then
       print*,'Reading in H2O absorption data'
       
C*************************************************************
C No need to edit below this line

       call datarchive(ipfile)

       open(12,file=ipfile,status='old')

       do 10 i=1,offset
        read(12,1)buffer
10     continue
1      format(a)

       do 20 i=1,nwave
        read(12,*)wh2ox(i),(kh2ox(i,k),k=ntemps,1,-1)
20     continue

       close(12)

      endif
      
      wres = wh2ox(2) - wh2ox(1)
      tres = temph2ox(2) - temph2ox(1)
      tmin = temph2ox(1) 
      tmax = temph2ox(1) 

C     Check that we're inside the wavelength range
      xl = 1e3*1e4/nu
      if(xl.lt.wh2ox(1).or.xl.gt.wh2ox(nwave))then
       interp_h2o_xsections=0.0
       return
      endif

      ioff = 1+int((xl-wh2ox(1))/wres)
      if(ioff.eq.nwave) ioff=nwave-1
      xioff=(xl-wh2ox(ioff))/wres

C     Check that we're inside the temperature range
      temp1=temp
      if(temp.lt.tmin)then
        temp1=tmin
      endif
      if(temp.gt.tmax)then
        temp1=tmax
      endif
      joff=1+int((temp1-tmin)/tres)
      if(joff.eq.size(temph2ox)) joff=1
      xjoff = (temp1-temph2ox(joff))/tres

C      print*,'joff,xjoff',joff,xjoff

      y1 = kh2ox(ioff,joff)
      y2 = kh2ox(ioff+1,joff)
      y3 = kh2ox(ioff+1,joff+1)
      y4 = kh2ox(ioff,joff+1)

      x = (1.0-xioff)*(1-xjoff)*y1 + xioff*(1.0-xjoff)*y2
      x = x + xioff*xjoff*y3 + (1.0-xioff)*xjoff*y4

      interp_h2o_xsections = x

      return
      
      end
