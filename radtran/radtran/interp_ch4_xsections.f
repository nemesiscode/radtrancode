    
      real function interp_ch4_xsections(nu,temp)
C*************************************************************
C     Function to read in and interpolate ch4 absorption data.
C
C     Input variables
C	  nu	real	Required wavenumber
C	  temp	real	Required temperature
C
C     Output variable
C	  interp_ch4_xsections	real x-section (cm2/molecule)
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
      parameter (nwave=23)
      parameter (ntemps=1)
      parameter (offset=9)
      real wch4x(nwave),kch4x(nwave,ntemps),tempch4x(ntemps)
      common /ch4xtable/wch4x,tempch4x,kch4x      
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      ipfile = 'uvXsect_CH4_ChenWu_2004.dat'

C     Define the temperatures of the x section columns
      tempch4x(1:1) = (/ 295 /)

C     Check if the array has been populated
      if(int(wch4x(1)).ne.120) then
       if(idiag.gt.0)print*,'Reading in CH4 absorption data'
       
C*************************************************************
C No need to edit below this line

       call datarchive(ipfile)

       open(12,file=ipfile,status='old')

       do 10 i=1,offset
        read(12,1)buffer
10     continue
1      format(a)

       do 20 i=1,nwave
        read(12,*)wch4x(i),(kch4x(i,k),k=ntemps,1,-1)
20     continue

       close(12)

      endif
      
      wres = wch4x(2) - wch4x(1)
      tres = tempch4x(2) - tempch4x(1)
      tmin = tempch4x(1) 
      tmax = tempch4x(1) 

C     Check that we're inside the wavelength range
      xl = 1e3*1e4/nu
      if(xl.lt.wch4x(1).or.xl.gt.wch4x(nwave))then
       interp_ch4_xsections=0.0
       return
      endif

      ioff = 1+int((xl-wch4x(1))/wres)
      if(ioff.eq.nwave) ioff=nwave-1
      xioff=(xl-wch4x(ioff))/wres

C     Check that we're inside the temperature range
      temp1=temp
      if(temp.lt.tmin)then
        temp1=tmin
      endif
      if(temp.gt.tmax)then
        temp1=tmax
      endif
      joff=1+int((temp1-tmin)/tres)
      if(joff.eq.size(tempch4x)) joff=1
      xjoff = (temp1-tempch4x(joff))/tres

C      print*,'joff,xjoff',joff,xjoff

      y1 = kch4x(ioff,joff)
      y2 = kch4x(ioff+1,joff)
      y3 = kch4x(ioff+1,joff+1)
      y4 = kch4x(ioff,joff+1)

      x = (1.0-xioff)*(1-xjoff)*y1 + xioff*(1.0-xjoff)*y2
      x = x + xioff*xjoff*y3 + (1.0-xioff)*xjoff*y4

      interp_ch4_xsections = x

      return
      
      end
