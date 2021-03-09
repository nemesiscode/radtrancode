    
      real function interp_c2h2_xsections(nu,temp)
C*************************************************************
C     Function to read in and interpolate c2h2 absorption data.
C
C     Input variables
C	  nu	real	Required wavenumber
C	  temp	real	Required temperature
C
C     Output variable
C	  interp_c2h2_xsections	real x-section (cm2/molecule)
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
      parameter (nwave=9001)
      parameter (ntemps=2)
      parameter (offset=11)
      real wc2h2x(nwave),kc2h2x(nwave,ntemps),tempc2h2x(ntemps)
      common /c2h2xtable/wc2h2x,tempc2h2x,kc2h2x      
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      ipfile = 'uvXsect_C2H2_smith_1991.dat'
      ipfile = 'uvXsect_C2H2_smith_1991_benilan_2000.dat'

C     Define the temperatures of the x section columns
      tempc2h2x(1:2) = (/195.0, 295.0/)

C     Check if the array has been populated
      if(int(wc2h2x(1)).ne.140) then
       if(idiag.gt.0)print*,'Reading in C2H2 absorption data'
       
C*************************************************************
C No need to edit below this line

       call datarchive(ipfile)

       open(12,file=ipfile,status='old')

       do 10 i=1,offset
        read(12,1)buffer
10     continue
1      format(a)

       do 20 i=1,nwave
        read(12,*)wc2h2x(i),(kc2h2x(i,k),k=ntemps,1,-1)
20     continue

       close(12)

      endif
      
      wres = wc2h2x(2) - wc2h2x(1)
      tres = tempc2h2x(2) - tempc2h2x(1)
      tmin = tempc2h2x(1) 
      tmax = tempc2h2x(1) 

C     Check that we're inside the wavelength range
      xl = 1e3*1e4/nu
      if(xl.lt.wc2h2x(1).or.xl.gt.wc2h2x(nwave))then
       interp_c2h2_xsections=0.0
       return
      endif

      ioff = 1+int((xl-wc2h2x(1))/wres)
      if(ioff.eq.nwave) ioff=nwave-1
      xioff=(xl-wc2h2x(ioff))/wres

C     Check that we're inside the temperature range
      temp1=temp
      if(temp.lt.tmin)then
        temp1=tmin
      endif
      if(temp.gt.tmax)then
        temp1=tmax
      endif
      joff=1+int((temp1-tmin)/tres)
      if(joff.eq.size(tempc2h2x)) joff=1
      xjoff = (temp1-tempc2h2x(joff))/tres

C      print*,'joff,xjoff',joff,xjoff

      y1 = kc2h2x(ioff,joff)
      y2 = kc2h2x(ioff+1,joff)
      y3 = kc2h2x(ioff+1,joff+1)
      y4 = kc2h2x(ioff,joff+1)

      x = (1.0-xioff)*(1-xjoff)*y1 + xioff*(1.0-xjoff)*y2
      x = x + xioff*xjoff*y3 + (1.0-xioff)*xjoff*y4

      interp_c2h2_xsections = x

      return
      
      end
