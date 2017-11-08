      subroutine extractcb(runname,iflagcb,nvar,varident,varparam,
     1 xn,ncloud,cbpbot,cbptop,ncblaycloud,cbodepth,cbfsh)
C     ************************************************************************
C     Subroutine to extract cloud definition parameters for Creme Brulee model
C
C     ASB	Creation				07.12.2016
C     ASB	Added extended haze layer option	24.05.2017
C
C     ************************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
 
      integer nvar,mcloud,ivar
      integer ix,np,npvar,iflagcb
      parameter (mcloud=6)
      integer ncblaycloud(mcloud),ncloud
      real cbpbot(mcloud),cbptop(mcloud),cbodepth(mcloud),
     1 cbfsh(mcloud),pcut
      real xn(mx),varparam(mvar,mparam)
      integer varident(mvar,3)
      character*100 runname
      integer npro,nvmr
      logical gasgiant
      common /srom223/pcut 

      call readrefhead(runname,npro,nvmr,gasgiant)

	print*, 'extractCB'

      ix=1
      do 10 ivar=1,nvar
       if(varident(ivar,1).ne.iflagcb)then
C       Skip to right point in xn array
        np=1
        if(varident(ivar,1).le.100)then
          np = npvar(varident(ivar,3),npro,varparam(ivar,1))
        endif
        if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))      
        if(varident(ivar,1).eq.444)np = 2+int(varparam(ivar,1))
        if(varident(ivar,1).eq.445)np = 3+int(varparam(ivar,1))
        if(varident(ivar,1).eq.222)np = 8
        if(varident(ivar,1).eq.223)np = 9
        if(varident(ivar,1).eq.224)np = 9
        if(varident(ivar,1).eq.225)np = 11
        if(varident(ivar,1).eq.227)np = 7
        ix=ix+np

       else

C       LTC
        cbpbot(1)=exp(xn(ix))
        ix=ix+1
        cbodepth(1)=exp(xn(ix))
        ix=ix+1
        cbfsh(1)=exp(xn(ix))
	ix=ix+1
        cbptop(1)=exp(xn(ix))
        ncblaycloud(1)=12

C        Upper tropospheric haze
         ix=ix+1
         cbodepth(2)=exp(xn(ix))
         cbpbot(2)=cbptop(1)
         if(varparam(ivar,1).lt.1e-4.or.varparam(ivar,1).gt.0.8)then
C         Creme Brulee model as specified by Sromovsky et al (2017)
          cbptop(2)=(0.9*cbpbot(2))
          ncblaycloud(2)=2
         else
C         Option of a more extended haze 
          cbptop(2)=(varparam(ivar,1)*cbpbot(2))
          ncblaycloud(2)=5  
         endif
         cbfsh(2)=1.

C        Stratospheric haze
         ix=ix+1
         cbpbot(3)=exp(xn(ix))/1.013
         cbptop(3)=exp(xn(ix))*0.93
	 ix=ix+1
         cbodepth(3)=exp(xn(ix))
         ncblaycloud(3)=1
         cbfsh(3)=1.

         ncloud=3

       endif

10    continue

      return

      end
