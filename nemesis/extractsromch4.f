      subroutine extractsromch4(runname,nvar,varident,varparam,nx,xn,
     1 ncloud,cpbot,cptop,nlaycloud,codepth,cfsh)
C     **************************************************************
C     Subroutine to extract cloud definition parameters for Sromovsky
C     cloud model (223)
C
C     Pat Irwin		19/6/14
C
C     **************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
 
      integer nvar,nx,mcloud,ivar
      integer ix,np,npvar
      parameter (mcloud=10)
      integer nlaycloud(mcloud),ncloud
      real cpbot(mcloud),cptop(mcloud),codepth(mcloud),
     1 cfsh(mcloud),pcut
      real xn(mx),varparam(mvar,mparam)
      integer varident(mvar,3)
      character*100 runname
      integer npro,nvmr
      logical gasgiant
      common /srom223/pcut

      call readrefhead(runname,npro,nvmr,gasgiant)

      ix=1
      do 10 ivar=1,nvar
       if(varident(ivar,1).ne.223)then
C       Skip to right point in xn array
        np=1
        if(varident(ivar,1).le.100)then
          np = npvar(varident(ivar,3),npro,varparam(ivar,1))
        endif
        if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))      
        if(varident(ivar,1).eq.887)np = int(varparam(ivar,1))      
        if(varident(ivar,1).eq.444)np = 2+int(varparam(ivar,1))
        if(varident(ivar,1).eq.445)np = 3+int(varparam(ivar,1))
        if(varident(ivar,1).eq.222)np = 8
        if(varident(ivar,1).eq.224)np = 9
        if(varident(ivar,1).eq.227)np = 7
        ix=ix+np

       else

C       LTC
        cpbot(1)=exp(xn(ix))/1.013
        cptop(1)=cpbot(1)*0.98
        nlaycloud(1)=1
        cfsh(1)=1.
        ix=ix+1
        codepth(1)=exp(xn(ix))
      
C       MTC
        ix=ix+1
        cpbot(2)=exp(xn(ix))/1.013
        cptop(2)=cpbot(2)*0.93
        nlaycloud(2)=1
        cfsh(2)=1.
        ix=ix+1
        codepth(2)=exp(xn(ix))

C       UTC
C        cpbot(3)=exp(xn(ix))/1.013
C       Set CH4 cloud at where pressure, pcut, where depleted CH4 profile 
C       meets SVP curve. pcut comes in via an ugly common block and has 
C       previously been calculated in subprofretg.

        cpbot(3)=pcut
        cptop(3)=cpbot(3)*0.93
        nlaycloud(3)=1
        cfsh(3)=1.
C       Skip pch4 where depletion starts and depletion factor to get 
C       odepth 
        ix=ix+3
        codepth(3)=exp(xn(ix))

C       Upper tropospheric haze
        ix=ix+1
        codepth(4)=exp(xn(ix))
        cpbot(4)=varparam(ivar,1)/1.013
        cptop(4)=varparam(ivar,2)/1.013
        nlaycloud(4)=int(varparam(ivar,4))
        cfsh(4)=1.

C       Stratospheric haze
        ix=ix+1
        codepth(5)=exp(xn(ix))
        cpbot(5)=varparam(ivar,2)/1.013
        cptop(5)=varparam(ivar,3)/1.013
        nlaycloud(5)=int(varparam(ivar,5))
        cfsh(5)=1.

        ncloud=5

       endif

10    continue

      return

      end
