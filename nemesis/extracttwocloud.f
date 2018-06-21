      subroutine extracttwocloud(runname,iflagsrom,nvar,varident,
     1 varparam,nx,xn,ncloud,cpbot,cptop,nlaycloud,codepth,cfsh)
C     **************************************************************
C     Subroutine to extract cloud definition parameters for two vertically
C     thin cloud model 226
C
C     Pat Irwin	19/6/14	Original
C     Pat Irwin	13/1/15	Revised.
C     Pat Irwin 14/2/17 Revised for two cloud
C
C     **************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
 
      integer nvar,nx,mcloud,ivar,i
      integer ix,np,npvar,iflagsrom
      parameter (mcloud=10)
      integer nlaycloud(mcloud),ncloud
      real cpbot(mcloud),cptop(mcloud),codepth(mcloud),
     1 cfsh(mcloud)
      real xn(mx),varparam(mvar,mparam)
      integer varident(mvar,3)
      character*100 runname
      integer npro,nvmr
      logical gasgiant

      call readrefhead(runname,npro,nvmr,gasgiant)

      ix=1
      do 10 ivar=1,nvar
       if(varident(ivar,1).ne.iflagsrom)then
C       Skip to right point in xn array
        np=1
        if(varident(ivar,1).le.100)then
          np = npvar(varident(ivar,3),npro,varparam(ivar,1))
        endif
        if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))      
        if(varident(ivar,1).eq.887)np = int(varparam(ivar,1))      
        if(varident(ivar,1).eq.444)np = 2+int(varparam(ivar,1))
        if(varident(ivar,1).eq.445)np = 3+(2*int(varparam(ivar,1)))
        if(varident(ivar,1).eq.223)np = 8
        if(varident(ivar,1).eq.223)np = 9
        if(varident(ivar,1).eq.224)np = 9
        if(varident(ivar,1).eq.225)np = 11
        if(varident(ivar,1).eq.226)np = 8
        ix=ix+np

       else

C       LTC
C        cpbot(1)=exp(xn(ix))/1.013
C        ix=ix+1
C        cptop(1)=exp(xn(ix))/1.013
C        nlaycloud(1)=5
C        ix=ix+1
C        cfsh(1)=exp(xn(ix))
C        ix=ix+1
C        codepth(1)=exp(xn(ix))

        cpbot(1)=exp(xn(ix))/1.013
        ix=ix+1
C        cptop(1)=exp(xn(ix))/1.013
        cptop(1)=cpbot(1)*0.93
C        nlaycloud(1)=5
        nlaycloud(1)=1
        ix=ix+1
C        cfsh(1)=exp(xn(ix))
        cfsh(1)=1.0
        ix=ix+1
        codepth(1)=exp(xn(ix))
      
C       MTC
        ix=ix+1
        cpbot(2)=exp(xn(ix))/1.013
        ix=ix+1
        cptop(2)=exp(xn(ix))/1.013
        nlaycloud(2)=5
        ix=ix+1
        cfsh(2)=exp(xn(ix))
        ix=ix+1
        codepth(2)=exp(xn(ix))

        ncloud=2

       endif

10    continue

      return

      end
