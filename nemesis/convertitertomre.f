      program convertitertomre
C     **************************************************************
C     Program to convert last succesful iteration in the .itr file to
C     a .mre file for cases where the retrieval has fallen over.
C
C     Pat Irwin		25/5/16
C
C     **************************************************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      integer nx,ny,kiter,iter,iform,ispace,lout,ispec
      integer npro,nvar,varident(mvar,3),ngeom
      integer ilbl,nspec,ioff,i,j,nvmr,lspec,nav(mgeom)
      integer igeom,nconv(mgeom),nconv1,vconv1(mconv)
      integer k,nwave1,nwave(mgeom),lpre,jsurf
      integer jalb,jxsc,jtan,lx(mx)
      real xlat,xlon,varparam(mvar,mparam),fwhm,sa(mx,mx)
      real angles(mgeom,mav,3),wgeom(mgeom,mav),flat(mgeom,mav)
      real vkstart,vkend,vkstep,rerr(mgeom,mconv)
      real flon(mgeom,mav),xerr,vwave(mgeom,mwave)
      real vwave1(mwave)
      real vconv(mgeom,mconv),jpre,jlogg,jrad,iscat,lin,jfrac
      real xn1(mx),xa(mx),y(my),se1(my),yn1(my),yn(my)
      real chisq,phi,kk(my,mx),err1(mx),woff,phlimit
      real ochisq
      character*100 runname,itname,buffer,ename
      logical gasgiant

      CALL prompt('Enter run name : ')
      READ(5,1)buffer
1     FORMAT(a)
      runname = buffer(1:36)

      CALL readrefhead(runname,npro,nvmr,gasgiant)
      if(npro.gt.maxpro)then
       print*,'Error in Nemesis. npro > maxpro : ',npro,maxpro
       stop
      endif


      CALL file(runname,runname,'inp')
      OPEN(32,file=runname,status='old')

C     Read in whether to calculate with wavenumbers(0) or wavelength(1)
C     Also read in whether scattering is required (iscat)
C     Also read in whether lbl calculation is required (ilbl)
      READ(32,*)ispace,iscat,ilbl

      if(ilbl.gt.0) then
       print*,'Nemesis - LBL calculation'
      endif
      if(ilbl.eq.0)CALL readkkhead(runname,vkstart,vkend,vkstep)


C     Read any wavenumber offset to add to measured spectra
      READ(32,*)woff

C     Read in name of forward modelling error file
      READ(32,1)ename

C     Read in number of iterations
      READ(32,*)kiter

C     Read limiting % change in phi
      READ(32,*)phlimit

C     Read in total number of spectra to fit and starting offset
      READ(32,*)nspec,ioff

      READ(32,*)lin
      iform=0
      READ(32,*,END=999)iform
999   continue
      CLOSE(32)

      print*,'iform = ',iform


      lspec=36
      CALL file(runname,runname,'spx')
      open(lspec,file=runname,status='old')

C     skip first ioff-1 spectra
      do ispec=1,ioff-1

       call readnextspavX(lspec,iform,woff,xlat,xlon,ngeom,nav,ny,y,se1,
     & fwhm,nconv,vconv,angles,wgeom,flat,flon)
       print*,'A iform = ',iform
      enddo
    

      do 2999 ispec=ioff,ioff-1+nspec

C     Read in measurement vector, obs. geometry and covariances
      call readnextspavX(lspec,iform,woff,xlat,xlon,ngeom,nav,ny,y,
     1  se1,fwhm,nconv,vconv,angles,wgeom,flat,flon)

      print*,'B iform = ',iform

C     Read in forward modelling errors
      call forwarderr(ename,ngeom,nconv,vconv,woff,rerr)

C     Add forward errors to measurement covariances
      k=0
      DO i=1,ngeom
       do j=1,nconv(i)
        k = k+1
        xerr=rerr(i,j)
        if(iform.eq.3)xerr=xerr*1e-18
        se1(k)=se1(k)+xerr**2
       enddo
      ENDDO

C     If we are doing an SCR cell, then we only need half the
C     wavelengths in the spx file as we're calculating both WB and SB
C     outputs

      if(ilbl.eq.0)then
C     Calculate the tabulated wavelengths of c-k look up tables
      do igeom=1,ngeom
       nconv1=nconv(igeom)
       do j=1,nconv1
        vconv1(j)=vconv(igeom,j)
       enddo
       CALL wavesetb(runname,vkstart,vkend,vkstep,nconv1,vconv1,
     1  fwhm,nwave1,vwave1)
       do j=1,nwave1
        vwave(igeom,j)=vwave1(j)
       enddo
       nwave(igeom)=nwave1
      enddo
      endif

      print*,'C iform = ',iform

C     set up a priori of x and its covariance
      CALL readapriori(runname,lin,lpre,xlat,npro,nvar,varident,
     1  varparam,jsurf,jalb,jxsc,jtan,jpre,jrad,jlogg,jfrac,nx,xa,
     2  sa,lx)


      print*,'D iform = ',iform

C     Calculate retrieval errors.
C     Simple errors, set to sqrt of diagonal of ST
      do i=1,nx
       err1(i)=sqrt(abs(sa(i,i)))
      enddo

      call file(runname,itname,'itr')
      open(37,file=itname,status='unknown')
      read(37,*)nx,ny,kiter   

      ochisq=1e20

      do 401 iter=1,kiter
       call readnextiter(nx,ny,xn1,xa,y,se1,yn1,yn,
     &		chisq,phi)


       print*,iter,chisq,ochisq
       if(chisq.lt.ochisq)then
        lout=38
        CALL file(runname,runname,'mre')
        open(lout,file=runname,status='unknown')
        write(lout,*)nspec,' ! Total number of retrievals'

        print*,'Writing out mre :',iform
        CALL writeout(iform,runname,ispace,lout,ispec,xlat,xlon,npro,
     1   nvar,varident,varparam,nx,ny,y,yn,se1,xa,sa,xn1,err1,ngeom,
     2   nconv,vconv,gasgiant,jpre,jrad,jlogg,jfrac,iscat,lin)

        close(lout)
        ochisq=chisq

       endif

401   continue

2999  continue

      end
      
