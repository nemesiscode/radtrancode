      subroutine computeFOVA(iread,ichan,ipixA,ipixB,thcentre,
     1 thbore,thetrot,nfov,thfov,rfov)
      implicit none
      integer ichan,ipixA,ipixB,mfov,nfov,iread,ipix
      integer k1,k2,k3,jangle,iangle,j,k,i1
      parameter (mfov=200)
      real thfov(mfov),rfov(mfov),x,f,fvcen
      real th1(mfov),th2(mfov),r1(mfov),r2(mfov)
      real thint(mfov),rint1(mfov),rint2(mfov)
      double precision fovs(11,9,21,2,mfov)
      real thcentre,thbore,thetrot,fovcentre(11,9,21)
      character*100 buffer,aname
      common /mcsfov/fovs,fovcentre

      print*,'ComputeFOVA, iread,ichan,ipixA,ipixB,thcentre,thbore'
      print*, iread,ichan,ipixA,ipixB,thcentre,thbore
      nfov=50
      if(iread.eq.1)then
       aname='fov_50.dat'
       call datarchive(aname)
       open(12,file=aname,status='old')
       do iangle=1,11
        read(12,1)buffer
1       format(a)
        do j=1,9
         do k=1,21
          read(12,*)x,fovcentre(iangle,j,k),k1,k2,k3
          do i1=1,nfov
           read(12,*)fovs(iangle,j,k,1,i1),fovs(iangle,j,k,2,i1)
          enddo
         enddo
        enddo
       enddo
      endif

      if(ichan.le.6) then
       ipix = ipixA
      else
       ipix = ipixB
      endif

      if(thetrot.ge.1.0) then
       print*,'Cant extrapolate rotation'
       jangle=10
       f=1.0
      elseif(thetrot.lt.-1.0) then
       print*,'Cant extrapolate rotation'
       jangle=1
       f=0.0
      else
       jangle = int(5*(thetrot+1.0))
       f = 5*(thetrot+1.0)-jangle
       jangle=jangle+1
      endif

      print*,'computeFOVA thetrot, jangle, f  :',thetrot,jangle,f
      fvcen = (1-f)*fovcentre(jangle,ichan,ipix)+
     1  f*fovcentre(jangle+1,ichan,ipix)

      do i1=1,nfov 
       th1(i1)= fovs(jangle,ichan,ipix,1,i1)
       th2(i1)= fovs(jangle+1,ichan,ipix,1,i1)
       r1(i1) = fovs(jangle,ichan,ipix,2,i1)
       r2(i1) = fovs(jangle+1,ichan,ipix,2,i1)

C      angle arrays may not be identical so we need to interpolate one 
C      FOV array on to the other. Set angles to array with greatest
C      weight.
       if(f.gt.0.5) then
        thint(i1)=th2(i1)
       else
        thint(i1)=th1(i1)
       endif
      enddo

C     Now interpolate FOVs in to same angle array
      call fovinterp(nfov,thint,th1,r1,rint1)
      call fovinterp(nfov,thint,th2,r2,rint2)

      do i1=1,nfov
       rfov(i1)=(1.0-f)*rint1(i1) + f*rint2(i1)
C      reverse angles to Nemesis assumption and add in boresight angle 
       thfov(i1)=thbore-thint(i1)
      enddo

      return

      end
      

      subroutine fovinterp(nfov,thint,th1,r1,rint1)
      integer nfov,mfov
      parameter(mfov=200)
      real thint(mfov),th1(mfov),r1(mfov),rint1(mfov)

      do i=1,nfov
       call verint(th1,r1,nfov,rint1(i),thint(i))
      enddo

      return

      end
