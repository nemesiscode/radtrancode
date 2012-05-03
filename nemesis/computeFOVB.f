      subroutine computeFOVB(ichan,ipixA,ipixB,thcentre,
     1 thbore,thetrot,nfov,thfov,rfov)
      implicit none
      integer ichan,ipixA,ipixB,mfov,nfov
      integer j,k,i1,ipix
      parameter (mfov=200)
      real thfov(mfov),rfov(mfov),x,f,fvcen,thcentre,thbore,thetrot
      character*100 buffer

        ! Set array dimensions
        integer ndet,narr,iy
        parameter (ndet = 21, narr = 9, iy = 2**10)


        integer iv, ih, w, g, ny
        integer Lh(ndet,narr), Lv(ndet,narr),iangle
        real angle, vgrid(iy), vout(iy,ndet,narr), xout(mfov),yout(mfov)
        real VDet(ndet,narr), HDet(ndet,narr), tmpx(1024), tmpy(1024)
        real x1,x2
        real, allocatable :: Yfov(:,:,:), Vfov(:,:,:)
        real, allocatable :: Xfov(:,:,:), Hfov(:,:,:)
        integer i

      print*,'ComputeFOVB,ichan,ipixA,ipixB,thcentre,thbore'
      print*, ichan,ipixA,ipixB,thcentre,thbore
      print*,'computeFOVB thetrot  :',thetrot


        ! Call the first step of reading the FOV calibration data, so as to get
        !   the dimensions iv, & ih
        call FOV_ReadIn_Sizes(ndet,narr,iv,ih)
C        print*, '         iv = ', iv
C        print*, '         ih = ', ih

        ! Now that we have iv & ih, allocate array space for the FOV variables
        allocate( Yfov(iv,ndet,narr), Vfov(iv,ndet,narr) )
        allocate( Xfov(ih,ndet,narr), Hfov(ih,ndet,narr) )

        ! Now call the second step of reading the FOV calibration data
        call FOV_ReadIn(ndet,narr,iv,ih,Yfov,Xfov,Vfov,Hfov,
     1    VDet,HDet,Lv,Lh)

        angle=thetrot
        ! Call FOV_Rotation routine
         call FOV_Rotation( Yfov, Xfov, Vfov, Hfov, VDet, HDet, Lv, 
     1    Lh,iv, ih, ndet, narr, angle, iy, vgrid, vout )


         if(ichan.le.6) then
          ipix = ipixA
         else
          ipix = ipixB
         endif

         ny=0
         do i1=1,iy
              if(vout(i1,ipix,ichan).gt.0)ny=ny+1
         enddo
         k=0
         do i1=1,iy
              if(vout(i1,ipix,ichan).gt.0)then
                k=k+1
                tmpx(k)=vgrid(i1)
                tmpy(k)=vout(i1,ipix,ichan)
              endif
         enddo
         x1 = tmpx(1)
         x2 = tmpx(ny)

         nfov=mfov

         do i1=1,nfov
             xout(i1) = x1 + (i1-1)*(x2-x1)/(1.0*nfov-1.0)
             call verint(tmpx,tmpy,ny,yout(i1),xout(i1))
             thfov(i1)=thbore-xout(i1)
             rfov(i1)=yout(i1)
C             print*,thfov(i1),xout(i1),rfov(i1)
         enddo

         fvcen = thbore-Vdet(ipix,ichan)
         print*,'thbore, fvcen : ',thbore,fvcen

      return

      end
      
