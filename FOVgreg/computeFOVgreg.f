subroutine computeFOVgreg(ichan,ipixA,ipixB,thcentre,angle,nfov,thfov,rfov)
  implicit none
  integer, intent(in) :: ichan,ipixA,ipixB
  integer, intent(out) :: nfov

  integer mfov
  parameter (mfov=1024)
      real thfov(mfov),rfov(mfov)
      real thcentre,thetrot,fovcentre

        ! Set array dimensions
        integer ndet,narr,iy
        parameter (ndet = 21, narr = 9, iy = 2**10)

        integer k, iv, ih, w, g
        integer Lh(ndet,narr), Lv(ndet,narr)
        real angle, vgrid(iy), vout(iy,ndet,narr)
        real VDet(ndet,narr), HDet(ndet,narr)
        real, allocatable :: Yfov(:,:,:), Vfov(:,:,:)
        real, allocatable :: Xfov(:,:,:), Hfov(:,:,:)
        integer i,j
        ! Call the first step of reading the FOV calibration data, so as to get
        !   the dimensions iv, & ih
        call FOV_ReadIn_Sizes(ndet,narr,iv,ih)
        print*, '         iv = ', iv
        print*, '         ih = ', ih

        ! Now that we have iv & ih, allocate array space for the FOV variables
        allocate( Yfov(iv,ndet,narr), Vfov(iv,ndet,narr) )
        allocate( Xfov(ih,ndet,narr), Hfov(ih,ndet,narr) )

        ! Now call the second step of reading the FOV calibration data
        call FOV_ReadIn(ndet,narr,iv,ih,Yfov,Xfov,Vfov,Hfov,VDet,HDet,Lv,Lh)
      
        ! Call FOV_Rotation routine
        call FOV_Rotation( Yfov, Xfov, Vfov, Hfov, VDet, HDet, Lv, Lh,iv, ih, ndet, narr, angle, iy, vgrid, vout )


        ! Find centre of required pixel wrt to the boresight
        fovcentre = vdet(ichan,ipixA)
       
        nfov = iy
   
        print*,'ichan,ipixA = ',ichan,ipixA

        do i=1,iy
         thfov(i)=vgrid(i) - fovcentre + thcentre
         rfov(i)=vout(i,ichan,ipixA)
        enddo

end subroutine computeFOVgreg
      
