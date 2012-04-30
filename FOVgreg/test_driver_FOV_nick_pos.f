! This is an example driver program that calls subroutines to read in the 
!   FOV calibration data and then rotate the data by some angle.

	program test_driver_FOV_nick_pos

  	implicit none

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

  	! As a sanity check, print some of the variables to make sure there has
  	!  been no obvious I/O problems
  	print*, 'Yfov(1,1,1) = ', Yfov(1,1,1)
  	print*, 'Vfov(1,1,1) = ', Vfov(1,1,1)
  	print*, 'Xfov(1,1,1) = ', Xfov(1,1,1)
  	print*, 'Hfov(1,1,1) = ', Hfov(1,1,1)
  	print*, '  VDet(1,1) = ', VDet(1,1)
  	print*, '  HDet(1,1) = ', HDet(1,1)
  	print*, '    Lv(1,1) = ', Lv(1,1)
  	print*, '    Lh(1,1) = ', Lh(1,1)

  	! Set desired angle
  	angle = 5.0

  	! Call FOV_Rotation routine
  	call FOV_Rotation( Yfov, Xfov, Vfov, Hfov, VDet, HDet, Lv, Lh,iv, ih, ndet, narr, angle, iy, vgrid, vout )

  	open(unit=100,file='det.dat',status='unknown')
      do j=1,narr
	  do i=1,ndet
          write(100,*) HDet(i,j),VDet(i,j),j,i
        enddo
      enddo
  	close(100)

  	! Write out the rotated final effective vertical FOV
  	w = 1;  g = 1
  	open(unit=100,file='vf.dat',status='unknown')
  	do i=1,iy
  	  write(100,'(2e30.15)') vgrid(i),vout(i,w,g)
  	end do
  	close(100)

	end


