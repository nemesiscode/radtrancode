! This is an example driver program that calls subroutines to read in the 
!   FOV calibration data and then rotate the data by some angle.

	program outputfov

  	implicit none

  	! Set array dimensions
  	integer ndet,narr,iy
  	parameter (ndet = 21, narr = 9, iy = 2**10)

  	integer k, iv, ih, w, g, ny, mfov,jc
  	integer Lh(ndet,narr), Lv(ndet,narr),iangle,i1,j1
  	real angle, vgrid(iy), vout(iy,ndet,narr), xout(200),yout(200)
  	real VDet(ndet,narr), HDet(ndet,narr), tmpx(1024), tmpy(1024)
        real x1,x2,dx(1024),sum,xfac,XMIN,xgrid(1024),xcen
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

  	open(unit=100,file='fov.dat',status='unknown')

        do iangle = 1,11
        
         angle = -1.0+0.2*(iangle-1)        
         print*,'iangle ',iangle
         write(100,*)iangle,angle

  	! Call FOV_Rotation routine
  	 call FOV_Rotation( Yfov, Xfov, Vfov, Hfov, VDet, HDet, Lv, Lh,iv, ih, ndet, narr, angle, iy, vgrid, vout )

         mfov = 50
         XMIN = 1e-5
         do j=1,narr
          print*,'jarr = ',j
 	  do i=1,ndet
            print*,'idet = ',i
            ny=0
            do i1=1,iy
     	      if(vout(i1,i,j).gt.XMIN)ny=ny+1
            enddo
            write(100,*)HDet(i,j),VDet(i,j),j,i,mfov

            k=0
            do i1=1,iy
              if(vout(i1,i,j).gt.XMIN)then
                k=k+1
                tmpx(k)=vgrid(i1)
                tmpy(k)=vout(i1,i,j)
              endif
            enddo

            sum=-200.0
            do i1=1,ny
             dx(i1)=1-alog(tmpy(i1))
             if(tmpy(i1).gt.sum)then
               sum=tmpy(i1)
               j1=i1
             endif
            enddo

            x1 = tmpx(1)
            x2 = tmpx(ny)
            xcen = tmpx(j1)
            print*,x1,xcen,x2,j1,ny
            xfac=1.0

            jc=0
202         xgrid(j1)=xcen
            do i1=j1+1,ny
             xgrid(i1)=xgrid(i1-1)+dx(i1-1)*xfac
            enddo
            do i1=j1-1,1,-1
             xgrid(i1)=xgrid(i1+1)-dx(i1+1)*xfac
            enddo
 
            k=0
            do i1=1,ny
             if(xgrid(i1).ge.x1.and.xgrid(i1).le.x2)k=k+1
            enddo
          
            if(k.ne.mfov)then
             xfac = xfac*k/(1.0*mfov)
           !  print*,k,xfac
              jc=jc+1
             goto 202 
            else
             print*,'niter = ',jc
             k=1
             do i1=1,ny
              if(xgrid(i1).ge.x1.and.xgrid(i1).le.x2)then
               xout(k)=xgrid(i1)
               k=k+1
              endif
             enddo
            endif

            do i1=1,mfov
             call verint(tmpx,tmpy,ny,yout(i1),xout(i1))
       	     write(100,'(2e30.15)') xout(i1),yout(i1)
            enddo

          !  do i1=1,iy
          !    if(vout(i1,i,j).gt.0)then
       	  !      write(100,'(2e30.15)') vgrid(i1),vout(i1,i,j)
          !    endif
          !  enddo
          enddo
         enddo

        enddo

  	close(100)

	end


