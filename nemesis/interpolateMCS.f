      subroutine interpolateMCS(nview,thview,
     1 nconv,iconv,vconv,nlayer,nx,esurf,jsurf,ispace,tsurf,
     2 calcoutL,gradientsL,nfov,thfov,rfov,radmean,gradmean)

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      integer mview,mfov
      parameter (mview=100,mfov=200)
      integer nview,nconv,iconv,nlayer,nx,jsurf,nfov
      integer ipath,i,j,ioff1,ioff2,ipathA
      integer ioff1A,ioff2A,ifov,ispace
      real sum1,sum2,radview(mview),radnow,grad(mview)
      real gradview(mview,mx),tsurf,planck_wave,planckg_wave
      real hview(mview),thview(mview),tv1(mview)
      real esurf,calcoutL(maxout3)
      real gradientsL(maxout4),vconv(mconv),vv
      real thfov(mfov),rfov(mfov),radmean,gradmean(mx)
C     First set up radiance array
C
      vv = vconv(iconv)
C     First limb paths
      do i=1,nlayer
       ipath=i
       ioff1=nconv*(ipath-1)+iconv
       radview(1+nlayer-i)=calcoutL(ioff1)
       tv1(1+nlayer-i)=thview(i)
       do j=1,nx
        ioff2 = nconv*nx*(ipath-1)+(j-1)*nconv + iconv
        gradview(1+nlayer-i,j)=gradientsL(ioff2)
       enddo
      enddo

C     Then ground intersecting paths
      do i = nlayer+1,nview
       ipath=nlayer+1 + (i-nlayer-1)*2
       ipathA = ipath+1
       ioff1=nconv*(ipath-1)+iconv
       ioff1A=nconv*(ipathA-1)+iconv

       radview(i)=calcoutL(ioff1) + calcoutL(ioff1A)*
     1   planck_wave(ispace,vv,tsurf)*esurf
       tv1(i)=thview(i)

C     1 planck_wave(ispace,vv,tsurf),tsurf,esurf
       do j=1,nx
        ioff2 = nconv*nx*(ipath-1)+(j-1)*nconv + iconv
        ioff2A = nconv*nx*(ipathA-1)+(j-1)*nconv + iconv
        if(j.eq.jsurf) then
         gradview(i,j)=calcoutL(ioff1A)*
     1    planckg_wave(ispace,vv,tsurf)*esurf
        else
         gradview(i,j)=gradientsL(ioff2) + gradientsL(ioff2A)*
     1          planck_wave(ispace,vv,tsurf)*esurf
        endif
       enddo
      enddo

C      open(12,file = 'interp.dat',status='unknown') 
C      write(12,*)'interpolateMCS: Radiance profile to interpolate'
C      write(12,*)nview,nx
C      do i=1,nview
C       write(12,*)tv1(i),radview(i)
C      enddo     

C      do i=1,nview
C       write(12,*)tv1(i),(gradview(i,j),j=1,nx)
C      enddo     

C      write(12,*)nfov
C      do i=1,nfov
C       write(12,*)thfov(i),rfov(i)
C      enddo

C      close(12)

C     Now convolve with FOV array

      sum1=0.0
      sum2=0.0
      do ifov=1,nfov
       call verint(tv1,radview,nview,radnow,thfov(ifov))
       sum1 = sum1 + rfov(ifov)*radnow
       sum2 = sum2 + rfov(ifov)
      enddo

      radmean = sum1/sum2


      do i=1,nx
       do j=1,nview
        grad(j) = gradview(j,i)
       enddo

       sum1=0.0
       sum2=0.0
       do ifov=1,nfov
        call verint(tv1,grad,nview,radnow,thfov(ifov))
        sum1 = sum1 + rfov(ifov)*radnow
        sum2 = sum2 + rfov(ifov)
       enddo

       gradmean(i)=sum1/sum2

      enddo

      return

      end

