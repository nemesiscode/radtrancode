      subroutine interpolateMCSnog(nview,npath,thview,
     1 nconv,iconv,vconv,nlayer,nx,esurf,jsurf,ispace,tsurf,
     2 calcoutL,nfov,thfov,rfov,radmean)

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      integer mview,mfov,npath
      parameter (mview=100,mfov=200)
      integer nview,nconv,iconv,nlayer,nx,jsurf,nfov
      integer ipath,i,j,ioff1,ioff2,ipathA
      integer ioff1A,ioff2A,ifov,ispace
      real sum1,sum2,radview(mview),radnow,grad(mview)
      real tsurf,planck_wave,planckg_wave
      real hview(mview),thview(mview),tv1(mview)
      real esurf,calcoutL(maxout3)
      real vconv(mconv)
      real thfov(mfov),rfov(mfov),radmean,gradmean(mx)
C     First set up radiance array
C
C     First limb paths
      do i=1,nlayer
       ipath=i
       ioff1=ipath+(iconv-1)*npath
       radview(1+nlayer-i)=calcoutL(ioff1)
       tv1(1+nlayer-i)=thview(i)
      enddo

C     Then ground intersecting paths
      do i = nlayer+1,nview
       ipath=nlayer+1 + (i-nlayer-1)*2
       ipathA = ipath+1
       ioff1=ipath + (iconv-1)*npath
       ioff1A=ipathA + (iconv-1)*npath
       radview(i)=calcoutL(ioff1) + calcoutL(ioff1A)*
     1   planck_wave(ispace,vconv(iconv),tsurf)*esurf
       tv1(i)=thview(i)
      enddo
 

C     Now convolve with FOV arrays

      sum1=0.0
      sum2=0.0
      do ifov=1,nfov
       call verint(tv1,radview,nview,radnow,thfov(ifov))
       sum1 = sum1 + rfov(ifov)*radnow
       sum2 = sum2 + rfov(ifov)
      enddo

      radmean = sum1/sum2


      return

      end

