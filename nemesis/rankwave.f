      subroutine rankwave(ngeom,nwave,vwave,nconv,vconv,nwaveT,vwaveT,
     1 nconvT,vconvT)
C     ************************ VARIABLES *******************************
      implicit none
C     Set measurement vector and source vector lengths here.
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      integer ngeom, nwave(mgeom), nconv(mgeom)
      real vwave(mgeom,mwave),vconv(mgeom,mconv)
      real vwaveT(mwave),vconvT(mconv)
      integer nwaveT,nconvT,igeom,iconv,iwave,icheck,i

      nwaveT=0
      nconvT=0

      do igeom=1,ngeom

       do iwave=1,nwave(igeom)
        icheck=0
        do i=1,nwaveT
         if(vwave(igeom,iwave).eq.vwaveT(i))icheck=1
        enddo
        if(icheck.eq.0)then
         nwaveT=nwaveT+1
         vwaveT(nwaveT)=vwave(igeom,iwave)
        endif
       enddo

       do iconv=1,nconv(igeom)
        icheck=0
        do i=1,nconvT
         if(vconv(igeom,iconv).eq.vconvT(i))icheck=1
        enddo
        if(icheck.eq.0)then
         nconvT=nconvT+1
         vconvT(nconvT)=vconv(igeom,iconv)
        endif        
       enddo

      enddo


      if(nwaveT.gt.1)call sort(nwaveT,vwaveT)
      if(nconvT.gt.1)call sort(nconvT,vconvT)

      return

      end 
