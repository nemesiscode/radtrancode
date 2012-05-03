      subroutine impfluxsol(LUNIS,ioff,nlays,nmu,nf,umift,uplft,
     1  x,nwavef,vwavef,Ig,ng)

      implicit none
      INCLUDE '../includes/arrdef.f'

      integer nlays,nmu,iwave,Ig,ng,joff,i,j,k,LUNIS,IOFF
      integer nwavef,irec,nf,i1,ksol
      real umift(maxmu,maxmu,maxscatlay,maxf)
      real uplft(maxmu,maxmu,maxscatlay,maxf)
      real vwavef(maxbin),PREC,dx,dxmin,x
      parameter(PREC = 1e-6)
      iwave=0
C      print*,'impflux, nwavef = ',nwavef
      dxmin = 1e10
      do i=1,nwavef 
       dx = abs(x-vwavef(i))
       if(dx.lt.dxmin)then
        iwave=i
        dxmin = dx
       endif
      enddo
C      print*,'impflux,x,iwave',x,iwave,vwavef(iwave)
      if(dxmin.gt.PREC)then
       print*,'Error in impflux. Match could not be found'
       print*,'within the precision set'
       print*,'wavelength : ',x
       print*,'Wavelengths in table are : '
       print*,nwavef
       do i=1,nwavef
        print*,i,vwavef(i)
       enddo
       print*,'Closest match at ',vwavef(iwave)
       print*,'Precision set = ',PREC
       stop
      endif
  
C      print*,ioff,iwave,ng,ig,nlays,nmu

      irec = ioff+2*((iwave-1)*ng + (ig-1))*nlays*nmu*nmu*(nf+1)

      do j=1,nlays
       do ksol=1,nmu
        do k=1,nmu
         do i1=1,nf+1
          read(LUNIS,REC=IREC)umift(ksol,k,j,i1)
          irec=irec+1
          read(LUNIS,REC=IREC)uplft(ksol,k,j,i1)
          irec=irec+1
         enddo
        enddo
       enddo
      enddo

      if(ig.eq.5)then
        open(48,file='testA.dat',status='unknown')
        write(48,*)nlays
        write(48,*)nf
        do j=1,nlays
         do i1 = 1,nf+1
          do ksol=1,nmu
           write(48,*)(umift(ksol,k,j,i1),k=1,nmu)
          enddo
         enddo
        enddo
        do j=1,nlays
         do i1 = 1,nf+1
          do ksol=1,nmu
           write(48,*)(uplft(ksol,k,j,i1),k=1,nmu)
          enddo
         enddo
        enddo
        close(48)
      endif

      return

      end
