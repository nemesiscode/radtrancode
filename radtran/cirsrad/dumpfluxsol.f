      subroutine dumpfluxsol(LUNIS,ioff,nlays,nmu,nf,radg,
     1  umift,uplft,iwave,x,nwave,Ig,ng)

      implicit none
      INCLUDE '../includes/arrdef.f'
      integer LUNIS,ioff,nwave,irec,ksol
      integer nlays,nmu,iwave,Ig,ng,joff,i,j,k,i1,nf
      real umift(maxmu,maxmu,maxscatlay,maxf)
      real uplft(maxmu,maxmu,maxscatlay,maxf)
      real x, radg(maxmu),tmpzero

      irec = ioff+2*((iwave-1)*ng + (ig-1))*nlays*nmu*nmu*(nf+1)

      do j=1,nlays
       do ksol=1,nmu
        do k=1,nmu
         do i1=1,nf+1
          write(LUNIS,REC=IREC)umift(ksol,k,j,i1)
          irec=irec+1
          write(LUNIS,REC=IREC)uplft(ksol,k,j,i1)
          irec=irec+1
         enddo
         enddo
        enddo
      enddo


      return

      end
