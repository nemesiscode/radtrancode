      subroutine dumpflux(LUNIS,ioff,nlays,nmu,nf,radg,umif,uplf,iwave,
     1  x,nwave,Ig,ng,j0,xmu,solar,taus)

      implicit none
      INCLUDE '../includes/arrdef.f'
      integer LUNIS,ioff,nwave,irec,j0
      integer nlays,nmu,iwave,Ig,ng,joff,i,j,k,i1,nf
      real x,umif(maxmu,maxscatlay,maxf),uplf(maxmu,maxscatlay,maxf)
      real radg(maxmu),tmpzero,taus(maxscatlay),xmu,solar,tmp,sol1

      irec = ioff+2*((iwave-1)*ng + (ig-1))*nlays*nmu*(nf+1)

C     Add in bottom radiance to complete radiance field at base of
C     atmosphere. NB. Scattering flux defines uplf to be radiance going 
C     downwards out of bottom of layer and umif to be the radiance going
C     upwards out of the top. We change the meaning here to be the 
C     total radiance field (up and down) at the bottom of each layer.  
C
C     In addition, the order of layers of the umif,uplf arrays is from 
C     top to bottom. We here reverse the order to be from bottom to top
C     so that the alignment with the baseh array is correct.

C     Also deduct direct solar part here as it screws up the fourier
C     decomposition

C     First strip out direct solar component from uplf
      tmp=0.0
C      print*,xmu,j0
      do j=1,nlays
       tmp = tmp+taus(j)/xmu
       sol1 = solar*exp(-tmp)
C       if(j.eq.15) then 
C         print*,j,solar,taus(j),tmp,sol1
C         do i1=1,nf-1 
C          print*,uplf(j0,j,i1),uplf(j0,j,i1)-sol1
C         enddo
C         stop
C       endif


       do i1=1,nf+1
        uplf(j0,j,i1)=uplf(j0,j,i1)-sol1
       enddo


      enddo

      tmpzero=0.0
C     Reverse order of layers here. Scattering code has 1 as the top
C     layer. Rest of code have 1 as the bottom layer.
 
      do k=1,nmu
        do i1=1,nf+1
          if(i1.eq.1) then
            write(LUNIS,REC=IREC)radg(k)
          else
            write(LUNIS,REC=IREC)tmpzero
          endif
          irec=irec+1
          write(LUNIS,REC=IREC)uplf(k,nlays,i1)
          irec=irec+1
        enddo
      enddo
      
C      print*,'umif'
      do j=nlays-1,1,-1
       do k=1,nmu
        do i1=1,nf+1
         write(LUNIS,REC=IREC)umif(k,j+1,i1)
         irec=irec+1
         write(LUNIS,REC=IREC)uplf(k,j,i1)
         irec=irec+1
        enddo
       enddo
C       print*,(umif(k,j+1,1),k=1,nmu)
      enddo


C      do j=nlays,1,-1
C       write(LUNIS,REC=IREC)taus(j)
C      enddo


      return

      end
