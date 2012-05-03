      subroutine combineflux(nsol,tsurf,nem,vem,emissivity,ispace)
C     ******************************************************************
C     Subroutine for combining internal radiation fields calculated
C     for different solar zenith angles into a single file.
C     ******************************************************************

      implicit none
      INCLUDE '../includes/arrdef.f'
      integer nsol,nlays,nmu,iwave,Ig,ng,joff,i,j,k,LUNIS,IOFF
      integer nwavef,irec,nf,i1,i2,ipop,nft,LUNOUT
      integer nlayer,isol,nwave
      real x,umif(maxmu,maxscatlay,maxf),uplf(maxmu,maxscatlay,maxf)
      real umift(maxmu,maxmu,maxscatlay,maxf),radg(maxmu)
      real uplft(maxmu,maxmu,maxscatlay,maxf)
      real vwavef(maxbin),dx,dxmin,tsurf,esurf,xrad
      real basep(maxscatlay),baseh(maxscatlay),planck_wave
      integer tabnf(maxf),imu,imu0,ilayer,nem,ispace
      real vem(maxsec),emissivity(maxsec),interpem
      character*100 fintname,tmpname,outname
      ipop=0
      nft=0      
      LUNIS=60
      LUNOUT=100
      do 10 isol=1,nsol
C       print*,'combineflux: isol = ',isol
       fintname = 'internal**.fld'
       i1=int(isol/10)
       i2=isol-i1*10
       fintname(9:9)=char(i1+48)
       fintname(10:10)=char(i2+48)
        
       call openflux(fintname,ipop,LUNIS+isol,nlayer,nmu,
     1  nwave,ng,nf,ioff,vwavef,basep,baseh)
  
       tabnf(isol)=nf
       if(nf.gt.nft)nft=nf

10    continue
      ipop=1
      fintname='internalsol.fld'

      call openflux(fintname,ipop,LUNOUT,nlayer,nmu,
     1  nwave,ng,nft,ioff,vwavef,basep,baseh)

      tmpname = 'individ******.unf'
      outname = 'combine****.unf'

      do 200 iwave=1,nwave
       x=vwavef(iwave)
       esurf = interpem(nem,vem,emissivity,x)
       xrad = esurf*planck_wave(ispace,x,tsurf)

       i1 = int(iwave/10)
       i2 = iwave-i1*10
       tmpname(8:8)=char(i1+48)
       tmpname(9:9)=char(i2+48)
       outname(8:8)=char(i1+48)
       outname(9:9)=char(i2+48)

       do 300 ig=1,ng

        i1 = int(ig/10)
        i2 = ig-i1*10
        tmpname(10:10)=char(i1+48)
        tmpname(11:11)=char(i2+48)
        outname(10:10)=char(i1+48)
        outname(11:11)=char(i2+48)

        do imu0=1,nmu
         radg(imu0)=xrad
         do imu=1,nmu
          do ilayer=1,nlayer
           do i1=1,nft
            umift(imu0,imu,ilayer,i1)=0.0
            uplft(imu0,imu,ilayer,i1)=0.0
           enddo     
          enddo
         enddo
        enddo


        do 400 isol=1,nsol

         i1=int(isol/10)
         i2=isol-i1*10
         tmpname(12:12)=char(i1+48)
         tmpname(13:13)=char(i2+48)

         nf=tabnf(isol)
C         print*,'isol, nf = ',isol,nf
         call impflux(LUNIS+isol,ioff,nlayer,nmu,nf,umif,uplf,
     1  x,nwave,vwavef,Ig,ng)

         if(Ig.eq.5.and.iwave.eq.1)then
          open(48,file=tmpname,status='unknown',form='unformatted')
           write(48)umif
           write(48)uplf
          close(48)
         endif

         do imu=1,nmu
          do ilayer=1,nlayer
           do i1=1,nf
            umift(isol,imu,ilayer,i1)= umif(imu,ilayer,i1)
            uplft(isol,imu,ilayer,i1)= uplf(imu,ilayer,i1)
           enddo     
          enddo
         enddo

400     continue

        if(Ig.eq.5.and.iwave.eq.1)then
         open(48,file=outname,status='unknown',form='unformatted')
          write(48)umift
          write(48)uplft
         close(48)
        endif

        call dumpfluxsol(LUNOUT,ioff,nlayer,nmu,nft,radg,umift,uplft,
     1  iwave,x,nwave,Ig,ng)

300    continue
200   continue

      close(lunout)
      do isol=1,nsol
       close(lunis+isol)
      enddo

      return

      end
