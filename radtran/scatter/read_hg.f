      subroutine read_hg(vwave,icont,ncont,f,g1,g2)
      implicit none
      integer mcon,mphas
      parameter(mcon=10,mphas=300)
      real vwave,xwave(mphas),xf(mcon,mphas),xg1(mcon,mphas)
      real xg2(mcon,mphas)
      real twave,tico,x,y,z,frac,xw,tnco,xtest
      double precision f,g1,g2
      integer iunit,ico,i,ncont,icont,nco,j
      character*100 ipfile

      common /hgphas/xwave,xf,xg1,xg2,tnco,twave,frac,tico

      iunit=11
      if(xwave(1).lt.0.0)then

         if(ncont.gt.mcon)then
          print*,'Error in read_hg. Too many particle types'
          print*,ncont,mcon
          stop
         endif

         do 50 i=1,ncont
          ipfile='hgphase*.dat'
          ipfile(8:8)=char(i+48)
          print*,'Reading : ',ipfile
          ico = 0
          open(iunit,file=ipfile,status='old')
10        read(iunit,*,end=99)xw,x,y,z
          ico=ico+1
          if(ico.gt.mphas)then
           print*,'Error in read_hg: Too many wavelengths in'
           print*,'file : ',ipfile
           stop
          endif
          xwave(ico)=xw
          xf(i,ico) =x
          xg1(i,ico)=y
          xg2(i,ico)=z
          goto 10
99        continue
          close(iunit)
          nco=ico
          tnco = nco+0.1
50       continue

      else   
          nco = int(tnco)
      endif

      if(vwave.ne.twave)then
        ico=0
        do 60 i=1,nco-1
         if(xwave(i).le.vwave.and.xwave(i+1).gt.vwave)ico=i
60      continue
        xtest = abs(vwave - xwave(nco))
        if(xtest.lt.0.1)ico=nco-1
        if(ico.eq.0)then
         print*,'read_hg: wavenumber out of range : ',
     1			(xwave(j),j=1,nco),vwave
         print*,'icont,ncont : ',icont,ncont
         do i=1,nco
           print*,xwave(i),xf(icont,i),xg1(icont,i),xg2(icont,i)
         end do

         stop
        endif
        frac = (vwave - xwave(ico))/(xwave(ico+1)-xwave(ico))
        tico = ico + 0.1

        twave = vwave

      else

        ico=int(tico)

      endif
  
      f  = (1.0-frac)*xf(icont,ico)  + frac*xf(icont,ico+1)
      g1 = (1.0-frac)*xg1(icont,ico) + frac*xg1(icont,ico+1)
      g2 = (1.0-frac)*xg2(icont,ico) + frac*xg2(icont,ico+1)


      return
      end
         
     
