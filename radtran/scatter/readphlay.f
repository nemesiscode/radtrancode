      subroutine readphlay(vwave,nlay,ff,gg1,gg2)
C     $Id: readphlay.f,v 1.2 2011-06-17 15:57:54 irwin Exp $
C     ***************************************************************

C     ***************************************************************
      include '../includes/arrdef.f'
      integer nphas,mphas
      parameter (mphas=60)
      real ff(MAXSCATLAY),gg1(MAXSCATLAY),gg2(MAXSCATLAY),vwave,nlay1
      real wphas(mphas),layphas(mphas,MAXSCATLAY,3)

      common /phlay/iread,nphas,wphas,layphas

      if(iread.ne.-1) then
       open(54,file='layphas.dat',status='old',
     &		form='unformatted' )
        read(54)nphas
        read(54)wphas
        read(54)layphas
       close(54)
      endif

      xl = 1e4/vwave

      k=-1
      do i=1,nphas-1
       if(xl.ge.wphas(i).and.xl.lt.wphas(i+1))k=i
      enddo
      if(xl.eq.wphas(nphas))k=nphas-1

      if(k.lt.0)then
       print*,'Error in readphlay'
       print*,'Wavelength not covered : ',xl
       stop
      endif
      frac = (xl-wphas(k))/(wphas(k+1)-wphas(k)) 

      do ilay=1,nlay
       ff(ilay)  = (1-frac)*layphas(k,ilay,1) + frac*layphas(k+1,ilay,1)
       gg1(ilay) = (1-frac)*layphas(k,ilay,2) + frac*layphas(k+1,ilay,2)
       gg2(ilay) = (1-frac)*layphas(k,ilay,3) + frac*layphas(k+1,ilay,3)
      enddo

      return

      end
