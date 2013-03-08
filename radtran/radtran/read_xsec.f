      subroutine read_xsec(xscfil)
C     $Id: read_xsec.f,v 1.5 2011-06-17 15:40:27 irwin Exp $
C     ***********************************************************************
C 
C     Modified routine for reading in dust cross-section files. Reads in
C     wavelengths as well so that cross-sections can be interpolated
C
C     Pat Irwin		18/10/94
C
C     ***********************************************************************
      include '../includes/arrdef.f'
      include '../includes/pathcom.f'
      integer i,j,ncont1, iunit
      real dummy,omega(10)
      character*100 xscfil,buffer

      iunit=19

      call file(xscfil,xscfil,'xsc')
      print*,'Reading_xsec: X-section file : ',xscfil

      open(iunit,file=xscfil,status='old')
1     format(a)
54    read(iunit,1)buffer
      if(buffer(1:1).eq.'#')goto 54
      read(buffer,*)ncont1
      if(ncont1.ne.ncont)then
       print*,'Error in read_xsec.'
       print*,'ncont in .xsc file is not consistent with driver file'
       print*,'ncont, ncont1 = ',ncont,ncont1
       stop
      end if
      if(ncont.gt.maxcon)then
       print*,'Error in read_xsec.'
       print*,'Too many particle types in xsc file'
       print*,'ncont, maxcon = ',ncont,maxcon
       stop
      end if

      j=0
105   j=j+1
      if(j.gt.maxsec)then
       print*,'Error in read_xsec: Too many wavelengths in xsc file'
       print*,'j, maxsec : ',j,maxsec
       stop
      endif
      read(iunit,*,end=106)vsec(j),(xsec(1,i,j),i=1,ncont)
      read(iunit,*,end=106)(omega(i),i=1,ncont)
      do 344 i=1,ncont
       xsec(2,i,j)=omega(i)*xsec(1,i,j)
344   continue
      goto 105
106   continue
      nsec=j-1

      print*,'read_xsec'
      print*,nsec,ncont
      do i=1,nsec
       print*,i,vsec(i),(xsec(1,j,i),j=1,ncont)
       print*,(xsec(2,j,i),j=1,ncont)
      enddo
      
      close(iunit)
      print*,'Ncont is', ncont


      return

      end
