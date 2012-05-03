      subroutine read_xsc(xscfil,ncont,nsec,vsec,xsec)
C     $Id: read_xsc.f,v 1.2 2011-06-17 15:57:54 irwin Exp $
C     ***********************************************************************
C 
C     Modified routine for reading in dust cross-section files. Reads in
C     wavelengths as well so that cross-sections can be interpolated
C
C     Pat Irwin		18/10/94
C
C     ***********************************************************************
      integer i,j,ncont1, iunit, nsec
      include '../includes/arrdef.f'

      real dummy,omega(maxcon),vsec(maxsec),xsec(2,maxsec)
      character*100 xscfil,buffer

      call file(xscfil,xscfil,'xsc')
      iunit=19

      open(iunit,file=xscfil,status='old')
1     format(a)
54    read(iunit,1)buffer
      if(buffer(1:1).eq.'#')goto 54
      read(buffer,*)ncont1
      if(ncont1.ne.ncont)then
       print*,'Error in read_xsc.'
       print*,'ncont in .xsc file is not consistent'
       print*,'ncont, ncont1 = ',ncont,ncont1
       stop
      end if

      ncont=ncont1

      j=0
105   j=j+1
      read(iunit,*,end=106)vsec(j),xsec(1,j)
      read(iunit,*,end=106)xsec(2,j)
      goto 105
106   continue
      nsec=j-1

      close(iunit)

      return

      end
