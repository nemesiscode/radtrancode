      program normxsc
C     $Id: normxsc.f,v 1.2 2011-06-17 15:57:53 irwin Exp $
C
***********************************************************************
C 
C     Routine for normalising aerosol x-section spectra to OD at a 
C     specific wavenumber
C
C     Pat Irwin         28/8/00
C
C
***********************************************************************
      integer j,k,ncont, iunit, nsec
      real vsec(1000),xsec(2,1000,10),ods(10)
      character*100 xscfil,buffer

      call prompt('Enter filename : ')
      read(5,1)xscfil
1     format(a)

      call file(xscfil,xscfil,'xsc')
      iunit=19

      open(iunit,file=xscfil,status='old')
54    read(iunit,1)buffer
      if(buffer(1:1).eq.'#')goto 54
      read(buffer,*)ncont

      write(*,*)'Wavenumber ordinates are : '
      j=0
105   j=j+1
      read(iunit,*,end=106)vsec(j),(xsec(1,j,k),k=1,ncont)
      read(iunit,*,end=106)(xsec(2,j,k),k=1,ncont)
      write(*,*)j,vsec(j)
      goto 105
106   continue
      nsec=j-1

      close(iunit)

      call prompt('Enter wavelength ordinate and desired OD : ')
      read*,iwave,xod

      do k=1,ncont
       ods(k)=xsec(1,iwave,k)
      enddo

      do j=1,nsec
       do k=1,ncont
        xsec(1,j,k) = xod*xsec(1,j,k)/ods(k)
       enddo
      enddo


      open(iunit,file=xscfil,status='unknown')
      write(iunit,*)ncont

      do j=1,nsec
       write(iunit,*)vsec(j),(xsec(1,j,k),k=1,ncont)
       write(iunit,*)(xsec(2,j,k),k=1,ncont)
      enddo

      close(iunit)

      end
