      program write_xsec
C     $Id: write_xsec.f,v 1.2 2011-06-17 15:53:03 irwin Exp $
C     ***********************************************************************
C
C     Simple program to write out a dust cross section file.
C
C     Pat Irwin		17/10/94
C
C     ***********************************************************************
      include '../includes/arrdef.f'
      include '../includes/pathcom.f'

      integer i,j
      character*100 ipfile

      call prompt('Enter dust x-section filename : ')
      read(5,1)ipfile
1     format(a60)

      call file(ipfile,opfile,'xsc')

      open(12,file=opfile,status='new')
      call prompt('Enter ncont : ')
      read*,ncont
      write(12,*)ncont

      call prompt('Enter number of spectral points : ')
      read*,nsec

      do 20 i=1,nsec
       print*,'Spectral point : ',i
       call prompt('Enter wavenumber : ')
       read*,vsec(i)

       do 10 j=1,ncont
        print*,'Aerosol type : ',j
        call prompt('Enter extinction xsection (cm2) : ')
        read*,x
        xsec(1,j,i)=x
        call prompt('Enter single scattering albedo : ')
        read*,x
        xsec(2,j,i)=x
10     continue

       write(6,*)vsec(i),(xsec(1,j,i),j=1,ncont)
       write(12,*)vsec(i),(xsec(1,j,i),j=1,ncont)
       write(6,*)vsec(i),(xsec(2,j,i),j=1,ncont)
       write(12,*)vsec(i),(xsec(2,j,i),j=1,ncont)

20    continue

      close(12)

      end
