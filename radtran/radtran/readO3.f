      subroutine readO3(icheck)
C     ***************************************************************
C     Subroutine to read in ozone absorption coefficients supplied
C     by Don Grainger
C
C     Pat Irwin  	7/6/99
C
C     Updated to Serdyuchenko 2012 ozone cross sections, 213-1100 nm
C
C     Dan Dawson    9/4/13
C
C     ***************************************************************

      implicit none
      integer i,j,iout,k,iread,icheck
      real o3k(5,88668),xdum,ydum
      character*100 ipfile

      common /o3table/o3k,iread

      ipfile = 'xsection_o3_serd12.dat'
      call datarchive(ipfile)

      open(12,file=ipfile,status='old')

c      do 20 i=1,88668
c        read(12,*)(o3k(k,iout),k=1,5),xdum,ydum
c20    continue
c      iread=-1
c      close(12)

c   note that 88668 = 108 * 821
      
      do 20 i=1,108
10     continue
       do 15 j=1,821
        iout = j+(i-1)*821
        read(12,*)(o3k(k,iout),k=1,5),xdum,ydum
15     continue
20    continue
      iread=-1
      close(12)
       
      icheck=1

      return

      end

