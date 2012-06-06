      subroutine readbass1(icheck)
C     ***************************************************************
C     Subroutine to read in ozone absorption coefficients supplied
C     by Don Grainger
C
C     Pat Irwin  	7/6/99
C
C     ***************************************************************

      implicit none
      integer i,j,iout,k,iread,icheck
      real bassk(5,1900),xdum,ydum
      character*60 buffer
      character*80 ipfile

      common /basstable/bassk,iread

      ipfile = 'bass_conv.dat'
      call datarchive(ipfile)

      open(12,file=ipfile,status='old')

      do 20 i=1,38
10     continue
       do 15 j=1,50
        iout = j+(i-1)*50
        read(12,*)(bassk(k,iout),k=1,5),xdum,ydum
15     continue
20    continue
      iread=-1
      close(12)

      icheck=1

      return

      end

