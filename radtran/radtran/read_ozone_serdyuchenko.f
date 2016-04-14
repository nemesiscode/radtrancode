      subroutine read_ozone_serdyuchenko(icheck)
C     ***************************************************************
C     Subroutine to read in ozone absorption coefficients of
C     Serdyuchenko et al.
C
C     Pat Irwin  	15/4/16
C
C     ***************************************************************

      implicit none
      integer i,j,iout,k,iread,icheck,nozone
      parameter (nozone=88668)
      real wozone(nozone),kozone(nozone,11),tempozone(11)
      character*100 ipfile,buffer

      common /serdozonetable/wozone,tempozone,kozone

      ipfile = 'serdyuchenkogorshelev5digits.dat'

      call datarchive(ipfile)

      open(12,file=ipfile,status='old')

      do 10 i=1,45
       read(12,1)buffer
10    continue
1     format(a)

      do 20 i=1,nozone
        read(12,*)wozone(i),(kozone(i,k),k=11,1,-1)
20    continue

      close(12)

      do 30 i=1,11
       tempozone(i)=193+(i-1)*10.
30    continue
       
      icheck=1

      return

      end

