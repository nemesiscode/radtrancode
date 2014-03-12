      subroutine get_xsecA(xscfil,nmode,nwave,wave,xsec)
C     ***********************************************************************
C 
C     Modified routine for reading in dust cross-section files. 
C
C     Pat Irwin		18/10/94 	Original read_xsec,f
C     Pat Irwin		25/2/14		Revised
C
C     ***********************************************************************
      implicit none
      integer max_mode, max_wave,nmode,nwave,j
      parameter (max_mode = 10)
      parameter (max_wave = 1000)
      integer i,iunit
      real wave(max_wave),xsec(max_mode,max_wave,2)
      character*100 xscfil,buffer

      iunit=19

      call file(xscfil,xscfil,'xsc')
      print*,'Reading_xsec: X-section file : ',xscfil

      open(iunit,file=xscfil,status='old')
1     format(a)
54    read(iunit,1)buffer
      if(buffer(1:1).eq.'#')goto 54
      read(buffer,*)nmode
      if(nmode.gt.max_mode)then
       print*,'Error in get_xsecA: Too many particle types in xsc file'
       print*,'nmode, max_mode : ',nmode,max_mode
       stop
      endif

      j=0
105   j=j+1
      if(j.gt.max_wave)then
       print*,'Error in get_xsecA: Too many wavelengths in xsc file'
       print*,'j, max_wave : ',j,max_wave
       stop
      endif
      read(iunit,*,end=106)wave(j),(xsec(i,j,1),i=1,nmode)
      read(iunit,*,end=106)(xsec(i,j,2),i=1,nmode)
      goto 105
106   continue
      nwave=j-1

      close(iunit)

      return

      end
