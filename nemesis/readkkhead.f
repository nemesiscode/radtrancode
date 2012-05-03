      subroutine readkkhead(runname,vkstart,vkend,vkstep)
C     ****************************************************************
C     Subroutine to read in the header of the first .kta file to determine
C     the range and step of the k-tables.
C
C     Pat Irwin		10/2/04
C  
C     ****************************************************************
      implicit none
      character*100 runname
      character*100 kname
      integer lun0,npoint,idgas,isogas,np,nt,ng,irec0,mpoint
      parameter(mpoint=12000)
      real vcen(mpoint)
      real vmin,delv,fwhm,press(20),temp(20),g_ord(21),del_g(21)
      real vkstart,vkend,vkstep

      call file(runname,runname,'kls')
1     format(a100)

C      print*,'readkkhead: reading kta information from : '
C      write(6,1)runname

      open(12,file=runname,status='old')
        read(12,1)kname
      close(12)


C      print*,'Reading header information from : '
C      write(6,1)kname

      lun0 = 100

      call read_khead(kname,lun0,npoint,vmin,delv,fwhm,vcen,
     1  idgas,isogas,press,temp,np,nt,g_ord,del_g,ng,irec0)
      close(100)


      vkstep = delv
      vkstart = vmin
      vkend = vkstart + delv*(npoint-1)

      return
 
      end
