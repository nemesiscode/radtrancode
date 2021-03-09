      subroutine readkkhead(runname,vkstart,vkend,vkstep)
C     ****************************************************************
C     Subroutine to read in the header of the first .kta file to determine
C     the range and step of the k-tables.
C
C     Pat Irwin		10/2/04
C  
C     ****************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      character*100 runname
      character*200 kname
      integer lun0,npoint,idgas,isogas,np,nt,ng,irec0
      real vcen(MAXBIN)
      real vmin,delv,fwhm,press(MAXK),temp(MAXK),g_ord(MAXG),del_g(MAXG)
      real vkstart,vkend,vkstep,temp2(MAXK,MAXK)
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      call file(runname,runname,'kls')
1     format(a200)


      open(12,file=runname,status='old')
        read(12,1)kname
      close(12)


C      if(idiag.gt.0)print*,'Reading header information from : '
C      if(idiag.gt.0)write(6,1)kname

      lun0 = 100

      call read_khead(kname,lun0,npoint,vmin,delv,fwhm,vcen,
     1  idgas,isogas,press,temp,temp2,np,nt,g_ord,del_g,ng,
     2  irec0)
      close(100)


      vkstep = delv
      vkstart = vmin
      vkend = vkstart + delv*(npoint-1)

      return
 
      end
