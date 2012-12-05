      program merge
C     ******************************************************************
C     Program to merge together two sequential access line data files
C
C     Pat Irwin		8/9/95
C     ******************************************************************
      implicit none
      integer iold,iextra,inew,irecl,ichoice
      parameter (iold=12,iextra=13,inew=14)           
      character*256 buffer1,buffer2
      character*100 ipfile1,ipfile2,opfile
      real vv1,vv2

2     format(a)
      call prompt('Enter record length : ')
      read*,irecl
      call prompt('Enter first line data file : ')
      read(5,2)ipfile1
      call prompt('Enter second line data file : ')
      read(5,2)ipfile2
      call prompt('Enter output line data file : ')
      read(5,2)opfile

      open(iold,file=ipfile1,status='old')
      open(iextra,file=ipfile2,status='old')
      open(inew,file=opfile,status='unknown')

C     skip headers

23    read(unit=iold,fmt=80,end=999)buffer1(1:irecl)
      print*,buffer1(1:80)
      if(buffer1(2:2).eq.'#')then
       goto 23
      end if

24    read(unit=iextra,fmt=80,end=998)buffer2(1:irecl)
      print*,buffer2(1:80)
      if(buffer2(2:2).eq.'#')then
       goto 24
      end if


C     begin merge

20    continue
      if(irecl.eq.80.or.irecl.eq.82)then
       read(buffer1,10)vv1
       read(buffer2,10)vv2
10     format(f10.3)
      else
       read(buffer1(4:15),11)vv1
       read(buffer2(4:15),11)vv2
11     format(f12.6)
      endif      

80    format(a)

      if(vv1.lt.vv2)then
       write(unit=inew,fmt=80)buffer1(1:irecl)
       read(unit=iold,fmt=80,end=999)buffer1(1:irecl)
      else if(vv1.eq.vv2)then

C       print*,'Wavelengths are the same. Enter which line to use'
C       print*,'1'
C       write(6,1)buffer1
C       print*,'2'
C       write(6,1)buffer2
C       print*,'3 = both '
C       read*,ichoice

       ichoice=3
       if(buffer1(1:irecl).eq.buffer2(1:irecl))ichoice=1
       if(ichoice.eq.1)then
        write(unit=inew,fmt=80)buffer1(1:irecl)
        read(unit=iold,fmt=80,end=999)buffer1(1:irecl)
       else if(ichoice.eq.2)then
        write(unit=inew,fmt=80)buffer2(1:irecl)
        read(unit=iextra,fmt=80,end=998)buffer2(1:irecl)
       else
        write(unit=inew,fmt=80)buffer1(1:irecl)
        read(unit=iold,fmt=80,end=999)buffer1(1:irecl)
        write(unit=inew,fmt=80)buffer2(1:irecl)
        read(unit=iextra,fmt=80,end=998)buffer2(1:irecl)
       end if
      else
       write(unit=inew,fmt=80)buffer2(1:irecl)
       read(unit=iextra,fmt=80,end=998)buffer2(1:irecl)
      end if
      goto 20

998   close(iextra)
887   write(unit=inew,fmt=80)buffer1(1:irecl)
      read(unit=iold,fmt=80,end=888)buffer1(1:irecl)
      goto 887
888   close(iold)
      close(inew)
      stop

999   close(iold) 
777   write(unit=inew,fmt=80)buffer2(1:irecl)
      read(unit=iextra,fmt=80,end=778)buffer2(1:irecl)
      goto 777
778   close(iextra)
      close(inew)
      stop
      end











