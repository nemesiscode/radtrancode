      program sortlines
C     ******************************************************************
C     Program to sort a sequential-access linedata file
C
C     Pat Irwin		8/9/95
C     ******************************************************************
      implicit none
      integer iold,iextra,inew,irecl,ichoice,iswap
      parameter (iold=12,iextra=13,inew=14)           
      character*256 buffer1,buffer2
      character*100 ipfile1,opfile,tmpfile,swapfile
      real vv1,vv2

2     format(a)
      call prompt('Enter record length : ')
      read*,irecl
      call prompt('Enter line data file : ')
      read(5,2)ipfile1
      tmpfile='temp.aaa'
      call prompt('Enter output line data file : ')
      read(5,2)opfile

17    continue

      open(iold,file=ipfile1,status='old')
      open(inew,file=tmpfile,status='unknown')

      read(unit=iold,fmt=80,end=999)buffer1(1:irecl)
      read(unit=iold,fmt=80,end=999)buffer2(1:irecl)
80    format(a)
      iswap = 0
            
      if(irecl.eq.80)then
       read(buffer1,10)vv1
       read(buffer2,10)vv2
10     format(f10.3)
      else
       read(buffer1(4:15),11)vv1
       read(buffer2(4:15),11)vv2
11     format(f12.6)
      endif      


20    if(vv1.le.vv2)then
       write(unit=inew,fmt=80)buffer1(1:irecl)
       buffer1 = buffer2
       vv1 = vv2
      else
       write(unit=inew,fmt=80)buffer2(1:irecl)
       iswap = 1
      end if

      read(unit=iold,fmt=80,end=999)buffer2(1:irecl)
      if(irecl.eq.80)then
       read(buffer2,10)vv2
      else
       read(buffer2(4:15),11)vv2
      endif      

      goto 20

999   continue
      close(iold) 
      write(unit=inew,fmt=80)buffer1(1:irecl)
      close(inew)

      if(iswap.gt.0)then
       swapfile = ipfile1
       ipfile1 = tmpfile
       tmpfile = swapfile
       goto 17
      else
       open(iold,file=tmpfile,status='old')
       open(inew,file=opfile,status='unknown')
18     read(unit=iold,fmt=80,end=888)buffer1(1:irecl)
       write(unit=inew,fmt=80)buffer1(1:irecl)
       goto 18

888    close(inew)
       close(iold) 

      endif

      end











