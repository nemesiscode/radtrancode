      program merge_multi_ascending
C     ******************************************************************
C     Program to merge together nfile no. of sequential line database files
C     Files must not overlap in wavelength/wavenumber space, but be contiguous
C     Purpose is to get a single file with the whole range (all lines)
C     The files must be input in ascending wavenumber.
C     Change ipfile dimension to allow more than 30 files (i.e. default ipfile(30))
C     Only works if tables don't overlap at all  
C     Pat Irwin		8/9/95
C     Mahmuda Afrin Badhan 18/02/2015 Merge nfiles together in ascending wavenumber
C     ******************************************************************
      implicit none
      integer nfile, I
      integer iold,inew,irecl
      character*256 buffer
      character*100 filename, opfile
      character*100 ipfile(30)
      
2     format(a)
      call prompt('Enter record length : ')
      read*,irecl
      call prompt('Enter number of files to be merged: ')
      read*,nfile
      
      call prompt('Enter output line data file : ')
      read(5,2)opfile
      open(inew,file=opfile,status='unknown')
      
      call prompt('Enter data file names (ascending wavenumber): ')
      do I=1,nfile
       READ(5,2)ipfile(I)
       print*,ipfile(I)
       write(filename,2)ipfile(I)
       print*,filename
       open(I+6,file=filename,status='old')
       iold = I+6
23     read(unit=iold,fmt=80,end=999)buffer(1:irecl)
       print*,buffer(1:80)
       if(buffer(2:2).eq.'#')then
        goto 23
       end if
887    write(unit=inew,fmt=80)buffer(1:irecl)
       read(unit=iold,fmt=80,end=999)buffer(1:irecl)
       goto 887
999    close(iold)
       if(I.eq.nfile)GOTO 998
      ENDDO
      
80    format(a)
      
998   close(inew)
      stop
      end
      
