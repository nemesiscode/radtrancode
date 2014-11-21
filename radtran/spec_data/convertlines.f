      program convertlines
C     *****************************************************************
C     Program to convert sequential access linedata output by Select to
C     another line data format. 
C
C     Output from Select should be for a single gas only as here you have
C     to define what the new isotope format is for the new format.
C
C     Pat Irwin		Original	11/6/12
C
C     *****************************************************************
      implicit none
      character*100 apfile,bpfile
      character*256 buffer
      integer nform,newrecl,isocode(20),niso,i,j,gasout
      INCLUDE '../includes/dbcom.f'

C      data isocode /26,36,28,27,38,37/

      
      print*,'Enter input filename : '
      print*,'NB: This file should contain lines from one gas only'

      READ(5,1)apfile
1     format(a)

      call prompt('Enter dbform and dbrecl for input file : ')
      read*,dbform,dbrecl

      call prompt('Enter number of isotopes for this gas (NISO): ')
      read*,niso

      call prompt('Enter new isotopes IDs for the NISO isotopes : ')
      read*,(isocode(i),i=1,niso)

      call prompt('Enter converted Gas ID : ')
      read*,GASOUT

      call prompt('Enter output filename : ')
      READ(5,1)bpfile

      call prompt('Enter dbform and dbrecl for output file : ')
      read*,nform,newrecl

      open(12,file=apfile,status='old')
      open(13,file=bpfile,status='unknown')

      
10    continue
      read(12,1,end=999)buffer
C     Skip past header, if any
      if(buffer(2:2).eq.'#') goto 10

      call rdline(buffer)
      print*,'LNISO',LNISO
      do i=1,NISO
       if(lniso.eq.isocode(i))j=i
      enddo
      print*,'LNISO',LNISO,j
      lniso = j
      lnid = gasout
      call printline(buffer,nform,newrecl)
      write(13,1)buffer(1:newrecl)

      goto 10

999   close(12)
      close(13)
      end      
