      subroutine init_partf
C     ****************************************************************
C     Subroutine to initialise and read in the partition functions for
C     gases where we require tabulated partition function arrays rather
C     than cubic interpolations as is done by default.
C
C     Pat Irwin		27/10/14
C 
C     ****************************************************************
      character*100 partfil,buffer
      include '../includes/partfcom.f'

      partfil='partfextra.dat'
      call datarchive(partfil)
      open(12,file=partfil,status='old')
1     format(a)
10    read(12,1)buffer
      if(buffer(1:1).eq.'#'.or.buffer(2:2).eq.'#') GOTO 10
      read(buffer,*)NPARTEXTRA
      if(NPARTEXTRA.gt.MPARTEXTRA)then
       print*,'Error in init_partf.f: NPARTEXTRA > MPARTEXTRA'
       print*,NPARTEXTRA,MPARTEXTRA
       stop
      endif
      print*,'NPARTEXTRA = ',NPARTEXTRA
      do 20 i=1,NPARTEXTRA
       read(12,*)IDEXTRA(i),ISOEXTRA(i)
20    continue
      do 30 i=1,NPARTEXTRA
25     read(12,1)buffer
       if(buffer(1:1).eq.'#'.or.buffer(2:2).eq.'#') GOTO 25
       read(buffer,*)NTABEXTRA(i)

       if(NTABEXTRA(i).gt.MTABEXTRA)then
        print*,'Error in init_partf.f: NTABEXTRA(I) > MTABEXTRA'
        print*,I,NTABEXTRA(i),MTABEXTRA
        stop
       endif

       do 40 j=1,NTABEXTRA(i)
        read(12,*)TEMPEXTRA(i,j),PARTFEXTRA(i,j)
40     continue

30    continue

      close(12)
  
      ireadextra=-1
      return

      end
