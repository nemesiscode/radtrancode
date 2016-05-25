      program calc_chisq
C     ************************************************************************
C     Simple program to print out fitted chisq/ny calculated from the final
C     retrieved .mre file.
C
C     Pat Irwin  27/4/16
C
C     ************************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'
      real ydat(7,mconv),sum
      character*100 ipfile,buffer
      integer idum(5),i,j,ny
  
      call prompt('Enter .mre filename : ')
      read(5,1)ipfile
1     format(a)


      open(12,file=ipfile,status='old')
      read(12,1)buffer
      read(12,*)(idum(j),j=1,5)
      ny=idum(3)
      read(12,1)buffer
      read(12,1)buffer
      read(12,1)buffer
      do 10 i=1,ny
       read(12,*)(ydat(j,i),j=1,7)
10    continue
      close(12)

      sum=0.

      do 20 i=1,ny
       sum=sum+((ydat(3,i)-ydat(6,i))/ydat(4,i))**2
20    continue

      sum=sum/float(ny)

      print*,'Chisq/ny = ',sum

      end
