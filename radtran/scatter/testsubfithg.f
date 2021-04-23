      program testsubfithg
C     ******************************************************************
      implicit none
      real f,g1,g2,pi,x,x1,x2,calpha,sum1,sum2,dtr
      real du,henyey,alpha
      integer max_thet,nphase,i,j,itemp,nover,nc
      parameter (max_thet=100,pi=3.1415927)
      real theta(max_thet),phase(max_thet),rms
      character*100 ipfile,opfile
     
      call prompt('Enter name of test phase function file ')
      read(5,1)ipfile
1     format(a)

      open(12,file=ipfile,status='old')
      read(12,*)nphase
      do i=1,nphase
       read(12,*)theta(i),phase(i)
      enddo
      close(12)

C     Subfithgm is expecting phase functions normalised
C     to 4pi, rather than 1, so pass phase1 instead of
C     phase.
      call subfithgm(nphase,theta,4*pi*phase,f,g1,g2,rms)

      call prompt('Enter output file name : ')
      read(5,1)opfile

      open(13,file=opfile,status='unknown')
      write(13,*)nphase

      dtr=pi/180.0
      sum1=0.0
      sum2=0.0
      do i=1,nphase-1
       du=cos(theta(i)*dtr)-cos(theta(i+1)*dtr)
       alpha = cos(theta(i)*dtr)
       x1 = 0.25*henyey(alpha,f,g1,g2)/pi
       alpha = cos(theta(i+1)*dtr)
       x2 = 0.25*henyey(alpha,f,g1,g2)/pi
       sum1=sum1+0.5*(phase(i)+phase(i+1))*du
       sum2=sum2+0.5*(x1+x2)*du
      enddo
      sum1=sum1*2*pi
      sum2=sum2*2*pi
      write(13,*)sum1,sum2
      write(13,*)f,g1,g2
      do i=1,nphase
       alpha = cos(theta(i)*dtr)
       x = 0.25*henyey(alpha,f,g1,g2)/pi
       write(13,*)theta(i),phase(i),x
      enddo

      close(13)

      end

