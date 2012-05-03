      subroutine read_scatter1(radfile,nmu,mu1,wt1,isol,dist,lowbc,galb,
     1ncont,liscat,lnorm,lncons,lcons,scatfile,
     2sol_ang,emiss_ang,aphi,nf)
C     $Id: read_scatter1.f,v 1.3 2011-06-17 15:40:27 irwin Exp $
C     *******************************************************************
C     Subroutine to read in the new format .sca files with angles
C     substituted for IMU0, IMU and the defunct KJ removed.
C
C     Pat Irwin 	19/5/97
C
C     *******************************************************************
      integer maxcon,ulog
      parameter(maxcon=10,ulog=17)
      integer lowbc,liscat(maxcon),lncons(maxcon),lnorm(maxcon)
      integer maxmu,isol,ncont,nf
      parameter (maxmu=21)
      double precision mu1(maxmu), wt1(maxmu), galb
      real dist,lcons(maxcon,10),aphi,sol_ang,emiss_ang
      character*30 header
      character*30 scatfile(maxcon)
      character*100 radfile

      write(ulog,223)radfile
223   format(1X,'Reading scattering information from ',A)

      open(49,file=radfile,status='old')
      read(49,1)header
      write(ulog,1)header
1     format(a)
      read(49,*)nmu
      write(ulog,501)nmu
501   format(1X,'nmu = ',i3)
      do 10 i=1,nmu
       read(49,*)mu1(i),wt1(i)
C       print*,mu1(i),wt1(i)
       write(ulog,*)mu1(i),wt1(i)
10    continue
      read(49,*)isol
      write(ulog,*)'isol = ',isol
      read(49,*)dist
      write(ulog,*)'dist = ',dist
      read(49,*)lowbc
      write(ulog,*)'lowbc = ',lowbc
      read(49,*)galb
      write(ulog,*)'galb = ',galb
      read(49,*)sol_ang,emiss_ang
      write(ulog,*)'sol_ang, emiss_ang = ',sol_ang,emiss_ang
      read(49,*)aphi
      write(ulog,*)'aphi = ',aphi
      read(49,*)nf
      write(ulog,*)'nf = ',nf
      read(49,*)ncont
      write(ulog,*)'ncont = ',ncont
      do 100 i=1,ncont
       read(49,*)liscat(i)
       write(ulog,*)'liscat(i) = ',liscat(i)
       read(49,*)lnorm(i)
       write(ulog,*)'lnorm(i) = ',lnorm(i)
       if(liscat(i).lt.4)then
        read(49,*)lncons(i)
        write(ulog,*)'lncons(i) = ',lncons(i)
        do 30 j=1,lncons(i)
         read(49,*)lcons(i,j)
         write(ulog,*)'lcons(i,j) = ',lcons(i,j)
30      continue
       else
        read(49,3)scatfile(i)
        write(ulog,*)'scatfile = ',scatfile(i)
       end if
3      format(1X,a30)
100   continue

      close(49)

      return
      end
 



