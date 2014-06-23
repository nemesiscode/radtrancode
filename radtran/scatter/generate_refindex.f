      program generate_refindex
C     ********************************************************
C     Program to generates real part of RI given only imaginary part
C     imaginary part.
C
C     Pat Irwin		23/6/14
C     ********************************************************

      implicit none
      integer mpoint,npoints,iwave,i
      real lambda,r1,r2,xm,vm,nm
      parameter (mpoint=2000)
      real x1,x2,dx,n(mpoint),k(mpoint),x(mpoint),v(mpoint),nr(mpoint)
      character*100 ipfile
 
      call prompt('Enter name of imaginary ref. index file : ')
      read(5,1)ipfile
1     format(a)

      open(12,file=ipfile,status='old')
       read(12,*)npoints
       read(12,*)iwave
C            NB: iwave=1=wavelength, iwave=2=wavenumber 
       read(12,*)xm,nm
       do i=1,npoints
        read(12,*)x(i),k(i)
        v(i)=x(i)
        if(iwave.eq.1)v(i)=1e4/x(i)      
       enddo
      close(12)
      vm=xm
      if(iwave.eq.1)vm=1e4/vm


      call kk_new_sub(npoints, v, k, vm, nm, n)

      open(12,file='genrefindex.dat',status='unknown')
      write(12,*)npoints
      write(12,*)iwave
      do i=1,npoints
       write(12,*)x(i),n(i),k(i)
      enddo
      close(12)

      end
