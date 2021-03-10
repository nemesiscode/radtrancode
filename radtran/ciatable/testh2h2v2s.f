      program testh2h2v2s
      implicit double precision (a-h,o-z)
      dimension f(601),alf(601),xout(3,601)
      integer normal
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet


      idiag=1
 
      fnumin=7000.
      fnumax=13000.
      dnu=10.

      print*,'Enter temperature : '
      read*,temp

      normal=0
      print*,normal,temp,fnumin,fnumax,dnu
      call h2h2_v2s(normal,temp,fnumin,fnumax,dnu,nf,f,alf)

      print*,nf
      do i=1,nf 
       xout(1,i)=f(i)
       xout(2,i)=alf(i)
       print*,f(i),alf(i)
      enddo

      normal=1
      call h2h2_v2s(normal,temp,fnumin,fnumax,dnu,nf,f,alf)

      do i=1,nf 
       xout(3,i)=alf(i)
      enddo

      open(12,file='test.dat',status='unknown')
      write(12,*)nf

      do i=1,nf 
       write(12,*)(xout(k,i),k=1,3)
      enddo

      close(12)

      end
