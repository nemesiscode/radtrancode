      program testbrown
      implicit none
      include '../radtran/includes/arrdef.f'
      real press(maxpro),temp(maxpro),dtempdx(maxpro,5)
      integer npro,i,atmp(4),height(maxpro)
      real T0,Teff,tau0,n,alpha
      character*100 ipfile,buffer
      print*,'Enter name of .prf file to read in : '
      read(5,1)ipfile
1     format(a)

      open(12,file=ipfile,status='old')
    
      read(12,1)buffer
      read(12,*)(atmp(i),i=1,4)
      npro=atmp(3) 
      do i=1,3
       read(12,1)buffer
      enddo
      do i=1,npro
       read(12,*)height(i),press(i),temp(i)
      enddo

      print*,'Enter values for tau0,n,Teff,alpha,T0'
      read*,tau0,n,Teff,alpha,T0

      call tbrownrc(npro,press,T0,Teff,tau0,n,alpha,temp,
     &  dtempdx)


      end
