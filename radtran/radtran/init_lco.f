      subroutine init_lco(lcofil)
      character*100 lcofil,buffer
      double precision percen
      integer idgas,i
      include '../includes/lcocom.f'
      ijlco=1
      open(12,file=lcofil,status='old')
1     format(a)      
      print*,'init_lco: reading in extra line continuum data'
      read(12,1)buffer
      write(6,1)buffer
      read(buffer(14:),*)vminlco,vmaxlco

      read(12,1)buffer
      write(6,1)buffer
      read(buffer(17:),*)lcobinsize,percen

      read(12,1)buffer
      write(6,1)buffer
      read(buffer(21:),*)tcalclco,iptflco

      read(12,1)buffer
      write(6,1)buffer
      read(buffer(11:),*)idlco,isolco

      read(12,1)buffer
      write(6,1)buffer
      read(buffer(8:),*)nbinlco

      print*,nbinlco

      do 10 i=1,nbinlco
       read(12,*)vlco(i),slco(i),lclse(i),lcwida(i),lcwids(i),
     &   lctdep(i)
       print*,vlco(i),slco(i),lclse(i),lcwida(i),lcwids(i),
     &   lctdep(i)
10    continue
      close(12)

      return

      end
