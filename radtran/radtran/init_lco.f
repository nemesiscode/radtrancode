      subroutine init_lco(lcofil)
      character*100 lcofil,buffer
      double precision percen
      integer idgas,i
      include '../includes/lcocom.f'
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      ijlco=1
      open(12,file=lcofil,status='old')
1     format(a)      
      if(idiag.gt.0)then
       print*,'init_lco: reading in extra line continuum data'
       print*,lcofil
      endif
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

      if(idiag.gt.0)print*,nbinlco

C      read(12,1)buffer
C      write(6,1)buffer
      do 10 i=1,nbinlco
       read(12,*)vlco(i),slco(i),lclse(i),lcwida(i),lcwids(i),
     &   lctdep(i)
C       print*,vlco(i),slco(i),lclse(i),lcwida(i),lcwids(i),
C     &   lctdep(i)
10    continue
      close(12)

      return

      end
