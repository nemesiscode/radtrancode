      subroutine readsurfem(xscfil,nem,vem,emissivity)
C     $Id:
C     ***********************************************************************
C     Routine for reading in surfce emissivity spectrum file.
C
C     Input varaibles
C	xscfil	character*100	Run-name
C
C     Output variables
C	nem	integer	Number of wavenumber points in table
C	vem(maxsec) real	Table wavenumbers
C	emissivity(maxsec) real Table emissivities
C
C     Pat Irwin		25/11/03
C
C     ***********************************************************************
      integer nem
      include '../includes/arrdef.f'
      real emissivity(maxsec),vem(maxsec)
      character*100 xscfil,buffer
      integer idiag,iquiet
      common/diagnostic/idiag,iquiet

      iunit=19

      call file(xscfil,xscfil,'sur')
      if(idiag.gt.0)then
       print*,'readsurfem: surface emissivity file : ',xscfil
      endif
      open(iunit,file=xscfil,status='old')
1     format(a)
54    read(iunit,1)buffer
      if(buffer(1:1).eq.'#')goto 54
      read(buffer,*)nem
      if(nem.gt.maxsec)then
       print*,'nem > maxsec ',nem,maxsec
       stop
      endif
      do 344 i=1,nem
        read(iunit,*)vem(i),emissivity(i)
344   continue

      close(iunit)

      return

      end
