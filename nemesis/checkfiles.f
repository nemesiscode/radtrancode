      subroutine checkfiles(runname)
C     $Id:
C     *****************************************************************
C     Subroutine to ensure input files to Nemesis run are consistent 
C     with each other.
C
C     Input files
C	runname		character*100	Run name
C
C     Pat Irwin		21/10/03
C
C     *****************************************************************
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      character*100 runname,buffer
      integer npro,nvmr,npro1,ncont
      logical gasgiant

      if(my.gt.idim)then
       print*,'Error on checkfiles'
       print*,'my > idim', my,idim
       stop
      endif

      call readrefhead(runname,npro,nvmr,gasgiant)

      print*,'Runname : ',runname
      print*,'NPRO, NVMR : ',npro,nvmr
1     format(a)

      open(12,file='aerosol.ref',status='old')
C      First skip header
54     read(12,1)buffer
       if(buffer(1:1).EQ.'#') goto 54
       read(buffer,*)npro1,ncont
      close(12)

      if(npro1.ne.npro)then
       print*,'Error on checkfiles'
       print*,'npro in aerosol.ref not the same as in runname.ref'
       print*,'npro in runname.ref = ',npro
       print*,'npro in aerosol.ref = ',npro1
       stop
      endif

      call file(runname,runname,'xsc')
      open(12,file=runname,status='old')
       read(12,*)ncont1
      close(12)
   
      if(ncont1.ne.ncont)then
       print*,'Error on checkfiles'
       print*,'ncont in runname.xsc not the same as in aerosol.ref'
       stop
      endif

      return

      end
