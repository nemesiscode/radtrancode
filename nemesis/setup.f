      subroutine setup(runname,gasgiant,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf,nphi,layht,tsurf,nlayer,laytyp,layint)
C     $Id:
C     ******************************************************************
C
C     Subroutine to read in the default atmospheric layering 
C     and scattering details from the '.set' setup file.
C
C     Input variables 
C       runname         character*100    Root run name.
C	gasgiant	logical		Indicates if planet is a gas giant
C
C     Output variables
C	nmu		integer		Number of zenith ordinates
C	mu(maxmu)	double prec.	Cos(zenith) points
C	wtmu(maxmu)	double prec.	Quadrature weights
C	isol		integer		Sunlight on/off
C	dist		real		Solar distance (AU)
C	lowbc		integer		lower boundary condition
C	galb		real		ground albedo
C	nf		integer		Required number of Fourier components
C	nphi		integer		Number of azimuth angles
C       layht           real            Base height of lowest layer
C	tsurf		real		Surface temperature (if planet is
C					  not a gas giant)
C	nlayer		integer		Number if vertical levels to split
C					  atmosphere into
C	laytup		integer		How layering is performed (Radtran)
C	layint		integer		How layer amounts are calculated
C					   (Radtran)	
C
C     Pat Irwin	23/5/97		Original
C     Pat Irwin 17/10/03	Tidied for Nemesis
C
C     ******************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'

      integer isol,lowbc,nf,intbuf,nlayer,laytyp
      integer nmu,i,layint,nphi,lun
      real dist,galb,tsurf
      double precision mu(maxmu),wtmu(maxmu)
      real realbuf,layht
      character*100 runname,setfile
      character*80 buffer
      logical gasgiant

      lun=73

      call file(runname,setfile,'set')

      open(lun,file=setfile,status='old')

      read(lun,1)buffer			! Header
1     format(a)

      read(lun,1)buffer
      nmu = intbuf(buffer)		! Number of zenith points
      do i=1,nmu
       read(lun,*)mu(i),wtmu(i)		! Zenith points
      enddo

      read(lun,1)buffer
      nf = intbuf(buffer)               ! Number of fourier components

      read(lun,1)buffer
      nphi = intbuf(buffer)             ! Azimuth angles for Fourier analysis

      read(lun,1)buffer
      isol = intbuf(buffer)             ! Sunlight switch

      read(lun,1)buffer
      dist = realbuf(buffer)            ! Solar distance (AU)

      read(lun,1)buffer
      lowbc = intbuf(buffer)            ! Lower boundary cond.

      read(lun,1)buffer
      galb = realbuf(buffer)            ! ground albedo

      read(lun,1)buffer
      tsurf = realbuf(buffer)		! Surface temperature

      if(tsurf.le.0.0)then
       print*,'Error in nemesis/setup.f. Surface temperature must'
       print*,'be greater than 0.0'
       stop
      endif

      read(lun,1)buffer			! Header

      read(lun,1)buffer
      layht = realbuf(buffer)		! Lower altitude

      read(lun,1)buffer
      nlayer = intbuf(buffer)		! Number of layers

      read(lun,1)buffer
      laytyp = intbuf(buffer)		! Layer type

      read(lun,1)buffer
      layint = intbuf(buffer)		! Layint

      read(lun,1)buffer			!header
   
      close(lun)

C      print*,'Setup Read OK'

C     rewrite .set file to keep correctly formatted
      open(lun,file=setfile,status='unknown')

      buffer='*********************************************************'
      write(lun,1)buffer
      write(lun,102)nmu
102   format(1x,'Number of zenith angles : ',i2)

      do i=1,nmu
       write(lun,*)mu(i),wtmu(i)		! Zenith points
      enddo

      write(lun,103)nf
103   format(1x,'Number of fourier components : ',i2)

      write(lun,104)nphi
104   format(1x,'Number of azimuth angles for fourier analysis : ',i3)

      write(lun,105)isol
105   format(1x,'Sunlight on(1) or off(0) : ',i2)

      write(lun,106)dist
106   format(1x,'Distance from Sun (AU) : ',f7.3)

      write(lun,107)lowbc
107   format(1x,'Lower boundary cond. Thermal(0) Lambert(1) : ',i2)

      write(lun,108)galb
108   format(1x,'Ground albedo : ',f7.3)

      write(lun,113)tsurf
113   format(1x,'Surface temperature : ',f8.3)

      buffer='*********************************************************'
      write(lun,1)buffer

      write(lun,109)layht
109   format(1x,'Alt. at base of bot.layer (not limb) : ',f8.3)

      write(lun,111)nlayer
111   format(1x,'Number of atm layers : ',i3)

      write(lun,121)laytyp
121   format(1x,'Layer type : ',i2)

      write(lun,112)layint
112   format(1x,'Layer integration : ',i2)


      buffer='*********************************************************'
      write(lun,1)buffer

      close(lun)
      
      return
      end



      real function realbuf(buffer)
      implicit none
      character*80 buffer,temp          
      integer icopy,i,j
      real x
      
      icopy = 0
      j = 0
      do i=1,80
       if (icopy.eq.1) then
        j=j+1
        temp(j:j) = buffer(i:i)
       end if
       if (buffer(i:i).eq.':') icopy = 1
       buffer(i:i)=' '
      end do

      read(temp(1:j),*)x

      realbuf = x
      
      return
      end
        
      integer function intbuf(buffer)
      implicit none
      character*80 buffer,temp
      integer icopy,j,i,ix
      
      icopy = 0
      j = 0
      do i=1,80
       if (icopy.eq.1) then
        j=j+1  
        temp(j:j) = buffer(i:i)
       end if  
       if (buffer(i:i).eq.':') icopy = 1
       buffer(i:i)=' '
      end do   
      
      read(temp(1:j),*)ix
 
      intbuf = ix
 
      return 
      end

