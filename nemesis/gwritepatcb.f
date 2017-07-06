      SUBROUTINE GWRITEPATCB(RUNNAME,GASGIANT,ISCAT,SOL_ANG,EMISS_ANG,
     1 NCONV,VCONV,FWHM,LAYHT,NLAYER,LAYTYP,LAYINT,FLAGH2P,NCLOUD,
     2 cbpbot,cbptop,ncblaycloud,cbodepth,cbfsh,IFLAGCB)

C     $Id:
C     *******************************************************************
C     Subroutine to write out the .pat file needed for a CIRSradg run.
C
C     Input variables:
C       RUNNAME         character*100    Run root name
C	GASGIANT	logical		Indicates if planet is a Gas Giant
C       ISCAT           integer         0=non-scattering
C					1=plane-parallel scattering
C					2=limb/near-limb scattering
C					3=single-scattering
C       SOL_ANG         real            Solar zenith angle
C       EMISS_ANG       real            Emission zenith angle
C       NCONV           integer         Number of convolution wavelengths
C       VCONV(MCONV)   	real 	        Convolution wavelengths
C	FWHM		real		FWHM of convolved spectrum
C       LAYHT           real            Base height of bottom layer 
C					 (if nadir)
C       NLAYER          integer         Number of layers
C       LAYTYP          integer         Type of layering
C       LAYINT          integer         Type of layer integration
C	FLAGH2P		integer		Equals 1 if para-H2 is variable
C	IFLAGCB		integer		Equals 1 if cloud layering scheme follows creme brulee model
C
C     Output variables:
C	NCLOUD		integer		Number of cloud layers (=3 for creme brulee scheme)
C	CBPBOT		real		Pressure level of bottom of each cloud layer
C	CBPTOP		real		Pressure level of top of each cloud layer
C	NCBLAYCLOUD	integer		Number of layers to split each cloud layer into
C	CBODEPTH	real		Optical depth of each cloud layer
C	CBFSH		real		Fractional scale height of each cloud layer
C
C     Ashwin Braude	08.12.2016
C
C     *******************************************************************

      IMPLICIT NONE
      include '../radtran/includes/arrdef.f'
      INCLUDE 'arraylen.f'
      CHARACTER*100 RUNNAME
      REAL LAYHT,EMISS_ANG,VCONV(MCONV),LAYANG,LAYHT1,DELV,FWHM
      REAL SOL_ANG,E1
      INTEGER ISCAT,LAYTYP,LAYINT,NLAYER,NCONV,LAYBOT,FLAGH2P,I
      CHARACTER*80 TEXT
      LOGICAL GASGIANT
      INTEGER IFLAGCB,KEXT
c  ** variable for reflected atmos
      CHARACTER*100 rflfile,dummy
      REAL angle_inc,angle_rfl
      INTEGER layer_rfl
      LOGICAL fexist
c  ** variables for solar reflected cloud **
      real refl_cloud_albedo
      logical reflecting_atmos
      integer mcloud,icloud
      parameter (mcloud=10)
      real cbpbot(mcloud),cbptop(mcloud),cbodepth(mcloud),cbfsh(mcloud)
      integer ncblaycloud(mcloud),nlayg,nlaybot,nlaytop,ncloud

      common /refl_cloud_params/refl_cloud_albedo,reflecting_atmos

      integer cellngas,cellid(maxgas),celliso(maxgas),icread
      real cellength,cellpress,celltemp,cellvmr(maxgas)
      common/celldat/icread,cellngas,cellid,celliso,cellvmr,cellength,
     1  cellpress,celltemp
      INTEGER IPZEN

      
      CALL FILE(RUNNAME,RUNNAME,'pat')

      OPEN(31,FILE=RUNNAME,STATUS='UNKNOWN')

      WRITE(31,1)' '
      TEXT='interval'
      WRITE(31,1)TEXT
1     FORMAT(A)
      IF(NCONV.GT.1)THEN
       DELV = (VCONV(NCONV)-VCONV(1))/FLOAT(NCONV-1)
      ELSE
       DELV=FWHM
      ENDIF 
      WRITE(31,*)VCONV(1),VCONV(NCONV),DELV,FWHM
      WRITE(31,1)'   24  0  0'
      WRITE(31,1)' '
 
      TEXT(1:10)='spec data '
      CALL FILE(RUNNAME,RUNNAME,'kls')
      TEXT(11:60)=RUNNAME(1:50)
      WRITE(31,1)TEXT
      WRITE(31,1)' '


      CALL FILE(RUNNAME,RUNNAME,'prf')
      TEXT=' '
      TEXT(1:6)='model '
      TEXT(7:60)=RUNNAME(1:54)
      WRITE(31,1)TEXT
      WRITE(31,*)' '

      TEXT = 'dust model aerosol.prf'
      WRITE(31,1)TEXT
      WRITE(31,*)' '


      CALL FILE(RUNNAME,RUNNAME,'xsc')
      TEXT=' '
      TEXT(1:13)='dust spectra '
      TEXT(14:60)=RUNNAME(1:47)
      WRITE(31,1)TEXT
      WRITE(31,*)' '

      IF(FLAGH2P.EQ.1)THEN
       TEXT = 'fparah2 model parah2.prf'
       WRITE(31,1)TEXT
       WRITE(31,*)' '
      ENDIF

      IF(ISCAT.GT.0)THEN
       TEXT = 'fcloud model fcloud.prf'
       WRITE(31,1)TEXT
       WRITE(31,*)' '
      ENDIF
          
2     format(a,i3)
3     format(a,f8.2)

      LAYHT1 = LAYHT
      LAYANG = 0.0
C     Is observation at limb? (coded with -ve emission angle where
C     sol_ang is then the tangent altitude)
      IF(EMISS_ANG.LT.0.0)THEN
       LAYHT1 = SOL_ANG
       LAYANG = 90.0
      ENDIF       

      print*,'gwritepatcb: IFLAG',IFLAGCB,NCLOUD
      WRITE(31,1)'ncomplayer'!this may need changing
      if(ncloud.eq.3)then
       nlayg=11-ncblaycloud(2)
       nlaybot=7
       nlaytop=8
      else
       nlayg=1
       nlaybot=7
       nlaytop=8
      endif
      write(31,2)'nlayg ',nlayg
      write(31,2)'nlaybot ',nlaybot
      write(31,2)'nlaytop ',nlaytop
      write(31,2)'ncloud ',ncloud
      do icloud=1,ncloud
       write(31,*)cbpbot(icloud),cbptop(icloud),ncblaycloud(icloud),
     &  cbodepth(icloud),cbfsh(icloud)
      enddo
      WRITE(31,3)'layang ',layang
      WRITE(31,2)'layint ',layint

      WRITE(31,1)' '

      CALL FILE(RUNNAME,RUNNAME,'zen')
      INQUIRE(FILE=RUNNAME,EXIST=FEXIST)
      IPZEN=0
      IF(FEXIST)THEN
       OPEN(12,FILE=RUNNAME,STATUS='OLD')
       READ(12,*)IPZEN
      ENDIF

      WRITE(31,1)'atm'
      IF(EMISS_ANG.GE.0.0) THEN
       TEXT='nadir'
       LAYBOT=1
       E1 = EMISS_ANG
       IF(ISCAT.EQ.1)E1 = 0.0
       WRITE(31,4)TEXT,E1,LAYBOT,IPZEN
4      FORMAT(A6,F7.2,I4,I4)
      ELSE
       TEXT='limb'
       LAYBOT=1
       WRITE(31,5)TEXT,LAYBOT
5      FORMAT(A6,I4)
      ENDIF
    
      IF(ISCAT.EQ.0)THEN
       WRITE(31,1)'therm'
       WRITE(31,1)'noscatter'
      ELSE
       IF(ISCAT.EQ.1.OR.ISCAT.EQ.2)THEN
        WRITE(31,1)'notherm'
        WRITE(31,1)'scatter'
        IF(ISCAT.EQ.2)write(31,1)'nearlimb'
       ELSE
        WRITE(31,1)'notherm'
        WRITE(31,1)'single'
       ENDIF
      ENDIF
      WRITE(31,1)'nowf'
      WRITE(31,1)'nocg'
      WRITE(31,1)'noabsorb'
      WRITE(31,1)'binbb'
      WRITE(31,1)'nobroad'
      WRITE(31,1)' '

      
c  ** if rflfile exists then insert stuff for reflecting path **
      CALL FILE(RUNNAME,rflfile,'rfl')
      inquire(file=rflfile,exist=fexist)
C      print*,'file=',rflfile
C      print*,'fexist=',fexist
      if ( fexist ) then
c	** read in parameters from .rfl file **
         open(67,file=rflfile,status='old')
         read(67,*) dummy
         read(67,*) angle_inc
         read(67,*) angle_rfl
         read(67,*) layer_rfl
         read(67,*) refl_cloud_albedo
         close(67)
c	** replace angle_inc=sol_ang and angle_rfl=emiss_ang **
	   angle_inc = SOL_ANG
	   angle_rfl = EMISS_ANG
c	** write reflecting atmosphere to pat file **
         WRITE(31,1)'reflatm'
         WRITE(31,6) angle_inc,angle_rfl,layer_rfl
6	   FORMAT('angles ',f8.3,x,f8.3,i3)
         WRITE(31,1)'nocg'
         WRITE(31,1)'noabsorb'
         WRITE(31,1)' '
         reflecting_atmos = .TRUE.
      else 
         reflecting_atmos = .FALSE.
      endif

      if(icread.eq.1)then
       write(31,1)' '
       write(31,1)'cell'
       WRITE(31,7) cellngas
7	   FORMAT('gases ',i3)
       do i=1,cellngas
        write(31,8)cellid(i),celliso(i),cellvmr(i)
8       format(i3,i3,f8.5)
       enddo
       write(31,9)cellength
9      format('length ',f9.4)
       write(31,1)'scr'
       write(31,*)cellpress,celltemp
       write(31,1)' '
      endif      
      
      WRITE(31,1)'clrlay'
      if(icread.ne.1)then
       WRITE(31,1)'nocombine'
      else
       WRITE(31,1)'combine'
      endif

      CALL FILE(RUNNAME,RUNNAME,'pra')

      OPEN(12,FILE=RUNNAME,STATUS='OLD',ERR=222)
221    CONTINUE
       READ(12,1,END=222)DUMMY
       WRITE(31,1)' '
       WRITE(31,1)DUMMY             
       GOTO 221
222   CLOSE(12)

      CLOSE(31)

      RETURN
      END



