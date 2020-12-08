      SUBROUTINE GWRITEPATV(RUNNAME,GASGIANT,ISCAT,SOL_ANG,EMISS_ANG,
     1 NCONV,VCONV,FWHM,LAYHT,NLAYER,LAYTYP,LAYINT,FLAGH2P)
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
C					4=single-scattering in spherical atm.
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
C
C     Pat Irwin	29/7/96		Original
C     Pat Irwin 17/10/03	Tidied for Nemesis
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
c  ** variable for reflected atmos
      CHARACTER*100 rflfile,dummy
      REAL angle_inc,angle_rfl
      INTEGER layer_rfl
      LOGICAL fexist
c  ** variables for solar reflected cloud **
      real refl_cloud_albedo
      logical reflecting_atmos
      common /refl_cloud_params/refl_cloud_albedo,reflecting_atmos

      integer cellngas,cellid(maxgas),celliso(maxgas),icread
      real cellength,cellpress,celltemp,cellvmr(maxgas)
      common/celldat/icread,cellngas,cellid,celliso,cellvmr,cellength,
     1  cellpress,celltemp
      INTEGER IPZEN

      integer iread,nspt,iform
      real solrad,swave(maxbin),srad(maxbin)
      common /solardat/iread,iform,solrad,swave,srad,nspt

      
      CALL FILE(RUNNAME,RUNNAME,'pat')

      OPEN(31,FILE=RUNNAME,STATUS='UNKNOWN')


      print*,'gwritepatV: iscat=',iscat
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

      print*,'AA ISCAT = ',ISCAT
      IF(ISCAT.GT.0)THEN
       TEXT = 'fcloud model fcloud.prf'
       WRITE(31,1)TEXT
       WRITE(31,*)' '
      ENDIF
          
2     format(a,i3)
3     format(a,f8.2)

      LAYHT1 = LAYHT
      LAYANG = 0.0

      WRITE(31,1)'layer'
      WRITE(31,2)'nlay ',nlayer
      WRITE(31,3)'layht ',layht1
      WRITE(31,3)'layang ',layang
      WRITE(31,2)'layint ',layint
      WRITE(31,2)'laytyp ',laytyp

      WRITE(31,1)' '

      CALL FILE(RUNNAME,RUNNAME,'zen')
      INQUIRE(FILE=RUNNAME,EXIST=FEXIST)
      IPZEN=0
      IF(FEXIST)THEN
       OPEN(12,FILE=RUNNAME,STATUS='OLD')
       READ(12,*)IPZEN
      ENDIF

      WRITE(31,1)'atm'
      TEXT='nadir'
      LAYBOT=1
      E1 = 0.0
      WRITE(31,4)TEXT,E1,LAYBOT,IPZEN
4     FORMAT(A6,F7.2,I4,I4)
    
      IF(ISCAT.EQ.5)THEN
        WRITE(31,1)'notherm'
        WRITE(31,1)'scatter'
        WRITE(31,1)'netflux'
      ELSE
        print*,'Error in gwritepatV - ISCAT <> 5'
        stop
      ENDIF
      WRITE(31,1)'nowf'
      WRITE(31,1)'nocg'
      WRITE(31,1)'noabsorb'
      WRITE(31,1)'binbb'
      WRITE(31,1)'nobroad'
      WRITE(31,1)' '
      
      reflecting_atmos = .FALSE.

      WRITE(31,1)'clrlay'
      WRITE(31,1)'nocombine'

      CALL FILE(RUNNAME,RUNNAME,'pra')

      OPEN(12,FILE=RUNNAME,STATUS='OLD',ERR=222)
221    CONTINUE
       READ(12,1,END=222)DUMMY
       WRITE(31,1)' '
       WRITE(31,1)'process '//DUMMY             
       GOTO 221
222   CLOSE(12)

      CLOSE(31)

      RETURN
      END



