      subroutine readflags(runname,inormal,iray,ih2o,ich4,io3,
     &  iptf)
C     ********************************************************************
C     Subroutine to read in the runname.fla file, which contains a number
C     of spectral processing flags.
C
C     Input variables
C     	runname	character*100	calculation root file name
C
C     Output variables
C	inormal	integer	Hydrogen ortho/para ratio. 0=eqm, 1=normal
C	iray	integer	Rayleigh optical depth off(0) or on(1)
C	ih2o	integer	Include H2O continuum (1)
C	ich4	integer	Include Karkoschka methane continuum (1)
C	io3	integer	Include Bass UV ozone absorption (1)
C	iptf	integer	Use modified CH4 high-T part. func. (1)
C
C     Pat Irwin	21/2/12
C
C     ********************************************************************
      character*(*) runname
      integer inormal,iray,ih2o,ich4,io3,iptf

      call file(runname,runname,'fla')

      open(12,file=runname,status='old')
      read(12,*)inormal
      read(12,*)iray
      read(12,*)ih2o
      read(12,*)ich4
      read(12,*)io3
      read(12,*)iptf
      close(12)

      return

      end
