      subroutine readflags(runname,inormal,iray,ih2o,ich4,io3,
     &  inh3,iptf,imie, iuvscat)
C     ********************************************************************
C     Subroutine to read in the runname.fla file, which contains a number
C     of spectral processing flags.
C
C     Input variables
C      runname	character*100	calculation root file name
C
C     Output variables
C	 inormal	integer	Hydrogen ortho/para ratio. 0=eqm, 1=normal
C	 iray		integer	Rayleigh optical depth off(0) or on(1)
C	 ih2o		integer	Include H2O continuum (1)
C	 ich4		integer	Include Karkoschka methane continuum (1)
C	 io3		integer	Include Bass UV ozone absorption (1)
C	 inh3		integer	Include Lutz and Owen ammonia continuum
C	 iptf		integer	Use modified CH4 high-T part. func. (1)
C      imie		integer	Use PHASEN.DAT phase file (1) or hgphaseN.dat phase
C					file otherwise
C	 iuvscat 	integer 	Include additional UV opacity (Henrik Melin)
C
C	Output variables passed via a common block
C      inlte_flag real		Type of ntle approximation to use (0=none)
C					see radtran/radtran/calc_nlte_t for details
C
C     ********************************************************************
C     Pat Irwin	21/2/12	Original version
C	Nick Teanby	25/10/23	Added inlte_flag and common block to avoid
C				changing loads of function calls (this is only needed
C				by cirsradg_wave which is called from multiple places)
C	Nick Teanby	14/4/25	non-LTE stuff moved to an optional .lte file,  
C                             makes more sense than having it here
C					nLTE flags now handled by:
C                             cirsrtfg_wave.f, cirsradg_wave.f, calc_nlte_t.f
C
C     ********************************************************************
      character*(*) runname
      integer inormal,iray,ih2o,ich4,io3,iptf,imie,inh3,iuvscat

      call file(runname,runname,'fla')

      open(12,file=runname,status='old')
      read(12,*)inormal
      read(12,*)iray
      read(12,*)ih2o
      read(12,*)ich4
      read(12,*)io3
      read(12,*)inh3
      read(12,*)iptf
      read(12,*)imie
      read(12,*)iuvscat
      close(12)

      return

      end
