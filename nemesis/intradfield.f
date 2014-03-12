      subroutine intradfield(runname,ispace,xlat,nwave,vwave,
     1 nconv,vconv,gasgiant,lin,nvar,varident,
     2 varparam,jsurf,jalb,jtan,jpre,nx,xn)
C
C     $Id:
C     **************************************************************
C     Subroutine to calculate the internal radiation field of a planetary 
C     atmosphere.
C
C     Input variables:
C       runname(60)   character Name of run.
C       ispace           integer Indicates if wavelengths in vconv and
C                               vwave are in wavenumbers(0) or
C                               wavelengths (1)
C       iscat           integer 0=non-scattering, 1=scattering
C       fwhm            real    Desired FWHM of final spectrum
C       nwave           integer Total number of calculation wavelengths
C       vwave(mwave)    real    Calculation wavelengths
C       nconv           integer Total number of convolution wavelengths
C       vconv(mconv)    real    Convolution wavelengths
C       gasgiant        logical Indicates if planet is a gas giant
C       lin             integer indicates role of previous retrieval (if any)
C       nvar    integer Number of variable profiles (gas,T,aerosol)
C       varident(nvar,3) integer identity of constituent to retrieved and
C                                       parameterisation method
C       varparam(nvar,mparam) real Additional arameters constraining profile.
C       jsurf           integer Position of surface temperature element in
C                               xn (if included)
C       jalb            integer position of first surface albedo element in
C                               xn (if included)
C       jtan            integer position of tangent ht. correction element in
C                               xn (if included)
C       jpre            integer position of tangent pressure element in
C                               xn (if included)
C       nx              integer Number of elements in state vector
C       xn(mx)          real    State vector
C
C     Output variables
C       none.
C
C     Pat Irwin	3/7/07		Original, modified from forwardnogX.f
C
C     **************************************************************

      implicit none
      integer i,j,ispace,ulog
      parameter (ulog=17)
      integer ioff,lin,iscat
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/gascom.f'
      include 'arraylen.f'
      real xlat,xref,dx
      integer layint,inormal,itype,nlayer,laytyp
      integer iray,iptf,imie,imie1
      integer nwave,ix,ix1,iav,nwave1
      integer nsol,iloop,i1,i2
      real vwave(mwave),interpem
      real calcout(maxout3),fwhm,planck_wave,output(maxout3)
      integer nx,nconv,npath,ioff1,ioff2,nconv1
      real vconv(mconv)
      real layht,tsurf
      real xn(mx),dtr
      integer jsurf,jalb,jtan,jpre,nem
      integer nphi,ipath,fintrad
      character*100 fintname
      integer nmu,isol,lowbc,nf,nx2,kiter
      real dist,galb,sol_ang,emiss_ang,z_ang,aphi
      double precision mu(maxmu),wtmu(maxmu)
      character*100 runname,logname,intname
      real xmap(maxv,maxgas+2+maxcon,maxpro),xrad
      common /imiescat/imie1

      integer nvar,varident(mvar,3)
      real varparam(mvar,mparam)
      logical gasgiant
      real vem(maxsec),emissivity(maxsec)

      common/intrad/fintrad,fintname


      iscat=1
      fwhm=0.0
      dtr = 3.1415927/180.0
      intname = 'intrad'

      call setup(intname,gasgiant,nmu,mu,wtmu,isol,
     1 dist,lowbc,galb,nf,nphi,layht,tsurf,nlayer,laytyp,layint)

      if(isol.eq.0)then
       nsol=1
       sol_ang = 0.0
       emiss_ang=0.0
       aphi = 0.0
       fintrad=0
      else
       nsol=nmu
       aphi = 0.0
       emiss_ang=0.0
       fintrad=1
      endif

      call file(intname,logname,'log')

      open(ulog,file=logname,status='unknown')

      print*,'intradfield: gasgiant',gasgiant
C     If planet is not a gas giant then we need to read in the surface 
C      emissivity spectrum
      if(.not.gasgiant)then
          call readsurfem(runname,nem,vem,emissivity)
      else
          nem=2
          vem(1)=-100.0
          vem(2)=1e7
          emissivity(1)=1.0
          emissivity(2)=1.0
      endif
         
      CALL READFLAGS(runname,INORMAL,IRAY,IH2O,ICH4,IO3,IPTF,IMIE)
      IMIE1=IMIE

      itype=13			        ! scloud13flux

C     Need to run this over the different solar zenith angles required 
C     and then concatenate the internal.fld files together so that they 
C     can be properly accessed by cirsrad. Needs a major shake-up to do 
C     this. Have started by making code in cirsrad deal with multiple 
C     fourier components, but routines such as impflux and dumpflux need 
C     a lot more modification
 
      do 1000 iloop=1,nsol
       if(fintrad.eq.1)then
        sol_ang = acos(sngl(mu(iloop)))/dtr
        fintname = 'internal**.fld'
        i1=int(iloop/10)
        i2=iloop-i1*10
        fintname(9:9)=char(i1+48)
        fintname(10:10)=char(i2+48)
C        nf = int(20*sol_ang/90.0)
         nf=20
       else
        fintname='internal.fld'
        nf=0
       endif
       print*,'fintrad,fintname',fintrad,fintname
C      Set up all files for a direct cirsrad run
       print*,'Calling gsetrad - gasgiant = ',gasgiant
       call gsetrad(intname,iscat,nmu,mu,wtmu,isol,dist,
     1     lowbc,galb,nf,nconv,vconv,fwhm,ispace,gasgiant,
     2     layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,lin,
     3     nvar,varident,varparam,nx,xn,jalb,jtan,jpre,tsurf,xmap)

       print*,'Calling cirsrtf_wave'
       call CIRSrtf_wave(intname, dist, inormal, iray, fwhm, ispace,
     1     vwave,nwave,npath, output, vconv, nconv, itype, nem, vem, 
     2     emissivity, tsurf, calcout)

1000  continue


C     Now internal field calculated for all solar zenith angles (if isol 
C     eq 1). Combine into a single lookup file if required.

      if(fintrad.eq.1)call combineflux(nsol,tsurf,nem,vem,
     1 emissivity,ispace)

      close(ulog)

      return

      end
