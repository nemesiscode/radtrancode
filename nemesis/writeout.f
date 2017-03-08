      subroutine writeout(iform,runname,ispace,lout,ispec,xlat,xlon,
     1  npro,nvar,varident,varparam,nx,ny,y,yn,se,xa,sa,xn,err,ngeom,
     2  nconv,vconv,gasgiant,jpre,jrad,jlogg,iscat,lin)
C     $Id:
C     ***********************************************************************
C     Output the results of retrieval code
C
C     Input variables
C	runname		character*100	run name
C	iform		integer		Output format: 0 = radiance
C						       1 = F_plan/F_star
C						       2 = 100*A_plan/A_star
C						       3 = planet spectral flux
C							    i.e. F_plan
C	ispace		integer		0=cm-1,1=microns
C	lout		integer		Output unit number
C	ispec		integer		Spectrum ID
C	xlat		real		Latitude
C	xlon		real		Longitude 
C	npro		integer		Number of levels in .ref file
C	nvar		integer		Number of variable profiles
C	varident(mvar,3) integer	Identity of profiles and 
C						parameterisation
C	varparam(mvar,mparam) real 	Extra parameters as required
C	nx		integer		Number of elements in state vector
C	ny		integer		Number of elements in measurement
C						vector
C	y(my)		real		Measured spectrum
C	yn(my)		real		Best calculated spectrum
C	se(my)		real		Variances of measured spectrum
C	xa(mx)		real		A priori measurement vector
C	sa(mx,mx)	real		A priori covariance matrix
C	xn(mx)		real		Retrieved vector x
C	err(mx)		real		Retrieved errors
C	ngeom		integer		Number of observation geometries
C	nconv(mgeom)	integer		Number of convolution wavenumbers
C					 at each observation angle
C	vconv(mgeom,mconv) real		Convolution wavenumbers
C	gasgiant	logical		Gas giant flag
C	jpre		integer		Indicates if pressure retrieval 
C					performed.
C	jrad		integer		Indicates if radiusretrieval 
C					performed.
C	jlogg		integer		Indicates if surface gravity retrieval 
C					performed.
C       iscat		integer		Flag to indicate scattering calc.
C	lin		integer		Previous retrieval flag
C
C     Pat Irwin		29/7/96
C     Steve Smith       10/12/96 Added extra space to xn op format
C     Steve Smith       17/3/97 New version for anneal_tng
C     Pat Irwin		17/10/03 Revised and recommented for Nemesis
C     Pat Irwin		11/5/12	 Updated for variable formats and new Nemesis
C     Pat Irwin		30/1/13	Updated to write out last retrieved .prf
C				 files.
C
C     ***********************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      include '../radtran/includes/planrad.f'
      include 'arraylen.f'

      integer ny,nx,i,ngeom,igeom,lout,ispec,nsubspec,iform
      real y(my),xa(mx),xn(mx),err(mx),se(my),yn(my)
      real err1,xerr1,sa(mx,mx),xfac,xlat,xlon
      integer nconv(mgeom),j,ioff,varident(mvar,3),nvar,npro
      integer nxtemp,ivar,np,ix,iflag,ispace,npvar
      integer logflag,xflag,jpara,flagh2p,jpre,ncont,npro1
      integer iscat,lin,jrad,jlogg,iplanet
      real xa1,ea1,xn1,en1,iav,xdnu,RADIUS,Grav
      parameter (Grav=6.672E-11)
      real vconv(mgeom,mconv),varparam(mvar,mparam)
      real relerr
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      real xmapx(maxv,maxgas+2+maxcon,maxpro)
      real xlatx,varidentx(mvar,3),varparamx(mvar,mparam)
      real stx(mx,mx)
      integer nxx,xnx(mx),nvarx,nprox,jtanx,jprex,jradx
      integer icread,jloggx,ierr,ierrx
      character*100 runname,aname,buffer,cellfile
      logical gasgiant,cellexist

      integer iwave,imode,nx1,nvmr,np1,nmode,nwave,max_wave
      integer max_mode,inorm
      parameter (max_wave = 1000,max_mode = 10)
      real r0,v0,clen,vm,nm,lambda0
      real wave(max_wave),xsec(max_mode,max_wave,2),nimag(max_wave)
      real v1(max_wave),k1(max_wave),vm1,n1(max_wave)
      real nreal(max_wave),minlam
      real srefind(max_wave,2),parm(3),rs(3)



C     See if a cell file is present. If so we have two outputs per
C     wavelength: SB and WB
      call file(runname,cellfile,'cel')
      inquire(file=cellfile,exist=cellexist)
      icread=0
      if(cellexist)icread=1
     
1     FORMAT(A)

C     Output ny instead of nconv to keep format of mre file the same
      write(lout,901) ispec,ngeom,ny,nx,ny,
     & '   ! ispec,ngeom,ny,nx,ny'
901   format(1x,i4,i3,i5,i4,i5,a)
      write(lout,*)xlat,xlon,'Latitude, Longitude'

      if(ispace.eq.0) then
C      Wavenumber space

C      Default format
       if(iform.eq.0)then
        write(lout,*)'Radiances expressed as nW cm-2 sr-1 cm'
        xfac=1e9

C      F_plan/F_star format
       elseif(iform.eq.1)then
         write(lout,*)'F_plan/F_star Ratio of planet'
         xfac=1.0

C      Spectral power format
       elseif(iform.eq.3)then
      write(lout,*)'Spectral Radiation of planet: W (cm-1)-1'
         xfac=1e18

C      NemesisPT format
       elseif(iform.eq.2) then
        write(lout,*)'Transit depth: 100*Planet_area/Stellar_area'
        xfac=1.

C      Transmission*solar_flux
       elseif(iform.eq.4) then
        write(lout,*)'Solar flux: W cm-2 (cm-1)-1'
        xfac=1.

C      Default
       else
        print*,'Error in writeout - iform not defined. Default=0'
        write(lout,*)'Radiances expressed as nW cm-2 sr-1 cm'
        xfac = 1e9
       endif

      else
C      Wavelength space

C      Default format
       if(iform.eq.0)then
        write(lout,*)'Radiances expressed as uW cm-2 sr-1 um-1'
        xfac = 1e6

C      F_plan/F_star format
       elseif(iform.eq.1)then
         write(lout,*)'F_plan/F_star Ratio of planet'
         xfac=1.0

C      Spectral irradiance format
       elseif(iform.eq.3)then
      write(lout,*)'Spectral Radiation of planet: W um-1'
         xfac=1e18

C      NemesisPT format
       elseif(iform.eq.2)then
        write(lout,*)'Transit depth: 100*Planet_area/Stellar_area'
        xfac=1.

C      Transmission*solar_flux
       elseif(iform.eq.4) then
        write(lout,*)'Solar flux: W cm-2 um-1'
        xfac=1.

C      Default format
       else
        print*,'Error in writeout - iform not defined. Default=0'
        write(lout,*)'Radiances expressed as uW cm-2 sr-1 um-1'
        xfac=1e6
       endif

      endif


      print*,'Writeout: ispace,iform,xfac : ',ispace,iform,xfac

      write(lout,*)
     1  '  i  lambda  R_meas     error   %err  R_fit     Diff%'
      ioff = 0
      do igeom=1,ngeom
       do j=1,nconv(igeom)
        i = ioff+j
        err1 = sqrt(se(i))
        if(y(i).ne.0)then
         xerr1 = abs(100.0*err1/y(i))
         relerr = abs(100.0*(y(i)-yn(i))/y(i))
        else
         xerr1 = -1.0
         relerr = -1.0
        endif
        if(xerr1.gt.100.0)xerr1=100.0
        if(relerr.gt.100.0)relerr=100.0

        if(iform.eq.0)then
         write(lout,1000)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr

        elseif(iform.eq.1)then
         write(lout,1010)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr

        elseif(iform.eq.2)then
         write(lout,1020)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr

        elseif(iform.eq.3)then
         write(lout,1030)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr

        else
C        Going back to default
         write(lout,1000)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr
        endif

       end do
       ioff = ioff+nconv(igeom)

       if(icread.eq.1)then
        do j=1,nconv(igeom)
         i = ioff+j
         err1 = sqrt(se(i))
         if(y(i).ne.0)then
          xerr1 = abs(100.0*err1/y(i))
          relerr = abs(100.0*(y(i)-yn(i))/y(i))
         else
          xerr1 = -1.0
          relerr = -1.0
         endif
         if(xerr1.gt.100.0)xerr1=100.0
         if(relerr.gt.100.0)relerr=100.0

         if(iform.eq.0)then
          write(lout,1000)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr

         elseif(iform.eq.1)then
          write(lout,1010)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr

         elseif(iform.eq.2)then
          write(lout,1020)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr

         elseif(iform.eq.3)then
          write(lout,1030)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr

         else
C         Going back to default
          write(lout,1000)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr
         endif

        end do
        ioff = ioff+nconv(igeom)

       endif

      enddo

C1000  format(1x,i4,1x,f10.4,1x,e15.8,1x,e15.8,1x,f7.2,1x,e15.8,1x,f9.5)
1000  format(1x,i4,1x,f14.8,1x,e15.8,1x,e15.8,1x,f7.2,1x,e15.8,1x,f9.5)
1010  format(1x,i4,1x,f10.4,1x,e15.8,1x,e15.8,1x,f7.2,1x,e15.8,1x,f9.5)
1020  format(1x,i4,1x,f9.4,1x,e12.6,1x,e12.6,1x,f6.2,1x,e12.6,1x,f6.2)
1030  format(1x,i4,1x,f10.4,1x,e15.8,1x,e15.8,1x,f7.2,1x,e15.8,1x,f9.5)


      write(lout,*)' '
      write(lout,*)'nvar = ',nvar

      nxtemp=0

      do 299 ivar=1,nvar
       write(lout,*)'Variable ',ivar
       write(lout,*)(varident(ivar,j),j=1,3)
       write(lout,*)(varparam(ivar,j),j=1,mparam)
       np=1
       if(varident(ivar,1).le.100)then
         np = npvar(varident(ivar,3),npro)
       endif
       if(varident(ivar,1).eq.888)np = int(varparam(ivar,1))
       if(varident(ivar,1).eq.444)np = 2+int(varparam(ivar,1))
       if(varident(ivar,1).eq.222)np = 8
       if(varident(ivar,1).eq.223)np = 9
       if(varident(ivar,1).eq.224)np = 9
       if(varident(ivar,1).eq.225)np = 11
       if(varident(ivar,1).eq.226)np = 8

       write(lout,*)
     &  '   i, ix, xa          sa_err       xn          xn_err'
       do i = 1,np
        ix = nxtemp+i

        xa1 = xa(ix)
        ea1 = sqrt(abs(sa(ix,ix)))
        xn1 = xn(ix)
        en1 = err(ix)

        iflag = logflag(varident(ivar,1),varident(ivar,3),i)

        if(iflag.eq.1)then
          xa1 = exp(xa1)
          ea1 = xa1*ea1
          xn1 = exp(xn1)
          en1 = xn1*en1
        endif

        write(lout,1015)i,ix,xa1,ea1,xn1,en1

       enddo

       nxtemp = nxtemp+np

299   continue


1015  format(1x,i4,i4,' ',e12.5,e12.5,' ',e12.5,e12.5)

      return

      end


