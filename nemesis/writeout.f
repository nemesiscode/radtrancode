      subroutine writeout(iform,runname,ispace,lout,ispec,xlat,xlon,
     1  npro,nvar,varident,varparam,nx,ny,y,yn,se,xa,sa,xn,err,ngeom,
     2  nconv,vconv)
C     $Id:
C     ***********************************************************************
C     Output the results of retrieval code
C
C     Input variables
C	runname		character*100	run name
C	iform		integer		Output format: 0 = regular
C						       1 = nemesisdisc
C						       2 = nemesisPT
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
C
C     Pat Irwin		29/7/96
C     Steve Smith       10/12/96 Added extra space to xn op format
C     Steve Smith       17/3/97 New version for anneal_tng
C     Pat Irwin		17/10/03 Revised and recommented for Nemesis
C     Pat Irwin		11/5/12	 Updated for variable formats and new Nemesis
C
C     ***********************************************************************
      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'

      integer ny,nx,i,ngeom,igeom,lout,ispec,nsubspec,iform
      real y(my),xa(mx),xn(mx),err(mx),se(my),yn(my)
      real err1,xerr1,sa(mx,mx),xfac,xlat,xlon
      integer nconv(mgeom),j,ioff,varident(mvar,3),nvar,npro
      integer nxtemp,ivar,np,ix,iflag,ispace,npvar
      integer logflag
      real xa1,ea1,xn1,en1,iav
      real vconv(mgeom,mconv),varparam(mvar,mparam)
      real relerr
      character*100 runname
      logical fexist



      if(iform.eq.1)then
C      Output from Nemesisdisc - see if solar file is present
       CALL file(runname,runname,'sol')  
       inquire(file=runname,exist=fexist)
      endif

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

C      Default Nemesisdisc format
       elseif(iform.eq.1)then
        if(fexist)then
         write(lout,*)'Units are F_plan/F_star'
         xfac=1.0
        else
         write(lout,*)'Units are nW cm-2 (cm-1)-1'
         xfac=1e9
        endif

C      NemesisPT format
       elseif(iform.eq.2) then
        write(lout,*)'Units are 100*Planet_area/Stellar_area (i.e. %)'
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

C      Nemesisdisc format
       elseif(iform.eq.1)then
        if(fexist)then
         write(lout,*)'Units are F_plan/F_star'
         xfac=1.0
        else
         write(lout,*)'Units uW cm-2 um-1'
         xfac=1e6
        endif

C      NemesisPT format
       elseif(iform.eq.2)then
        write(lout,*)'Units are 100*Planet_area/Stellar_area (i.e. %)'
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

        else
C        Going back to default
         write(lout,1000)i,vconv(igeom,j),y(i)*xfac,err1*xfac,
     & xerr1,yn(i)*xfac,relerr
        endif

       end do
       ioff = ioff+nconv(igeom)
      enddo

C1000  format(1x,i4,1x,f10.4,1x,f10.4,1x,f10.4,1x,f7.2,1x,f10.4,1x,f7.2)
1000  format(1x,i4,1x,f10.4,1x,e15.8,1x,e15.8,1x,f7.2,1x,e15.8,1x,f9.5)
1010  format(1x,i4,1x,f10.4,1x,e15.8,1x,e15.8,1x,f7.2,1x,e15.8,1x,f9.5)
1020  format(1x,i4,1x,f9.4,1x,e12.6,1x,f9.4,1x,f6.2,1x,e12.6,1x,f6.2)


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
       if(varident(ivar,1).eq.888)np = varparam(ivar,1)

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


