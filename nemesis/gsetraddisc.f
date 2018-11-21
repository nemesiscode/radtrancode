      subroutine gsetraddisc(runname,iscat,nmu,mu,wtmu,isol,dist,
     1 lowbc,galb,nf,nconv,vconv,fwhm,ispace,gasgiant,
     2 layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,xlon,lin,
     3 nvar,varident,varparam,nx,xn,jalb,jxsc,jtan,jpre,tsurf,xmap)
C     $Id:
C     ************************************************************************
C     Subroutine to write out the .pat, .prf, .xsc and .sca and aerosol 
C     files needed for a CIRSradg run. Routine also returns xmap
C     (calculated by subprofretg) which relates the functional
C     derivatives calculated by CIRSRADG to the state vector elements
C
C     Input variables
C       runname         character*100    Root run name.
C	iscat		integer		0=thermal emission	
C					1=plane parallel scattering 
C					2= limb/near-limb scattering 
C	nmu		integer		Number of zenith ordinates
C	mu(maxmu)	double pr.	Cos(zenith angles)
C	wtmu(maxmu)	double pr.	Quadrature weights
C	isol		integer		Sunlight on(1)/off(0)
C       dist            real            Solar distance (AU)
C       lowbc           integer         lower boundary condition
C       galb            real            ground albedo
C       nf              integer         Required number of Fourier components
C       nconv           integer         Number of calculation wavelengths
C       vconv(mconv)   	real            Calculation wavelength array
C	fwhm		real		FWHM of convolved spectrum
C	ispace		integer		0=cm-1, 1= wavelengths
C	gasgiant	logical		Indicates if planet is a gas giant
C	layht		real		Altitude of base layer (if non-limb,
C					this is set to -emiss_ang for limb
C					calculations)
C	nlayer		integer		Number of layers
C	laytyp		integer		How layers are separated
C	layint		integer		How layer amounts are integrated
C	sol_ang		real		Solar zenith angle
C	emiss_ang	real		Thermal emission angle
C	aphi		real		azimuth angle
C	xlat		real		latitude of observation
C	xlon		real		longitude of observation
C	lin		integer		Integer to indicate role of previous
C                                       retrieval (if any)
C       nvar    	integer 	Number of variable profiles 
C					  (e.g. gas,T,aerosol)
C       varident(mvar,3) integer 	identity of constituent to retrieved
C					 and parameterisation
C       varparam(mvar,mparam) real 	Additional parameters constraining
C					  profile.
C	nx		integer		Number of elements in state vector
C	xn(mx)		real		state vector
C	jalb		integer		Position of surface albedo spectrum
C					in xn
C	jxsc		integer		Position of x-section spectrum
C					in xn
C	jtan		integer		Position of tangent height correction
C					in xn
C	jpre		integer		Position of tangent pressure
C					in xn
C	tsurf		real		Surface temperature
C
C    Output variables      
C	xmap(maxv,maxgas+2+maxcon,maxpro) real Mapping array relate
C		functional derivatives calculated by CIRSRADG to the
C		state vector elements
C	tsurf		real		Updated surface temperature (if
C					using previously retrieved value)
C
C     Pat Irwin	17/10/03	Modified from oxcirsg for Nemesis
C
C     ************************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'

      integer nconv,lin,ispace,xflag
      real xlat,fwhm,xlatx,tsurf,xlon,xlonx
      integer nlayer,laytyp,iscat,nx,nxx,ncont
      integer layint,ierr,ierrx,nxsc,ncont1,icont
      real layht
      real vconv(mconv)
      integer nmu,isol,lowbc,nf,flagh2p,jalb,jxsc,jtan,jpre
      integer jsurfx,jalbx,jtanx,jprex,nprox,jradx,jpara
      integer jloggx,jxscx
      double precision mu(maxmu),wtmu(maxmu)
      real dist,galb,xn(mx),xnx(mx),aphi,emiss_ang,sol_ang
      real stx(mx,mx),xdnu,xtest
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      real xmapx(maxv,maxgas+2+maxcon,maxpro)
      character*100 runname,aname,buffer

      integer nvar,ivar,varident(mvar,3),i,j,nalb,nalb1
      real varparam(mvar,mparam),alb(maxsec),valb(maxsec)
      integer nvarx,varidentx(mvar,3),ivarx
      real varparamx(mvar,mparam),xsc(maxsec,maxgas)
      real ssa(maxsec,maxgas)

      logical gasgiant

      print*,'gsetraddisc, lin = ',lin

1     format(a)

C     Look to see if the CIA file refined has variable para-H2 or not.
      call file(runname,runname,'cia')

      open(12,FILE=runname,STATUS='OLD')
       read(12,1) aname
       read(12,*) xdnu
       read(12,*) jpara
      close(12)

      flagh2p=0
      if(jpara.ne.0)then
       flagh2p=1
      endif
 
      xflag=0
      call subprofretg(xflag,runname,ispace,iscat,gasgiant,xlat,xlon,
     1  nvar,varident,varparam,nx,xn,jpre,ncont,flagh2p,xmap,ierr)


      if(lin.eq.1.or.lin.eq.3)then

       call readxtmp(runname,xlatx,xlonx,nvarx,varidentx,varparamx,
     1  nprox,nxx,xnx,stx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,jloggx)

       call stripvar(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  varparam,nxx,xnx)

       xflag=1
       call subprofretg(xflag,runname,ispace,iscat,gasgiant,xlat,xlon,
     1  nvarx,varidentx,varparamx,nxx,xnx,jprex,ncont,flagh2p,xmapx,
     2  ierrx)

       do ivarx = 1,nvarx

        if(varidentx(ivarx,1).eq.887)then

C        ********* Adjust cloud cross-section spectrum  **********
         nxsc = int(varparamx(ivarx,1))
         icont = int(varparamx(ivarx,2))
         call file(runname,runname,'rxsc')
         open(9,file=runname,status='old')
51       read(9,1)buffer
         if(buffer(1:1).eq.'#')goto 51
         read(buffer,*)ncont1
         if(ncont1.ne.ncont)then
          print*,'Error in gsetrad ncont1 <> ncont'
          print*,ncont1,ncont
          print*,'file : ',runname
         endif

         do i=1,nxsc
          read(9,*)valb(i),(xsc(i,j),j=1,ncont)
          read(9,*)(ssa(i,j),j=1,ncont)
         enddo
         close(9)

         call file(runname,runname,'xsc')
         open(9,file=runname,status='unknown')
         write(9,*)ncont
         do i=1,nxsc
          xsc(i,icont)=exp(xnx(jxscx+i-1))
          write(9,*)valb(i),(xsc(i,j),j=1,ncont)
          write(9,*)(ssa(i,j),j=1,ncont)
         enddo

        endif


        if(varident(ivarx,1).eq.888)then

C        ********* reset surface albedo spectrum  **********
         nalb = int(varparamx(ivarx,1))
         call file(runname,runname,'sur')
         open(9,file=runname,status='old')
54       read(9,1)buffer
         if(buffer(1:1).eq.'#')goto 54
         read(buffer,*)nalb1
         if(nalb1.ne.nalb)then
          print*,'Error in gsetraddisc nalbx <> nalb1'
          print*,nalb,nalb1
          print*,'file : ',runname
         endif

         do i=1,nalb
          read(9,*)valb(i),alb(i)
         enddo
         close(9)

         open(9,file=runname,status='unknown')
         write(9,*)nalb
         do i=1,nalb
          xtest = 1.-exp(xnx(jalbx+i-1))
C         Stop emissivity going negative
          if(xtest.lt.0.0)xtest=0.0
          write(9,*)valb(i),xtest
         enddo
         close(9)

        endif


       if (varident(ivar,1).eq.889)then
C       ********* surface albedo spectrum multiplier retrieval **********
        call file(runname,runname,'rsu')
        open(9,file=runname,status='old')
57      read(9,1)buffer
        if(buffer(1:1).eq.'#')goto 57
        read(buffer,*)nalb
        do i=1,nalb
         read(9,*)valb(i),alb(i)
        enddo
        close(9)

        call file(runname,runname,'sur')
        open(9,file=runname,status='unknown')
        write(9,*)nalb
        do i=1,nalb
         xtest = 1.0 - (1.-alb(i))*exp(xnx(jalbx))
C        Stop emissivity going negative
         if(xtest.lt.0.0)xtest=0.0
         write(9,*)valb(i),xtest
        enddo
        close(9)

       endif



        if(varidentx(ivarx,1).eq.777)then
C       ************ reset tangent heights   ***********
         if(emiss_ang.gt.0)then
          print*,'Can not do tangent ht correction for non-limb case'
         else 
          sol_ang = sol_ang+xnx(jtanx)
         endif

        endif

        if(varidentx(ivar,1).eq.999)then
C       ***************** Surface temperature correction ***********
         tsurf = xnx(jsurfx)
        endif

       enddo

      endif


      do ivar=1,nvar

C       print*,ivar

       if(varident(ivar,1).eq.887)then

C        ********* Adjust cloud cross-section spectrum  **********
         nxsc = int(varparam(ivar,1))
         icont = int(varparam(ivar,2))
         call file(runname,runname,'rxsc')
         open(9,file=runname,status='old')
52       read(9,1)buffer
         if(buffer(1:1).eq.'#')goto 52
         read(buffer,*)ncont1
         if(ncont1.ne.ncont)then
          print*,'Error in gsetrad ncont1 <> ncont'
          print*,ncont1,ncont
          print*,'file : ',runname
         endif

         do i=1,nxsc
          read(9,*)valb(i),(xsc(i,j),j=1,ncont)
          read(9,*)(ssa(i,j),j=1,ncont)
         enddo
         close(9)

         call file(runname,runname,'xsc')
         open(9,file=runname,status='unknown')
         write(9,*)ncont
         do i=1,nxsc
          xsc(i,icont)=exp(xn(jxsc+i-1))
          write(9,*)valb(i),(xsc(i,j),j=1,ncont)
          write(9,*)(ssa(i,j),j=1,ncont)
         enddo

       endif


       if (varident(ivar,1).eq.888)then
C       ********* surface albedo spectrum retrieval **********
        nalb = int(varparam(ivar,1))
        call file(runname,runname,'sur')
        open(9,file=runname,status='old')
56      read(9,1)buffer
        if(buffer(1:1).eq.'#')goto 56
        read(buffer,*)nalb1
        if(nalb1.ne.nalb)then
         print*,'Error in gsetraddisc nalb <> nalb1'
         print*,nalb,nalb1
         print*,'file : ',runname
        endif

        do i=1,nalb
         read(9,*)valb(i),alb(i)
        enddo
        close(9)

        open(9,file=runname,status='unknown')
        write(9,*)nalb
        do i=1,nalb
         xtest = 1.-exp(xnx(jalb+i-1))
C        Stop emissivity going negative
         if(xtest.lt.0.0)xtest=0.0
         write(9,*)valb(i),xtest
        enddo
        close(9)

       endif


       if (varident(ivar,1).eq.889)then
C       ********* surface albedo spectrum multiplier retrieval **********
        call file(runname,runname,'rsu')
        open(9,file=runname,status='old')
58      read(9,1)buffer
        if(buffer(1:1).eq.'#')goto 58
        read(buffer,*)nalb
        do i=1,nalb
         read(9,*)valb(i),alb(i)
        enddo
        close(9)

        open(9,file=runname,status='unknown')
        write(9,*)nalb
        do i=1,nalb
          xtest = 1.0 - (1.-alb(i))*exp(xn(jalb))
C         Stop emissivity going negative
          if(xtest.lt.0.0)xtest=0.0
         write(9,*)valb(i),xtest
        enddo
        close(9)

       endif

       if(varident(ivar,1).eq.777)then
C      ***************** Tangent height correction ***********
        if(emiss_ang.gt.0)then
         print*,'Can not do tangent ht correction for non-limb case'
        else
         sol_ang = sol_ang+xn(jtan)
        endif

       endif

      enddo

      call gwritepatdisc(runname,gasgiant,iscat,sol_ang,
     1 emiss_ang,nconv,vconv,fwhm,layht,nlayer,
     2 laytyp,layint,flagh2p)

      call mod_scatter(runname,ncont,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf,sol_ang,emiss_ang,aphi)


      return

      end
    
