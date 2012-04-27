      subroutine gsetrad(runname,iscat,nmu,mu,wtmu,isol,dist,
     1 lowbc,galb,nf,nconv,vconv,fwhm,ispace,gasgiant,
     2 layht,nlayer,laytyp,layint,sol_ang,emiss_ang,aphi,xlat,lin,
     3 nvar,varident,varparam,nx,xn,jalb,jtan,jpre,tsurf,xmap)
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

      integer nconv,lin,ispace,ncont1
      real xlat,fwhm,xlatx,tsurf
      integer nlayer,laytyp,iscat,nx,nxx,ncont,nx1
      integer layint
      real layht,xod(maxcon),xscal(maxcon)
      real vconv(mconv)
      integer nmu,isol,lowbc,nf,flagh2p,jalb,jtan,jpre
      integer jsurfx,jalbx,jtanx,jprex,nprox,icheck,icont
      double precision mu(maxmu),wtmu(maxmu)
      real dist,galb,xn(mx),xnx(mx),aphi,emiss_ang,sol_ang
      real stx(mx,mx),xdnu
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      integer jpara
      character*100 runname,buffer,aname

      integer nvar,ivar,varident(mvar,3),i,j,nalb,nalb1
      real varparam(mvar,mparam),alb(maxsec),valb(maxsec)
      integer nvarx,varidentx(mvar,3),ivarx
      real varparamx(mvar,mparam),od1
      logical gasgiant
      integer NN,NDUST
      REAL DUSTH(MAXLAY),DUST(MAXCON,MAXLAY)

      print*,'gsetrad, lin = ',lin

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

      call subprofretg(runname,ispace,iscat,gasgiant,xlat,nvar,
     1  varident,varparam,nx,xn,jpre,ncont,flagh2p,xmap)

      print*,'subprofretg OK'
     
C     Add in check to make sure we don't use variable para-H2
C     for wavenumbers greater than 1500cm-1 
      if(ispace.eq.0)then
       do j=1,nconv
        if(vconv(j).gt.1500.)flagh2p=0
       enddo
      endif

      if(lin.eq.1.or.lin.eq.3)then

       call readxtmp(runname,xlatx,nvarx,varidentx,varparamx,nprox,
     1  nxx,xnx,stx,jsurfx,jalbx,jtanx,jprex)

       call stripvar(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  nxx,xnx)

       call subprofretgx(runname,ispace,iscat,gasgiant,nvarx,varidentx,
     1  varparamx,nxx,xnx,jprex,flagh2p)

       do ivarx = 1,nvarx
        if(varident(ivarx,1).eq.888)then

C        ********* reset surface albedo spectrum  **********
         nalb = varparamx(ivarx,1)
         call file(runname,runname,'alb')
         open(9,file=runname,status='old')
         read(9,*)nalb1
         if(nalb1.ne.nalb)then
          print*,'Error in gsetrad nalbx <> nalb1'
          print*,nalb,nalb1
         endif

         do i=1,nalb
          read(9,*)valb(i),alb(i)
         enddo
         close(9)

         open(9,file=runname,status='unknown')
         write(9,*)nalb
         do i=1,nalb
          write(9,*)valb(i),exp(xnx(jalbx+i-1))
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
       if (varident(ivar,1).eq.888)then
C       ********* surface albedo spectrum retrieval **********
        nalb = varparam(ivar,1)
        call file(runname,runname,'alb')
        open(9,file=runname,status='old')
        read(9,*)nalb1
        if(nalb1.ne.nalb)then
         print*,'Error in gsetrad nalb <> nalb1'
         print*,nalb,nalb1
        endif

        do i=1,nalb
         read(9,*)valb(i),alb(i)
        enddo
        close(9)

        open(9,file=runname,status='unknown')
        write(9,*)nalb
        do i=1,nalb
         write(9,*)valb(i),exp(xn(jalb+i-1))
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

      call gwritepat(runname,gasgiant,iscat,sol_ang,
     1 emiss_ang,nconv,vconv,fwhm,layht,nlayer,
     2 laytyp,layint,flagh2p)

      call mod_scatter(runname,ncont,nmu,mu,wtmu,isol,dist,lowbc,
     1 galb,nf,sol_ang,emiss_ang,aphi)


      icheck=0
      do ivar=1,nvar
       if(varident(ivar,3).eq.9.or.varident(ivar,3).eq.10)icheck=1
      enddo

C     Compute the drv file to get the aerosol optical depths
      if(icheck.eq.1) then

        print*,'Check',runname,runname
        call subpath(runname)

        call readdustod(runname,ncont1,xod)
       
        do i=1,ncont1
         xscal(i)=1.0
        enddo

        open(12,file='aerosol.prf',status='old')
54       READ(12,1)BUFFER     
1        FORMAT(A)
         IF(BUFFER(1:1).EQ.'#') GOTO 54
         READ(BUFFER,*)NN,NDUST
         DO 105 J=1,NN
           READ(12,*)DUSTH(J),(DUST(I,J),I=1,NDUST)
105      CONTINUE
        close(12)

        nx1=0

C        do i=1,nx
C         print*,'gsetrad',i,xn(i),exp(xn(i))
C        enddo

        do ivar=1,nvar

         if(varident(ivar,1).le.100)then

          if(varident(ivar,3).eq.0)nx1=nx1+NN
          if(varident(ivar,3).eq.1)nx1=nx1+2
          if(varident(ivar,3).eq.2)nx1=nx1+1
          if(varident(ivar,3).eq.3)nx1=nx1+1
          if(varident(ivar,3).eq.4)nx1=nx1+3
          if(varident(ivar,3).eq.5)nx1=nx1+NN
          if(varident(ivar,3).eq.6)nx1=nx1+2
          if(varident(ivar,3).eq.7)nx1=nx1+2
          if(varident(ivar,3).eq.8)nx1=nx1+3
          if(varident(ivar,3).eq.9)then
              icont=abs(varident(ivar,1))
              od1=exp(xn(nx1+1))
              xscal(icont)=xod(icont)/od1

C              print*,'gsetrad',icont,nx1+1
C              print*,'gsetrad',od1,xod(icont),xscal(icont)

              do j=1,NN
               dust(icont,j)=dust(icont,j)/xscal(icont)
              enddo
              nx1=nx1+3
          endif
          if(varident(ivar,3).eq.10)then
              icont=int(varparam(ivar,1))
              od1=exp(xn(nx1+3))
              xscal(icont)=xod(icont)/od1

              do j=1,NN
               dust(icont,j)=dust(icont,j)/xscal(icont)
              enddo
              nx1=nx1+4
          endif

         else 

          nx1=nx1+1

         endif

        enddo

        open(12,file='aerosol.prf',status='unknown')
         BUFFER='# simple.prf'
         WRITE(12,1)BUFFER     
         WRITE(12,*)NN,NDUST
         DO 106 J=1,NN
           WRITE(12,*)DUSTH(J),(DUST(I,J),I=1,NDUST)
C           print*,DUSTH(J),(DUST(I,J),I=1,NDUST)
106      CONTINUE
        close(12)

C       check that rescaling has happened correctly
        print*,'Checking'
        call subpath(runname)

      endif


      return

      end
    
