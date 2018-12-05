      subroutine gsetradPT(runname,nconv,vconv,fwhm,ispace,iscat,
     1 gasgiant,layht,nlayer,laytyp,layint,xlat,lin,
     2 nvar,varident,varparam,nx,xn,xmap)
C     $Id:
C     ************************************************************************
C     Subroutine to write out the .pat, .prf, .xsc and .sca and aerosol 
C     files needed for a CIRSradg run. Routine also returns xmap
C     (calculated by subprofretg) which relates the functional
C     derivatives calculated by CIRSRADG to the state vector elements
C
C     Input variables
C       runname         character*100    Root run name.
C       nconv           integer         Number of calculation wavelengths
C       vconv(mconv)   	real            Calculation wavelength array
C	fwhm		real		FWHM of convolved spectrum
C	ispace		integer		0=cm-1, 1= wavelengths
C	iscat		integer		Scattering ID
C	gasgiant	logical		Indicates if planet is a gas giant
C	layht		real		Altitude of base layer
C	nlayer		integer		Number of layers
C	laytyp		integer		How layers are separated
C	layint		integer		How layer amounts are integrated
C	xlat		real		latitude of observation
C	lin		integer		Unit number of previous retrieval (if any)
C       nvar    	integer 	Number of variable profiles 
C					  (e.g. gas,T,aerosol)
C       varident(mvar,3) integer 	identity of constituent to retrieved
C					 and parameterisation
C       varparam(mvar,mparam) real 	Additional parameters constraining
C					  profile.
C	nx		integer		Number of elements in state vector
C	xn(mx)		real		state vector
C
C    Output variables      
C	xmap(maxv,maxgas+2+maxcon,maxpro) real Mapping array relate
C		functional derivatives calculated by CIRSRADG to the
C		state vector elements
C
C     Pat Irwin	17/10/03	Modified from oxcirsg for Nemesis
C
C     ************************************************************************

      implicit none
      include '../radtran/includes/arrdef.f'
      include 'arraylen.f'

      integer nconv,lin,ispace,iscat,xflag
      real xlat,fwhm,xlatx,hcorrx,tsurf,radius
      integer nlayer,laytyp,nx,nxx,ncont,jpara
      integer np,imode,inorm,iwave,nx1,np1
      integer layint,jsurfx,jalbx,jtanx,jprex,jradx,nprox
      integer jloggx,ierr,ierrx,jxscx
      real layht,xdnu,r0,v0
      real vconv(mconv),minlam,lambda0
      integer flagh2p, nmode, nwave, max_mode, max_wave
      double precision mu(maxmu),wtmu(maxmu)
      real xn(mx),xnx(mx),stx(mx,mx)
      real xmap(maxv,maxgas+2+maxcon,maxpro)
      real xmapx(maxv,maxgas+2+maxcon,maxpro)
      character*100 runname,aname,buffer

      integer nvar,varident(mvar,3),i,j,ivar,ivarx
      real varparam(mvar,mparam)
      integer nvarx,varidentx(mvar,3),jpre,icont
      real varparamx(mvar,mparam),xsc(maxsec,maxgas)
      real ssa(maxsec,maxgas)
      logical gasgiant

      parameter (max_wave = 1000)
      parameter (max_mode = 10)
      real wave(max_wave),xsec(max_mode,max_wave,2)
      real srefind(max_wave,2),parm(3),rs(3),vm,nm
      real nimag(max_wave),nreal(max_wave)

      jpre=-1
      
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
      call subprofretg(xflag,runname,ispace,iscat,gasgiant,xlat,
     1  nvar,varident,varparam,nx,xn,jpre,ncont,flagh2p,xmap,ierr)

      hcorrx=0.0

      if(lin.eq.1.or.lin.eq.3)then

       call readxtmp(runname,xlatx,nvarx,varidentx,varparamx,nprox,
     1 nxx,xnx,stx,jsurfx,jalbx,jxscx,jtanx,jprex,jradx,jloggx)

       call stripvar(nvarx,varidentx,varparamx,nprox,nvar,varident,
     1  varparam,nxx,xnx)

       print*,'gsetradPT - variables to be updated from .pre'
       print*,'xlatx = ',xlatx
       print*,'nvarx = ',nvarx
       do i=1,nvarx
        print*,'varidentx',(varidentx(i,j),j=1,3)
        print*,'varparamx',(varparamx(i,j),j=1,mparam)
       enddo

       xflag=1
       call subprofretg(xflag,runname,ispace,iscat,gasgiant,xlat,
     1  nvarx,varidentx,varparamx,nxx,xnx,jprex,ncont,flagh2p,xmapx,
     2  ierrx)


       do ivarx=1,nvarx
        if(varidentx(ivarx,1).eq.777)hcorrx=xnx(jtanx)
        if(varidentx(ivarx,1).eq.999)tsurf=xnx(jsurfx)
        if(varidentx(ivarx,1).eq.443)then
C         ********* Power law cross-section spectrum  **********

		 icont=1
         call get_xsecA(runname,nmode,nwave,wave,xsec)
         
         call file(runname,runname,'xsc')
         open(9,file=runname,status='unknown')
         write(9,*)ncont
         xsc(1,icont)=1.0
         ssa(1,icont)=1.0
         write(9,*)wave(1),xsc(1,icont)
         write(9,*)ssa(1,icont)
         do i=2,nwave
          xsc(i,icont)=(wave(i)/wave(1))
     1     **(varparamx(ivarx,1))
          do j=1,ncont
          	ssa(i,j)=1.0
          enddo
          write(9,*)wave(i),(xsc(i,j),j=1,ncont)
          write(9,*)(ssa(i,j),j=1,ncont)        
         enddo
         close(9)

      endif
             
        if(varidentx(ivarx,1).eq.442)then
C         ********* Power law cross-section spectrum  **********

		 icont=1
         call get_xsecA(runname,nmode,nwave,wave,xsec)
         
         call file(runname,runname,'xsc')
         open(9,file=runname,status='unknown')
         write(9,*)ncont
         xsc(1,icont)=1.0
         ssa(1,icont)=1.0
         write(9,*)wave(1),xsc(1,icont)
         write(9,*)ssa(1,icont)
         do i=2,nwave
          xsc(i,icont)=(wave(i)/wave(1))
     1     **(varparamx(ivarx,1))
          do j=1,ncont
          	ssa(i,j)=1.0
          enddo
          write(9,*)wave(i),(xsc(i,j),j=1,ncont)
          write(9,*)(ssa(i,j),j=1,ncont)        
         enddo
         close(9)

      endif
       if(varident(ivar,1).eq.441)then
C         ********* Power law cross-section spectrum  **********
         
		 icont=1
		 
         call get_xsecA(runname,nmode,nwave,wave,xsc)
         
         call file(runname,runname,'xsc')
         open(9,file=runname,status='unknown')
         write(9,*)ncont
         xsc(1,icont)=1.0
         ssa(1,icont)=1.0
         write(9,*)wave(1),xsc(1,icont)
         write(9,*)ssa(1,icont)
         do i=2,nwave
          xsc(i,icont)=(wave(i)/wave(1))
     1     **(varparam(ivar,1))
          do j=1,ncont
          	ssa(i,j)=1.0
          enddo
          write(9,*)wave(i),(xsc(i,j),j=1,ncont)
          write(9,*)(ssa(i,j),j=1,ncont)        
         enddo
         close(9)

      endif

      
      if(varident(ivar,1).eq.440)then
C       *****Calculate cross sections from refindex files****        

           iwave=ispace
           if(iwave.eq.0)iwave=2

           imode=1

           r0 = 10**varparam(ivar,1)
           v0 = varparam(ivar,2)  
           buffer='refindex.dat' 

           open(15,file=buffer,status='old')

           inorm=1
           iscat=2

           read(15,*)np1
           read(15,*)iwave
           read(15,*)lambda0
           do i=1,np1
            read(15,*),wave(i),nreal(i),nimag(i)
            srefind(i,1)=nreal(i)
            srefind(i,2)=nimag(i)
           enddo
           close(15)

           minlam=wave(1)

           parm(1)=r0
           parm(2)=v0
           parm(3)=(1. - 3 * parm(2))/parm(2)

           rs(1)=0.015*minlam
           rs(2)=0.
           rs(3)=rs(1)

           call modmakephase(iwave,imode,inorm,iscat,
     1   parm,rs,srefind,runname,lambda0)


       endif
       enddo

      endif
      
      
      do ivar=1,nvar
      if(varident(ivar,1).eq.443)then
C         ********* Power law cross-section spectrum  **********

		 icont=1
		 
         call get_xsecA(runname,nmode,nwave,wave,xsc)
         
         call file(runname,runname,'xsc')
         open(9,file=runname,status='unknown')
         write(9,*)ncont
         xsc(1,icont)=1.0
         ssa(1,icont)=1.0
         write(9,*)wave(1),xsc(1,icont)
         write(9,*)ssa(1,icont)
         do i=2,nwave
          xsc(i,icont)=(wave(i)/wave(1))
     1     **(varparam(ivar,1))
          do j=1,ncont
          	ssa(i,j)=1.0
          enddo
          write(9,*)wave(i),(xsc(i,j),j=1,ncont)
          write(9,*)(ssa(i,j),j=1,ncont)        
         enddo
         close(9)

      endif
        if(varident(ivar,1).eq.442)then
C         ********* Power law cross-section spectrum  **********
        

		 icont=1
		 
         call get_xsecA(runname,nmode,nwave,wave,xsc)
         
         call file(runname,runname,'xsc')
         open(9,file=runname,status='unknown')
         write(9,*)ncont
         xsc(1,icont)=1.0
         ssa(1,icont)=1.0
         write(9,*)wave(1),xsc(1,icont)
         write(9,*)ssa(1,icont)
         do i=2,nwave
          xsc(i,icont)=(wave(i)/wave(1))
     1     **(varparam(ivar,1))
          do j=1,ncont
          	ssa(i,j)=1.0
          enddo
          write(9,*)wave(i),(xsc(i,j),j=1,ncont)
          write(9,*)(ssa(i,j),j=1,ncont)        
         enddo
         close(9)

      endif
              if(varident(ivar,1).eq.441)then
C         ********* Power law cross-section spectrum  **********

		 icont=1
		 
         call get_xsecA(runname,nmode,nwave,wave,xsc)
         
         call file(runname,runname,'xsc')
         open(9,file=runname,status='unknown')
         write(9,*)ncont
         xsc(1,icont)=1.0
         ssa(1,icont)=1.0
         write(9,*)wave(1),xsc(1,icont)
         write(9,*)ssa(1,icont)
         do i=2,nwave
          xsc(i,icont)=(wave(i)/wave(1))
     1     **(varparam(ivar,1))
          do j=1,ncont
          	ssa(i,j)=1.0
          enddo
          write(9,*)wave(i),(xsc(i,j),j=1,ncont)
          write(9,*)(ssa(i,j),j=1,ncont)        
         enddo
         close(9)

      endif
      
      if(varident(ivar,1).eq.440)then
C       *****Calculate cross sections from refindex files****            

           iwave=ispace
           if(iwave.eq.0)iwave=2

           imode=1

           r0 = 10**varparam(ivar,1)
           v0 = varparam(ivar,2)     
           buffer='refindex.dat' 
           open(15,file=buffer,status='old')



           
           inorm=1
           iscat=2
           read(15,*)np1
           read(15,*)iwave
           read(15,*)lambda0
           do i=1,np1
            read(15,*),wave(i),nreal(i),nimag(i)
            srefind(i,1)=nreal(i)
            srefind(i,2)=nimag(i)
           enddo
           close(15)

           minlam=wave(1)
           
           parm(1)=r0
           parm(2)=v0
           parm(3)=(1. - 3 * parm(2))/parm(2)

           rs(1)=0.015*minlam
           rs(2)=0.
           rs(3)=rs(1)

           print*,'Calling modmakephase'

           call modmakephase(iwave,imode,inorm,iscat,
     1   parm,rs,srefind,runname,lambda0)

C           np=2

       endif
      enddo
      

      call gwritepatPT(runname,iscat,nconv,vconv,fwhm,layht,
     2 nlayer,laytyp,layint,flagh2p)

      return

      end
    
